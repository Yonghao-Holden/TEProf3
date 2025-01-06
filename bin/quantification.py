# TEProf3: TE-derived Promoter Finder 3
# Author: Holden Liang
# de novo assembly script
# -----------------------------------------------------------------------------
# Functions related to process gtf files post assembly
import argparse
import time
import subprocess # https://geekflare.com/python-run-bash/
import os
os.environ['NUMEXPR_MAX_THREADS'] = '4'
os.environ['NUMEXPR_NUM_THREADS'] = '2'
import numexpr as ne 
import multiprocess as mp
from multiprocessing import Pool
from itertools import repeat
import sys
import pandas as pd
from tqdm import tqdm
import glob
import numpy as np
import ast
from natsort import natsort_keygen
import concurrent.futures
from functools import partial
import pysam
import json

import misc
import filter_transcripts

def quantification(input_dataset, flags):
	quantification_mode = flags.quanmode
	filter_mode = flags.filtermode
	quantification_thread = flags.quansamplenumber

	## prepare gene annotation for quantification for each sample
	## TE.gtf, gencode (transcripts with the exact same intron structure as the TE.gtf will be removed)
	misc.print_time("start quantification, preparing gene annotation for quantification (remove duplciated transcripts in gencode) and combine genecode with TE trasncript consensus")

	# mode 1: quantification with only short reads
	if quantification_mode == "1":	
		# step 0. prepare quantification
		prepare_quantification_gene_annotation(flags)
		input_bam_SR_files=[]
		input_gtf_files=[]
		input_SJ_files = []
		samples=[]
		for sample in input_dataset:
			samples.append(sample)
			input_bam_SR_files.append(input_dataset[sample]['bam_SR'])
			input_SJ_files.append(input_dataset[sample]['SJ'])
			if flags.guided == "no":
				input_gtf_files.append(input_dataset[sample]['gtf'])
			else:
				input_gtf_files.append(input_dataset[sample]['bam_SR'].replace("bam/","assembled/").replace("bam","gtf"))

		# step 1. run stringtie
		misc.print_time("start stringtie for quantification")
		commands = []
		for input_bam_SR_file, input_gtf_file in zip(input_bam_SR_files, input_gtf_files):
			commands.append("samtools view -q 255 -h "+input_bam_SR_file+"| stringtie - -o "+input_gtf_file.replace("gtf","quantification.gtf")+" -e -b "+input_gtf_file.replace("gtf","stats")+" -p 2 -m 100 -c 1 -G assembled/TE_transcript_consensus.filtered.with_gene_annotation.gtf")

		with open("assembled/quantification_command.txt","w") as output_file:
			for command in commands:
				output_file.write(command+"\n")

		with mp.Pool(quantification_thread) as pool:
			useless = pool.map(misc.run_bash_command, commands)

		if flags.keepintermediate == False:
			misc.run_bash_command("rm assembled/TE_transcript_consensus.filtered.with_gene_annotation.gtf")

		# step 2. collect data after stringtie quantification
		misc.print_time("collect stringtie quantification output and concatenate result")
		replacement_dict = {}
		for sample, input_gtf_file in zip(samples, input_gtf_files):
			replacement_dict[input_gtf_file.replace("assembled/","").replace("gtf","stats")] = sample
		final_panda = SR_quantification_collect_data(flags)
		final_panda['sample'] = final_panda['stringtie_sample'].replace(replacement_dict)
		final_panda.drop_duplicates(subset=['transcript_id', 'stringtie_sample'], inplace=True)

		# step 3. generate raw transcript count
		if flags.quannoprepde == False:
			misc.print_time("generate raw transcript count using prepDE.py3")
			with open("assembled/sample_list_for_prepDE.txt", "w") as sample_list_file:
				for sample, input_gtf_file in zip(samples, input_gtf_files):
					sample_list_file.write('\t'.join([sample, input_gtf_file.replace("gtf","quantification.gtf")])+"\n")
			misc.run_bash_command("prepDE.py3 -i assembled/sample_list_for_prepDE.txt --length "+flags.quanreadlength)
			misc.run_bash_command("mv gene_count_matrix.csv assembled/; mv transcript_count_matrix.csv assembled/")
			prepDE_pandas = pd.read_table("assembled/transcript_count_matrix.csv", header=0, sep=",")
			prepDE_pandas = pd.melt(prepDE_pandas, id_vars='transcript_id', var_name='sample', value_name='raw_count')

			quantification_pandas = pd.merge(final_panda, prepDE_pandas, how="left", on=["transcript_id", "sample"])
		else:
			quantification_pandas = final_panda

		# step 4. save the quantification result
		if flags.quannogene == False:
			quantification_pandas.to_csv("assembled/teprof3_output_quantification.tsv", sep="\t", index=False)

		# step 5. check SJ support for novel splice junctions
		information_pandas = pd.read_table("assembled/teprof3_output_filter_transcript_TE_transcript_consensus.tsv")
		quantification_pandas = quantification_pandas[quantification_pandas["transcript_id"].isin(information_pandas["transcript_id"])]
		with open("assembled/TE_transcript_consensus.filtered.filter_tss.txt","w") as tss_file:
			useless = information_pandas.apply(filter_transcripts.SR_filter_step2_prepare, tss_file=tss_file, axis=1)

		misc.print_time("Check SJ support in each sample")
		with mp.Pool(quantification_thread) as pool:
			perfect_chimeric_mate_count_panda = pd.concat(pool.starmap(intersect_each_SJ_with_TACO, zip(samples, input_SJ_files, repeat(flags), repeat(information_pandas))), ignore_index=True)
			
		quantification_pandas = pd.merge(quantification_pandas, perfect_chimeric_mate_count_panda, how="left", on=["transcript_id","sample"])
		quantification_pandas = quantification_pandas.fillna(0)
		quantification_pandas.to_csv("assembled/teprof3_output_quantification.TE.tsv", sep="\t", index=False)

		if flags.keepintermediate == False:
			misc.run_bash_command("rm assembled/TE_transcript_consensus.filtered.filter_tss.txt ")	

		# add intron chain support info from LR to the quantification table
		if filter_mode == "2":
			misc.print_time("collect LR intron chain support for each transcript in each sample")
			quantification_pandas = pd.merge(quantification_pandas, information_pandas[["transcript_id","transcript_chr","transcript_coordinates","transcript_total_num_of_exon","transcript_start","transcript_stop","transcript_strand"]],on="transcript_id")

			split_dfs = np.array_split(quantification_pandas, quantification_thread)
			with mp.Pool(quantification_thread) as pool:
				quantification_pandas = pd.concat(pool.starmap(parallel_find_intron_chain_support, zip(split_dfs, repeat(flags), repeat(input_dataset))), ignore_index=True)

			quantification_pandas = quantification_pandas.drop(columns=["transcript_chr","transcript_coordinates","transcript_total_num_of_exon","transcript_start","transcript_stop","transcript_strand"])
			quantification_pandas.to_csv("assembled/teprof3_output_quantification.TE.tsv", sep="\t", index=False)

	# mode 2: quantification with only SJ files
	elif quantification_mode == "2":
		# step 0. prepare quantification
		input_SJ_files = []
		samples=[]
		for sample in input_dataset:
			samples.append(sample)
			input_SJ_files.append(input_dataset[sample]['SJ'])

		# step 1. check SJ support for novel splice junctions
		information_pandas = pd.read_table("assembled/teprof3_output_filter_transcript_TE_transcript_consensus.tsv")

		with open("assembled/TE_transcript_consensus.filtered.filter_tss.txt","w") as tss_file:
			useless = information_pandas.apply(filter_transcripts.SR_filter_step2_prepare, tss_file=tss_file, axis=1)

		misc.print_time("Check SJ support in each sample")
		with mp.Pool(quantification_thread) as pool:
			perfect_chimeric_mate_count_panda = pd.concat(pool.starmap(intersect_each_SJ_with_TACO, zip(samples, input_SJ_files, repeat(flags), repeat(information_pandas))), ignore_index=True)
			
		quantification_pandas = perfect_chimeric_mate_count_panda
		quantification_pandas = quantification_pandas.fillna(0)
		quantification_pandas.to_csv("assembled/teprof3_output_quantification.TE.tsv", sep="\t", index=False)

		if flags.keepintermediate == False:
			misc.run_bash_command("rm assembled/TE_transcript_consensus.filtered.filter_tss.txt ")
	else:
		misc.print_time("Quantification mode other than 1 and 2 is not available for now.")
		exit()

	return(input_dataset)
	
def parallel_find_intron_chain_support(input_pandas, flags, input_dataset):
	input_pandas["LR_intron_chain_support"] = input_pandas.apply(find_intron_chain_support, axis=1, flags=flags, input_dataset=input_dataset)
	return(input_pandas)

def find_intron_chain_support(input_row, flags, input_dataset):
	longread_sj_tolerance = flags.longreadsjtolerance
	longread_mono_tolerance = flags.longreadmtolerance

	input_bam_file = pysam.AlignmentFile(input_dataset[str(input_row["sample"])]['bam_LR'], "rb")
	transcript_chr = input_row["transcript_chr"]
	transcript_coordinates = ast.literal_eval(input_row["transcript_coordinates"])
	good_alignment_count = 0
	## 1. for mono-exonic transcripts
	if int(input_row["transcript_total_num_of_exon"]) == 1:
		for alignment in input_bam_file.fetch(transcript_chr,int(transcript_coordinates[0][0]),int(transcript_coordinates[0][1])):
			if alignment.mapping_quality>0:
				alignment_start = alignment.reference_start
				alignment_stop = alignment.reference_end
				if abs(input_row["transcript_start"] - alignment_start) <= longread_mono_tolerance & abs(input_row["transcript_stop"] - alignment_stop) <= longread_mono_tolerance:
					good_alignment_count+=1
		return(good_alignment_count)
	## 2. for multi-exonic transcripts
	elif int(input_row["transcript_total_num_of_exon"]) > 1:
		transcript_intron_coordinates = []
		for i in range(len(transcript_coordinates)-1):
			transcript_intron_coordinates.append((transcript_coordinates[i][1],transcript_coordinates[i+1][0]))
		if input_row["transcript_strand"] == "-":
			transcript_coordinates = transcript_coordinates[::-1]
		for alignment in input_bam_file.fetch(transcript_chr,int(transcript_coordinates[0][0]),int(transcript_coordinates[0][1])):
			if alignment.mapping_quality>0:
				alignment_intron_coordinates = filter_transcripts.get_intron_coordinates(alignment)
				if len(alignment_intron_coordinates) == len(transcript_intron_coordinates):
					matched_intron_count = 0
					for alignment_intron_coordinate, transcript_intron_coordinate in zip(alignment_intron_coordinates, transcript_intron_coordinates):
						transcript_intron_start = transcript_intron_coordinate[0]
						transcript_intron_stop = transcript_intron_coordinate[1]
						alignment_intron_start = alignment_intron_coordinate[0]
						alignment_intron_stop = alignment_intron_coordinate[1]
						if abs(transcript_intron_start-alignment_intron_start)<=longread_sj_tolerance and abs(transcript_intron_stop-alignment_intron_stop)<=longread_sj_tolerance:
							matched_intron_count+=1
					if matched_intron_count == len(transcript_intron_coordinates):
						good_alignment_count += 1
		return(good_alignment_count)


def SR_quantification_collect_data(flags):
	quantification_thread_for_concatenating = flags.quansamplenumbercon
	input_data_files = glob.glob("assembled/*stats/t_data.ctab")
	information_panda = pd.read_table("assembled/teprof3_output_filter_transcript_TE_transcript_consensus.tsv")
	gene_name_info_panda = information_panda[["transcript_id", "gene_name"]]

	all_input_panda_list = []
	for i in range(0, len(input_data_files), quantification_thread_for_concatenating):
		if i+quantification_thread_for_concatenating <= len(input_data_files):
			subset_input_data_files = input_data_files[i:i+quantification_thread_for_concatenating]
		else:
			subset_input_data_files = input_data_files[i:]
		input_panda_list = []
		for input_data_file in subset_input_data_files:
			misc.print_time("process "+input_data_file)
			sample_name = input_data_file.replace("assembled/","").replace("/t_data.ctab","")
			input_panda_temp = pd.read_table(input_data_file, header=0)
			input_panda_temp["sample"] = sample_name
			input_panda_temp["tpm"] = input_panda_temp['FPKM']/(input_panda_temp['FPKM'].sum())*1000000
			if flags.quannogene == True:
				input_panda_temp = input_panda_temp[input_panda_temp["t_name"].isin(information_panda["transcript_id"])]
			input_panda_list.append(input_panda_temp)
		input_panda = pd.concat(input_panda_list, ignore_index=True)
		input_panda_list=[]
		misc.print_time("collect all stringtie output")
		input_panda.columns = ["stringtie_index","stringtie_chr","stringtie_strand","stringtie_start","stringtie_stop","transcript_id","stringtie_number_of_exon","stringtie_length","stringtie_gene_id","stringtie_gene_name","stringtie_cov","stringtie_fpkm","stringtie_sample","stringtie_tpm"]
		all_input_panda = pd.merge(input_panda, gene_name_info_panda, how="left", on="transcript_id")
		input_panda=[]
		misc.print_time("get gene name info")
		all_input_panda["final_gene_name"] = all_input_panda.apply(find_gene_name, axis=1)
		misc.print_time("find gene name info")
		all_input_panda = all_input_panda[["transcript_id","stringtie_tpm","stringtie_sample", "final_gene_name"]]
		all_input_panda.columns = ["transcript_id","stringtie_tpm","stringtie_sample", "gene_name"]
	
		sum_of_expression = all_input_panda[["stringtie_sample", "gene_name", "stringtie_tpm"]].groupby(["stringtie_sample","gene_name"]).sum()
		misc.print_time("groupby to get total expression of each gene")
		sum_of_expression = sum_of_expression.reset_index()
		sum_of_expression.columns = ["stringtie_sample","gene_name","sum_tpm"]
		all_input_panda = pd.merge(all_input_panda, sum_of_expression, how="left", on=["stringtie_sample", "gene_name"])
		all_input_panda["percentage_of_expression"] = all_input_panda.apply(get_percentage_of_expression, axis=1)
		misc.print_time("get percentage of gene expression")
	
		all_input_panda = all_input_panda.sort_values(by=['transcript_id', "stringtie_sample"], ascending = True, key=natsort_keygen())
		all_input_panda = all_input_panda.reset_index(drop=True)

		all_input_panda_list.append(all_input_panda)

	return(pd.concat(all_input_panda_list, ignore_index=True))


def get_percentage_of_expression(input_row):
	if input_row["stringtie_tpm"] != 0:
		return(round(input_row["stringtie_tpm"]/input_row["sum_tpm"]*100, 2))
	elif input_row["stringtie_tpm"] == 0:
		return(0)

def find_gene_name(input_row):
	transcript_id = input_row["transcript_id"]
	if "ENST" in transcript_id:
		return(input_row["stringtie_gene_name"])
	elif "TU" in transcript_id:
		if input_row["gene_name"] != ".":
			return(input_row["gene_name"])
		elif input_row["gene_name"] == ".":
			return(transcript_id)

def prepare_quantification_gene_annotation(flags):
	## remove duplicated transcripts in gene annotation
	gene_annotation_file = os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/gene_annotation.gtf"
	with open(os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/gene_intron_coordinates_annotation.json", "r") as json_file:
		gene_annotation_intron_dict = json.load(json_file)

	taco_intron_chain = []
	with open("assembled/TE_transcript_consensus.filtered.intron_chain.txt", "r") as input_file:
		for line in input_file:
			taco_intron_chain.append(line.strip("\n"))

	removed_gene_annotation_transcript_count = 0
	with open("assembled/duplicated_transcript_ids_in_gene_annotation.txt","w") as output_file:
		for gene_transcript_id in gene_annotation_intron_dict:
			if gene_annotation_intron_dict[gene_transcript_id][0] + gene_annotation_intron_dict[gene_transcript_id][1] in taco_intron_chain:
				removed_gene_annotation_transcript_count += 1
				output_file.write("transcript_id \""+gene_transcript_id+"\"\n")

	misc.run_bash_command("grep -F -v -f assembled/duplicated_transcript_ids_in_gene_annotation.txt "+gene_annotation_file+"> assembled/deduplicated_transcript_ids_in_gene_annotation.gtf")
	misc.run_bash_command("cat assembled/deduplicated_transcript_ids_in_gene_annotation.gtf assembled/TE_transcript_consensus.filtered.gtf > assembled/TE_transcript_consensus.filtered.with_gene_annotation.gtf")
	if not flags.keepintermediate:
		misc.run_bash_command("rm assembled/TE_transcript_consensus.filtered.intron_chain.txt")
		misc.run_bash_command("rm assembled/duplicated_transcript_ids_in_gene_annotation.txt")
		misc.run_bash_command("rm assembled/deduplicated_transcript_ids_in_gene_annotation.gtf")

def intersect_each_SJ_with_TACO(sample, input_SJ_file, flags, information_pandas):
	misc.run_bash_command("grep -v chrM "+input_SJ_file+"|grep -v KI | grep -v GL | grep -v Un | grep -v H | grep -v CMV > "+input_SJ_file+".cleanup")
	misc.run_bash_command("bedtools intersect -a assembled/TE_transcript_consensus.filtered.filter_tss.txt -b "+input_SJ_file+".cleanup"+" -wa -wb >"+input_SJ_file+".intersected_with_taco.txt")
	chimeric_mate_count_panda, perfect_chimeric_mate_count_panda = filter_transcripts.SR_filter_step2_count_chimeric_mate(input_SJ_file+".intersected_with_taco.txt", flags)
	perfect_chimeric_mate_count_panda = pd.merge(perfect_chimeric_mate_count_panda, chimeric_mate_count_panda[["transcript_id", "SJ_uniqlymapped_read_upstream"]], how="left", on=["transcript_id"])
	if flags.keepintermediate == False:
		misc.run_bash_command("rm "+input_SJ_file+".cleanup")
	perfect_chimeric_mate_count_panda = perfect_chimeric_mate_count_panda[["transcript_id","perfect_SJ_uniqlymapped_read_downstream","SJ_uniqlymapped_read_upstream"]]
	perfect_chimeric_mate_count_panda = pd.merge(information_pandas.drop_duplicates(subset='transcript_id')["transcript_id"], perfect_chimeric_mate_count_panda, how="left", on=["transcript_id"])
	perfect_chimeric_mate_count_panda = perfect_chimeric_mate_count_panda.fillna(0)
	perfect_chimeric_mate_count_panda["sample"] = sample
	perfect_chimeric_mate_count_panda = perfect_chimeric_mate_count_panda[["transcript_id","sample","perfect_SJ_uniqlymapped_read_downstream","SJ_uniqlymapped_read_upstream"]]
	return(perfect_chimeric_mate_count_panda)
	
	








































