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
from concurrent.futures import ProcessPoolExecutor as Pool
import pysam

import process_assemble
import misc

def filter_transcripts_with_read_info(input_dataset, flags):
	''' filter them using read information before running mega assembly with TACO '''
	filter_mode = flags.filtermode
	filter_samplenumber = flags.filtersamplenumber

	## part 1. it's all short read here [filter_mode=1]
	if filter_mode == "1":
		misc.print_time("Filter mode is 1, filtering with short read information will be performed")
		misc.print_time("Filter transcripts with read information only using short read information")
	
		input_bam_SR_files = []
		input_gtf_files = []
		input_SJ_files = []
		for sample in input_dataset:
			input_bam_SR_files.append(input_dataset[sample]['bam_SR'])
			input_gtf_files.append(input_dataset[sample]['gtf'])
			input_SJ_files.append(input_dataset[sample]['SJ'])
	
		with mp.Pool(filter_samplenumber) as pool: # number of samples to be processed together
			pool.starmap(SR_filter_transcripts, zip(input_bam_SR_files, input_gtf_files, input_SJ_files, repeat(flags)))

	## part 2. if LR data is provided, use both short and long read data for filtering [filter_mode=2]
	elif filter_mode == "2":
		misc.print_time("Filter mode is 2, filtering with short read information and finding SJ support from LR data will be performed")
		
		input_bam_LR_files = []
		for sample in input_dataset:
			input_bam_LR_files.append(input_dataset[sample]['bam_LR'])
		if len(input_bam_LR_files) == 0:
			misc.print_time("No LR bam file is provided. Please check your input manifest file")
			exit()

		input_bam_SR_files = []
		input_gtf_files = []
		input_SJ_files = []
		for sample in input_dataset:
			input_bam_SR_files.append(input_dataset[sample]['bam_SR'])
			input_gtf_files.append(input_dataset[sample]['gtf'])
			input_SJ_files.append(input_dataset[sample]['SJ'])
		with mp.Pool(filter_samplenumber) as pool: # number of samples to be processed together
			pool.starmap(SR_filter_transcripts, zip(input_bam_SR_files, input_gtf_files, input_SJ_files, repeat(flags)))

		with mp.Pool(filter_samplenumber) as pool: # number of samples to be processed together
			pool.starmap(LR_filter_transcripts, zip(input_bam_LR_files, input_gtf_files, repeat(flags)))

	## part 4. intersect with HERV annotation
	#if os.path.isfile(os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/herv_annotation.txt") is True:
	#	misc.print_time("Intersect with transcripts with HERV annotation (using annotation from Telescope)")
	#	with mp.Pool(filter_samplenumber) as pool: # number of samples to be processed together
	#		pool.starmap(intersect_with_herv, zip(input_gtf_files, repeat(flags)))
#
	### part 5. intersect with full length LINE1 annotation
	#if os.path.isfile(os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/LINE1.bed") is True:
	#	misc.print_time("Intersect with transcripts with full length LINE1 annotation (using annotation from L1Base 2)")
	#	with mp.Pool(filter_samplenumber) as pool: # number of samples to be processed together
	#		pool.starmap(intersect_with_line1, zip(input_gtf_files, repeat(flags)))

	return(input_dataset)

def SR_filter_transcripts(input_bam_SR_file, input_gtf_file, input_SJ_file, flags):
	filter_thread = flags.filterthread # filter_thread = 20
	intron_retention = flags.filterintronretention # intron_retention=3
	filter_annotated = flags.filterannotated # filter_annotated = True
	filter_monoexon = flags.filtermonoexon # filter_monoexon = False
	filter_monoexon_tpm = flags.filtermonoexontpm # filter_monoexon = 2

	misc.print_time("Start filtering for "+input_gtf_file+", remove annotated transcripts: "+str(filter_annotated)+", remove mono-exonic transcripts: "+str(filter_monoexon))
	information_file = "assembled/teprof3_output_transcript_info_"+input_gtf_file.replace("assembled/","").replace("gtf","tsv") # input_gtf_file="assembled/Trimmed_mRNA_JHU006_NEG_BRep2_R1.fastqAligned.sortedByCoord.out.gtf" information_file="assembled/teprof3_output_process_assemble_Trimmed_mRNA_JHU006_NEG_BRep2_R1.fastqAligned.sortedByCoord.out.tsv"
	transcript_TE_panda = pd.read_table(information_file) # 8206

	keep_intermediate_files = flags.keepintermediate
	if not keep_intermediate_files:
		misc.run_bash_command("rm "+information_file)

	## step 1.1 mono-exonic speicific-filterings
	if filter_monoexon == True:
		transcript_TE_panda = transcript_TE_panda[transcript_TE_panda["transcript_total_num_of_exon"]>1]
	elif filter_monoexon == False:
		transcript_TE_panda = transcript_TE_panda[~((transcript_TE_panda["transcript_total_num_of_exon"]==1)&(transcript_TE_panda["transcript_TPM"]<filter_monoexon_tpm))]
	with open(input_gtf_file+".stats.txt", "a") as stats_file:	
		stats_file.write("\t".join([input_gtf_file, "# of TE-derived transcripts passing filter step 1 (remove lowly expressed mono-exonic)", str(len(transcript_TE_panda.index))])+"\n")	

	## step 2. chimeric mate filter
	## check if there're mate splicing from exon 1 to exon 2 of TE-derived transcript and less mate splicing into exon 1 of TE-derived transcripts (chimeric mate using SJ.out.tab file)
	## stringtie sometimes splits a transcript into 2 and call one transcript as 2 artifically. This is specially real when we have a TE in the exon. 
	## Because of low mappability in the TE region, less reads are mapped to the middle of the exon where the TE is
	## so stringtie would thought this is not an exon and break the transcript into half.
	## If the latter half of the transcript happens to have the fake tss in TE, 
	## it doesn't mean this transcript starts from TE, but simply because we have an exonic TE and the transcripts happens got split there by stringtie.
	misc.print_time("Start filtering (chimeric mate) for "+input_gtf_file)
	with open(input_gtf_file.replace("gtf","filter_tss.txt"),"w") as tss_file:
		useless = transcript_TE_panda.apply(SR_filter_step2_prepare, tss_file=tss_file, axis=1)

	## clean up input_SJ_file
	misc.run_bash_command("grep -v chrM "+input_SJ_file+"|grep -v KI | grep -v GL | grep -v Un | grep -v H | grep -v CMV > "+input_SJ_file+".cleanup")
	## run bedtools intersect
	useless = misc.run_bash_command("bedtools intersect -a "+input_gtf_file.replace("gtf","filter_tss.txt")+" -b "+input_SJ_file+".cleanup"+" -wa -wb >"+input_gtf_file.replace("gtf","filter_tss_SJ.txt"))
	chimeric_mate_count_panda, perfect_chimeric_mate_count_panda = SR_filter_step2_count_chimeric_mate(input_gtf_file.replace("gtf","filter_tss_SJ.txt"), flags)
	if not keep_intermediate_files:
		misc.run_bash_command("rm "+input_gtf_file.replace("gtf","filter_tss.txt"))	

	transcript_TE_panda = pd.merge(transcript_TE_panda, chimeric_mate_count_panda, how="left", on="transcript_id")
	transcript_TE_panda = pd.merge(transcript_TE_panda, perfect_chimeric_mate_count_panda, how="left", on="transcript_id")
	transcript_TE_panda = transcript_TE_panda.fillna(0)

	## step 3. count intron retention
	misc.print_time("Start filtering (intron retention) for "+input_gtf_file)
	exon_gene_file = input_gtf_file+".exon.gene_exon.txt"
	intron_retention_count_panda = SR_filter_step4_find_intron_retention_transcripts(exon_gene_file, transcript_TE_panda, flags)
	transcript_TE_panda = pd.merge(transcript_TE_panda, intron_retention_count_panda, how="left", on="transcript_id")
	transcript_TE_panda = transcript_TE_panda.fillna(0)

	## step 4. start filter transcripts based on information collected above
	misc.print_time("Start filtering using collected information for "+input_gtf_file)
	with open(input_gtf_file+".stats.txt", "a") as stats_file:
		## step 5.1 chimeric mate filter
		transcript_TE_panda["chimeric_mate_filter"] = transcript_TE_panda.apply(SR_filter_step5_check_chimeric_mate, flags=flags, axis=1)
		stats_file.write("\t".join([input_gtf_file, "# of TE-derived transcripts passing filter step 2 (chimeric mate)", str(len(transcript_TE_panda[transcript_TE_panda["chimeric_mate_filter"]=="keep"].index))])+"\n")
		## step 5.2 chimeric read filter
		#transcript_TE_panda["chimeric_read_filter"] = transcript_TE_panda.apply(SR_filter_step5_check_chimeric_read, flags=flags, axis=1)
		## step 5.3 intron retention filter
		transcript_TE_panda["intron_retention_filter"] = transcript_TE_panda["number_of_exons_it_overlaps_with"].apply(lambda x: "keep" if x<intron_retention else "no")
		stats_file.write("\t".join([input_gtf_file, "# of TE-derived transcripts passing filter step 3 (intron retention)", str(len(transcript_TE_panda[(transcript_TE_panda["chimeric_mate_filter"]=="keep") & (transcript_TE_panda["intron_retention_filter"]=="keep")].index))])+"\n")
		## step 5.4 check if this transcript is annotated (1. tss is in exon1 of an annotated transcript; 2. exon1 of TE-derived transcript overlaps with exon1 of an annotated transcript)
		transcript_TE_panda['annotation_exon1'] = transcript_TE_panda.apply(SR_filter_step5_check_annotated_or_not, axis=1)
		stats_file.write("\t".join([input_gtf_file, "# of TE-derived transcripts passing filter step 4 (annotation filter)", str(len(transcript_TE_panda[(transcript_TE_panda["chimeric_mate_filter"]=="keep") & (transcript_TE_panda["intron_retention_filter"]=="keep") & (transcript_TE_panda["annotation_exon1"]!="yes")].index))])+"\n")
		## step 5.5 combine all filter information
		transcript_TE_panda['filter_result'] = transcript_TE_panda.apply(SR_filter_step5_combine_all_info, flags=flags, axis=1) 
		stats_file.write("\t".join([input_gtf_file, "# of TE-derived transcripts passing all filters", str(len(transcript_TE_panda[transcript_TE_panda["filter_result"]=="keep"].index))])+"\n")
		## step 5.7 output transcripts information	
		transcript_TE_panda.to_csv("assembled/teprof3_output_filter_transcript_"+input_gtf_file.replace("assembled/","").replace("gtf","tsv"), sep="\t", index=False)
		## step 5.8 check how many of them are annotated (with the exact same intron chain)
		stats_file.write("\t".join([input_gtf_file, "# of annotated TE-derived transcripts (same intron chain)", str(len(transcript_TE_panda[(transcript_TE_panda["filter_result"]=="keep")&(transcript_TE_panda["annotation_whole_intron_chain"]=="yes")].index))])+"\n")
		stats_file.write("\t".join([input_gtf_file, "# of non-annotated TE-derived transcripts (same intron chain)", str(len(transcript_TE_panda[(transcript_TE_panda["filter_result"]=="keep")&(transcript_TE_panda["annotation_whole_intron_chain"]=="no")].index))])+"\n")
		stats_file.write("\t".join([input_gtf_file, "# of mono-exonic TE-derived transcripts", str(len(transcript_TE_panda[(transcript_TE_panda["filter_result"]=="keep")&(transcript_TE_panda["annotation_whole_intron_chain"]=="mono_exonic")].index))])+"\n")

	## step 5. output gtf files after filtering for TACO mega assembly
	#for transcript_class in ["TE transcript","TE coding gene transcript","TE noncoding gene transcript","TE no gene transcript"]:
	#	with open(input_gtf_file.replace("gtf",transcript_class.replace(" ","_")+".gtf"), "w") as output_file:
	#		useless = transcript_TE_panda[transcript_TE_panda["transcript_class"]==transcript_class].apply(extract_filtered_info_for_output_gtf, output_file=output_file, axis=1)
	transcript_TE_panda=transcript_TE_panda[transcript_TE_panda["filter_result"]=="keep"]
	with open(input_gtf_file.replace("gtf","TE.gtf"), "w") as output_file:
		useless = transcript_TE_panda.apply(extract_filtered_info_for_output_gtf, output_file=output_file, axis=1)

def SR_filter_step2_prepare(input_row, tss_file):
	transcript_coordinates = ast.literal_eval(input_row["transcript_coordinates"])
	if input_row["transcript_strand"] == "-":
		transcript_coordinates = transcript_coordinates[::-1]
	if input_row["transcript_total_num_of_exon"] > 1:
		## 1. for TE transcripts with tss in exonic-TE:
		if "exon" in input_row["tss_genic_info"] and input_row["transcript_strand"] == "+": # input_row = transcript_TE_panda[transcript_TE_panda["transcript_id"]=="STRG.173.1"]
			start = min(int(input_row["tss_intronexon_start"]), int(transcript_coordinates[0][0]))
			stop = int(transcript_coordinates[0][1])
			tss_file.write('\t'.join([input_row["transcript_chr"], str(start-20), str(start), "upstream", input_row["transcript_strand"], input_row["transcript_id"], str(input_row["tss_intronexon_start"]), str(input_row["tss_intronexon_stop"]), str(transcript_coordinates[0][0]), str(transcript_coordinates[0][1]), ".", ".",str(input_row["transcript_total_num_of_exon"])])+"\n")
			tss_file.write('\t'.join([input_row["transcript_chr"], str(stop), str(stop+20), "downstream", input_row["transcript_strand"], input_row["transcript_id"], str(input_row["tss_intronexon_start"]), str(input_row["tss_intronexon_stop"]), str(transcript_coordinates[0][0]), str(transcript_coordinates[0][1]), str(transcript_coordinates[1][0]), str(transcript_coordinates[1][1]),str(input_row["transcript_total_num_of_exon"])])+"\n")
		elif "exon" in input_row["tss_genic_info"] and input_row["transcript_strand"] == "-":
			start = int(transcript_coordinates[0][0])
			stop = max(int(input_row["tss_intronexon_stop"]), int(transcript_coordinates[0][1]))
			tss_file.write('\t'.join([input_row["transcript_chr"], str(stop), str(stop+20), "upstream", input_row["transcript_strand"], input_row["transcript_id"], str(input_row["tss_intronexon_start"]), str(input_row["tss_intronexon_stop"]), str(transcript_coordinates[0][0]), str(transcript_coordinates[0][1]), ".", ".",str(input_row["transcript_total_num_of_exon"])])+"\n")
			tss_file.write('\t'.join([input_row["transcript_chr"], str(start-20), str(start), "downstream", input_row["transcript_strand"], input_row["transcript_id"], str(input_row["tss_intronexon_start"]), str(input_row["tss_intronexon_stop"]), str(transcript_coordinates[0][0]), str(transcript_coordinates[0][1]), str(transcript_coordinates[1][0]), str(transcript_coordinates[1][1]),str(input_row["transcript_total_num_of_exon"])])+"\n")	
		## 2. for TE transcripts with tss in intronic-TE and intergenic-TE:
		elif "exon" not in input_row["tss_genic_info"] and input_row["transcript_strand"] == "+":
			tss_file.write('\t'.join([input_row["transcript_chr"], str(int(transcript_coordinates[0][0])-20), str(transcript_coordinates[0][0]), "upstream", input_row["transcript_strand"], input_row["transcript_id"], ".", ".", str(transcript_coordinates[0][0]), str(transcript_coordinates[0][1]), ".", ".",str(input_row["transcript_total_num_of_exon"])])+"\n")
			tss_file.write('\t'.join([input_row["transcript_chr"], str(transcript_coordinates[0][1]), str(int(transcript_coordinates[0][1])+20), "downstream", input_row["transcript_strand"], input_row["transcript_id"], ".", ".", str(transcript_coordinates[0][0]), str(transcript_coordinates[0][1]), str(transcript_coordinates[1][0]), str(transcript_coordinates[1][1]),str(input_row["transcript_total_num_of_exon"])])+"\n")
		elif "exon" not in input_row["tss_genic_info"] and input_row["transcript_strand"] == "-":
			tss_file.write('\t'.join([input_row["transcript_chr"], str(transcript_coordinates[0][1]), str(int(transcript_coordinates[0][1])+20), "upstream", input_row["transcript_strand"], input_row["transcript_id"], ".", ".", str(transcript_coordinates[0][0]), str(transcript_coordinates[0][1]), ".", ".",str(input_row["transcript_total_num_of_exon"])])+"\n")
			tss_file.write('\t'.join([input_row["transcript_chr"], str(int(transcript_coordinates[0][0])-20), str(transcript_coordinates[0][0]), "downstream", input_row["transcript_strand"], input_row["transcript_id"], ".", ".", str(transcript_coordinates[0][0]), str(transcript_coordinates[0][1]), str(transcript_coordinates[1][0]), str(transcript_coordinates[1][1]),str(input_row["transcript_total_num_of_exon"])])+"\n")
	elif input_row["transcript_total_num_of_exon"] == 1:
		## 1. for TE transcripts with tss in exonic-TE:
		if "exon" in input_row["tss_genic_info"] and input_row["transcript_strand"] == "+": # input_row = transcript_TE_panda[transcript_TE_panda["transcript_id"]=="STRG.173.1"]
			start = min(int(input_row["tss_intronexon_start"]), int(transcript_coordinates[0][0]))
			stop = int(transcript_coordinates[0][1])
			tss_file.write('\t'.join([input_row["transcript_chr"], str(start-20), str(start), "upstream", input_row["transcript_strand"], input_row["transcript_id"], str(input_row["tss_intronexon_start"]), str(input_row["tss_intronexon_stop"]), str(transcript_coordinates[0][0]), str(transcript_coordinates[0][1]), ".", ".",str(input_row["transcript_total_num_of_exon"])])+"\n")
			tss_file.write('\t'.join([input_row["transcript_chr"], str(stop), str(stop+20), "downstream", input_row["transcript_strand"], input_row["transcript_id"], str(input_row["tss_intronexon_start"]), str(input_row["tss_intronexon_stop"]), str(transcript_coordinates[0][0]), str(transcript_coordinates[0][1]), str(stop),"999999999999999",str(input_row["transcript_total_num_of_exon"])])+"\n")
		elif "exon" in input_row["tss_genic_info"] and input_row["transcript_strand"] == "-":
			start = int(transcript_coordinates[0][0])
			stop = max(int(input_row["tss_intronexon_stop"]), int(transcript_coordinates[0][1]))
			tss_file.write('\t'.join([input_row["transcript_chr"], str(stop), str(stop+20), "upstream", input_row["transcript_strand"], input_row["transcript_id"], str(input_row["tss_intronexon_start"]), str(input_row["tss_intronexon_stop"]), str(transcript_coordinates[0][0]), str(transcript_coordinates[0][1]), ".", ".",str(input_row["transcript_total_num_of_exon"])])+"\n")
			tss_file.write('\t'.join([input_row["transcript_chr"], str(start-20), str(start), "downstream", input_row["transcript_strand"], input_row["transcript_id"], str(input_row["tss_intronexon_start"]), str(input_row["tss_intronexon_stop"]), str(transcript_coordinates[0][0]), str(transcript_coordinates[0][1]),"0",str(start),str(input_row["transcript_total_num_of_exon"])])+"\n")	
		## 2. for TE transcripts with tss in intronic-TE and intergenic-TE:
		elif "exon" not in input_row["tss_genic_info"] and input_row["transcript_strand"] == "+":
			tss_file.write('\t'.join([input_row["transcript_chr"], str(int(transcript_coordinates[0][0])-20), str(transcript_coordinates[0][0]), "upstream", input_row["transcript_strand"], input_row["transcript_id"], ".", ".", str(transcript_coordinates[0][0]), str(transcript_coordinates[0][1]), ".", ".",str(input_row["transcript_total_num_of_exon"])])+"\n")
			tss_file.write('\t'.join([input_row["transcript_chr"], str(transcript_coordinates[0][1]), str(int(transcript_coordinates[0][1])+20), "downstream", input_row["transcript_strand"], input_row["transcript_id"], ".", ".", str(transcript_coordinates[0][0]), str(transcript_coordinates[0][1]), str(transcript_coordinates[0][1]),"999999999999999",str(input_row["transcript_total_num_of_exon"])])+"\n")
		elif "exon" not in input_row["tss_genic_info"] and input_row["transcript_strand"] == "-":
			tss_file.write('\t'.join([input_row["transcript_chr"], str(transcript_coordinates[0][1]), str(int(transcript_coordinates[0][1])+20), "upstream", input_row["transcript_strand"], input_row["transcript_id"], ".", ".", str(transcript_coordinates[0][0]), str(transcript_coordinates[0][1]), ".", ".",str(input_row["transcript_total_num_of_exon"])])+"\n")
			tss_file.write('\t'.join([input_row["transcript_chr"], str(int(transcript_coordinates[0][0])-20), str(transcript_coordinates[0][0]), "downstream", input_row["transcript_strand"], input_row["transcript_id"], ".", ".", str(transcript_coordinates[0][0]), str(transcript_coordinates[0][1]), "0", str(transcript_coordinates[0][0]),str(input_row["transcript_total_num_of_exon"])])+"\n")		


def SR_filter_step2_count_chimeric_mate(input_file, flags): # input_file = "Trimmed_mRNA_JHU006_NEG_BRep2_R1.fastqAligned.sortedByCoord.out.filter_tss_SJ.txt"
	input_panda = pd.read_table(input_file, names=["transcript_chr", "start", "stop", "direction", "strand", "transcript_id", "gene_exon_start", "gene_exon_stop", "transcript_exon_start", "transcript_exon_stop", "transcript_next_exon_start", "transcript_next_exon_stop", "transcript_total_num_of_exon", "SJ_chr", "SJ_start", "SJ_stop", "SJ_strand", "SJ_intron_motif", "SJ_annotation", "SJ_uniqlymapped_read", "SJ_multimapped_read", "SJ_overhang"])
	if flags.keepintermediate == False:
		misc.run_bash_command("rm "+input_file)
	## find all splice junctions within the exon
	input_panda_all = input_panda.copy()
	input_panda_all["goodbad"] = input_panda_all.apply(SR_filter_step2_find_overlapped_splice_junctions,axis=1)
	input_panda_all = input_panda_all[input_panda_all["goodbad"]=="good"]
	input_panda_all = input_panda_all.groupby(["transcript_id","direction"]).agg({"SJ_uniqlymapped_read":'sum', "SJ_multimapped_read":'sum'}).reset_index()
	input_panda_all = input_panda_all.pivot(index='transcript_id', columns='direction', values=['SJ_uniqlymapped_read', 'SJ_multimapped_read'])
	input_panda_all.columns = input_panda_all.columns.map('_'.join)
	if len(input_panda_all.columns) != 5:
		if "SJ_uniqlymapped_read_upstream" not in input_panda_all.columns:
			input_panda_all["SJ_uniqlymapped_read_upstream"] = 0
			with open("teprof3_warning_please_read.txt", "a") as warning_file:
				warning_file.write(input_file+" has no SJ_uniqlymapped_read_upstream detected. This is not a good sign. It usually means there's something wrong with this sample. Please check."+"\n")
		if "SJ_multimapped_read_upstream" not in input_panda_all.columns:
			input_panda_all["SJ_multimapped_read_upstream"] = 0
			with open("teprof3_warning_please_read.txt", "a") as warning_file:
				warning_file.write(input_file+" has no SJ_multimapped_read_upstream detected. This is not a good sign. It usually means there's something wrong with this sample. Please check."+"\n")
	input_panda_all = input_panda_all.reset_index()
	input_panda_all = input_panda_all.fillna(0)
	input_panda_all = input_panda_all[["transcript_id","SJ_uniqlymapped_read_upstream","SJ_uniqlymapped_read_downstream","SJ_multimapped_read_upstream", "SJ_multimapped_read_downstream"]]
	## find perfect splice junctions within the exon
	input_panda_perfect = input_panda.copy()
	input_panda_perfect["goodbad"] = input_panda_perfect.apply(SR_filter_step2_find_perfectly_overlapped_splice_junctions,axis=1)
	input_panda_perfect = input_panda_perfect[input_panda_perfect["goodbad"]=="good"]
	input_panda_perfect = input_panda_perfect.groupby(["transcript_id","direction"]).agg({"SJ_uniqlymapped_read":'sum', "SJ_multimapped_read":'sum'}).reset_index()
	input_panda_perfect = input_panda_perfect.pivot(index='transcript_id', columns='direction', values=['SJ_uniqlymapped_read', 'SJ_multimapped_read'])
	input_panda_perfect.columns = input_panda_perfect.columns.map('_'.join)
	input_panda_perfect.columns = ["perfect_"+x for x in input_panda_perfect.columns]
	if len(input_panda_perfect.columns) != 5:
		if "perfect_SJ_uniqlymapped_read_upstream" not in input_panda_perfect.columns:
			input_panda_perfect["perfect_SJ_uniqlymapped_read_upstream"] = 0
			with open("teprof3_warning_please_read.txt", "a") as warning_file:
				warning_file.write(input_file+" has no perfect_SJ_uniqlymapped_read_upstream detected. This is not a good sign. It usually means there's something wrong with this sample. Please check."+"\n")
		if "perfect_SJ_multimapped_read_upstream" not in input_panda_perfect.columns:
			input_panda_perfect["perfect_SJ_multimapped_read_upstream"] = 0
			with open("teprof3_warning_please_read.txt", "a") as warning_file:
				warning_file.write(input_file+" has no perfect_SJ_multimapped_read_upstream detected. This is not a good sign. It usually means there's something wrong with this sample. Please check."+"\n")
	input_panda_perfect = input_panda_perfect.reset_index()
	input_panda_perfect = input_panda_perfect.fillna(0)
	input_panda_perfect = input_panda_perfect[["transcript_id","perfect_SJ_uniqlymapped_read_upstream","perfect_SJ_uniqlymapped_read_downstream","perfect_SJ_multimapped_read_upstream", "perfect_SJ_multimapped_read_downstream"]]
	return(input_panda_all, input_panda_perfect)

def SR_filter_step2_find_overlapped_splice_junctions(input_row):
	SJ_start = int(input_row["SJ_start"]) - int(input_row["SJ_overhang"])
	SJ_stop = int(input_row["SJ_stop"]) + int(input_row["SJ_overhang"])
	if input_row["gene_exon_start"] != ".":
		if input_row["strand"] == "+":
			exon_start = min(int(input_row["transcript_exon_start"]), int(input_row["gene_exon_start"]))
			exon_stop = int(input_row["transcript_exon_stop"])
		elif input_row["strand"] == "-":
			exon_start = int(input_row["transcript_exon_start"])
			exon_stop = max(int(input_row["transcript_exon_stop"]), int(input_row["gene_exon_stop"]))
	elif input_row["gene_exon_start"] == ".":
		exon_stop = int(input_row["transcript_exon_stop"])
		exon_start = int(input_row["transcript_exon_start"])
	if input_row["direction"] == "downstream":
		next_exon_start = int(input_row["transcript_next_exon_start"])
		next_exon_stop = int(input_row["transcript_next_exon_stop"])
	if input_row["direction"] == "upstream" and input_row["strand"] == "+" and SJ_start < exon_start and exon_start <= SJ_stop and SJ_stop <= exon_stop:
		return("good")
	elif input_row["direction"] == "downstream" and input_row["strand"] == "+" and exon_start <= SJ_start and SJ_start <= exon_stop and next_exon_start <= SJ_stop and SJ_stop <= next_exon_stop:
		return("good")
	elif input_row["direction"] == "upstream" and input_row["strand"] == "-" and exon_start <= SJ_start and SJ_start <= exon_stop and exon_stop < SJ_stop:
		return("good")
	elif input_row["direction"] == "downstream" and input_row["strand"] == "-" and next_exon_start <= SJ_start and SJ_start <= next_exon_stop and exon_start <= SJ_stop and SJ_stop <= exon_stop:
		return("good")
	else:
		return("bad")

def SR_filter_step2_find_perfectly_overlapped_splice_junctions(input_row): # input_row=input_panda.iloc[0]
	SJ_start = int(input_row["SJ_start"])-1
	SJ_stop = int(input_row["SJ_stop"])+1
	if input_row["gene_exon_start"] != ".":
		if input_row["strand"] == "+":
			exon_start = min(int(input_row["transcript_exon_start"]), int(input_row["gene_exon_start"]))
			exon_stop = int(input_row["transcript_exon_stop"])
		elif input_row["strand"] == "-":
			exon_start = int(input_row["transcript_exon_start"])
			exon_stop = max(int(input_row["transcript_exon_stop"]), int(input_row["gene_exon_stop"]))
	elif input_row["gene_exon_start"] == ".":
		exon_stop = int(input_row["transcript_exon_stop"])
		exon_start = int(input_row["transcript_exon_start"])
	if input_row["direction"] == "downstream":
		next_exon_start = int(input_row["transcript_next_exon_start"])
		next_exon_stop = int(input_row["transcript_next_exon_stop"])
	if input_row["direction"] == "upstream" and input_row["strand"] == "+" and SJ_stop == exon_start:
		return("good")
	elif input_row["direction"] == "downstream" and input_row["strand"] == "+" and SJ_start == exon_stop and next_exon_start == SJ_stop:
		return("good")
	elif input_row["direction"] == "upstream" and input_row["strand"] == "-" and SJ_start == exon_stop:
		return("good")
	elif input_row["direction"] == "downstream" and input_row["strand"] == "-" and SJ_start == next_exon_stop and exon_start == SJ_stop:
		return("good")
	else:
		return("bad")

def SR_filter_step4_find_intron_retention_transcripts(exon_gene_file, transcript_TE_panda, flags):
	exon_gene_panda = pd.read_table(exon_gene_file, names=['exon_chr','exon_start','exon_stop','transcript_id','transcript_cov','transcript_strand','transcript_id_2','exon_number','gene_chr','gene_start','gene_stop','gene_name','useless1','gene_strand','gene_transcript_id','gene_exon_number','gene_id','gene_type','gene_transcript_type','gene_transcript_name','useless2'])
	
	keep_intermediate_files = flags.keepintermediate
	if not keep_intermediate_files:	
		subprocess.run("rm "+exon_gene_file,shell=True, executable='/bin/bash')

	exon_gene_panda = exon_gene_panda[['transcript_id','gene_transcript_id','exon_number','gene_exon_number']]
	exon_gene_panda = exon_gene_panda[exon_gene_panda['transcript_id'].isin(list(transcript_TE_panda.transcript_id))]
	exon_gene_panda = exon_gene_panda.groupby(['transcript_id', 'exon_number', 'gene_transcript_id']).size().reset_index(name='count')
	exon_gene_panda = exon_gene_panda[exon_gene_panda["count"]!=1]
	exon_gene_panda = exon_gene_panda[exon_gene_panda.groupby(['transcript_id'])['count'].transform(max) == exon_gene_panda['count']]
	exon_gene_panda = exon_gene_panda[exon_gene_panda['transcript_id'].duplicated() == False]
	exon_gene_panda = exon_gene_panda.drop("gene_transcript_id", axis=1)

	exon_gene_panda.columns = ["transcript_id","exon_number_with_intron_retention","number_of_exons_it_overlaps_with"]
	return(exon_gene_panda)


def SR_filter_step5_check_chimeric_mate(input_row, flags):
	ratio_of_upstream_and_downstream_mates = flags.filterratio
	number_of_downstream_mate = flags.filterdownstreammate
	filternochimericfilter = flags.filternochimericfilter
	if filternochimericfilter == True:
		return("keep")
	if input_row["transcript_total_num_of_exon"] == 1:
		return("keep")
	elif ((input_row["SJ_uniqlymapped_read_upstream"]+input_row["SJ_multimapped_read_upstream"]) < input_row["perfect_SJ_uniqlymapped_read_downstream"] * ratio_of_upstream_and_downstream_mates) and (input_row["perfect_SJ_uniqlymapped_read_downstream"]>=number_of_downstream_mate):
		return("keep")
	else:
		return("no")

def SR_filter_step5_check_annotated_or_not(input_row):
	if input_row["exon1_gene_exon_number"]=="exon_1":
		if str(input_row["transcript_strand"]) == "+":	 	
			if int(input_row["TE_start"]) <= input_row["exon1_gene_start"] <= input_row["TE_stop"]:
				return("yes, but both are TE-derived")
			else:
				return("yes")
		elif str(input_row["transcript_strand"]) == "-":
			if int(input_row["TE_start"]) <= input_row["exon1_gene_stop"] <= input_row["TE_stop"]:
				return("yes, but both are TE-derived")
			else:
				return("yes")
	else:
		return("no")

def SR_filter_step5_check_chimeric_read(input_row, flags):
	filter_thread = flags.filterthread # filter_thread = 20
	return("keep")

def SR_filter_step5_combine_all_info(input_row, flags):
	filter_annotated = flags.filterannotated # filter_annotated=True
	keep_flag = 0
	if (input_row["chimeric_mate_filter"]=="keep") and (input_row["intron_retention_filter"]=="keep"):
	#if (input_row["chimeric_mate_filter"]=="keep") and (input_row["chimeric_read_filter"]=="keep") and (input_row["intron_retention_filter"]=="keep"):
		keep_flag = 1
	if filter_annotated:
		if keep_flag==1 and (input_row['annotation_exon1']!="yes"):
			keep_flag = 1
		else:
			keep_flag = 0

	if keep_flag==1:
		return("keep")
	elif keep_flag==0:
		return("no")

def extract_filtered_info_for_output_gtf(input_row, output_file):
	chromosome = input_row["transcript_chr"]
	start = str(input_row["transcript_start"])
	stop = str(input_row["transcript_stop"])
	strand = input_row["transcript_strand"]
	gene_id = "gene_id \""+input_row["transcript_gene_id"]+"\""
	transcript_id = "transcript_id \""+input_row["transcript_id"]+"\""
	gene_name = "gene_name \""+input_row["gene_name"]+"\""
	if gene_name==".":
		gene_name = transcript_id
	coverage = "cov \""+str(input_row["transcript_cov"])+"\""
	fpkm = "FPKM \""+str(input_row["transcript_FPKM"])+"\""
	tpm = "TPM \""+str(input_row["transcript_TPM"])+"\""	
	info = [chromosome,"StringTie","transcript",start,stop,"1000",strand,".","; ".join([gene_id,transcript_id,coverage,fpkm,tpm])]
	output_file.write("\t".join(info)+";\n")
	transcript_coordinates = ast.literal_eval(input_row["transcript_coordinates"])
	for exon in range(len(transcript_coordinates)):
		start = str(transcript_coordinates[exon][0])
		stop = str(transcript_coordinates[exon][1])
		exon_number = "exon_number \""+str(exon+1)+"\""
		info = [chromosome,"StringTie","exon",start,stop,"1000",strand,".","; ".join([gene_id,transcript_id,exon_number])]
		output_file.write("\t".join(info)+";\n")
	return(None)



def LR_filter_transcripts(input_bam_LR_file, input_gtf_file, flags):
	# it will find splicing junction support from LR data (PS. only check the SJ between TE and gene, ie things that aren't annotated)
	misc.print_time("Start finding splicing junction support for "+input_gtf_file+" using long read data")
	information_file = "assembled/teprof3_output_filter_transcript_"+input_gtf_file.replace("assembled/","").replace("gtf","tsv") # input_gtf_file="assembled/Trimmed_mRNA_JHU006_NEG_BRep2_R1.fastqAligned.sortedByCoord.out.gtf" information_file="assembled/teprof3_output_process_assemble_Trimmed_mRNA_JHU006_NEG_BRep2_R1.fastqAligned.sortedByCoord.out.tsv"
	transcript_TE_panda = pd.read_table(information_file) # 8206
	input_bam_file = pysam.AlignmentFile(input_bam_LR_file, "rb") # input_bam_LR_file="/bar/xqu/scratch/p53/data/A549_pilot/long_read/processed/All_Samples/aligned/aln_nanofilt_headcrop_A549KO1.sortedByCoord.out.bam"
	
	## step 1. find support from LR that support the novel splice junctions that were not annotated in gene annotation
	transcript_TE_panda["LR_novel_SJ_support"] = transcript_TE_panda.apply(LR_filter_step1_find_SJ, flags=flags, input_bam_file=input_bam_file, axis=1)
	transcript_TE_panda["LR_novel_SJ_support_summary"] = transcript_TE_panda["LR_novel_SJ_support"].apply(lambda x: "yes" if "\'y\'" in str(x) else "no" if "\'n\'" in str(x) else "na")

	## step 2. find support from LR that support the whole transcript (intron structure)
	transcript_TE_panda["LR_intron_chain_support"] = transcript_TE_panda.apply(LR_filter_step2_find_transcript, flags=flags, input_bam_file=input_bam_file, axis=1)
	transcript_TE_panda.to_csv("assembled/teprof3_output_filter_transcript_"+input_gtf_file.replace("assembled/","").replace("gtf","tsv"), sep="\t", index=False)
	with open(input_gtf_file+".stats.txt", "a") as stats_file:
		stats_file.write("\t".join([input_gtf_file, "# of TE-derived transcripts with novel splice junction", str(len(transcript_TE_panda[(transcript_TE_panda["filter_result"]=="keep")&(transcript_TE_panda["LR_novel_SJ_support"]!="it's a mono-exonic transcript")&(transcript_TE_panda["LR_novel_SJ_support"]!="all introns are in gene annotation")].index))])+"\n")
		stats_file.write("\t".join([input_gtf_file, "# of TE-derived transcripts with LR support (novel splice junction)", str(len(transcript_TE_panda[(transcript_TE_panda["filter_result"]=="keep")&(transcript_TE_panda["LR_novel_SJ_support_summary"]=="yes")].index))])+"\n")
		stats_file.write("\t".join([input_gtf_file, "# of TE-derived transcripts with LR support (intron chain)", str(len(transcript_TE_panda[(transcript_TE_panda["filter_result"]=="keep")&(transcript_TE_panda["LR_intron_chain_support"]>0)].index))])+"\n")
	return(None)

def LR_filter_step1_find_SJ(input_row, flags, input_bam_file):
	longread_sj_tolerance = flags.longreadsjtolerance
	if str(input_row["transcript_total_num_of_exon"]) == "1":
		return("it's a mono-exonic transcript")
	if input_row["transcript_exon_number_overlap_gene"] == "exon_1" and input_row["gene_exon_number"] != ".":
		return("all introns are in gene annotation")
	## splice junctions detected by LR
	transcript_chr = input_row["transcript_chr"]
	transcript_coordinates = ast.literal_eval(input_row["transcript_coordinates"])
	if input_row["transcript_strand"] == "+":
		first_exon = transcript_coordinates[0]
	elif input_row["transcript_strand"] == "-":
		first_exon = transcript_coordinates[-1]
	introns = input_bam_file.find_introns(alignment for alignment in input_bam_file.fetch(transcript_chr,int(first_exon[0]),int(first_exon[1])) if alignment.mapping_quality>0) # introns = input_bam_file.find_introns(read for read in input_bam_file.fetch("chr7",27167131,27178233) if read.mapping_quality>0)
	## find SJ that needs to be validated
	if input_row["gene_name"] != ".":
		exon_that_splice_into_something = int(input_row["transcript_exon_number_overlap_gene"].split("_")[1])
	elif input_row["gene_name"] == ".":
		exon_that_splice_into_something = int(input_row["transcript_total_num_of_exon"])
	if input_row["transcript_strand"] == "+":
		interested_transcript_coordinates = transcript_coordinates[0:exon_that_splice_into_something]
	elif input_row["transcript_strand"] == "-":
		interested_transcript_coordinates = transcript_coordinates[::-1][0:exon_that_splice_into_something][::-1]
	interested_SJs = []
	for i in range(len(interested_transcript_coordinates)-1):
		interested_SJs += [int(interested_transcript_coordinates[i][1]),int(interested_transcript_coordinates[i+1][0])]
	interested_SJs = sorted(interested_SJs)
	## start checking each SJ
	output = []
	for i in range(0,len(interested_SJs),2):
		interested_SJ_start = interested_SJs[i]
		interested_SJ_stop = interested_SJs[i+1]
		found_LR_support_flag=0
		for intron in introns:
			intron_start = int(intron[0])
			intron_stop = int(intron[1])
			if abs(intron_start-interested_SJ_start)<=longread_sj_tolerance and abs(intron_stop-interested_SJ_stop)<=longread_sj_tolerance:
				output.append([interested_SJ_start,interested_SJ_stop,introns[intron],"y"])
				found_LR_support_flag=1
				break
		if found_LR_support_flag==0:
			output.append([interested_SJ_start,interested_SJ_stop,0,"n"])
	return(output)


def LR_filter_step2_find_transcript(input_row, flags, input_bam_file):
	longread_sj_tolerance = flags.longreadsjtolerance
	longread_mono_tolerance = flags.longreadmtolerance

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
				alignment_intron_coordinates = get_intron_coordinates(alignment)
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

def get_intron_coordinates(alignment):
	position = alignment.reference_start
	intron_coordinates = []
	for op, length in alignment.cigartuples:
		if op == 3:  # CREF_SKIP represents intron
			intron_start = position
			intron_end = position + length - 1
			intron_coordinates.append((intron_start, intron_end))
		# Update position based on CIGAR operation
		if op in (0, 2, 3, 7, 8):
			position += length
	return(intron_coordinates)


def intersect_with_herv(input_gtf_file, flags):
	# it will find splicing junction support from LR data (PS. only check the SJ between TE and gene, ie things that aren't annotated)
	misc.print_time("Start intersecting HERV annotation for "+input_gtf_file)
	information_file = "assembled/teprof3_output_filter_transcript_"+input_gtf_file.replace("assembled/","").replace("gtf","tsv") # input_gtf_file="assembled/Trimmed_mRNA_JHU006_NEG_BRep2_R1.fastqAligned.sortedByCoord.out.gtf" information_file="assembled/teprof3_output_process_assemble_Trimmed_mRNA_JHU006_NEG_BRep2_R1.fastqAligned.sortedByCoord.out.tsv"
	misc.run_bash_command("tail -n +2 "+information_file+">"+information_file+".temp")
	misc.run_bash_command("bedtools intersect -a "+information_file+".temp -b "+os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/herv_annotation.txt"+" -wa -wb -loj > "+information_file+".herv")
	misc.run_bash_command("header=`head -1 "+information_file+"`; echo -e \"$header\therv_chr\therv_start\therv_stop\therv_locus\therv_repFamily\therv_strand\" > "+information_file+".header")
	misc.run_bash_command("cat "+information_file+".header "+information_file+".herv > "+information_file)
	misc.run_bash_command("rm "+information_file+".temp")
	misc.run_bash_command("rm "+information_file+".herv")
	misc.run_bash_command("rm "+information_file+".header")
	return(None)


def intersect_with_line1(input_gtf_file, flags):
	# it will find splicing junction support from LR data (PS. only check the SJ between TE and gene, ie things that aren't annotated)
	misc.print_time("Start intersecting LINE1 annotation for "+input_gtf_file)
	information_file = "assembled/teprof3_output_filter_transcript_"+input_gtf_file.replace("assembled/","").replace("gtf","tsv") # input_gtf_file="assembled/Trimmed_mRNA_JHU006_NEG_BRep2_R1.fastqAligned.sortedByCoord.out.gtf" information_file="assembled/teprof3_output_process_assemble_Trimmed_mRNA_JHU006_NEG_BRep2_R1.fastqAligned.sortedByCoord.out.tsv"
	misc.run_bash_command("tail -n +2 "+information_file+">"+information_file+".temp")
	misc.run_bash_command("bedtools intersect -a "+information_file+".temp -b "+os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/LINE1.bed"+" -wa -wb -loj > "+information_file+".line1")
	misc.run_bash_command("header=`head -1 "+information_file+"`; echo -e \"$header\tline1_chr\tline1_start\tline1_stop\tline1_id\tline1_score1\tline1_strand\tline1_pos1\tline1_pos2\tline1_score2\" > "+information_file+".header")
	misc.run_bash_command("cat "+information_file+".header "+information_file+".line1 > "+information_file)
	misc.run_bash_command("rm "+information_file+".temp")
	misc.run_bash_command("rm "+information_file+".line1")
	misc.run_bash_command("rm "+information_file+".header")
	return(None)












