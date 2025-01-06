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
import json

import process_assemble
import filter_transcripts
import misc

def run_taco(input_dataset, flags):
	taco_thread = flags.tacothread # taco_thread="10"
	tpm_cutoff = min([flags.processtpm, flags.filtermonoexontpm]) # flags.processtpm=0.5, flags.filtermonoexontpm=1
	if tpm_cutoff == 0:
		tpm_cutoff = 0.0000001
	assemble_length = flags.assemblelength # 100
	misc.print_time("Start running TACO")
	ref_gtf_file = os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/gene_annotation.gtf"
	ref_genome_file = os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/ref_genome.fa"
	useless = subprocess.run("ls assembled/*TE.gtf > assembled/gtf_to_merge.txt", shell=True, executable='/bin/bash')
	taco_command = "taco_run --filter-min-length "+str(assemble_length)+" --verbose --gtf-expr-attr TPM --filter-min-expr "+str(tpm_cutoff)+" --ref-genome-fasta "+ref_genome_file+" --ref-gtf "+ref_gtf_file+" -o assembled/TACO_output -p "+taco_thread+" assembled/gtf_to_merge.txt &> debug/taco.log 2> debug/taco.err"

	with open("./assembled/taco_command.txt","w") as file:
		file.write(taco_command+"\n")
	misc.run_bash_command(taco_command)

	misc.run_bash_command("ln -s TACO_output/assembly.gtf assembled/TE_transcript_consensus.gtf")

	misc.print_time("Finished TACO, start process TACO data")

	return(input_dataset)

def process_taco(input_dataset, flags):
	intron_retention = flags.filterintronretention # intron_retention=3
	filter_annotated = flags.filterannotated # filter_annotated = True
	filter_monoexon = flags.filtermonoexon # filter_monoexon = False
	filter_monoexon_tpm = flags.filtermonoexontpm # filter_monoexon = 3

	## step 0. prepare gene_annotation_intron_dict
	misc.print_time("Prepare intron annotation")
	with open(os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/gene_intron_coordinates_annotation.json", "r") as json_file:
		gene_annotation_intron_dict = json.load(json_file)

	## step 1. convert each gtf files into tsv for easier handling
	if flags.guided != "no":
		misc.print_time("You are using guided mode. TEProf3 skipped steps before TACO and directly use the gtf file you provided as TACO output")
		if flags.guided != "teprof3_output_TE_transcript_consensus.gtf":
			misc.run_bash_command("mv "+flags.guided+" assembled/"+flags.guided)
			misc.run_bash_command("ln -s "+flags.guided+" assembled/TE_transcript_consensus.gtf")
		elif flags.guided == "teprof3_output_TE_transcript_consensus.gtf":
			misc.run_bash_command("mv "+flags.guided+" assembled/"+flags.guided+"_user_provided_for_guided_mode")
			misc.run_bash_command("ln -s "+flags.guided+"_user_provided_for_guided_mode"+" assembled/TE_transcript_consensus.gtf")

	misc.print_time("Extract exon and transcript information from taco gtf file")
	input_gtf_file = "assembled/TE_transcript_consensus.gtf"
	convert_taco_gtf_to_pandas(input_gtf_file)

	## step 2. intersect transcripts with repeatmasker and gene annotation (gencode) to classify transcripts
	input_gtf_exon_file = "assembled/TE_transcript_consensus.gtf.exon.txt"
	input_gtf_transcript_file = "assembled/TE_transcript_consensus.gtf.transcript.txt"
	input_gtf_tss_file = "assembled/TE_transcript_consensus.gtf.tss.txt"

	misc.print_time("Intersect transcript with TE annotation => to check if it's TE-derived")
	process_assemble.intersect_with_TE(input_gtf_transcript_file)

	misc.print_time("Intersect TSS with gene annotation => to check genic feature of TSS")
	process_assemble.intersect_with_gene_exon_intron(input_gtf_tss_file)

	misc.print_time("Intersect exons with gene annotation => to find the assocaited gene")
	process_assemble.intersect_with_gene_exon(input_gtf_exon_file)

	## step 3. based on the intersection result, classify transcripts into 4 categories
	misc.print_time("Classify transcripts into four categories")
	final_output_panda = categorize_taco_transcripts(input_gtf_file, flags, gene_annotation_intron_dict)

	## step 4. filter transcripts
	## step 4.1 count intron retention
	misc.print_time("Start filtering for "+input_gtf_file+", remove annotated transcripts: "+str(filter_annotated)+", remove mono-exonic transcripts: "+str(filter_monoexon))
	exon_gene_file = input_gtf_file+".exon.gene_exon.txt"
	intron_retention_count_panda = filter_transcripts.SR_filter_step4_find_intron_retention_transcripts(exon_gene_file, final_output_panda, flags)
	final_output_panda = pd.merge(final_output_panda, intron_retention_count_panda, how="left", on="transcript_id")
	final_output_panda = final_output_panda.fillna(0)
	## step 4.2 intron retention filter
	final_output_panda["intron_retention_filter"] = final_output_panda["number_of_exons_it_overlaps_with"].apply(lambda x: "keep" if x<intron_retention else "no")
	## step 4.3 check if this transcript is annotated (1. tss is in exon1 of an annotated transcript; 2. exon1 of TE-derived transcript overlaps with exon1 of an annotated transcript)
	final_output_panda['annotation_exon1'] = final_output_panda.apply(filter_transcripts.SR_filter_step5_check_annotated_or_not, axis=1)
	## step 4.4 combine all filter information
	final_output_panda['filter_result'] = final_output_panda.apply(combine_all_filter_info, flags=flags, axis=1) 
	## step 4.5 check for duplicated transcripts with the same intron chain and same TSS, but only different transcription end sites
	final_output_panda = label_duplicated_transcripts(final_output_panda)
	## step 4.6 save information
	final_output_panda.to_csv("assembled/teprof3_output_filter_transcript_"+input_gtf_file.replace("assembled/","").replace("gtf","tsv"), sep="\t", index=False)
	## step 4.7 intersect with HERV annotation
	#if os.path.isfile(os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/herv_annotation.txt") is True:
	#	misc.print_time("Intersect with transcripts with HERV annotation (using annotation from Telescope)")
	#	filter_transcripts.intersect_with_herv(input_gtf_file, flags)
	### step 4.8 intersect with LINE1 annotation
	#if os.path.isfile(os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/LINE1.bed") is True:
	#	misc.print_time("Intersect with transcripts with full length LINE1 annotation (using annotation from L1Base 2)")
	#	filter_transcripts.intersect_with_line1(input_gtf_file, flags)
	## step 4.9 keep good and unique candidates
	final_output_panda=final_output_panda[final_output_panda["filter_result"]=="keep"]
	final_output_panda=final_output_panda[final_output_panda['transcript_id'].duplicated()==False] # because sometime one transcript has multiple HERVs overlap with it.
	final_output_panda=final_output_panda[final_output_panda["duplicated_transcripts"]!="duplicated_transcript_nokeep"]

	## generate final gtf file
	with open(input_gtf_file.replace("gtf","filtered.gtf"), "w") as output_file:
		useless = final_output_panda.apply(extract_filtered_taco_info_for_output_gtf, output_file=output_file, axis=1)
	
	## save intron chain information for quantification:
	with open("assembled/TE_transcript_consensus.filtered.intron_chain.txt","w") as output_file:
		final_output_panda.apply(print_intron_coordinates, output_file=output_file, axis=1)

	return(input_dataset)


def convert_taco_gtf_to_pandas(input_gtf_file):
	""" convert gtf file to pandas dataframe for faster and easier downstream processing """
	exon_count = {}
	with open(input_gtf_file, "r") as gtf_file, open(input_gtf_file+".transcript.txt", "w") as transcript_file, open(input_gtf_file+".exon.txt", "w") as exon_file, open(input_gtf_file+".tss.txt","w") as tss_file:
		for line in tqdm(gtf_file, desc="Processing "+input_gtf_file):
			if line[0] != "#":
				entry = line.strip("\n").split("\t")
				detail_info = entry[8]
				strand = entry[6]
				gene_id = detail_info.split("gene_id \"")[1].split("\";")[0]
				#tss_id = detail_info.split("tss_id \"")[1].split("\";")[0]
				#locus_id = detail_info.split("locus_id \"")[1].split("\";")[0]
				transcript_id = detail_info.split("transcript_id \"")[1].split("\";")[0]
				#transcript_id = transcript_id+"."+locus_id+"."+tss_id
				if entry[2] == "transcript":
					transcript_file.write("\t".join([entry[0], entry[3], entry[4], transcript_id, ".", strand, gene_id])+'\n') # chr, start, stop, transcript_id, ".", strand, gene_id
					if strand == "+":
						tss_file.write("\t".join([entry[0], entry[3], str(int(entry[3])+1), transcript_id, ".", strand, gene_id])+'\n') # chr, start, stop, transcript_id, ".", strand, gene_id
					if strand == "-":
						tss_file.write("\t".join([entry[0], str(int(entry[4])-1), entry[4], transcript_id, ".", strand, gene_id])+'\n') # chr, start, stop, transcript_id, ".", strand, gene_id
					exon_count[transcript_id] = 1
				if entry[2] == "exon":
					if strand == "+" or strand == "-":
						exon_number = "exon_"+str(exon_count[transcript_id])
						exon_count[transcript_id] += 1
						exon_file.write("\t".join([entry[0], entry[3], entry[4], transcript_id, ".", strand, gene_id, exon_number])+'\n') # chr, start, stop, transcript_id, ".", strand, gene_id, exon_number			


def categorize_taco_transcripts(input_gtf_file, flags, gene_annotation_intron_dict):
	transcript_TE_file = input_gtf_file+".transcript.TE.txt" # input_gtf_file="guided_merged_taco.gtf"
	tss_gene_file = input_gtf_file+".tss.gene_exon_intron.txt"
	exon_gene_file = input_gtf_file+".exon.gene_exon.txt"

	misc.print_time("Read in intersected informations for "+input_gtf_file)
	transcript_TE_panda = pd.read_table(transcript_TE_file, names=['transcript_chr','transcript_start','transcript_stop','transcript_id','uesless1','transcript_strand','transcript_gene_id','TE_chr','TE_start','TE_stop','TE_subfamily','TE_sub','TE_strand','TE_family','TE_class'])
	tss_gene_panda = pd.read_table(tss_gene_file, names=['tss_chr','tss_start','tss_stop','tss_transcript_id','uesless1','tss_strand','tss_gene_id','gene_chr','gene_start','gene_stop','gene_gene_name','useless1','gene_strand','gene_transcript_id','gene_exon_intron_number','gene_gene_id','gene_gene_type','gene_transcript_type','gene_transcript_name','.'])
	exon_gene_panda = pd.read_table(exon_gene_file, names=['exon_chr','exon_start','exon_stop','transcript_id','uesless1','transcript_strand','transcript_gene_id','exon_number','gene_chr','gene_start','gene_stop','gene_name','useless1','gene_strand','gene_transcript_id','gene_exon_number','gene_id','gene_type','gene_transcript_type','gene_transcript_name','useless2'])

	keep_intermediate_files = flags.keepintermediate
	if not keep_intermediate_files:	
		subprocess.run("rm "+transcript_TE_file+" "+tss_gene_file, shell=True)

	##################
	## Part 1. process transcript_TE_panda to extract TE-derived transcripts
	##################
	misc.print_time("Process transcript_TE_panda to extract TE-derived transcripts for "+input_gtf_file)
	transcript_TE_panda['start_from_TE'] = transcript_TE_panda.apply(process_assemble.start_from_TE, axis=1)
	transcript_TE_panda['all_in_TE'] = transcript_TE_panda.apply(lambda x: "yes" if (x['TE_start'] <= x['transcript_start'] and x['transcript_stop'] <= x['TE_stop']) else "no", axis=1)
	transcript_TE_panda = transcript_TE_panda[transcript_TE_panda['start_from_TE']=="yes"]
	transcript_TE_panda = transcript_TE_panda[transcript_TE_panda['transcript_id'].duplicated() == False]
	transcript_TE_panda = transcript_TE_panda.reset_index(drop=True)

	TE_transcript_ids = set(list(transcript_TE_panda["transcript_id"]))
	
	##################
	## Part 2. process tss_gene_panda to get genic information of TSS
	##################
	misc.print_time("Process tss_gene_panda to get genic information of TSS for "+input_gtf_file)
	tss_gene_panda = tss_gene_panda[tss_gene_panda['tss_transcript_id'].isin(TE_transcript_ids)]
	tss_gene_panda = tss_gene_panda[tss_gene_panda["tss_strand"]!="."]
	tss_gene_panda["genic_info"] = tss_gene_panda["gene_exon_intron_number"].apply(lambda x: "intergenic" if (x == ".") else x)
	tss_gene_panda["genic_info_text"] = tss_gene_panda["genic_info"].apply(lambda x: x.split("_")[0] if ("_" in x) else x)
	tss_gene_panda["genic_info_text_rank"] = tss_gene_panda["genic_info"].apply(lambda x: 1 if ("exon" in x) else 2 if "intron" in x else 3)
	tss_gene_panda["genic_info_number"] = tss_gene_panda["genic_info"].apply(lambda x: int(x.split("_")[1]) if ("_" in x) else 0)
	tss_gene_panda["gene_transcript_type_number"] = tss_gene_panda["gene_transcript_type"].apply(lambda x: 1 if x=="coding" else 2 if x=="nonCoding" else 3 if x=="pseudo" else 4 if x=="problem" else 5)
	tss_gene_panda["distance_between_stringtie_tss_and_annotated_tss"] = tss_gene_panda.apply(process_assemble.get_distance_between_stringtie_tss_and_annotated_tss_overlap_tssandtranscript,axis=1)
	## pick coding first, and then exon over intron (eg exon1 over intron1), and then smallest exon_intron_number (eg exon1 over exon2), and then pick the annoated transcript that has the cloest TSS
	tss_gene_panda = tss_gene_panda.sort_values(by=['tss_transcript_id','gene_transcript_type_number','genic_info_text_rank','genic_info_number','distance_between_stringtie_tss_and_annotated_tss'])
	# tss_gene_panda.to_csv("/scratch/yliang/TEProf3/data/rna-seq_HNSCC/assembled/tss_gene_panda.tsv", sep="\t", index=False)
	tss_gene_panda = tss_gene_panda[['tss_transcript_id','gene_gene_name','gene_strand','gene_start','gene_stop','gene_transcript_id','genic_info','gene_gene_id','gene_gene_type','gene_transcript_type','gene_transcript_name']]
	tss_gene_panda.columns = ["transcript_id","tss_gene_name","tss_gene_strand",'tss_intronexon_start','tss_intronexon_stop',"tss_transcript_id","tss_genic_info","tss_gene_id","tss_gene_type","tss_transcript_type","tss_transcript_name"]
	# By default, for each set of duplicated values, the first occurrence is set on False and all others on True. so if a transcript overlaps with both exon 1 and exon 2 of an annotated gene, it will take save exon 1 instead of 2 as it's associated gene for tss.
	tss_gene_panda = tss_gene_panda[tss_gene_panda['transcript_id'].duplicated() == False]
	tss_gene_panda = tss_gene_panda.reset_index(drop=True)

	exon_coordinates_file = os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/gene_exon_coordinates_annotation.txt"
	exon_coordinates_panda = pd.read_table(exon_coordinates_file, names=["tss_transcript_id", "uesless", "tss_transcript_coordinates", "tss_transcript_start", "tss_transcript_stop"])
	exon_coordinates_panda = exon_coordinates_panda.drop("uesless", axis=1)
	tss_gene_panda = pd.merge(tss_gene_panda, exon_coordinates_panda, how="left", on="tss_transcript_id")

	##################
	## Part 3. process exon_gene_panda to get the spliced gene information (compare intron and intron chain) and annotation filter (compare the exon 1 of TE-derived transcript) 
	##################
	misc.print_time("Process exon_gene_panda to get splicing-gene information for "+input_gtf_file)
	exon_gene_panda = exon_gene_panda.drop(['useless1','useless2'], axis=1)
	exon_gene_panda = exon_gene_panda[exon_gene_panda['transcript_id'].isin(TE_transcript_ids)]

	exon_gene_panda['exon_number_clean'] = exon_gene_panda['exon_number'].apply(lambda x: int(x.split("_")[1]) if "_" in x else x)
	exon_gene_panda['gene_exon_number_clean'] = exon_gene_panda['gene_exon_number'].apply(lambda x: int(x.split("_")[1]) if "_" in x else 999)
	exon_gene_panda["gene_transcript_type_number"] = exon_gene_panda["gene_transcript_type"].apply(lambda x: 1 if x=="coding" else 2 if x=="nonCoding" else 3 if x=="pseudo" else 4 if x=="problem" else 5)

	exon_gene_panda_cleanup_list = []
	for transcript_id in tqdm(list(exon_gene_panda.transcript_id.unique()), desc="Collect splice gene target for "+input_gtf_file):
		##################
		## Part 3.1 correct exon numbers
		##################
		exon_gene_panda_temp = exon_gene_panda[exon_gene_panda['transcript_id']==transcript_id] # exon_gene_panda_temp = exon_gene_panda[exon_gene_panda['transcript_id']=="STRG.20.1"]
		exon_gene_panda_temp = exon_gene_panda_temp.reset_index(drop=True)
		exon_gene_panda_temp['transcript_total_num_of_exon'] = max(exon_gene_panda_temp['exon_number_clean'])
		if exon_gene_panda_temp.iloc[0]['transcript_strand'] == "+":
			exon_gene_panda_temp['corrected_exon_number_clean'] = exon_gene_panda_temp['exon_number_clean']
		elif exon_gene_panda_temp.iloc[0]['transcript_strand'] == "-":
			exon_gene_panda_temp['corrected_exon_number_clean'] = max(exon_gene_panda_temp['exon_number_clean']) - exon_gene_panda_temp['exon_number_clean'] + 1

		##################
		## Part 3.2 collect transcript coordinates
		##################
		## save the coordinate information for this TE-derived transcript
		## !! in the future, if you want to do exon coordiinate correction, do this per gene_transcript_id. Because different gene_transcript_id have different exon coordinate, 
		## !! this correction is performed corresponds to each gene_transcript_id. So, different gene_transcript_id will have different corrected_exon_coodinates
		## !! I actually don't recommend doing this as it's too complicated and increase the usage of memory and time signficantly. What this correction tries to solve is sequencing error, especially for long read
		## !! but this issue should not be our concern. 
		exon_gene_panda_temp = exon_gene_panda_temp.sort_values(by=['corrected_exon_number_clean'])
		exon_gene_panda_temp = exon_gene_panda_temp.reset_index(drop=True)
		transcript_coordinates = []
		exon_number = 1
		for index in exon_gene_panda_temp.index:
			if exon_gene_panda_temp.iloc[index]["corrected_exon_number_clean"] == exon_number:
				transcript_coordinates.append((exon_gene_panda_temp.iloc[index]["exon_start"], exon_gene_panda_temp.iloc[index]["exon_stop"]))
				exon_number+=1
		transcript_coordinates = sorted(transcript_coordinates, key=lambda x: x[0])

		##################
		## Part 3.3 find exon_1 overlapped exons to prepare for annotation filter
		##################
		## for annotation filter, find if exon_1 of TE transcript overlaps with any annotated exon_1
		## this step should be separated from finding the splice gene target
		exon_1_panda = exon_gene_panda_temp.copy()
		exon_1_panda = exon_1_panda[exon_1_panda["corrected_exon_number_clean"]==1]
		exon_1_panda['annotated_exon1_in_the_same_TE'] = exon_1_panda.apply(process_assemble.find_good_annotated_exon1, transcript_TE_panda=transcript_TE_panda[transcript_TE_panda["transcript_id"]==transcript_id], axis=1)
		exon_1_panda['distance_between_stringtie_exonedge_and_annotated_exonedge'] = exon_1_panda.apply(process_assemble.get_distance_between_stringtie_exonedge_and_annotated_exonedge_overlap_exonandgene, axis=1)
		exon_1_panda['overlap_between_stringtie_exon_and_annotated_exon'] = exon_1_panda.apply(process_assemble.get_overlap_between_stringtie_exon_and_annotated_exon, axis=1)
		exon_1_panda = exon_1_panda.sort_values(by=["annotated_exon1_in_the_same_TE", "gene_exon_number_clean", "distance_between_stringtie_exonedge_and_annotated_exonedge", "overlap_between_stringtie_exon_and_annotated_exon", "gene_transcript_type_number"], ascending=[False, True,True,False,True])
		exon_1_panda = exon_1_panda.reset_index(drop=True)
		exon_1_panda = exon_1_panda[["transcript_id","gene_transcript_id","gene_exon_number","gene_start","gene_stop","gene_transcript_type"]]
		exon_1_panda.columns = ["transcript_id","exon1_gene_transcript_id","exon1_gene_exon_number","exon1_gene_start","exon1_gene_stop","exon1_gene_transcript_type"]

		##################
		## Part 3.4 find splice gene target by comparing the intron chain and picking the transcript with the same intron chain
		##################
		compare_intron_chain_result = {}
		intron_coordinates = []
		if len(transcript_coordinates) == 1:
			exon_gene_panda_temp["intron_coordinates"] = "mono_exonic"
			exon_gene_panda_temp["position_of_the_first_overlap_intron_in_TE_transcript"] = exon_gene_panda_temp["gene_transcript_id"].apply(lambda x: 0 if x == exon_1_panda.iloc[0]["exon1_gene_transcript_id"] else 999)
			exon_gene_panda_temp["number_of_overlap_intron"] = exon_gene_panda_temp["gene_transcript_id"].apply(lambda x: 999 if x == exon_1_panda.iloc[0]["exon1_gene_transcript_id"] else 0)
			exon_gene_panda_temp["overlap_intron_versus_gene_all_intron"] = exon_gene_panda_temp["gene_transcript_id"].apply(lambda x: 999 if x == exon_1_panda.iloc[0]["exon1_gene_transcript_id"] else 0)
		elif len(transcript_coordinates) >= 2:
			for i in range(len(transcript_coordinates)-1):
				intron_coordinates.append((transcript_coordinates[i][1], transcript_coordinates[i+1][0]))
			exon_gene_panda_temp["intron_coordinates"] = [intron_coordinates for i in exon_gene_panda_temp.index]
			if exon_gene_panda_temp.iloc[0]['transcript_strand'] == "-":
				intron_coordinates = intron_coordinates[::-1]
			for gene_transcript_id in list(exon_gene_panda_temp.gene_transcript_id.unique()):
				if gene_transcript_id == ".":
					compare_intron_chain_result[gene_transcript_id] = [999, 0, 1]
				elif gene_transcript_id != ".":
					a_intron_overlap_flag = 0
					position_of_the_first_overlap_intron_in_TE_transcript = 999
					number_of_overlap_intron_count = 0
					if exon_gene_panda_temp.iloc[0]['transcript_strand'] == "+":
						gene_intron_coordinates = ast.literal_eval(gene_annotation_intron_dict[gene_transcript_id][1])
					elif exon_gene_panda_temp.iloc[0]['transcript_strand'] == "-":
						gene_intron_coordinates = ast.literal_eval(gene_annotation_intron_dict[gene_transcript_id][1])[::-1]
					for intron in intron_coordinates:
						if a_intron_overlap_flag == 0:
							if intron in gene_intron_coordinates:
								a_intron_overlap_flag = 1
								number_of_overlap_intron_count += 1
								position_of_overlap_intron_in_gene_transcript = gene_intron_coordinates.index(intron)
								position_of_the_first_overlap_intron_in_TE_transcript = intron_coordinates.index(intron)
						elif a_intron_overlap_flag == 1:
							if position_of_overlap_intron_in_gene_transcript == len(gene_intron_coordinates)-1:
								break
							elif intron == gene_intron_coordinates[position_of_overlap_intron_in_gene_transcript+1]:
								position_of_overlap_intron_in_gene_transcript+=1
								number_of_overlap_intron_count += 1	
							else:
								break
					compare_intron_chain_result[gene_transcript_id] = [position_of_the_first_overlap_intron_in_TE_transcript, number_of_overlap_intron_count, len(gene_intron_coordinates)]
			exon_gene_panda_temp["position_of_the_first_overlap_intron_in_TE_transcript"] = exon_gene_panda_temp["gene_transcript_id"].apply(lambda x: compare_intron_chain_result[x][0])
			exon_gene_panda_temp["number_of_overlap_intron"] = exon_gene_panda_temp["gene_transcript_id"].apply(lambda x: compare_intron_chain_result[x][1])
			exon_gene_panda_temp["overlap_intron_versus_gene_all_intron"] = exon_gene_panda_temp["gene_transcript_id"].apply(lambda x: compare_intron_chain_result[x][1]/compare_intron_chain_result[x][2] if compare_intron_chain_result[x][2]!= 0 else 0)

		exon_gene_panda_temp = exon_gene_panda_temp.sort_values(by=["position_of_the_first_overlap_intron_in_TE_transcript","number_of_overlap_intron","overlap_intron_versus_gene_all_intron","gene_transcript_type_number"], ascending=[True,False,False,True])
		exon_gene_panda_temp = exon_gene_panda_temp.reset_index(drop=True)
		exon_gene_panda_temp = exon_gene_panda_temp[["transcript_id","transcript_total_num_of_exon","intron_coordinates","corrected_exon_number_clean","gene_chr","gene_start","gene_stop","gene_name","gene_strand","gene_transcript_id","gene_exon_number","gene_id","gene_type","gene_transcript_type","gene_transcript_name"]]
		exon_gene_panda_temp["corrected_exon_number_clean"] = exon_gene_panda_temp["corrected_exon_number_clean"].apply(lambda x: "exon_"+str(int(x)))
		exon_gene_panda_temp.columns = ["transcript_id","transcript_total_num_of_exon","transcript_intron_coordinates","transcript_exon_number_overlap_gene","gene_chr","gene_start","gene_stop","gene_name","gene_strand","gene_transcript_id","gene_exon_number","gene_id","gene_type","gene_transcript_type","gene_transcript_name"]

		##################
		## Part 3.5 combine all collected information into one dataframe for one transcript
		##################
		exon_gene_panda_temp = pd.merge(pd.DataFrame(exon_gene_panda_temp.iloc[0]).transpose(), pd.DataFrame(exon_1_panda.iloc[0]).transpose(), on="transcript_id")
		exon_gene_panda_temp["transcript_coordinates"] = [transcript_coordinates for i in exon_gene_panda_temp.index]
		exon_gene_panda_temp = exon_gene_panda_temp[["transcript_id","exon1_gene_transcript_id","exon1_gene_exon_number","exon1_gene_start","exon1_gene_stop","exon1_gene_transcript_type","transcript_coordinates","transcript_total_num_of_exon","transcript_intron_coordinates","transcript_exon_number_overlap_gene","gene_chr","gene_start","gene_stop","gene_name","gene_strand","gene_transcript_id","gene_exon_number","gene_id","gene_type","gene_transcript_type","gene_transcript_name"]]
		exon_gene_panda_cleanup_list.append(exon_gene_panda_temp)

	exon_gene_panda_cleanup = pd.concat(exon_gene_panda_cleanup_list, ignore_index=True)
	exon_coordinates_panda.columns = ["gene_transcript_id", "gene_transcript_coordinates", "gene_transcript_start", "gene_transcript_stop"]
	exon_gene_panda_cleanup = pd.merge(exon_gene_panda_cleanup, exon_coordinates_panda, how="left", on="gene_transcript_id")

	##################
	## Part 4. combine all informations together for each sample about TE-derived transcripts
	##################
	final_output_panda = pd.merge(transcript_TE_panda, tss_gene_panda, how="left", on="transcript_id")
	final_output_panda = pd.merge(final_output_panda, exon_gene_panda_cleanup, how="left", on="transcript_id")

	## get annotation of different type of TE-derived transcripts
	final_output_panda["transcript_class"] = final_output_panda.apply(process_assemble.tell_the_class,axis=1)
	## get gene name
	final_output_panda['teprof3_gene_name'] = final_output_panda.apply(lambda x: x["TE_subfamily"]+"-"+x["gene_name"] if x["gene_name"]!="." else x["TE_subfamily"]+"-"+x["transcript_id"],axis=1)
	## save the information as a table
	#final_output_panda.to_csv("assembled/teprof3_output_transcript_info_"+input_gtf_file.replace("assembled/","").replace("gtf","tsv"), sep="\t", index=False)

	##################
	## Part 5. find out transcript id with the same intron chain in gene annotation compare to TE-derived transcripts
	##################
	misc.print_time("identify annotated TE-derived transcripts with same intron chain for "+input_gtf_file)
	final_output_panda["annotation_whole_intron_chain"] = final_output_panda.apply(process_assemble.check_overlap_with_whole_intron_chain_of_annotated_genes, gene_annotation_intron_dict=gene_annotation_intron_dict ,axis=1)
	final_output_panda["annotation_gene_intron_chain"] = final_output_panda.apply(process_assemble.check_overlap_with_intron_chain_of_annotated_genes, gene_annotation_intron_dict=gene_annotation_intron_dict ,axis=1)
	final_output_panda = final_output_panda.drop("transcript_intron_coordinates", axis=1)

	return(final_output_panda)


def extract_filtered_taco_info_for_output_gtf(input_row, output_file):
	chromosome = input_row["transcript_chr"]
	start = str(input_row["transcript_start"])
	stop = str(input_row["transcript_stop"])
	strand = input_row["transcript_strand"]
	transcript_id = "transcript_id \""+input_row["transcript_id"]+"\""
	gene_name = "gene_name \""+input_row["teprof3_gene_name"]+"\""
	gene_id = "gene_id \""+input_row["transcript_gene_id"]+"\""
	info = [chromosome,"TACO","transcript",start,stop,"1000",strand,".","; ".join([gene_id, gene_name,transcript_id])]
	output_file.write("\t".join(info)+";\n")
	transcript_coordinates = input_row["transcript_coordinates"]
	for exon in range(len(transcript_coordinates)):
		start = str(transcript_coordinates[exon][0])
		stop = str(transcript_coordinates[exon][1])
		exon_number = "exon_number \""+str(exon+1)+"\""
		info = [chromosome,"TACO","exon",start,stop,"1000",strand,".","; ".join([gene_id, gene_name,transcript_id,exon_number])]
		output_file.write("\t".join(info)+";\n")
	return(None)

def combine_all_filter_info(input_row, flags):
	filter_annotated = flags.filterannotated # filter_annotated=True

	keep_flag = 0
	if input_row["intron_retention_filter"]=="keep":
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


def print_intron_coordinates(input_row, output_file):
	if input_row["transcript_total_num_of_exon"] > 1:
		exon_coordinates = input_row["transcript_coordinates"]
		strand = input_row["transcript_strand"]
		intron_coordinates = []
		for i in range(len(exon_coordinates)-1):
			intron_coordinates.append((exon_coordinates[i][1], exon_coordinates[i+1][0]))
		output_file.write(str(strand)+str(intron_coordinates)+"\n")

def get_intron_coordinates(input_row):
	exon_coordinates = input_row["transcript_coordinates"]
	strand = input_row["transcript_strand"]
	intron_coordinates = []
	for i in range(len(exon_coordinates)-1):
		intron_coordinates.append((exon_coordinates[i][1], exon_coordinates[i+1][0]))
	return(str(strand)+str(intron_coordinates))

def label_duplicated_transcripts(input_dataframe):
	original_input_dataframe = input_dataframe
	## only look at transcripts that pass all filters
	input_dataframe=input_dataframe[input_dataframe["filter_result"]=="keep"]
	input_dataframe["intron_chain"] = input_dataframe.apply(get_intron_coordinates, axis=1)
	input_dataframe["TSS"] = np.where(input_dataframe["transcript_strand"] == "+", input_dataframe["transcript_start"], input_dataframe["transcript_stop"])
	## this function is removing duplicated transcripts with multiple exons. For mono-exonic transcripts, as we haven't seen such a case where there're duplicated mono-exonic transcripts, so i won't think about mono-exonic transcripts here.
	input_dataframe = input_dataframe[~((input_dataframe["intron_chain"] == "-[]") | (input_dataframe["intron_chain"] == "+[]"))]
	## mark multi-exonic transcripts with the same intron chian and TSS
	duplicate = input_dataframe[input_dataframe.duplicated(['intron_chain','TSS'], keep=False)]
	## nrow(input_dataframe) = 11508
	## nrow(duplicate) = 299
	## number of duplicated transcripts are not a lot. But we will still take care of them

	## keep the one with the longest last exon
	transcript_checked = []
	transcript_id_for_duplicated_transcripts_keep = []
	transcript_id_for_duplicated_transcripts_no_keep = []
	for intron_chain, TSS in zip(list(duplicate.intron_chain), list(duplicate.TSS)):
		if intron_chain+str(TSS) not in transcript_checked:
			transcript_checked.append(intron_chain+str(TSS))
			duplicate_temp = duplicate[(duplicate['intron_chain']==intron_chain) & (duplicate['TSS']==TSS)] # exon_gene_panda_temp = exon_gene_panda[exon_gene_panda['transcript_id']=="STRG.20.1"]			
			if intron_chain[0] == "+":
				duplicate_temp = duplicate_temp.sort_values(by=["transcript_stop"], ascending=[False])
			elif intron_chain[0] == "-":
				duplicate_temp = duplicate_temp.sort_values(by=["transcript_start"], ascending=[True])
			#print(intron_chain, TSS)
			#print(duplicate_temp)
			#print(list(duplicate_temp["transcript_id"])[0])
			#print(list(duplicate_temp["transcript_id"])[1:])
			transcript_id_for_duplicated_transcripts_keep.append(list(duplicate_temp["transcript_id"])[0])
			transcript_id_for_duplicated_transcripts_no_keep+=list(duplicate_temp["transcript_id"])[1:]
	original_input_dataframe["duplicated_transcripts"] = original_input_dataframe["transcript_id"].apply(lambda x: "duplicated_transcript_keep" if x in transcript_id_for_duplicated_transcripts_keep else ("duplicated_transcript_nokeep" if x in transcript_id_for_duplicated_transcripts_no_keep else "unique_transcript"))
	return(original_input_dataframe)

























