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
import json
import ast

import misc

def classify_transcripts(input_dataset, flags):
	''' classify transcripts into different types (TE-gene, TE-non-coding, TE transcript) '''
	processsamplenumber = flags.processsamplenumber # processsamplenumber = 10
	tpmcutoff = flags.processtpm # tpmcutoff = 0.5
	filter_mono_exon_tpm = flags.filtermonoexontpm # filter_mono_exon_tpm = 1

	## step 0. prepare gene_annotation_intron_dict
	misc.print_time("Prepare intron annotation")
	with open(os.path.dirname(os.path.realpath(sys.argv[0]))+"/../reference/gene_intron_coordinates_annotation.json", "r") as json_file:
		gene_annotation_intron_dict = json.load(json_file)

	## step 1. convert each gtf files into tsv for easier handling
	input_gtf_files = []
	input_samples = []
	for sample in input_dataset:
		input_samples.append(sample)
		input_gtf_files.append(input_dataset[sample]['gtf'])

	misc.print_time("Extract exon and transcript information from gtf files")
	with mp.Pool(processsamplenumber) as pool:
		uesless = pool.map(convert_gtf_to_pandas, input_gtf_files)

	## step 2. intersect transcripts with repeatmasker and gene annotation (gencode) to classify transcripts
	input_gtf_exon_files = glob.glob("./assembled/*exon.txt")
	input_gtf_transcript_files = glob.glob("./assembled/*transcript.txt")
	input_gtf_tss_files = glob.glob("./assembled/*tss.txt")

	misc.print_time("Intersect transcript with TE annotation => to check if it's TE-derived")
	with mp.Pool(processsamplenumber) as pool:
		uesless = pool.map(intersect_with_TE, input_gtf_transcript_files)

	misc.print_time("Intersect TSS with gene annotation => to check genic feature of TSS")
	with mp.Pool(processsamplenumber) as pool:
		uesless = pool.map(intersect_with_gene_exon_intron, input_gtf_tss_files)

	misc.print_time("Intersect exons with gene annotation => to find the assocaited gene")
	with mp.Pool(processsamplenumber) as pool:
		uesless = pool.map(intersect_with_gene_exon, input_gtf_exon_files)

	## step 3. based on the intersection result, classify transcripts into 4 categories, TE-coding gene, TE-noncoding gene, TE-no gene and TE transcripts
	misc.print_time("Classify transcripts into four categories with tpm filtering ("+str(tpmcutoff)+"TPM)")
	with mp.Pool(processsamplenumber) as pool:
		bad_samples = pool.starmap(categorize_transcripts, zip(input_samples, input_gtf_files, repeat(flags), repeat(tpmcutoff), repeat(filter_mono_exon_tpm), repeat(gene_annotation_intron_dict)))

	## step 4. remove bad samples from the manifest dictionary
	input_dataset_filtered = {}
	for sample in input_dataset:
		if sample not in bad_samples:
			input_dataset_filtered[sample] = input_dataset[sample]

	if len(input_dataset_filtered) == 0:
		misc.print_time("There's no sample that have > "+str(flags.processtranscriptnumber)+" TE-derived transcripts detected. Please check your sample quality. If quality looks good, change the --processtranscriptnumber flag to change the cutoff.")
		exit()

	return(input_dataset_filtered)

def convert_gtf_to_pandas(input_gtf_file):
	""" convert gtf file to pandas dataframe for faster and easier downstream processing """
	## check the first column and see if it's "chr1" or "1", bedtools is not compatible with "1", only with "chr1"
	with open(input_gtf_file, "r") as test_gtf_file:
		for i in range(10):  # Skip the first 10 lines
			useless = test_gtf_file.readline()
		test_line = test_gtf_file.readline()  # Read the 11th line
		if "chr" not in test_line.strip("\n").split("\t")[0]:
			misc.print_time(input_gtf_file+" was aligned to a reference fasta file with chromosome names as \"1\", instead of \"chr1\". Please re-run the alignment.")
			exit()
	with open(input_gtf_file, "r") as gtf_file, open(input_gtf_file+".transcript.txt", "w") as transcript_file, open(input_gtf_file+".exon.txt", "w") as exon_file, open(input_gtf_file+".tss.txt","w") as tss_file:
		for line in tqdm(gtf_file, desc="Processing "+input_gtf_file):
			if line[0] != "#":
				entry = line.strip("\n").split("\t")
				if "GL" not in entry[0] and "KI" not in entry[0] and "chrM" not in entry[0]:
					if entry[2] == "transcript":
						detail_info = entry[8]
						strand = entry[6]
						gene_id = detail_info.split("gene_id \"")[1].split("\";")[0]
						transcript_id = detail_info.split("transcript_id \"")[1].split("\";")[0]
						cov = detail_info.split("cov \"")[1].split("\";")[0]
						FPKM = detail_info.split("FPKM \"")[1].split("\";")[0]
						TPM = detail_info.split("TPM \"")[1].split("\";")[0]
						transcript_file.write("\t".join([entry[0], entry[3], entry[4], transcript_id, cov, strand, gene_id, FPKM, TPM])+'\n') # chr, start, stop, transcript_id, cov, strand, gene_id, FPKM, TPM
						if strand == "+":
							tss_file.write("\t".join([entry[0], entry[3], str(int(entry[3])+1), transcript_id, cov, strand, gene_id])+'\n') # chr, start, stop, transcript_id, cov, strand, gene_id, exon_number
						elif strand == "-":
							tss_file.write("\t".join([entry[0], str(int(entry[4])-1), entry[4], transcript_id, cov, strand, gene_id])+'\n') # chr, start, stop, transcript_id, cov, strand, gene_id, exon_number
					if entry[2] == "exon":
						detail_info = entry[8]
						strand = entry[6]
						gene_id = detail_info.split("gene_id \"")[1].split("\";")[0]
						transcript_id = detail_info.split("transcript_id \"")[1].split("\";")[0]
						exon_number = "exon_"+str(detail_info.split("exon_number \"")[1].split("\";")[0])
						cov = detail_info.split("cov \"")[1].split("\";")[0]
						exon_file.write("\t".join([entry[0], entry[3], entry[4], transcript_id, cov, strand, gene_id, exon_number])+'\n') # chr, start, stop, transcript_id, cov, strand, gene_id, exon_number			


def intersect_with_TE(input_file):
	repeatmasker_file = os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/repeatmasker_sorted.txt"
	subprocess.run("bedtools intersect -a "+input_file+" -b "+repeatmasker_file+" -wa -wb -loj >"+input_file.replace("txt","TE.txt"), shell=True, executable='/bin/bash')
	subprocess.run("rm "+input_file, shell=True)

def intersect_with_gene_exon_intron(input_file):
	gene_annotation_file = os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/gene_exon_intron_annotation_sorted.txt"
	subprocess.run("bedtools intersect -a "+input_file+" -b "+gene_annotation_file+" -wa -wb -loj -s >"+input_file.replace("txt","gene_exon_intron.txt"), shell=True, executable='/bin/bash')
	subprocess.run("rm "+input_file, shell=True)

def intersect_with_gene_exon(input_file):
	gene_annotation_file = os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/gene_exon_annotation_sorted.txt"
	subprocess.run("bedtools intersect -a "+input_file+" -b "+gene_annotation_file+" -wa -wb -loj -s >"+input_file.replace("txt","gene_exon.txt"), shell=True, executable='/bin/bash')
	subprocess.run("rm "+input_file, shell=True)

def categorize_transcripts(input_sample, input_gtf_file, flags, tpmcutoff, filter_mono_exon_tpm, gene_annotation_intron_dict):
	transcript_TE_file = input_gtf_file+".transcript.TE.txt" # input_gtf_file="Trimmed_mRNA_JHU006_NEG_BRep2_R1.fastqAligned.sortedByCoord.out.gtf"
	tss_gene_file = input_gtf_file+".tss.gene_exon_intron.txt"
	exon_gene_file = input_gtf_file+".exon.gene_exon.txt" # input_gtf_file="Trimmed_MPNST_PRC2_LOF_10.R1.fastqAligned.sortedByCoord.out.gtf"

	misc.print_time("Read in intersected informations for "+input_gtf_file)
	transcript_TE_panda = pd.read_table(transcript_TE_file, names=['transcript_chr','transcript_start','transcript_stop','transcript_id','transcript_cov','transcript_strand','transcript_gene_id','transcript_FPKM','transcript_TPM','TE_chr','TE_start','TE_stop','TE_subfamily','TE_sub','TE_strand','TE_family','TE_class'])
	tss_gene_panda = pd.read_table(tss_gene_file, names=['tss_chr','tss_start','tss_stop','tss_transcript_id','tss_transcript_cov','tss_strand','tss_gene_id','gene_chr','gene_start','gene_stop','gene_gene_name','useless1','gene_strand','gene_transcript_id','gene_exon_intron_number','gene_gene_id','gene_gene_type','gene_transcript_type','gene_transcript_name','.'])
	exon_gene_panda = pd.read_table(exon_gene_file, names=['exon_chr','exon_start','exon_stop','transcript_id','transcript_cov','transcript_strand','transcript_id_2','exon_number','gene_chr','gene_start','gene_stop','gene_name','useless1','gene_strand','gene_transcript_id','gene_exon_number','gene_id','gene_type','gene_transcript_type','gene_transcript_name','useless2'])

	keep_intermediate_files = flags.keepintermediate
	if not keep_intermediate_files:
		subprocess.run("rm "+transcript_TE_file+" "+tss_gene_file, shell=True)

	##################
	## Part 1. process transcript_TE_panda to extract TE-derived transcripts
	##################
	misc.print_time("Process transcript_TE_panda to extract TE-derived transcripts for "+input_gtf_file)
	transcript_TE_panda = transcript_TE_panda[transcript_TE_panda["transcript_strand"]!="."] # remove transcripts without strand information

	total_id_of_transcripts = set(transcript_TE_panda[transcript_TE_panda['transcript_TPM']>=tpmcutoff]["transcript_id"])
	total_number_of_transcripts = len(total_id_of_transcripts)

	with open(input_gtf_file+".stats.txt", "a") as stats_file:
		stats_file.write("\t".join([input_gtf_file, "# of transcripts", str(len(set(transcript_TE_panda["transcript_id"])))])+"\n")
		stats_file.write("\t".join([input_gtf_file, "# of transcripts that pass tpm cutoff (>="+str(tpmcutoff)+"TPM)", str(total_number_of_transcripts)])+"\n")

	transcript_TE_panda['start_from_TE'] = transcript_TE_panda.apply(start_from_TE, axis=1)
	transcript_TE_panda['all_in_TE'] = transcript_TE_panda.apply(lambda x: "yes" if (x['TE_start'] <= x['transcript_start'] and x['transcript_stop'] <= x['TE_stop']) else "no", axis=1)
	transcript_TE_panda = transcript_TE_panda[transcript_TE_panda['start_from_TE']=="yes"]
	transcript_TE_panda = transcript_TE_panda[transcript_TE_panda['transcript_id'].duplicated() == False]
	transcript_TE_panda = transcript_TE_panda[transcript_TE_panda['transcript_TPM']>=tpmcutoff]
	transcript_TE_panda = transcript_TE_panda.reset_index(drop=True)
	
	## remove transcripts less than tpm cutoff
	with open(input_gtf_file+".stats.txt", "a") as stats_file:
		stats_file.write("\t".join([input_gtf_file, "# of non-TE-derived transcripts that pass tpm cutoff (>="+str(tpmcutoff)+"TPM)", str(total_number_of_transcripts-len(transcript_TE_panda.index))])+"\n")
		stats_file.write("\t".join([input_gtf_file, "# of TE-derived transcripts that pass tpm cutoff (>="+str(tpmcutoff)+"TPM)", str(len(transcript_TE_panda.index))])+"\n")
	
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
	tss_gene_panda["distance_between_stringtie_tss_and_annotated_tss"] = tss_gene_panda.apply(get_distance_between_stringtie_tss_and_annotated_tss_overlap_tssandtranscript,axis=1)
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
	exon_gene_panda = exon_gene_panda.drop(['transcript_id_2','useless1','useless2'], axis=1)
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
		exon_1_panda['annotated_exon1_in_the_same_TE'] = exon_1_panda.apply(find_good_annotated_exon1, transcript_TE_panda=transcript_TE_panda[transcript_TE_panda["transcript_id"]==transcript_id], axis=1)
		exon_1_panda['distance_between_stringtie_exonedge_and_annotated_exonedge'] = exon_1_panda.apply(get_distance_between_stringtie_exonedge_and_annotated_exonedge_overlap_exonandgene, axis=1)
		exon_1_panda['overlap_between_stringtie_exon_and_annotated_exon'] = exon_1_panda.apply(get_overlap_between_stringtie_exon_and_annotated_exon, axis=1)
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
	final_output_panda = pd.merge(transcript_TE_panda, tss_gene_panda, how="inner", on="transcript_id")
	final_output_panda = pd.merge(final_output_panda, exon_gene_panda_cleanup, how="inner", on="transcript_id")

	## get annotation of different type of TE-derived transcripts
	final_output_panda["transcript_class"] = final_output_panda.apply(tell_the_class,axis=1)

	with open(input_gtf_file+".stats.txt", "a") as stats_file:
		stats_file.write("\t".join([input_gtf_file, "# of mono-exonic TE-derived transcripts that pass tpm cutoff (>="+str(tpmcutoff)+"TPM)", str(len(final_output_panda[(final_output_panda['transcript_total_num_of_exon']==1)&(final_output_panda['transcript_TPM']>=tpmcutoff)].index))])+"\n")
		stats_file.write("\t".join([input_gtf_file, "# of mono-exonic TE-derived transcripts that pass tpm cutoff (>="+str(filter_mono_exon_tpm)+"TPM)", str(len(final_output_panda[(final_output_panda['transcript_total_num_of_exon']==1)&(final_output_panda['transcript_TPM']>=filter_mono_exon_tpm)].index))])+"\n")

	##################
	## Part 5. find out transcript id with the same intron chain in gene annotation compare to TE-derived transcripts
	##################
	misc.print_time("identify annotated TE-derived transcripts with same intron chain for "+input_gtf_file)
	final_output_panda["annotation_whole_intron_chain"] = final_output_panda.apply(check_overlap_with_whole_intron_chain_of_annotated_genes, gene_annotation_intron_dict=gene_annotation_intron_dict ,axis=1)
	final_output_panda["annotation_gene_intron_chain"] = final_output_panda.apply(check_overlap_with_intron_chain_of_annotated_genes, gene_annotation_intron_dict=gene_annotation_intron_dict ,axis=1)
	final_output_panda = final_output_panda.drop("transcript_intron_coordinates", axis=1)

	##################
	## Part 8. save output files
	################## 
	## save the information as a table
	final_output_panda.to_csv("assembled/teprof3_output_transcript_info_"+input_gtf_file.replace("assembled/","").replace("gtf","tsv"), sep="\t", index=False)

	if len(final_output_panda) <= flags.processtranscriptnumber:
		with open("teprof3_warning_please_read.txt", "a") as warning_file:
			warning_file.write(input_sample+" has "+str(len(final_output_panda))+" (<="+str(flags.processtranscriptnumber)+") TE-derived transcripts identified, thus it's removed from further analysis. This usually indicates there's severe 3' bias in the sample due to RNA degradation or library preparation method. Change the --processtranscriptnumber flag to adjust the cutoff."+"\n")
		return(input_sample)
	else:
		return("good_sample")

def start_from_TE(input_row):
	if input_row['transcript_strand'] == "+":
		if input_row['TE_start'] <= input_row['transcript_start'] and input_row['transcript_start'] <= input_row['TE_stop']:
			return("yes")
		else:
			return("no")
	if input_row['transcript_strand'] == "-":
		if input_row['TE_start'] <= input_row['transcript_stop'] and input_row['transcript_stop'] <= input_row['TE_stop']:
			return("yes")
		else:
			return("no")

def extract_info_for_output_gtf(input_row, output_file):
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
	transcript_coordinates = input_row["transcript_coordinates"]
	if strand == "-":
		transcript_coordinates = transcript_coordinates[::-1]
	for exon in range(len(transcript_coordinates)):
		start = str(transcript_coordinates[exon][0])
		stop = str(transcript_coordinates[exon][1])
		exon_number = "exon_number \""+str(exon+1)+"\""
		info = [chromosome,"StringTie","exon",start,stop,"1000",strand,".","; ".join([gene_id,transcript_id,exon_number])]
		output_file.write("\t".join(info)+";\n")
	return(None)

def tell_the_class(input_row):
	## first class, the whole transcript is within TE
	if input_row["all_in_TE"] == "yes":
		return("TE transcript")
	## second class, TE-gene transcript
	elif input_row["gene_transcript_type"] == "coding":
		return("TE coding gene transcript")
	## third class, TE-non-coding transcript
	elif input_row["gene_transcript_type"] != ".":
		return("TE noncoding gene transcript")
	## fourth class, TE-nothing transcript
	elif input_row["gene_transcript_type"] == ".":
		return("TE no gene transcript")

def get_distance_between_stringtie_tss_and_annotated_tss_overlap_tssandtranscript(input_row):
	if str(input_row["gene_start"]) == "-1":
		return(-999)
	elif str(input_row["tss_strand"]) == "+":
		return(int(input_row["gene_start"]) - int(input_row["tss_start"]))
	elif str(input_row["tss_strand"]) == "-":
		return(int(input_row["tss_start"]) - int(input_row["gene_stop"]))

def get_distance_between_stringtie_exonedge_and_annotated_exonedge_overlap_exonandgene(input_row):
	if str(input_row["gene_start"]) == "-1":
		return(-1)
	elif str(input_row["transcript_strand"]) == "+":
		return(abs(int(input_row["exon_start"]) - int(input_row["gene_start"])))
	elif str(input_row["transcript_strand"]) == "-":
		return(abs(int(input_row["exon_stop"]) - int(input_row["gene_stop"])))

def get_overlap_between_stringtie_exon_and_annotated_exon(input_row):
	if str(input_row["gene_start"]) == "-1":
		return(-1)
	else:
		temp_list = [int(input_row["exon_start"]), int(input_row["exon_stop"]), int(input_row["gene_start"]), int(input_row["gene_stop"])]
		sorted_temp_list = sorted(temp_list)
		return(sorted_temp_list[2]-sorted_temp_list[1])

def check_overlap_with_whole_intron_chain_of_annotated_genes(input_row, gene_annotation_intron_dict):
	if input_row["transcript_total_num_of_exon"] == 1:
		return("mono_exonic")
	elif input_row["gene_transcript_id"] == ".":
		return("no") 
	elif str(input_row["transcript_intron_coordinates"]) == gene_annotation_intron_dict[input_row["gene_transcript_id"]][1]:
		return("yes")
	else:
		return("no")

def check_overlap_with_intron_chain_of_annotated_genes(input_row, gene_annotation_intron_dict): ## on top of compare the whole intron chain, also compare the canonical intron chain of TE-derived transcripts and see if they overlapped 100% as a QC
	if input_row["transcript_total_num_of_exon"] == 1:
		return("mono_exonic")
	elif input_row["gene_transcript_id"] == ".":
		return("no")
	elif str(input_row["transcript_intron_coordinates"][(int(input_row["transcript_exon_number_overlap_gene"].split("_")[1])-1):]) == gene_annotation_intron_dict[input_row["gene_transcript_id"]][1][(int(input_row["gene_exon_number"].split("_")[1])-1):]:
		return("yes")
	else:
		return("no")

def find_good_annotated_exon1(input_row, transcript_TE_panda):
	transcript_TE_panda = transcript_TE_panda.reset_index(drop=True)
	if input_row["gene_exon_number"] == "exon_1":
		if input_row["transcript_strand"] == "+":
			if transcript_TE_panda.loc[0,"TE_start"] <= input_row["gene_start"] <= transcript_TE_panda.loc[0,"TE_stop"]:
				return(1)
			else:
				return(0)
		elif input_row["transcript_strand"] == "-":
			if transcript_TE_panda.loc[0,"TE_start"] <= input_row["gene_stop"] <= transcript_TE_panda.loc[0,"TE_stop"]:
				return(1)
			else:
				return(0)
	else:
		return(0)








