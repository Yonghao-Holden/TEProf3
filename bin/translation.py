# TEProf3: TE-derived Promoter Finder 3
# Author: Holden Liang
# de novo assembly script
# -----------------------------------------------------------------------------
# Functions related to processing the bam files and perform de novo assembly
import argparse
import time
import subprocess # https://geekflare.com/python-run-bash/
import os
import multiprocess as mp
from multiprocessing import Pool
from itertools import repeat
import sys
import pandas as pd
from tqdm import tqdm
import glob
import numpy as np
import genomepy # genomepy.install_genome("hg38", annotation=True, provider="UCSC", genomes_dir="/bar/yliang/genomes/private/genomepy")
import math
import ast
import json

import misc
import process_meta_assemble

def start_translation(input_file, flags, version):
	""" in silico translation """
	misc.print_time("Start in silico translation")

	with open("command_used_for_teprof3_translation.txt", "w") as command_file:
		command_file.write(' '.join(sys.argv)+"\n")
		command_file.write("You are using TEProf3 version: "+version)

	## step 1. obtain essential information from previous step or teprof2
	if flags.teprof2 != "no":
		misc.print_time("You are using bonus function: directly taking teprof2 output for in silico translation")
		# use it like teprof3 --teprof2 /scratch/yliang/HNSCC/data/CPTAC/2_test2_actual_run_TEProf2/test2_database/annotatedcufftranscripts.tsv
		# input_panda = convert_teprof2_output_for_teprof3_translation("/scratch/yliang/HNSCC/data/CPTAC/2_test2_actual_run_TEProf2/test2_database/annotatedcufftranscripts.tsv")
		input_panda = convert_teprof2_output_for_teprof3_translation(flags.teprof2)
	else:
		misc.print_time("Clean up to prepare for in silico translation")
		misc.run_bash_command("mkdir teprof3_output_translation")
		misc.run_bash_command("mv "+input_file+" teprof3_output_translation/"+input_file)
		input_file = "teprof3_output_translation/"+input_file

		## check if need to re-annotate identified TE-derived transcripts
		if flags.translationmode == "1":
			input_file = re_annotate_transcripts(input_gtf_file=input_file, flags=flags)
		elif flags.translationmode == "2":
			pass
		else:
			misc.print_time("Please provide either 1 or 2 for flag --translationmode")
			exit()

		input_panda = pd.read_table(input_file)
		if flags.translationmode == "2":
			input_panda = input_panda[(input_panda["filter_result"]=="keep") & (input_panda["duplicated_transcripts"]!="duplicated_transcript_nokeep")]
		input_panda = input_panda.drop_duplicates(subset='transcript_id')
		input_panda = input_panda.apply(get_exon_coordinates, axis=1)

	## step 2. get RNA sequence and run in silico translation
	translation_length = flags.translationlength # translation_length=20

	genome_name = flags.translationgenome # genome_name = "hg38"
	genomepy_folder = os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/genomepy" # genomepy_folder = "/bar/yliang/genomes/private/genomepy/"
	genome = genomepy.Genome(genome_name, genomes_dir=genomepy_folder)

	input_panda_protein = translation(input_panda, translation_length, genome)

	## step 3. split the protein sequence to TE-derived and gene-derived
	input_panda_protein = input_panda_protein.apply(split_TE_and_gene_protein_sequence, axis=1) # input_panda_protein.to_csv("/scratch/yliang/HNSCC/data/CPTAC/2_test2_actual_run_TEProf2/test2_database/input_panda_protein.tsv", sep="\t", index=False)

	input_panda_protein.to_csv("teprof3_output_translation/teprof3_output_protein_information.tsv", sep="\t", index=False)

	## step 4. export sequences to fasta files
	get_fasta(input_panda_protein)

	## step 5. the complete tsv file is too big, drop the rna and protein sequence from the table to have a light weight table
	#input_panda_protein = input_panda_protein.drop(columns=["transcript_TE_exons_rna","transcript_gene_exons_rna"]) 
	input_panda_protein = input_panda_protein[["transcript_chr","transcript_start","transcript_stop","transcript_id","transcript_strand","transcript_coordinates","transcript_exon_number_overlap_gene","gene_name","gene_transcript_id","gene_exon_number","transcript_exons_rna","start_codons","start_codon","protein_sequence","transcript_TE_unannotated_exons_aa","transcript_nonTE_unannotated_exons_aa","transcript_gene_annotated_exons_aa"]]
	input_panda_protein.to_csv("teprof3_output_translation/teprof3_output_protein_information_light.tsv", sep="\t", index=False)

	misc.print_time("In silico translation using all ORFs is done. Please enjoy!")

	## step 6. blastp on TE-derived, gene-derived and whole protein sequence
	#run_blastp(input_panda_protein)

def re_annotate_transcripts(input_gtf_file, flags):
	## step 0. prepare gene_annotation_intron_dict
	misc.print_time("Prepare intron annotation")
	with open(os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/gene_intron_coordinates_annotation.json", "r") as json_file:
		gene_annotation_intron_dict = json.load(json_file)

	## step 1. convert the gtf file into tsv for easier handling
	misc.print_time("Extract exon and transcript information from the provided gtf file")
	uesless = process_meta_assemble.convert_taco_gtf_to_pandas(input_gtf_file)

	## step 2. intersect transcripts with repeatmasker and gene annotation (gencode) to classify transcripts
	input_gtf_exon_file = input_gtf_file+".exon.txt"
	input_gtf_transcript_file = input_gtf_file+".transcript.txt"
	input_gtf_tss_file = input_gtf_file+".tss.txt"

	misc.print_time("Intersect transcript with TE annotation => to check if it's TE-derived")
	uesless = intersect_with_TE(input_gtf_transcript_file)

	misc.print_time("Intersect TSS with gene annotation => to check genic feature of TSS")
	uesless = intersect_with_gene_exon_intron(input_gtf_tss_file)

	misc.print_time("Intersect exons with gene annotation => to find the assocaited gene")
	uesless = intersect_with_gene_exon(input_gtf_exon_file)

	## step 3. based on the intersection result, classify transcripts into 4 categories, TE-coding gene, TE-noncoding gene, TE-no gene and TE transcripts
	misc.print_time("Classify transcripts into four categories")
	final_output_panda = process_meta_assemble.categorize_taco_transcripts(input_gtf_file, flags, gene_annotation_intron_dict)

	final_output_panda.to_csv("teprof3_output_translation/teprof3_output_filter_transcript_TE_transcript_consensus.reannotated.tsv", sep="\t", index=False)

	return("teprof3_output_translation/teprof3_output_filter_transcript_TE_transcript_consensus.reannotated.tsv")


def intersect_with_TE(input_file):
	repeatmasker_file = os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/repeatmasker_sorted.txt"
	subprocess.run("bedtools intersect -a "+input_file+" -b "+repeatmasker_file+" -wa -wb -loj >"+input_file.replace("txt","TE.txt"), shell=True, executable='/bin/bash')
	subprocess.run("rm "+input_file, shell=True)

def intersect_with_gene_exon_intron(input_file):
	gene_annotation_file = os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/gene_exon_intron_annotation_sorted.for_translation.txt"
	subprocess.run("bedtools intersect -a "+input_file+" -b "+gene_annotation_file+" -wa -wb -loj -s >"+input_file.replace("txt","gene_exon_intron.txt"), shell=True, executable='/bin/bash')
	subprocess.run("rm "+input_file, shell=True)

def intersect_with_gene_exon(input_file):
	gene_annotation_file = os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/gene_exon_annotation_sorted.for_translation.txt"
	subprocess.run("bedtools intersect -a "+input_file+" -b "+gene_annotation_file+" -wa -wb -loj -s >"+input_file.replace("txt","gene_exon.txt"), shell=True, executable='/bin/bash')
	subprocess.run("rm "+input_file, shell=True)

def get_exon_coordinates(input_row):
	transcript_coordinates = ast.literal_eval(input_row["transcript_coordinates"])
	transcript_strand = input_row["transcript_strand"]

	## Part 1. get coordinates for annotated exons
	if input_row["gene_exon_number"] == ".":
		input_row["transcript_gene_annotated_exons_coord"] = []
		unannotated_transcript_coordinates = transcript_coordinates

	elif input_row["gene_exon_number"] != ".":
		gene_transcript_coordinates = ast.literal_eval(input_row["gene_transcript_coordinates"])
		transcript_exon_number_overlap_gene = int(input_row["transcript_exon_number_overlap_gene"].split("exon_")[1])
		splice_gene_target = int(input_row["gene_exon_number"].split("exon_")[1])

		## looking at the overlapped exon
		if transcript_strand == "+":
			unannotated_transcript_coordinates = transcript_coordinates[0:transcript_exon_number_overlap_gene-1]

			TE_exon_start = transcript_coordinates[transcript_exon_number_overlap_gene-1][0]
			TE_exon_stop = transcript_coordinates[transcript_exon_number_overlap_gene-1][1]
			gene_exon_start = gene_transcript_coordinates[splice_gene_target-1][0]
			gene_exon_stop = gene_transcript_coordinates[splice_gene_target-1][1]
		
			if TE_exon_start < gene_exon_start:
				unannotated_transcript_coordinates.append((TE_exon_start,gene_exon_start-1))
				if len(transcript_coordinates) > transcript_exon_number_overlap_gene:
					input_row["transcript_gene_annotated_exons_coord"] = [(gene_exon_start,TE_exon_stop)] + transcript_coordinates[transcript_exon_number_overlap_gene:]
				else:
					input_row["transcript_gene_annotated_exons_coord"] = [(gene_exon_start,TE_exon_stop)]
			elif gene_exon_start <= TE_exon_start:
				input_row["transcript_gene_annotated_exons_coord"] = transcript_coordinates[transcript_exon_number_overlap_gene-1:]

		elif input_row["transcript_strand"] == "-":
			if transcript_exon_number_overlap_gene != 1:
				unannotated_transcript_coordinates = transcript_coordinates[(-1)*(transcript_exon_number_overlap_gene-1):]
			elif transcript_exon_number_overlap_gene == 1:
				unannotated_transcript_coordinates = []

			TE_exon_start = transcript_coordinates[(-1)*transcript_exon_number_overlap_gene][0]
			TE_exon_stop = transcript_coordinates[(-1)*transcript_exon_number_overlap_gene][1]
			gene_exon_start = gene_transcript_coordinates[(-1)*splice_gene_target][0]
			gene_exon_stop = gene_transcript_coordinates[(-1)*splice_gene_target][1]

			if gene_exon_stop < TE_exon_stop:
				unannotated_transcript_coordinates.append((gene_exon_stop+1,TE_exon_stop))
				if len(transcript_coordinates) > transcript_exon_number_overlap_gene:
					input_row["transcript_gene_annotated_exons_coord"] = transcript_coordinates[:(-1)*transcript_exon_number_overlap_gene] + [(TE_exon_start,gene_exon_stop)]
				else:
					input_row["transcript_gene_annotated_exons_coord"] = [(TE_exon_start,gene_exon_stop)]
			elif TE_exon_stop <= gene_exon_stop:
				input_row["transcript_gene_annotated_exons_coord"] = transcript_coordinates[:(-1)*(transcript_exon_number_overlap_gene-1)]

		input_row["transcript_gene_annotated_exons_coord"] = sorted(input_row["transcript_gene_annotated_exons_coord"], key=lambda x: x[0])
		unannotated_transcript_coordinates = sorted(unannotated_transcript_coordinates, key=lambda x: x[0])

	## Part 2. get coordinates for unannotated exons
	TE_start = input_row["TE_start"]
	TE_stop = input_row["TE_stop"]
	input_row["transcript_TE_unannotated_exons_coord"] = []
	input_row["transcript_nonTE_unannotated_exons_coord"] = []
	for exon in unannotated_transcript_coordinates:
		if transcript_strand == "+":
			if exon[1] <= TE_stop:
				input_row["transcript_TE_unannotated_exons_coord"].append(exon)
			elif exon[0] < TE_stop and TE_stop < exon[1]:
				input_row["transcript_TE_unannotated_exons_coord"].append((exon[0],TE_stop))
				input_row["transcript_nonTE_unannotated_exons_coord"].append((TE_stop+1,exon[1]))
			elif TE_stop <= exon[0]:
				input_row["transcript_nonTE_unannotated_exons_coord"].append(exon)
		elif transcript_strand == "-":
			if exon[1] <= TE_start:
				input_row["transcript_nonTE_unannotated_exons_coord"].append(exon)
			elif exon[0] < TE_start and TE_start < exon[1]:
				input_row["transcript_nonTE_unannotated_exons_coord"].append((exon[0], TE_start-1))
				input_row["transcript_TE_unannotated_exons_coord"].append((TE_start, exon[1]))
			elif TE_start <= exon[0]:
				input_row["transcript_TE_unannotated_exons_coord"].append(exon)

	input_row["transcript_TE_unannotated_exons_coord"] = sorted(input_row["transcript_TE_unannotated_exons_coord"], key=lambda x: x[0])
	input_row["transcript_nonTE_unannotated_exons_coord"] = sorted(input_row["transcript_nonTE_unannotated_exons_coord"], key=lambda x: x[0])

	return(input_row)

def translation(input_panda, translation_length, genome):	
	input_panda = input_panda.apply(get_RNA_sequence, genome = genome, axis=1)
	# only consider start_codons within TE-derived (TE) sequence, so the resulted protein product will have TE-derived amino acid sequence to be detected in MS
	input_panda["start_codons"] = input_panda.apply(find_start_codon, axis=1)
	output_proteins_dict = {} 
	row_index = 0
	for transcript_id in tqdm(list(input_panda.transcript_id), desc="Translating each transcript"):	
		input_panda_temp = input_panda[input_panda['transcript_id']==transcript_id].reset_index(drop=True)
		rna_sequence = input_panda_temp.iloc[0]["transcript_exons_rna"]
		start_codons = input_panda_temp.iloc[0]["start_codons"]
		if len(start_codons)==0:
			output_proteins_dict[row_index] = [transcript_id,".","."]
			row_index+=1
		else:
			for start_codon in start_codons:
				rna_sequence_to_be_translated = rna_sequence[start_codon:]
				output_proteins_dict[row_index] = [transcript_id, start_codon, RNA_to_protein(rna_sequence_to_be_translated, translation_length)]
				row_index+=1
	output_panda = pd.DataFrame.from_dict(output_proteins_dict, orient='index',columns=['transcript_id', 'start_codon', 'protein_sequence']) # protein_type (in-frame truncated, etc.)
	output_panda = output_panda[output_panda["protein_sequence"]!="."]
	output_panda = input_panda.merge(output_panda, on="transcript_id")
	return(output_panda)


def get_RNA_sequence(input_row, genome):
	# input_row = input_panda[input_panda["transcript_id"]=="TCONS_00016150"]
	transcript_chr = input_row["transcript_chr"]
	transcript_strand = input_row["transcript_strand"]
	TE_exons = input_row["transcript_TE_unannotated_exons_coord"]
	nonTE_exons = input_row["transcript_nonTE_unannotated_exons_coord"]
	gene_exons = input_row["transcript_gene_annotated_exons_coord"]
	transcript_coordinates = ast.literal_eval(input_row["transcript_coordinates"])

	if len(TE_exons) != 0:
		if transcript_strand == "+":
			input_row["transcript_TE_unannotated_exons_rna"] = str(genome.get_spliced_seq(transcript_chr, TE_exons))
		elif transcript_strand == "-":
			input_row["transcript_TE_unannotated_exons_rna"] = str(genome.get_spliced_seq(transcript_chr, TE_exons, rc=True))
	else:
		input_row["transcript_TE_unannotated_exons_rna"] = ""

	if len(nonTE_exons) != 0:
		if transcript_strand == "+":
			input_row["transcript_nonTE_unannotated_exons_rna"] = str(genome.get_spliced_seq(transcript_chr, nonTE_exons))
		elif transcript_strand == "-":
			input_row["transcript_nonTE_unannotated_exons_rna"] = str(genome.get_spliced_seq(transcript_chr, nonTE_exons, rc=True))
	else:
		input_row["transcript_nonTE_unannotated_exons_rna"] = ""

	if len(gene_exons) != 0:
		if transcript_strand == "+":
			input_row["transcript_gene_annotated_exons_rna"] = str(genome.get_spliced_seq(transcript_chr, gene_exons))
		elif transcript_strand == "-":
			input_row["transcript_gene_annotated_exons_rna"] = str(genome.get_spliced_seq(transcript_chr, gene_exons, rc=True))	
	else:
		input_row["transcript_gene_annotated_exons_rna"] = ""

	if transcript_strand == "+":
		input_row["transcript_exons_rna"] = str(genome.get_spliced_seq(transcript_chr, transcript_coordinates))
	elif transcript_strand == "-":
		input_row["transcript_exons_rna"] = str(genome.get_spliced_seq(transcript_chr, transcript_coordinates, rc=True))

	return(input_row)

def find_start_codon(input_row):
	input_sequence = input_row["transcript_TE_unannotated_exons_rna"] + input_row["transcript_nonTE_unannotated_exons_rna"]
	input_sequence = input_sequence.upper()
	start_codons = []
	for i in range(len(input_sequence)-2):
		if input_sequence[i:i+3] == "ATG":
			start_codons.append(i)
	return(start_codons)

def RNA_to_protein(seq, length):
	seq = seq.upper()
	if "N" in seq:
		return(".")
	codon_table = {
		'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
		'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
		'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
		'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
		'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
		'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
		'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
		'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
		'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
		'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
		'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
		'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
		'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
		'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
		'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
		'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
	}
	protein =""
	for i in range(0, len(seq), 3): # len(seq)=10. 0, 3, 6, 9
		if i+3 <= len(seq):
			codon = seq[i:i + 3]
			if codon_table[codon] == "_":
				break
			else:
				protein+=codon_table[codon]
	if len(protein)>=length:
		return(protein)
	else:
		return(".")
	

def split_TE_and_gene_protein_sequence(input_row):
	start_codon = input_row["start_codon"]
	protein_sequence = input_row["protein_sequence"]

	transcript_TE_unannotated_exons_rna = input_row["transcript_TE_unannotated_exons_rna"]
	transcript_nonTE_unannotated_exons_rna = input_row["transcript_nonTE_unannotated_exons_rna"]
	transcript_gene_annotated_exons_rna = input_row["transcript_gene_annotated_exons_rna"]

	if start_codon <= len(transcript_TE_unannotated_exons_rna):
		input_row["transcript_TE_unannotated_exons_aa"]    = protein_sequence[0:math.ceil((len(transcript_TE_unannotated_exons_rna) - int(start_codon))/3)]
		input_row["transcript_nonTE_unannotated_exons_aa"] = protein_sequence[math.ceil((len(transcript_TE_unannotated_exons_rna) - int(start_codon))/3):math.ceil((len(transcript_TE_unannotated_exons_rna)+len(transcript_nonTE_unannotated_exons_rna)-int(start_codon))/3)]
		input_row["transcript_gene_annotated_exons_aa"] = protein_sequence[math.ceil((len(transcript_TE_unannotated_exons_rna)+len(transcript_nonTE_unannotated_exons_rna)-int(start_codon))/3):]
	elif len(transcript_TE_unannotated_exons_rna) < start_codon:
		input_row["transcript_TE_unannotated_exons_aa"]    = ""
		input_row["transcript_nonTE_unannotated_exons_aa"] = protein_sequence[0:math.ceil((len(transcript_TE_unannotated_exons_rna)+len(transcript_nonTE_unannotated_exons_rna)-int(start_codon))/3)]
		input_row["transcript_gene_annotated_exons_aa"] = protein_sequence[math.ceil((len(transcript_TE_unannotated_exons_rna)+len(transcript_nonTE_unannotated_exons_rna)-int(start_codon))/3):]

	return(input_row)


def get_fasta(input_panda):
	## creating 3 fasta files with all the protein sequences
	## reason why i pick minimum size for blastp-short as 6 amino acids: https://www.ebi.ac.uk/ipd/mhc/blast/#:~:text=The%20minimum%20length%20is%2011,cause%20the%20search%20to%20fail.
	## reason why i pick 30: https://www.ncbi.nlm.nih.gov/books/NBK279684/
	with open("teprof3_output_translation/teprof3_output_protein_sequence_gene_blastp.fa", "w") as gene_sequence_file_blastp, open("teprof3_output_translation/teprof3_output_protein_sequence_gene_blastpshort.fa", "w") as gene_sequence_file_blastpshort,\
	  open("teprof3_output_translation/teprof3_output_protein_sequence.fa", "w") as protein_sequence_file:
		for index in tqdm(input_panda.index, desc="Exporting protein sequences to fasta"):
			## for blast, gene-derived protein sequence
			if len(input_panda.iloc[index]["transcript_gene_annotated_exons_aa"])>=30:
				gene_sequence_file_blastp.write("> "+input_panda.iloc[index]["transcript_id"]+"_gene_"+str(input_panda.iloc[index]["start_codon"])+"\n"+input_panda.iloc[index]["transcript_gene_annotated_exons_aa"]+"\n")
			elif 6<=len(input_panda.iloc[index]["transcript_gene_annotated_exons_aa"])<30:
				gene_sequence_file_blastpshort.write("> "+input_panda.iloc[index]["transcript_id"]+"_gene_"+str(input_panda.iloc[index]["start_codon"])+"\n"+input_panda.iloc[index]["transcript_gene_annotated_exons_aa"]+"\n")
			## for MS search
			protein_sequence_file.write(">tr|"+input_panda.iloc[index]["transcript_id"]+"_whole_"+str(input_panda.iloc[index]["start_codon"])+"_"+input_panda.iloc[index]["TE_subfamily"]+"-"+input_panda.iloc[index]["gene_name"]+"|"+input_panda.iloc[index]["transcript_id"]+"_whole_"+str(input_panda.iloc[index]["start_codon"])+"_"+input_panda.iloc[index]["TE_subfamily"]+"-"+input_panda.iloc[index]["gene_name"]+" OS=Homo sapiens OX=9606 GN="+input_panda.iloc[index]["gene_name"]+" PE=5 SV=1"+"\n"+input_panda.iloc[index]["protein_sequence"]+"\n")			


def convert_teprof2_output_for_teprof3_translation(input_file):
	## open teprof2 output rdata and run the following commands:
	######################################
	## for HNSCC
	######################################
	## load("/scratch/yliang/HNSCC/analysis/HNSCC_TSTEA_candidates_v1/assembled/Step11_FINAL.RData")
	## write.table(annotatedcufftranscripts, file="/scratch/yliang/HNSCC/analysis/HNSCC_TSTEA_candidates_v1/assembled/annotatedcufftranscripts.tsv", sep="\t",col.names=TRUE, row.names=FALSE, quote = FALSE)
	## load the annotatedcufftranscripts.tsv file into this function
	## input_file = "/scratch/yliang/HNSCC/analysis/HNSCC_TSTEA_candidates_v1/assembled/annotatedcufftranscripts.tsv"
	######################################
	## for TCGA
	######################################
	## load("/scratch/nakul/TCGA33runALL/TCGAAnalysis.RData")
	## table2 = openxlsx::read.xlsx("/scratch/yliang/HNSCC/data/CPTAC/2_test2_actual_run_TEProf2/test2_database/Shah_Supplementary_Table2.xlsx", colNames=TRUE, startRow=2)
	## output_table = dplyr::filter(annotatedcufftranscripts, transcriptname %in% table2$Transcript.ID)
	## write.table(output_table, file="/scratch/yliang/HNSCC/data/CPTAC/2_test2_actual_run_TEProf2/test2_database/annotatedcufftranscripts.tsv", sep="\t",col.names=TRUE, row.names=FALSE, quote = FALSE)
	## input_file = "/scratch/yliang/HNSCC/data/CPTAC/2_test2_actual_run_TEProf2/test2_database/annotatedcufftranscripts.tsv"
	###################################### 
	input_panda = pd.read_table(input_file, header=0)
	input_panda = input_panda[["transcriptname","transcoord","subfamTE","gene2","elements2","id2","exonintron2","number2"]]
	input_panda.columns = ["transcript_id","transcript_coord","TE_subfamily","gene_name","gene_coord","gene_gene_id","gene_splice_target","gene_splice_target_number"]
	input_panda["strand"] = input_panda["transcript_coord"].apply(lambda x: x.split(",")[0])
	input_panda["transcript_chr"] = input_panda["transcript_coord"].apply(lambda x: x.split(",")[1])
	input_panda["gene_splice_target"] = input_panda["gene_splice_target"] + "_" + input_panda["gene_splice_target_number"].astype(str)
	input_panda = input_panda.apply(teprof2_cleanup_coord, axis=1)
	return(input_panda)

def teprof2_cleanup_coord(input_row):
	# input_row = input_panda.iloc[0,]
	# input_row = input_panda.iloc[1727,] "TCONS_00040931"
	# input_row = input_panda.iloc[1766,] "TCONS_00071273"
	strand = input_row['strand']
	gene_coordinates = [int(x) for x in input_row["gene_coord"].split(",")]
	transcript_coord = [int(x) for x in input_row['transcript_coord'].split(",")[2:]]
	output_gene_exons_coordinates = []
	output_transcript_exons_coordinates = []
	if strand == "+":
		gene_coordinates.sort()
		for i in range(0,len(gene_coordinates),4):
			output_gene_exons_coordinates.append((gene_coordinates[i], gene_coordinates[i+1])) # (left, right) of exons, from exon 1 to N
		transcript_coord.sort()
		for i in range(0,len(transcript_coord),2):
			output_transcript_exons_coordinates.append((transcript_coord[i], transcript_coord[i+1]))
	elif strand == "-":
		gene_coordinates.sort(reverse=True)
		for i in range(0,len(gene_coordinates),4):
			output_gene_exons_coordinates.append((gene_coordinates[i+1],gene_coordinates[i])) # (left, right) of exons, from exon 1 to N
		transcript_coord.sort(reverse=True)
		for i in range(0,len(transcript_coord),2):
			output_transcript_exons_coordinates.append((transcript_coord[i+1], transcript_coord[i]))
		# test if Nakul correct the exon coordinates, it turns out he didn't
		#>>> output_transcript_exons_coordinates
		#[(5619331, 5619413), (5584293, 5584500), (5583687, 5584103), (5582670, 5582846), (5581859, 5582062), (5558339, 5560043), (5553386, 5553556), (5541857, 5542027), (5539415, 5539585), (5536851, 5536940), (5533897, 5533988), (5533304, 5533384), (5532822, 5532984), (5530481, 5530704), (5521524, 5521786), (5520881, 5521012), (5517746, 5517887), (5515473, 5515517), (5514127, 5515073)]
		#>>> output_gene_exons_coordinates
		#[(5584293, 5584512), (5583687, 5584103), (5582670, 5582846), (5581859, 5582062), (5558339, 5560043), (5553386, 5553556), (5541857, 5542027), (5539415, 5539585), (5536851, 5536940), (5533897, 5533988), (5533304, 5533384), (5532822, 5532984), (5530481, 5530704), (5521524, 5521786), (5520881, 5521012), (5517746, 5517887), (5515473, 5515517), (5514754, 5515073)]
		#for i in range(len(output_gene_exons_coordinates)):
		#	if output_gene_exons_coordinates[i] in output_transcript_exons_coordinates:
		#		if i == int(input_row['gene_splice_target_number'])-1:
		#			print("yes")
		#			break
		#		else:
		#			print("no")
		#			print(input_row)
	start_of_splice_exon = output_gene_exons_coordinates[int(input_row['gene_splice_target_number'])-1][0]
	stop_of_splice_exon = output_gene_exons_coordinates[int(input_row['gene_splice_target_number'])-1][1]
	output_transcript_TE_exons_coordinates = []
	output_transcript_gene_exons_coordinates = []
	if strand == "+":
		for i in range(len(output_transcript_exons_coordinates)):
			start = output_transcript_exons_coordinates[i][0]
			stop = output_transcript_exons_coordinates[i][1]
			if stop <= start_of_splice_exon:
				output_transcript_TE_exons_coordinates.append((start, stop))
			elif start < start_of_splice_exon and start_of_splice_exon < stop:
				output_transcript_TE_exons_coordinates.append((start,start_of_splice_exon-1))
				output_transcript_gene_exons_coordinates.append((start_of_splice_exon,stop))
			elif start_of_splice_exon <= start:
				output_transcript_gene_exons_coordinates.append((start, stop))
	if strand == "-":
		for i in range(len(output_transcript_exons_coordinates)):
			start = output_transcript_exons_coordinates[i][0]
			stop = output_transcript_exons_coordinates[i][1]
			if stop_of_splice_exon <= start:
				output_transcript_TE_exons_coordinates.append((start,stop))
			elif start < stop_of_splice_exon and stop_of_splice_exon < stop:
				output_transcript_TE_exons_coordinates.append((stop_of_splice_exon+1,stop))
				output_transcript_gene_exons_coordinates.append((start,stop_of_splice_exon))
			elif stop <= stop_of_splice_exon:
				output_transcript_gene_exons_coordinates.append((start, stop))
	input_row['transcript_TE_exons_coord'] = output_transcript_TE_exons_coordinates
	input_row['transcript_gene_exons_coord'] = output_transcript_gene_exons_coordinates
	return(input_row)
































