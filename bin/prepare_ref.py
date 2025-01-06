# TEProf3: TE-derived Promoter Finder 3
# Author: Holden Liang
# this script is used to prepare reference files for TEProf3 run
# -----------------------------------------------------------------------------
import argparse
import time
import subprocess # https://geekflare.com/python-run-bash/
import os
import multiprocess as mp
from itertools import repeat
from tqdm import tqdm
import genomepy
import pandas as pd
import ast
import sys
import json

import misc

def prepare_repeat(input_file):
	## start generating a repeatmasker.txt file for downstream analysis
	misc.print_time("Generating reference file for transposable element annotation")

	## stop if there's already a repeatmasker.txt file in the reference folder
	if os.path.isfile('./repeatmasker.txt') is True:
		print("There's a repeatmasker.txt file in this folder, please rename/remove the repeatmasker.txt file before making a new one")
		exit()

	with open(input_file, "r") as file, open('./repeatmasker.txt', 'w') as output_file:
		for skip_header in range(3): # number of header rows
			file.readline()
		for line in file:
			entry = line.strip('\n').split()
			## get strand
			if entry[8] == "+":
				strand = "+"
			elif entry[8] == "C":
				strand = "-"
			## get family and class
			if "/" in entry[10]:
				TE_class = entry[10].split('/')[0]
				TE_family = entry[10].split('/')[1]
			elif "/" not in entry[10]:
				TE_class = entry[10]
				TE_family = entry[10]
			if TE_class in ["DNA", "LINE", "LTR", "Retroposon", "SINE"]:
				# format of repeatmasker file:
				#0 1306    = Smith-Waterman score of the match, usually complexity adjusted
				#       The SW scores are not always directly comparable. Sometimes
				#       the complexity adjustment has been turned off, and a variety of
				#       scoring-matrices are used.
				#1 15.6    = % substitutions in matching region compared to the consensus
				#2 6.2     = % of bases opposite a gap in the query sequence (deleted bp)
				#3 0.0     = % of bases opposite a gap in the repeat consensus (inserted bp)
				#4 HSU08988 = name of query sequence (chromosome)
				#5 6563    = starting position of match in query sequence
				#6 7714    = ending position of match in query sequence
				#7 (22462) = no. of bases in query sequence past the ending position of match
				#8 C       = match is with the Complement of the consensus sequence in the database
				#9 MER7A   = name of the matching interspersed repeat (subfamily)
				#10 DNA/MER2_type = the class of the repeat, in this case a DNA transposon 
				#           fossil of the MER2 group (see below for list and references)
				#11 (0)     = no. of bases in (complement of) the repeat consensus sequence 
				#           prior to beginning of the match (so 0 means that the match extended 
				#           all the way to the end of the repeat consensus sequence)
				#12 2418    = starting position of match in database sequence (using top-strand numbering)
				#13 1465    = ending position of match in database sequence					
				# chr, start, stop, subfamily, %substitutions, strand, family, class
				output_file.write('\t'.join([entry[4], entry[5], entry[6], entry[9], entry[1], strand, TE_family, TE_class]) + "\n")
	#subprocess.run("sort -k1,1 -k2,2n repeatmasker.txt | bgzip > repeatmasker_sorted.txt.gz; tabix -p bed repeatmasker_sorted.txt.gz", shell=True)
	subprocess.run("sort -k1,1 -k2,2n repeatmasker.txt > repeatmasker_sorted.txt", shell=True)
	subprocess.run("rm repeatmasker.txt", shell=True)
	misc.print_time("Reference file for repeatmasker is generated")

def prepare_gene_annotation(input_file):
	## start generating a bed file of gencode information for downstream analysis
	## each row is one exon
	misc.print_time("Generating reference file for gene annotation")

	## stop if there's already a gene_annotation files in the reference folder
	if os.path.isfile('./gene_annotation.txt') is True:
		print("There are gene_annotation files in this folder, please rename/remove them before making a new one")
		exit()

	## extract information from the input gtf file
	subprocess.run("ln -s "+input_file+" gene_annotation.gtf", shell=True)
	subprocess.run("awk '$3~/transcript/ {print $0}' "+input_file+" > "+input_file.replace("gtf","temp.gtf"), shell=True)
	subprocess.run("awk '$3~/start_codon/ {print $0}' "+input_file+" >> "+input_file.replace("gtf","temp.gtf"), shell=True)
	subprocess.run("awk '$3~/exon/ {print $0}' "+input_file+" >> "+input_file.replace("gtf","temp.gtf"), shell=True)

	count_line = subprocess.run("wc -l " + input_file.replace("gtf","temp.gtf"), shell = True, stdout=subprocess.PIPE, text=True)
	number_of_line = int(count_line.stdout.split()[0])

	with open(input_file.replace("gtf","temp.gtf"), "r") as file, open('./gene_exon_annotation.txt', 'w') as output_file_exon, open('./gene_transcript_annotation.txt','w') as output_file_transcript, open('./gene_start_codon_annotation.txt','w') as output_file_start_codon:
		intron_dict = {}
		exon_dict = {}
		for line in tqdm(file, total=number_of_line, desc="Preparing exon, transcript and start_codon information:"):
			entry = line.strip("\n").split("\t")
			chromosome = entry[0]
			start = entry[3]
			stop = entry[4]
			strand = entry[6]
			gene_id = entry[8].split("gene_id \"")[1].split("\";")[0]
			transcript_id = entry[8].split("transcript_id \"")[1].split("\";")[0]

			## ensembl annotation
			if entry[1] == "ensembl":
				gene_type = entry[8].split("gene_biotype \"")[1].split("\";")[0]
				gene_name = gene_id
				transcript_name = transcript_id
			## gencode annotation
			elif entry[1] == "HAVANA" or entry[1] == "ENSEMBL":
				gene_type = entry[8].split("gene_type \"")[1].split("\";")[0]
				gene_name = entry[8].split("gene_name \"")[1].split("\";")[0]
				transcript_name = entry[8].split("transcript_name \"")[1].split("\";")[0]
			else:
				misc.print_time("This is a new gene annotation format. Please report this issue on Github.")
				exit()
			transcript_type = get_simplified_transcript_type(gene_type)

			if entry[2] == "transcript":
				intron_dict[transcript_id] = [chromosome,[],[], gene_name, ".", strand, transcript_id, gene_id, gene_type, transcript_type, transcript_name] # save start and stop of introns
				exon_dict[transcript_id] = [strand, [], []] # save start and stop of exons
				output_file_transcript.write('\t'.join([chromosome, start, stop, gene_name, ".", strand, transcript_id, "transcript", gene_id, gene_type, transcript_type, transcript_name, "."])+"\n")
			if entry[2] == "start_codon":
				output_file_start_codon.write('\t'.join([transcript_id, start, stop])+"\n")
			if entry[2] == "exon":
				exon_number = entry[8].split("exon_number ")[1].split(";")[0]
				output_file_exon.write('\t'.join([chromosome, start, stop, gene_name, ".", strand, transcript_id, "exon_"+str(exon_number), gene_id, gene_type, transcript_type, transcript_name, "."])+"\n")
				intron_dict[transcript_id][1].append(int(start)) # start of exons
				intron_dict[transcript_id][2].append(int(stop)) # stop of exons
				exon_dict[transcript_id][1].append((int(start), int(stop)))
				exon_dict[transcript_id][2].append(int(start))
				exon_dict[transcript_id][2].append(int(stop))

	with open('./gene_intron_annotation.txt','w') as output_file_intron:
		for transcript_id in tqdm(intron_dict, total=len(intron_dict), desc="Preparing intron annotation"):
			information = intron_dict[transcript_id]
			chromosome = information[0]
			gene_name = information[3]
			strand = information[5]
			gene_id = information[7]
			gene_type = information[8]
			transcript_type = information[9]
			transcript_name = information[10]
			starts = sorted(information[1])
			stops = sorted(information[2])
			if len(starts) > 1:
				for index in range(len(starts)-1): # 4 exons => 0,1,2
					if strand == "+":
						intron_number = index+1 # 1,2,3
					elif strand == "-":
						intron_number = len(starts)-index-1 # 3,2,1
					start = str(stops[index]+1)
					stop = str(starts[index+1]-1)
					output_file_intron.write('\t'.join([chromosome, start, stop, gene_name, ".", strand, transcript_id, "intron_"+str(intron_number), gene_id, gene_type, transcript_type, transcript_name, "."])+"\n")		

	with open('./gene_exon_coordinates_annotation.txt', 'w') as output_file_exon, open('./gene_intron_coordinates_annotation.json', 'w') as output_file_intron:
		intron_coordinates_dict = {}
		for transcript_id in tqdm(exon_dict, total=len(exon_dict), desc="Preparing exon and intron coordinates"):
			information = exon_dict[transcript_id]
			strand = information[0]
			coordinates = information[1]
			transcript_start = min(information[2])
			transcript_stop = max(information[2])
			sorted_coordinates = sorted(coordinates, key=lambda x: x[0])

			output_file_exon.write(transcript_id+"\t"+strand+"\t"+str(sorted_coordinates)+"\t"+str(transcript_start)+"\t"+str(transcript_stop)+"\n")

			intron_coordinates = []
			for i in range(len(sorted_coordinates)-1):
				intron_coordinates.append((sorted_coordinates[i][1], sorted_coordinates[i+1][0]))
			intron_coordinates_dict[transcript_id] = [strand, str(intron_coordinates), str(transcript_start), str(transcript_stop)]

		json.dump(intron_coordinates_dict, output_file_intron)

	#subprocess.run("sort -k1,1 -k2,2n gene_annotation.txt | bgzip > gene_annotation_sorted.txt.gz; tabix -p bed gene_annotation_sorted.txt.gz", shell=True)
	subprocess.run("sort -k1,1 -k2,2n gene_exon_annotation.txt > gene_exon_annotation_sorted.txt", shell=True)
	subprocess.run("cat gene_intron_annotation.txt gene_exon_annotation.txt | sort -k1,1 -k2,2n > gene_exon_intron_annotation_sorted.txt", shell=True)
	subprocess.run("sort -k1,1 -k2,2n gene_transcript_annotation.txt > gene_transcript_annotation_sorted.txt", shell=True)
	subprocess.run("sort -k1,1 -k2,2n gene_start_codon_annotation.txt > gene_start_codon_annotation_sorted.txt", shell=True)

	subprocess.run("rm gene_exon_annotation.txt", shell=True)
	subprocess.run("rm gene_intron_annotation.txt", shell=True)
	subprocess.run("rm gene_transcript_annotation.txt", shell=True)
	subprocess.run("rm gene_start_codon_annotation.txt", shell=True)
	subprocess.run("rm "+input_file.replace("gtf","temp.gtf"), shell=True)
	misc.print_time("Reference file for gene annotation is generated")

def prepare_gene_annotation_for_translation():
	## this function removes TE-derived isoforms from the gene annotation
	if (os.path.isfile('gene_transcript_annotation_sorted.txt') is False) or (os.path.isfile('repeatmasker_sorted.txt') is False) or (os.path.isfile('gene_exon_intron_annotation_sorted.txt') is False) or (os.path.isfile('gene_exon_annotation_sorted.txt') is False):
		misc.print_time("Please prepare annotation for repeatmasker and gene annotation before running this command.")
		exit()

	misc.print_time("You are preparing gene annotation for tranlsation.")
	transcript_panda = pd.read_table("gene_transcript_annotation_sorted.txt", names=["chr","start","stop","gene_name","useless","strand","transcript_id","useless2","gene_id","transcript_type","transcript_class","gene_name_2","useless3"])
	transcript_panda[["tss","tss_next_base"]] = transcript_panda.apply(find_tss, axis=1)
	transcript_panda = transcript_panda[["chr","tss","tss_next_base","transcript_id"]]
	transcript_panda.to_csv("gene_transcript_annotation_sorted.temp.txt", sep="\t", index=False, header=False)

	misc.run_bash_command("bedtools intersect -a gene_transcript_annotation_sorted.temp.txt -b repeatmasker_sorted.txt -wa -wb > gene_transcript_annotation_sorted.temp.intersected_with_TE.txt")
	transcript_panda = pd.read_table("gene_transcript_annotation_sorted.temp.intersected_with_TE.txt", names=[f"col{i}" for i in range(1, 13)])
	transcript_id_of_TE_derived_isoform = transcript_panda["col4"].tolist()

	gene_exon_intron_annotation_sorted = pd.read_table("gene_exon_intron_annotation_sorted.txt", names=[f"col{i}" for i in range(1, 14)])
	gene_exon_intron_annotation_sorted = gene_exon_intron_annotation_sorted[~gene_exon_intron_annotation_sorted['col7'].isin(transcript_id_of_TE_derived_isoform)]
	gene_exon_intron_annotation_sorted.to_csv("gene_exon_intron_annotation_sorted.for_translation.txt", sep="\t", index=False, header=False)

	gene_exon_annotation_sorted = pd.read_table("gene_exon_annotation_sorted.txt", names=[f"col{i}" for i in range(1, 14)])
	gene_exon_annotation_sorted = gene_exon_annotation_sorted[~gene_exon_annotation_sorted['col7'].isin(transcript_id_of_TE_derived_isoform)]
	gene_exon_annotation_sorted.to_csv("gene_exon_annotation_sorted.for_translation.txt", sep="\t", index=False, header=False)

	subprocess.run("rm gene_transcript_annotation_sorted.temp.txt", shell=True)
	subprocess.run("rm gene_transcript_annotation_sorted.temp.intersected_with_TE.txt", shell=True)
	misc.print_time("Reference file for translation is generated")


def find_tss(input_row):
	if input_row["strand"] == "+":
		return pd.Series([int(input_row["start"]), int(input_row["start"]) + 1], index=["tss", "tss_next_base"])
	elif input_row["strand"] == "-":
		return pd.Series([int(input_row["stop"]) - 1, int(input_row["stop"])], index=["tss", "tss_next_base"])


def get_simplified_transcript_type(transcript_type):
	transcript_type_list_coding = ["IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene", "polymorphic_pseudogene", "protein_coding", "IG_LV_gene"]
	transcript_type_list_pseudo = ["rRNA_pseudogene", "translated_processed_pseudogene", "translated_unprocessed_pseudogene", "IG_D_pseudogene", "IG_pseudogene", "transcribed_processed_pseudogene", "unitary_pseudogene", "unprocessed_pseudogene", "transcribed_unprocessed_pseudogene", "transcribed_unitary_pseudogene", "IG_C_pseudogene", "IG_J_pseudogene", "IG_V_pseudogene", "TR_J_pseudogene", "TR_V_pseudogene", "pseudogene", "processed_pseudogene", "processed_pseudogene"]
	transcript_type_list_nonCoding = ["non_coding", "antisense_RNA", "ribozyme", "sRNA", "3prime_overlapping_ncRNA", "macro_lncRNA", "bidirectional_promoter_lncRNA", "scaRNA", "scRNA", "snoRNA", "snRNA", "sense_overlapping", "sense_intronic", "3prime_overlapping_ncrna", "Mt_rRNA", "Mt_tRNA", "antisense", "lincRNA", "miRNA", "misc_RNA", "processed_transcript","rRNA","lncRNA"]
	transcript_type_list_problem = ["TEC"]
	if transcript_type in transcript_type_list_coding:
		return("coding")
	elif transcript_type in transcript_type_list_pseudo:
		return("pseudo")
	elif transcript_type in transcript_type_list_nonCoding:
		return("nonCoding")
	elif transcript_type in transcript_type_list_problem:
		return("problem")
	else:
		return("other")
	
def prepare_herv_annotation(input_file):
	misc.print_time("prepare herv annotation file: "+input_file)
	with open(input_file, 'r') as input_file, open("herv_annotation.txt", "w") as output_file:
		for line in input_file:
			if line[0] == "#":
				gene_flag = 0
				exon_flag = 0
			elif line[0] != "#":
				entry=line.strip("\n").split("\t")
				if entry[2] == "gene" and gene_flag == 0:
					chromosome = entry[0]
					start = entry[3]
					stop  = entry[4]
					strand = entry[6]
					gene_flag = 1
				elif entry[2] == "exon" and exon_flag == 0:
					repFamily = entry[8].split("repFamily \"")[1].split("\";")[0]
					locus = entry[8].split("locus \"")[1].split("\";")[0]
					exon_flag = 1
				if exon_flag == 1 and gene_flag == 1:
					if "alt" not in chromosome and "random" not in chromosome and "Un" not in chromosome:
						output_file.write("\t".join([chromosome, start, stop, locus, repFamily, strand])+"\n")
						gene_flag = 0
						exon_flag = 0


def get_fasta_from_gene_annotation(flags):
	misc.print_time("start generate fasta file with protein sequence from gene annotation")
	# generate fasta file of protein sequences from allORFs translation from gene annotation file to identify all possible in-frame protein sequences 
	translation_length = flags.translationlength # translation_length=20
	genome_name = flags.translationgenome # genome_name = "hg38"

	genomepy_folder = os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/genomepy" # genomepy_folder = "/bar/yliang/genomes/private/genomepy/"
	genome = genomepy.Genome(genome_name, genomes_dir=genomepy_folder)

	transcript_panda = pd.read_table(os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/gene_transcript_annotation_sorted.txt", names=["chr","start","stop","gene_name","useless","strand","transcript_id","useless2","gene_id","transcript_type","transcript_class","gene_name_2","useless3"])
	transcript_panda.drop(['useless', 'useless2', 'useless3'], axis=1, inplace=True)
	exon_panda = pd.read_table(os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/gene_exon_coordinates_annotation.txt", names=["transcript_id","uesless","exon_coordinates","transcript_start","transcript_stop"])
	exon_panda = exon_panda.drop("uesless", axis=1)	

	input_panda = pd.merge(exon_panda, transcript_panda, on="transcript_id")
	input_panda["transcript_exons_rna"] = input_panda.apply(get_RNA_sequence_for_gene_annotation, genome=genome, axis=1)
	input_panda["start_codons"] = input_panda["transcript_exons_rna"].apply(find_start_codon)
	output_panda = translation(input_panda, translation_length)
	with open("gene_protein_sequence.fa", "w") as protein_sequence_file:
		output_panda.apply(get_fasta, protein_sequence_file=protein_sequence_file, axis=1)

def get_RNA_sequence_for_gene_annotation(input_row, genome):
	# input_row = input_panda[input_panda["transcript_id"]=="TCONS_00016150"]
	chromosome = input_row["chr"]
	strand = input_row["strand"]
	transcript_exons_coord = str(input_row["exon_coordinates"])
	if strand == "+":
		exons_coord = ast.literal_eval(transcript_exons_coord)
		transcript_exons_rna = str(genome.get_spliced_seq(chromosome, exons_coord))
	elif strand == "-":
		exons_coord = ast.literal_eval(transcript_exons_coord)
		transcript_exons_rna = str(genome.get_spliced_seq(chromosome, exons_coord, rc=True))
	return(transcript_exons_rna)

def find_start_codon(input_sequence):
	input_sequence = input_sequence.upper()
	start_codons = []
	for i in range(len(input_sequence)-2):
		if input_sequence[i:i+3] == "ATG":
			start_codons.append(i)
	return(start_codons)

def translation(input_panda, translation_length):
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

def get_fasta(input_row, protein_sequence_file):
	## create fasta file with protein sequences
	protein_sequence_file.write("> "+input_row["transcript_id"]+"_"+str(input_row["start_codon"])+"\n"+input_row["protein_sequence"]+"\n")
	return(None)


