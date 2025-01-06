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
pd.options.mode.chained_assignment = None
from tqdm import tqdm
import glob
import numpy as np
import genomepy # genomepy.install_genome("hg38", annotation=True, provider="UCSC", genomes_dir="/bar/yliang/genomes/private/genomepy")
import math
import ast
import matplotlib.pyplot as plt
import seaborn as sns

import misc

def run_split_on_sbatch_cluster(input_fasta_file): # input_fasta_file="protein_sequence_TE.fa"
	misc.print_time("Submitted split jobs to cluster")
	if not os.path.exists('subsampled_fasta'):
		misc.run_bash_command("mkdir subsampled_fasta")
	with open(input_fasta_file+"_split_commands.txt","w") as commands_file:
		commands_file.write('\n'.join(["jid=`sbatch <<- SPLIT | egrep -o -e \"\\b[0-9]+$\"","#!/bin/bash -l","#SBATCH -o "+input_fasta_file+".split.out",\
				"#SBATCH -e "+input_fasta_file+".split.err",\
				"#SBATCH -c 10","#SBATCH --ntasks=1",\
				"#SBATCH --mem=5G",\
				"#SBATCH -J split_"+input_fasta_file,\
				"date",\
				"split -l 100 "+input_fasta_file+" subsampled_fasta/subsampled_"+input_fasta_file,\
				"date",\
				"SPLIT`"])+"\n")
	misc.run_bash_command("bash "+input_fasta_file+"_split_commands.txt")

def run_blastp_on_sbatch_cluster(input_fasta_file, database, taxiddb, blastjobnum): # input_fasta_file="protein_sequence_TE.fa"
	## blast parameters: https://www.ncbi.nlm.nih.gov/books/NBK279684/
	misc.print_time("Submitting blastp jobs to cluster")
	subsampled_fasta_files = glob.glob("subsampled_fasta/subsampled_"+input_fasta_file+"??")
	number_of_jobs = str(len(subsampled_fasta_files))

	if database == "nr":
		blastp_database_folder = os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/blast_database_nr/nr"
	elif database == "gene":
		blastp_database_folder = os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/blast_database_gene_annotation/gene_protein_sequence"
	elif database == "uniprot":
		blastp_database_folder = os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/blast_database_uniprot/uniprot"
	elif database == "gencode":
		blastp_database_folder = os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/blast_database_gencode/gencode"
	elif database == "gencode_normal_tissue":
		blastp_database_folder = os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/blast_database_gencode_normal_tissue/gencode_normal_tissue"
	else:
		exit()

	if database == "nr":
		command = "blastp -query ${subsampled_fasta_files[$((SLURM_ARRAY_TASK_ID-1))]} -db "+blastp_database_folder+" -task 'blastp' -outfmt \"6 qseqid qlen slen sseqid salltitles sscinames qstart qend sstart send evalue length pident qcovs mismatch gapopen bitscore\" -num_threads 10 -word_size 5 -evalue 0.05 -max_target_seqs 5000 -taxids 9606 > ${subsampled_fasta_files[$((SLURM_ARRAY_TASK_ID-1))]}."+database+".blastp.out"
	elif database in ["gene","uniprot","gencode","gencode_normal_tissue"]:
		command = "blastp -query ${subsampled_fasta_files[$((SLURM_ARRAY_TASK_ID-1))]} -db "+blastp_database_folder+" -task 'blastp' -outfmt \"6 qseqid qlen slen sseqid salltitles sscinames qstart qend sstart send evalue length pident qcovs mismatch gapopen bitscore\" -num_threads 10 -word_size 5 -evalue 0.05 -max_target_seqs 5000 > ${subsampled_fasta_files[$((SLURM_ARRAY_TASK_ID-1))]}."+database+".blastp.out"

	with open(input_fasta_file+"_"+database+"_blastp_commands.txt","w") as command_file:
		command_file.write('\n'.join(["jid=$(sbatch <<'BLASTP' | egrep -o -e \"\\b[0-9]+$\"",\
			"#!/bin/bash -l",\
			"#SBATCH --array=1-"+number_of_jobs+"%"+blastjobnum,\
			"#SBATCH -o subsampled_fasta/"+input_fasta_file+".%A.%a."+database+".blastp.log",\
			"#SBATCH -e subsampled_fasta/"+input_fasta_file+".%A.%a."+database+".blastp.err",\
			"#SBATCH -c 10",\
			"#SBATCH --ntasks=1",\
			"#SBATCH --mem=5G",\
			"#SBATCH -J blastp_"+input_fasta_file+"_"+database,\
			"date",\
			"export BLASTDB=\""+taxiddb+"\"",\
			"subsampled_fasta_files=($(ls -1 subsampled_fasta/subsampled_"+input_fasta_file+"??))",\
			"echo ${subsampled_fasta_files[$((SLURM_ARRAY_TASK_ID-1))]}",\
			command,\
			"date",
			"BLASTP",\
			")"])+"\n")

	misc.run_bash_command("bash "+input_fasta_file+"_"+database+"_blastp_commands.txt")


def run_blastp_short_on_sbatch_cluster(input_fasta_file, database, taxiddb, blastjobnum): # input_fasta_file="protein_sequence_TE.fa"
	## blast parameters: https://www.ncbi.nlm.nih.gov/books/NBK279684/
	## qcovs means Query Coverage Per Subject (for all HSPs)
	## /scratch/yliang/HNSCC/analysis/CPTAC_TSTEA_candidates/2_test2_actual_run_TEProf2/blastp-short_setting_test
	misc.print_time("Submitting blastp jobs to cluster")
	subsampled_fasta_files = glob.glob("subsampled_fasta/subsampled_"+input_fasta_file+"??")
	number_of_jobs = str(len(subsampled_fasta_files))

	if database == "nr":
		blastp_database_folder = os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/blast_database_nr/nr"
	elif database == "gene":
		blastp_database_folder = os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/blast_database_gene_annotation/gene_protein_sequence"
	elif database == "uniprot":
		blastp_database_folder = os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/blast_database_uniprot/uniprot"
	elif database == "gencode":
		blastp_database_folder = os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/blast_database_gencode/gencode"
	elif database == "gencode_normal_tissue":
		blastp_database_folder = os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/blast_database_gencode_normal_tissue/gencode_normal_tissue"
	else:
		exit()

	if database == "nr":
		command = "blastp -query ${subsampled_fasta_files[$((SLURM_ARRAY_TASK_ID-1))]} -db "+blastp_database_folder+" -task 'blastp-short' -outfmt \"6 qseqid qlen slen sseqid salltitles sscinames qstart qend sstart send evalue length pident qcovs mismatch gapopen bitscore\" -num_threads 10 -threshold 11 -window_size 40 -evalue 200000 -max_target_seqs 5000 -taxids 9606 > ${subsampled_fasta_files[$((SLURM_ARRAY_TASK_ID-1))]}."+database+".blastp.out"
	elif database in ["gene","uniprot","gencode","gencode_normal_tissue"]:
		command = "blastp -query ${subsampled_fasta_files[$((SLURM_ARRAY_TASK_ID-1))]} -db "+blastp_database_folder+" -task 'blastp-short' -outfmt \"6 qseqid qlen slen sseqid salltitles sscinames qstart qend sstart send evalue length pident qcovs mismatch gapopen bitscore\" -num_threads 10 -threshold 11 -window_size 40 -evalue 200000 -max_target_seqs 5000 > ${subsampled_fasta_files[$((SLURM_ARRAY_TASK_ID-1))]}."+database+".blastp.out"

	with open(input_fasta_file+"_"+database+"_blastp_commands.txt","w") as command_file:
		command_file.write('\n'.join(["jid=$(sbatch <<'BLASTP' | egrep -o -e \"\\b[0-9]+$\"",\
			"#!/bin/bash -l",\
			"#SBATCH --array=1-"+number_of_jobs+"%"+blastjobnum,\
			"#SBATCH -o subsampled_fasta/"+input_fasta_file+".%A.%a."+database+".blastp.log",\
			"#SBATCH -e subsampled_fasta/"+input_fasta_file+".%A.%a."+database+".blastp.err",\
			"#SBATCH -c 10",\
			"#SBATCH --ntasks=1",\
			"#SBATCH --mem=5G",\
			"#SBATCH -J blastp_"+input_fasta_file+"_"+database,\
			"date",\
			"export BLASTDB=\""+taxiddb+"\"",\
			"subsampled_fasta_files=($(ls -1 subsampled_fasta/subsampled_"+input_fasta_file+"??))",\
			"echo ${subsampled_fasta_files[$((SLURM_ARRAY_TASK_ID-1))]}",\
			command,\
			"date",
			"BLASTP",\
			")"])+"\n")

	misc.run_bash_command("bash "+input_fasta_file+"_"+database+"_blastp_commands.txt")


def classify_protein_products(input_file):
	misc.print_time("Start classifying TE-derived proteins")
	misc.print_time("Concatenate blastp output")

	commands_list = ["cat subsampled_fasta/subsampled_protein_sequence_TE*nr.blast.out > subsampled_fasta/subsampled_protein_sequence_TE.nr.blast.out",\
	"cat subsampled_fasta/subsampled_protein_sequence_TE*gene.blast.out > subsampled_fasta/subsampled_protein_sequence_TE.gene.blast.out",\
	"cat subsampled_fasta/subsampled_protein_sequence_gene*nr.blast.out > subsampled_fasta/subsampled_protein_sequence_gene.nr.blast.out",\
	"cat subsampled_fasta/subsampled_protein_sequence_gene*gene.blast.out > subsampled_fasta/subsampled_protein_sequence_gene.gene.blast.out",\
	"cat subsampled_fasta/subsampled_protein_sequence_blast*nr.blast.out > subsampled_fasta/subsampled_protein_sequence.nr.blast.out"]
	with mp.Pool(5) as pool:
		uesless = pool.map(misc.run_bash_command, commands_list)

	commands_list = ["grep -P \"Homo sapiens\" subsampled_fasta/subsampled_protein_sequence_TE.nr.blast.out | cut -f 1,2,4,7,8,13,14 > subsampled_fasta/subsampled_protein_sequence.nr.blast.out.queryid_subjectid_pident_qcovs;\
	grep -P \"Homo sapiens\" subsampled_fasta/subsampled_protein_sequence_gene.nr.blast.out | cut -f 1,2,4,7,8,13,14 >> subsampled_fasta/subsampled_protein_sequence.nr.blast.out.queryid_subjectid_pident_qcovs",\
	"cut -f 1,2,4,7,8,13,14 subsampled_fasta/subsampled_protein_sequence_gene.gene.blast.out >  subsampled_fasta/subsampled_protein_sequence.gene.blast.out.queryid_subjectid_pident_qcovs;\
	cut -f 1,2,4,7,8,13,14 subsampled_fasta/subsampled_protein_sequence_TE.gene.blast.out   >> subsampled_fasta/subsampled_protein_sequence.gene.blast.out.queryid_subjectid_pident_qcovs",\
	"grep -P \"Homo sapiens\" subsampled_fasta/subsampled_protein_sequence.nr.blast.out  | cut -f 1,2,3,4,7,8,9,10,13,14 > subsampled_fasta/subsampled_protein_sequence_whole_protein.nr.blast.out.queryid_subjectid_pident_qcovs"]
	with mp.Pool(3) as pool:
		uesless = pool.map(misc.run_bash_command, commands_list)

	misc.print_time("Read in blastp result")
	blastp_nr = pd.read_table("subsampled_fasta/subsampled_protein_sequence.nr.blast.out.queryid_subjectid_pident_qcovs", names=["query_id","query_len","subject_id","query_start","query_end","pident","qcovs"])
	blastp_gene = pd.read_table("subsampled_fasta/subsampled_protein_sequence.gene.blast.out.queryid_subjectid_pident_qcovs", names=["query_id","query_len","subject_id","query_start","query_end","pident","qcovs"])

	## parse result from blastp against nr
	blastp_nr['blastp_result'] = blastp_nr['pident'] + blastp_nr['qcovs']
	blastp_nr = blastp_nr[blastp_nr['blastp_result'] == 200]
	blastp_nr['temp'] = blastp_nr["query_len"] - blastp_nr["query_end"]
	blastp_nr = blastp_nr[(blastp_nr["query_start"]==1) & (blastp_nr["temp"]==0)]
	blastp_nr[['col1','col2','protein_part','col4']] = blastp_nr['query_id'].str.split('_', expand=True)
	blastp_nr['protein_id'] = blastp_nr[['col1', 'col2', 'col4']].apply(lambda x: "_".join(x), axis=1)
	blastp_nr = blastp_nr[["protein_id","protein_part"]]
	protein_id_te_nr = list(set(blastp_nr[blastp_nr['protein_part'] == "TE"]["protein_id"]))
	protein_id_gene_nr = list(set(blastp_nr[blastp_nr['protein_part'] == "gene"]["protein_id"]))

	## parse result from blastp against gencode
	blastp_gene['blastp_result'] = blastp_gene['pident'] + blastp_gene['qcovs']
	blastp_gene = blastp_gene[blastp_gene['blastp_result'] == 200]
	blastp_gene['temp'] = blastp_gene["query_len"] - blastp_gene["query_end"]
	blastp_gene = blastp_gene[(blastp_gene["query_start"]==1) & (blastp_gene["temp"]==0)]
	blastp_gene[['col1','col2','protein_part','col4']] = blastp_gene['query_id'].str.split('_', expand=True)
	blastp_gene['protein_id'] = blastp_gene[['col1', 'col2', 'col4']].apply(lambda x: "_".join(x), axis=1)
	blastp_gene = blastp_gene[["protein_id","protein_part"]]
	protein_id_te_gene = list(set(blastp_gene[blastp_gene['protein_part'] == "TE"]["protein_id"]))
	protein_id_gene_gene = list(set(blastp_gene[blastp_gene['protein_part'] == "gene"]["protein_id"]))

	## parse database information
	protein_all_info = pd.read_table(input_file) # input_file="/scratch/yliang/HNSCC/data/CPTAC/2_test2_actual_run_TEProf2/test2_database/teprof3_protein_information.tsv"
	protein_all_info['start_codon'] = protein_all_info['start_codon'].astype(str)
	protein_all_info['protein_id'] = protein_all_info[['transcript_id','start_codon']].apply(lambda x: "_".join(x), axis=1)
	protein_info = protein_all_info[["protein_id","protein_sequence","protein_sequence_TE","protein_sequence_gene"]]

	## combine results and get protein type distribution of database and save it to a table
	protein_info["blastp_te_to_nr"] = protein_info["protein_id"].apply(lambda x: "yes" if(x in protein_id_te_nr) else "no")
	protein_info["blastp_gene_to_nr"] = protein_info["protein_id"].apply(lambda x: "yes" if(x in protein_id_gene_nr) else "no")
	protein_info["blastp_te_to_gene"] = protein_info["protein_id"].apply(lambda x: "yes" if(x in protein_id_te_gene) else "no")
	protein_info["blastp_gene_to_gene"] = protein_info["protein_id"].apply(lambda x: "yes" if(x in protein_id_gene_gene) else "no")
	protein_info['protein_sequence_gene'] = protein_info['protein_sequence_gene'].astype(str)
	protein_info["blastp_gene"] = protein_info.apply(get_blastp_gene_part_result, axis=1)
	protein_info["blastp_te"] = protein_info["blastp_te_to_nr"]
	protein_info["protein_type"] = protein_info.apply(get_protein_type_result, axis=1)
	protein_info = protein_info[["protein_id","protein_sequence","protein_sequence_TE","protein_sequence_gene","protein_type"]]
	protein_info.to_csv("./teprof3_database_classification.tsv", sep="\t", index=False)

	## visualize result as piechart
	protein_info['protein_type_count'] = protein_info.groupby('protein_type')['protein_type'].transform('count')
	protein_info = protein_info.drop_duplicates(subset='protein_type')
	protein_info = protein_info.sort_values(by='protein_type_count', ascending=False)
	plt.figure(figsize=(7, 4))
	plt.pie(protein_info['protein_type_count'], labels=None , colors = sns.color_palette("pastel"), autopct='%.0f%%')
	plt.legend(labels = protein_info['protein_type'], title="protein_type", loc="center left", bbox_to_anchor=(0.9, 0.5))
	plt.axis('equal')
	plt.title("Percentage of each protein_type of proteins in database (N= "+str(sum(protein_info['protein_type_count']))+")", loc='left', y=1)
	plt.tight_layout()
	plt.savefig("teprof3_percentage_of_database_protein_type_piechart.pdf", format='pdf', dpi=200)

	## find TE-derived proteins that have an exact match in nr
	blast_whole_protein = pd.read_table("subsampled_fasta/subsampled_protein_sequence_whole_protein.nr.blast.out.queryid_subjectid_pident_qcovs", names=["query_id","query_len","subject_len","subject_id","query_start","query_end","subject_start","subject_end","pident","qcovs"])
	blast_whole_protein['blastp_result'] = blast_whole_protein['pident'] + blast_whole_protein['qcovs']
	blast_whole_protein = blast_whole_protein[blast_whole_protein['blastp_result'] == 200]
	blast_whole_protein['temp'] = blast_whole_protein["query_len"] - blast_whole_protein["query_end"]
	blast_whole_protein = blast_whole_protein[(blast_whole_protein["query_start"]==1) & (blast_whole_protein["temp"]==0)]	
	blast_whole_protein['the_same_protein'] = blast_whole_protein.apply(lambda x: "yes" if (x["query_len"] == x["subject_len"]) else "no", axis=1)
	blast_whole_protein = blast_whole_protein[blast_whole_protein["the_same_protein"] == "yes"] 
	blast_whole_protein = blast_whole_protein[["query_id","subject_id"]]
	blast_whole_protein.to_csv("./teprof3_database_blast_to_the_same_protein.tsv", sep="\t", index=False)

	misc.print_time("Finished annotating database")

def get_blastp_gene_part_result(input_row):
	blastp_gene_to_nr = input_row["blastp_gene_to_nr"]
	blastp_gene_to_gene = input_row["blastp_gene_to_gene"]
	if input_row["protein_sequence_gene"]!="nan":
		if blastp_gene_to_nr == "yes":
			return("normal/truncated")
		elif blastp_gene_to_nr == "no" and blastp_gene_to_gene == "yes":
			return("unannotated")
		elif blastp_gene_to_nr == "no" and blastp_gene_to_gene == "no":
			return("frame shift")
	else:
		return("no_gene_protein")

def get_protein_type_result(input_row):
	blastp_gene = input_row["blastp_gene"]
	blastp_te = input_row["blastp_te"]
	if blastp_te=="no"  and blastp_gene=="no_gene_protein":
		return("solo unannotated TE")
	elif blastp_te=="no"  and blastp_gene=="normal/truncated":
		return("Chimeric uTE-normal/truncated")
	elif blastp_te=="no"  and blastp_gene=="unannotated":
		return("Chimeric uTE-unannotated")
	elif blastp_te=="no"  and blastp_gene=="frame shift":
		return("Chimeric uTE-frame shift")
	elif blastp_te=="yes" and blastp_gene=="no_gene_protein":
		return("solo annotated TE")
	elif blastp_te=="yes" and blastp_gene=="normal/truncated":
		return("Chimeric aTE-normal/truncated")
	elif blastp_te=="yes" and blastp_gene=="unannotated":
		return("Chimeric aTE-unannotated")
	elif blastp_te=="yes" and blastp_gene=="frame shift":
		return("Chimeric aTE-frame shift")

def detected_protein():
	## parse detected peptides result
	misc.print_time("Read in info of detected peptides")
	detected_peptides_dict = {}
	with open("detected_peptides_blastp_result/detected_peptides.fa", "r") as fasta_file:
		for line in fasta_file:
			entry = line.strip("\n")
			if entry[0] == ">":
				detected_peptide = entry.split(":")[1]
				if detected_peptide in detected_peptides_dict:
					continue
				else:
					info = entry.split(":")[0].split(" ")[1].split("_")
					protein_id = '_'.join([info[0],info[1],info[3]])
					gene_name = info[4]
					detected_peptides_dict[detected_peptide] = [protein_id, gene_name, detected_peptide]
	detected_peptides = pd.DataFrame.from_dict(detected_peptides_dict, orient="index", columns=["protein_id","gene_name","detected_peptide_seq"])
	detected_peptides = detected_peptides.reset_index(drop=True)

	## parse result from blastp against nr
	misc.print_time("Read in blastp result for detected peptides")
	misc.run_bash_command("cat detected_peptides_blastp_result/subsampled_fasta/*.out > detected_peptides_blastp_result/subsampled_fasta/detected_peptides.nr.blast.output")
	misc.run_bash_command("grep -P \"Homo sapiens\" detected_peptides_blastp_result/subsampled_fasta/detected_peptides.nr.blast.output | cut -f 1,2,4,7,8,13,14 > detected_peptides_blastp_result/subsampled_fasta/detected_peptides.nr.blast.output.queryid_subjectid_pident_qcovs")
	blastp = pd.read_table("detected_peptides_blastp_result/subsampled_fasta/detected_peptides.nr.blast.output.queryid_subjectid_pident_qcovs", names=["query_id","query_len","subject_id","query_start","query_end","pident","qcovs"])
	teprof3_database_blast_to_the_same_protein = pd.read_table("teprof3_database_blast_to_the_same_protein.tsv") # database_info_file_2 = "/scratch/yliang/HNSCC/analysis/CPTAC_TSTEA_candidates/2_test2_actual_run_TEProf2/annotate_database_use_teprof3/teprof3_database_blast_to_the_same_protein.tsv"
	
	blastp['blastp_result'] = blastp['pident'] + blastp['qcovs']
	blastp = blastp[blastp['blastp_result'] == 200]
	blastp['temp'] = blastp["query_len"] - blastp["query_end"]
	blastp = blastp[(blastp["query_start"]==1) & (blastp["temp"]==0)]
	blastp[['col1','col2','col3','col4','col5']] = blastp['query_id'].str.split('_', expand=True)
	blastp[['gene_name','detected_peptide_seq']] = blastp['col5'].str.split(':', expand=True)
	blastp['protein_id'] = blastp[['col1', 'col2', 'col4']].apply(lambda x: "_".join(x), axis=1)
	blastp = blastp[["query_id","subject_id","protein_id","gene_name","detected_peptide_seq"]]
	blastp = blastp.groupby(["query_id","protein_id","gene_name","detected_peptide_seq"])['subject_id'].agg(list).reset_index()
	blastp["subject_id"] = blastp["subject_id"].apply(lambda x: list(set(x)))
	blastp["peptide_blastp_to_one_protein_and_it's_same_as_TEP"] = blastp["subject_id"].apply(lambda x: "yes" if (len(x)==1 and x[0] in list(teprof3_database_blast_to_the_same_protein['subject_id'])) else "no")
	blastp = blastp[blastp["peptide_blastp_to_one_protein_and_it's_same_as_TEP"] == "no"]
	detected_peptides["detected_peptide_blastp_result"] = detected_peptides["detected_peptide_seq"].apply(lambda x: "hit" if (x in list(blastp["detected_peptide_seq"])) else "no_hit")

	## parse result from blat
	misc.print_time("Read in blat result for detected peptides")
	misc.run_bash_command("grep whole detected_peptides_blat_result/detected_peptides.blat.psl > detected_peptides_blat_result/detected_peptides.blat.tsv")
	blat = pd.read_table("detected_peptides_blat_result/detected_peptides.blat.tsv", names=["match","mismatch","repmatch","Ns","QgapCount","QgapBases","TgapCount","TgapBases","strand","Qname","Qsize","Qstart","Qend","Tname","Tsize","Tstart","Tend","blockCount","blockSizes","qStarts","tStarts"])
	blat['temp'] = blat['match'] - blat["Qsize"]
	blat = blat[blat['temp']==0]
	blat = blat[(blat['blockCount']==1) | (blat['blockCount']==2)]
	blat_blockcount_1 = blat[blat['blockCount']==1]
	blat_blockcount_1 = blat_blockcount_1.reset_index(drop=True)
	blat_blockcount_1_count = blat_blockcount_1.groupby('Qname').size().reset_index(name='detected_peptide_blat_count_1block')
	blat = pd.merge(blat, blat_blockcount_1_count, on="Qname")
	blat_blockcount_2 = blat[blat['blockCount']==2]
	blat_blockcount_2 = blat_blockcount_2.reset_index(drop=True)
	blat_blockcount_2_count = blat_blockcount_2.groupby('Qname').size().reset_index(name='detected_peptide_blat_count_2block')
	blat = pd.merge(blat, blat_blockcount_2_count, on="Qname")
	blat.to_csv("detected_peptides_blat_result/teprof3_blat_result.tsv", sep="\t", index=False)
	blat = blat[["Qname","detected_peptide_blat_count_1block",'detected_peptide_blat_count_2block']]
	blat = blat.drop_duplicates(subset='Qname')
	blat[["id","detected_peptide_seq"]] = blat['Qname'].str.split(':', expand=True)
	blat = blat[["detected_peptide_seq","detected_peptide_blat_count_1block",'detected_peptide_blat_count_2block']]

	detected_peptides = pd.merge(detected_peptides, blat, on="detected_peptide_seq", how="left")
	detected_peptides = detected_peptides.fillna(0)
	detected_peptides['detected_peptide_blat_count_1block'] = detected_peptides['detected_peptide_blat_count_1block'].astype(int)
	detected_peptides['detected_peptide_blat_count_2block'] = detected_peptides['detected_peptide_blat_count_2block'].astype(int)

	## parse information of database and check the position of detected peptides
	teprof3_database_classification = pd.read_table("teprof3_database_classification.tsv") # database_info_file_1 = "/scratch/yliang/HNSCC/analysis/CPTAC_TSTEA_candidates/2_test2_actual_run_TEProf2/annotate_database_use_teprof3/teprof3_database_classification.tsv"
	detected_peptides = pd.merge(detected_peptides, teprof3_database_classification, on="protein_id", how="inner")
	detected_peptides["detected_peptide_position"] = detected_peptides.apply(find_peptide_position, peptide_type="detected_peptide", axis=1)
	detected_peptides["detected_peptide_evidence_level"] = detected_peptides.apply(define_peptide_evidence_level, peptide_type="detected_peptide", axis=1)
	detected_peptides = detected_peptides.sort_values(by='detected_peptide_evidence_level', ascending=True)
	detected_peptides = detected_peptides.reset_index(drop=True)
	detected_peptides["detected_peptide_final_call"] = detected_peptides.apply(good_peptide_criteria, axis=1)
	detected_peptides["index"] = detected_peptides.index
	detected_peptides = detected_peptides[['index','protein_id','gene_name','protein_sequence','protein_sequence_TE','protein_sequence_gene','protein_type','detected_peptide_seq','detected_peptide_blastp_result','detected_peptide_blat_count_1block','detected_peptide_blat_count_2block',"detected_peptide_position",'detected_peptide_evidence_level','detected_peptide_final_call']]
	detected_peptides.to_csv("./teprof3_detected_proteins.tsv", sep="\t", index=False)

	with open("./teprof3_neoantigen_identification_statistic.tsv", "w+") as output_file:
		output_file.write("\t# of proteins\t# of peptides\n")
		output_file.write("detected by MassSpec\t"+str(len(set(detected_peptides["protein_id"])))+"\t"+str(len(set(detected_peptides["detected_peptide_seq"])))+"\n")
		detected_peptides = detected_peptides[detected_peptides["detected_peptide_final_call"]=="keep"]
		output_file.write("pass BLAST filter\t"+str(len(set(detected_peptides["protein_id"])))+"\t"+str(len(set(detected_peptides["detected_peptide_seq"])))+"\n")

	## output protein sequence with good peptide evidence support to a fasta file for netMHCpan processing
	with open("./teprof3_detected_proteins.fa", "w") as output_file:
		uesless = detected_peptides.apply(output_fasta_file_for_netmhcpan, output_file=output_file, axis=1)

	misc.print_time("Finish processing detected peptides, now you can run netMHCpan on the detected TE-derived proteins")

def parse_netmhcpan(netmhcpan_output):
	misc.print_time("Start parsing result from netMHCpan")
	## extract SB and WB neoantigen
	commands_list = ["grep -P \"<= SB\" "+netmhcpan_output+" > "+netmhcpan_output+".SB.txt",\
	"grep -P \"<= WB\" "+netmhcpan_output+" > "+netmhcpan_output+".WB.txt"] # netmhcpan_output = "teprof3_detected_proteins.fa_HLA-A02-01.txt"
	with mp.Pool(2) as pool:
		uesless = pool.map(misc.run_bash_command, commands_list)
	SB_peptides = pd.read_table(netmhcpan_output+".SB.txt", names=["position","MHC_allele","neoantigen_seq","neoantigen_core","Of","Gp","Gl","Ip","Il","Icore","protein_id","Socre_EL","%Rank_EL","Score_BA","%Rank_BA","neoantigen_affinity(nM)","useless1","neoantigen_bind_level"], delim_whitespace=True)
	WB_peptides = pd.read_table(netmhcpan_output+".WB.txt", names=["position","MHC_allele","neoantigen_seq","neoantigen_core","Of","Gp","Gl","Ip","Il","Icore","protein_id","Socre_EL","%Rank_EL","Score_BA","%Rank_BA","neoantigen_affinity(nM)","useless1","neoantigen_bind_level"], delim_whitespace=True)
	neoantigens = pd.concat([SB_peptides, WB_peptides], axis=0).reset_index(drop=True)

	## merge protein information with neoantigen table
	detected_peptides_table = pd.read_table("teprof3_detected_proteins.tsv")
	detected_peptides_table['index'] = detected_peptides_table['index'].astype(str)
	neoantigens[['index','useless2','useless3']] = neoantigens['protein_id'].str.split('_', expand=True)
	neoantigens = neoantigens[["MHC_allele","neoantigen_seq","neoantigen_core","neoantigen_affinity(nM)","neoantigen_bind_level","index"]]
	neoantigens['index'] = neoantigens['index'].astype(str)

	neoantigens = pd.merge(neoantigens, detected_peptides_table, on="index", how="inner")
	neoantigens = neoantigens[['protein_id','gene_name','protein_sequence','protein_sequence_TE','protein_sequence_gene','protein_type',\
	'detected_peptide_seq','detected_peptide_blastp_result','detected_peptide_blat_count_1block','detected_peptide_blat_count_2block',"detected_peptide_position",'detected_peptide_evidence_level','detected_peptide_final_call',\
	"MHC_allele","neoantigen_seq","neoantigen_core","neoantigen_affinity(nM)","neoantigen_bind_level"]]
	neoantigens["neoantigen_position"] = neoantigens.apply(find_peptide_position, peptide_type="neoantigen", axis=1)
	neoantigens["neoantigen_evidence_level"] = neoantigens.apply(define_peptide_evidence_level, peptide_type="neoantigen", axis=1)
	neoantigens = neoantigens.sort_values(by=['neoantigen_evidence_level','neoantigen_bind_level'], ascending=True)
	neoantigens = neoantigens.reset_index(drop=True)
	neoantigens.to_csv("./teprof3_neoantigen_candidates_before_blast.tsv", sep="\t", index=False)
	with open("./teprof3_neoantigen_candidates_before_blast.fa", "w") as output_file:
		uesless = neoantigens.apply(output_fasta_file_for_blast_after_netmhcpan, output_file=output_file, axis=1)

	misc.print_time("Finished parsing result from netMHCpan, please run BLAST and BLAT on the potential neoantigens")

def get_neoantigen_candidates():
	neoantigens = pd.read_table("./teprof3_neoantigen_candidates_before_blast.tsv")
	## parse result from blastp against nr
	misc.print_time("Read in blastp result for neoantigens")
	misc.run_bash_command("cat neoantigen_blastp_result/subsampled_fasta/*.out > neoantigen_blastp_result/subsampled_fasta/teprof3_neoantigen_candidates_before_blast.nr.blast.output")
	misc.run_bash_command("grep -P \"Homo sapiens\" neoantigen_blastp_result/subsampled_fasta/teprof3_neoantigen_candidates_before_blast.nr.blast.output | cut -f 1,2,4,7,8,13,14 > neoantigen_blastp_result/subsampled_fasta/teprof3_neoantigen_candidates_before_blast.nr.blast.output.queryid_subjectid_pident_qcovs")
	blastp = pd.read_table("neoantigen_blastp_result/subsampled_fasta/teprof3_neoantigen_candidates_before_blast.nr.blast.output.queryid_subjectid_pident_qcovs", names=["query_id","query_len","subject_id","query_start","query_end","pident","qcovs"])
	teprof3_database_blast_to_the_same_protein = pd.read_table("teprof3_database_blast_to_the_same_protein.tsv") # database_info_file_2 = "/scratch/yliang/HNSCC/analysis/CPTAC_TSTEA_candidates/2_test2_actual_run_TEProf2/annotate_database_use_teprof3/teprof3_database_blast_to_the_same_protein.tsv"
	
	blastp['blastp_result'] = blastp['pident'] + blastp['qcovs']
	blastp = blastp[blastp['blastp_result'] == 200]
	blastp['temp'] = blastp["query_len"] - blastp["query_end"]
	blastp = blastp[(blastp["query_start"]==1) & (blastp["temp"]==0)]
	blastp[['col1','col2','col3','neoantigen_seq']] = blastp['query_id'].str.split('_', expand=True)
	blastp['protein_id'] = blastp[['col1', 'col2', 'col3']].apply(lambda x: "_".join(x), axis=1)
	blastp = blastp[["query_id","subject_id","protein_id","neoantigen_seq"]]
	blastp = blastp.groupby(["query_id","protein_id","neoantigen_seq"])['subject_id'].agg(list).reset_index()
	blastp["subject_id"] = blastp["subject_id"].apply(lambda x: list(set(x)))
	blastp["peptide_blastp_to_one_protein_and_it's_same_as_TEP"] = blastp["subject_id"].apply(lambda x: "yes" if (len(x)==1 and x[0] in list(teprof3_database_blast_to_the_same_protein['subject_id'])) else "no")
	blastp = blastp[blastp["peptide_blastp_to_one_protein_and_it's_same_as_TEP"] == "no"]
	neoantigens["neoantigen_blastp_result"] = neoantigens["neoantigen_seq"].apply(lambda x: "hit" if (x in list(blastp["neoantigen_seq"])) else "no_hit")

	## parse result from blat
	misc.print_time("Read in blat result for neoantigens")
	misc.run_bash_command("grep neoantigen neoantigen_blat_result/teprof3_neoantigen_candidates_before_blast.blat.psl > neoantigen_blat_result/teprof3_neoantigen_candidates_before_blast.blat.tsv")
	blat = pd.read_table("neoantigen_blat_result/teprof3_neoantigen_candidates_before_blast.blat.tsv", names=["match","mismatch","repmatch","Ns","QgapCount","QgapBases","TgapCount","TgapBases","strand","Qname","Qsize","Qstart","Qend","Tname","Tsize","Tstart","Tend","blockCount","blockSizes","qStarts","tStarts"])
	blat['temp'] = blat['match'] - blat["Qsize"]
	blat = blat[blat['temp']==0]
	blat = blat[(blat['blockCount']==1) | (blat['blockCount']==2)]
	blat_blockcount_1 = blat[blat['blockCount']==1]
	blat_blockcount_1 = blat_blockcount_1.reset_index(drop=True)
	blat_blockcount_1_count = blat_blockcount_1.groupby('Qname').size().reset_index(name='neoantigen_blat_count_1block')
	blat = pd.merge(blat, blat_blockcount_1_count, on="Qname")
	blat_blockcount_2 = blat[blat['blockCount']==2]
	blat_blockcount_2 = blat_blockcount_2.reset_index(drop=True)
	blat_blockcount_2_count = blat_blockcount_2.groupby('Qname').size().reset_index(name='neoantigen_blat_count_2block')
	blat = pd.merge(blat, blat_blockcount_2_count, on="Qname")	
	blat.to_csv("neoantigen_blat_result/teprof3_blat_result.tsv", sep="\t", index=False)
	blat = blat[["Qname","neoantigen_blat_count_1block",'neoantigen_blat_count_2block']]
	blat = blat.drop_duplicates(subset='Qname')
	blat[["col1","col2","col3",'col4',"neoantigen_seq"]] = blat['Qname'].str.split('_', expand=True)
	blat = blat[["neoantigen_seq","neoantigen_blat_count_1block",'neoantigen_blat_count_2block']]
	neoantigens = pd.merge(neoantigens, blat, on="neoantigen_seq", how="left")
	neoantigens = neoantigens.fillna(0)
	neoantigens['neoantigen_blat_count_1block'] = neoantigens['neoantigen_blat_count_1block'].astype(int)
	neoantigens['neoantigen_blat_count_2block'] = neoantigens['neoantigen_blat_count_2block'].astype(int)
	
	## calculate peptide characteristics 
	## follow this paper: https://www.pnas.org/doi/epdf/10.1073/pnas.1500973112
	## Use data from Table S2: https://www.pnas.org/action/downloadSupplement?doi=10.1073%2Fpnas.1500973112&file=pnas.1500973112.sapp.pdf
	neoantigens['neoantigen_length'] = neoantigens['neoantigen_seq'].apply(lambda x: len(x))
	neoantigens[["neoantigen_MW", "neoantigen_hydrophobicity","neoantigen_polarity"]] = neoantigens["neoantigen_seq"].apply(get_peptide_characteristics)

	## save the final table
	neoantigens.to_csv("./teprof3_neoantigen_candidates_after_blast.tsv", sep="\t", index=False)

	## plot hydrophobicity and polarity distribution for good candidates
	neoantigens = neoantigens[neoantigens["neoantigen_blastp_result"]=="no_hit"]
	misc.run_bash_command("mkdir neoantigen_characteristics")
	useless = neoantigens.apply(make_neoantigen_figures, axis=1)

	with open("./teprof3_neoantigen_identification_statistic.tsv", "w+") as output_file:
		output_file.write("\t# of proteins\t# of peptides\n")
		output_file.write("neoantigens\t"+str(len(set(neoantigens["protein_id"])))+"\t"+str(len(set(neoantigens["neoantigen_seq"])))+"\n")

	misc.print_time("Finished generating list of neoantigens!")

def find_peptide_position(input_row, peptide_type):
	if peptide_type == "detected_peptide":
		peptide_seq = input_row['detected_peptide_seq']
	elif peptide_type == "neoantigen":
		peptide_seq = input_row['neoantigen_seq']
	protein_sequence = input_row["protein_sequence"]
	protein_sequence_TE = input_row["protein_sequence_TE"]
	protein_sequence_gene = input_row["protein_sequence_gene"]
	if peptide_seq in protein_sequence:
		if peptide_seq in protein_sequence_TE:
			return("TE")
		elif peptide_seq in protein_sequence_gene:
			return("gene")
		else:
			return("chimeric")
	else:
		return("wrong")

def define_peptide_evidence_level(input_row, peptide_type):
	if peptide_type == "detected_peptide":
		peptide_position = input_row['detected_peptide_position']
	elif peptide_type == "neoantigen":
		peptide_position = input_row['neoantigen_position']
	protein_type = input_row['protein_type']
	if(protein_type=="solo unannotated TE"):
		return("level 2b")
	elif(protein_type=="solo annotated TE"):
		return("level 3")
	elif(protein_type=="Chimeric uTE-normal/truncated" and peptide_position=="TE"):
		return("level 2b")
	elif(protein_type=="Chimeric uTE-normal/truncated" and peptide_position=="chimeric"):
		return("level 1a")
	elif(protein_type=="Chimeric uTE-normal/truncated" and peptide_position=="gene"):
		return("level 3")
	elif(protein_type=="Chimeric uTE-unannotated" and peptide_position=="TE"):
		return("level 2b")
	elif(protein_type=="Chimeric uTE-unannotated" and peptide_position=="chimeric"):
		return("level 1a")
	elif(protein_type=="Chimeric uTE-unannotated" and peptide_position=="gene"):
		return("level 2a")
	elif(protein_type=="Chimeric uTE-frame shift" and peptide_position=="TE"):
		return("level 2b")
	elif(protein_type=="Chimeric uTE-frame shift" and peptide_position=="chimeric"):
		return("level 1a")
	elif(protein_type=="Chimeric uTE-frame shift" and peptide_position=="gene"):
		return("level 1b")
	elif(protein_type=="Chimeric aTE-normal/truncated" and peptide_position=="TE"):
		return("level 3")
	elif(protein_type=="Chimeric aTE-normal/truncated" and peptide_position=="chimeric"):
		return("level 1a")
	elif(protein_type=="Chimeric aTE-normal/truncated" and peptide_position=="gene"):
		return("level 3")
	elif(protein_type=="Chimeric aTE-unannotated" and peptide_position=="TE"):
		return("level 3")
	elif(protein_type=="Chimeric aTE-unannotated" and peptide_position=="chimeric"):
		return("level 1a")
	elif(protein_type=="Chimeric aTE-unannotated" and peptide_position=="gene"):
		return("level 2a")
	elif(protein_type=="Chimeric aTE-frame shift" and peptide_position=="TE"):
		return("level 3")
	elif(protein_type=="Chimeric aTE-frame shift" and peptide_position=="chimeric"):
		return("level 1a")
	elif(protein_type=="Chimeric aTE-frame shift" and peptide_position=="gene"):
		return("level 1b")
	else:
		return("wrong")

def good_peptide_criteria(input_row):
	blastp = input_row["detected_peptide_blastp_result"]
	peptide_evidence_level = input_row["detected_peptide_evidence_level"]
	if blastp == "no_hit" and peptide_evidence_level in ["level 1a", "level 1b", "level 2a", "level 2b"]:
		return("keep")
	else:
		return("no")

def output_fasta_file_for_netmhcpan(input_row, output_file):
	output_file.write("> "+str(input_row['index'])+"_"+str(input_row["protein_id"])+"\n"+str(input_row["protein_sequence"])+"\n")

def output_fasta_file_for_blast_after_netmhcpan(input_row, output_file):
	output_file.write("> "+str(input_row['protein_id'])+"_neoantigen_"+str(input_row["neoantigen_seq"])+"\n"+str(input_row["neoantigen_seq"])+"\n")


def get_peptide_characteristics(peptide_seq):
	peptide_seq = peptide_seq.upper()
	amino_acid_mw_dict = {
		'A': 89.0935,   # Alanine
		'R': 174.2017,  # Arginine
		'N': 132.1184,  # Asparagine
		'D': 133.1032,  # Aspartic acid
		'C': 121.1590,  # Cysteine
		'E': 147.1299,  # Glutamic acid
		'Q': 146.1451,  # Glutamine
		'G': 75.0669,   # Glycine
		'H': 155.1552,  # Histidine
		'I': 131.1736,  # Isoleucine
		'L': 131.1736,  # Leucine
		'K': 146.1882,  # Lysine
		'M': 149.2124,  # Methionine
		'F': 165.1900,  # Phenylalanine
		'P': 115.1310,  # Proline
		'S': 105.0930,  # Serine
		'T': 119.1197,  # Threonine
		'W': 204.2262,  # Tryptophan
		'Y': 181.1894,  # Tyrosine
		'V': 117.1469}   # Valine
	hydrophobicity_dict = {"A": 1.8, # Alanine
		"C": 2.5, # Cysteine
		"D": -3.5, # Aspartic acid
		"E": -3.5, # Glutamic acid
		"F": 2.8, # Phenylalanine
		"G": -0.4, # Glycine
		"H": -3.2, # Histidine
		"I": 4.5, # Isoleucine
		"K": -3.9, # Lysine
		"L": 3.8, # Leucine
		"M": 1.9, # Methionine
		"N": -3.5, # Asparagine
		"P": -1.6, # Proline
		"Q": -3.5, # Glutamine
		"R": -4.5, # Arginine
		"S": -0.8, # Serine
		"T": -0.7, # Threonine
		"V": 4.2, # Valine
		"W": -0.9, # Tryptophan
		"Y": -1.3} # Tyrosine
	polarity_dict = {"A":8, # Alanine
		"C":5.5, # Cysteine
		"D":13, # Aspartic acid
		"E":12.3, # Glutamic acid
		"F":5.2, # Phenylalanine
		"G":9, # Glycine
		"H":10.4, # Histidine
		"I":5.2, # Isoleucine
		"K":11.3, # Lysine
		"L":4.9, # Leucine
		"M":5.7, # Methionine
		"N":11.6, # Asparagine
		"P":8, # Proline
		"Q":10.5, # Glutamine
		"R":10.5, # Arginine
		"S":9.2, # Serine
		"T":8.6, # Threonine
		"V":5.9, # Valine
		"W":5.4, # Tryptophan
		"Y":6.2} # Tyrosine
	amino_acid_mw = 0
	hydrophobicity = []
	polarity = []
	for aa in peptide_seq:
		amino_acid_mw += amino_acid_mw_dict[aa]
		hydrophobicity.append(hydrophobicity_dict[aa])
		polarity.append(polarity_dict[aa])
	return(pd.Series([amino_acid_mw, hydrophobicity, polarity], index=["neoantigen_MW", "neoantigen_hydrophobicity","polarity"]))
	#return amino_acid_mw, hydrophobicity, polarity

def make_neoantigen_figures(input_row):
	neoantigen_seq = str(input_row["neoantigen_seq"])
	hydrophobicity = pd.Series(input_row["neoantigen_hydrophobicity"])
	polarity = pd.Series(input_row["neoantigen_polarity"])
	hydrophobicity.index = hydrophobicity.index + 1
	polarity.index = polarity.index + 1
	## plot hydrophobicity
	sns.set(font_scale=0.6)
	sns.set_style("ticks")
	plt.figure(figsize=(3, 2))
	fig = sns.lineplot(data=hydrophobicity, lw=1)
	fig = sns.scatterplot(data=hydrophobicity, legend=False, s=15)
	fig.set_xlabel("Residue position")
	fig.set_ylabel("Hydrophobicity")
	fig.set_title(neoantigen_seq+":"+str(input_row["neoantigen_affinity(nM)"])+"nM,"+str(input_row["neoantigen_bind_level"])+","+input_row["protein_type"]+","+input_row["neoantigen_position"])
	#plt.title(neoantigen_seq+":"+str(input_row["neoantigen_affinity(nM)"])+"nM,"+str(input_row["neoantigen_bind_level"])+","+input_row["protein_type"]+","+input_row["neoantigen_position"], loc='left', y=1)
	plt.xticks(range(1,len(neoantigen_seq)+1))
	plt.tick_params(axis='x', which='both', length=4)
	plt.tick_params(axis='y', which='both', length=4)
	plt.tight_layout()
	plt.savefig("./neoantigen_characteristics/"+neoantigen_seq+"_hydrophobicity.pdf", format='pdf', dpi=100)
	## plot polarity
	sns.set(font_scale=0.6)
	sns.set_style("ticks")
	plt.figure(figsize=(3, 2))
	fig = sns.lineplot(data=polarity, lw=1)
	fig = sns.scatterplot(data=polarity, legend=False, s=15)
	fig.set_xlabel("Residue position")
	fig.set_ylabel("Polarity")
	fig.set_title(neoantigen_seq+":"+str(input_row["neoantigen_affinity(nM)"])+"nM,"+str(input_row["neoantigen_bind_level"])+","+input_row["protein_type"]+","+input_row["neoantigen_position"])
	#plt.title(neoantigen_seq+":"+str(input_row["neoantigen_affinity(nM)"])+"nM,"+str(input_row["neoantigen_bind_level"])+","+input_row["protein_type"]+","+input_row["neoantigen_position"], loc='left', y=1)
	plt.axhline(y=8, linewidth = 0.6, linestyle ="--", color="black")
	plt.xticks(range(1,len(neoantigen_seq)+1))
	plt.tick_params(axis='x', which='both', length=4)
	plt.tick_params(axis='y', which='both', length=4)
	plt.tight_layout()
	plt.savefig("./neoantigen_characteristics/"+neoantigen_seq+"_polarity.pdf", format='pdf', dpi=100)












