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
import tabix
import sys
import pandas as pd
import tqdm
import glob
import numpy as np

def print_time(input_string):
	""" Print message with time stamp """
	program_start_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
	print("[ "+program_start_time+" ] "+input_string )

def intersect(input_dataset, flags):
	''' classify transcripts into different types (TE-gene, TE-non-coding, TE transcript) '''
	''' filter them using read information before running mega assembly with TACO '''
	input_gtf_files = []
	for sample in input_dataset:
		input_gtf_files.append(input_dataset[sample]['gtf'])

	## step 1. convert each gtf files into tsv for easier handling
	print_time("Extract exon and transcript information from gtf files")
	with mp.Pool(10) as pool:
		pool.map(convert_gtf_to_pandas, input_gtf_files)

	## step 2. intersect transcripts with repeatmasker and gene annotation (gencode) to classify transcripts
	print_time("Intersect transcripts with TE and gene annotation")
	input_gtf_exon_files = glob.glob("./assemble/*exon.txt")
	input_gtf_transcript_files = glob.glob("./assemble/*transcript.txt")

	with mp.Pool(10) as pool:
		pool.starmap(classfy_transcripts, zip(input_gtf_exon_files, input_gtf_transcript_files, repeat(flags)))


def convert_gtf_to_pandas(input_gtf_file):
	""" convert gtf file to pandas dataframe for faster and easier downstream processing """
	allowed_chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
	with open(input_gtf_file, "r") as gtf_file, open(input_gtf_file+".transcript.txt", "w") as transcript_file, open(input_gtf_file+".exon.txt", "w") as exon_file:
		transcript_file.write("\t".join(["chr", "start", "stop", "strand", "gene_id", "transcript_id", "cov", "FPKM", "TPM"])+'\n')
		#exon_file.write("\t".join(["chr", "start", "stop", "strand", "gene_id", "transcript_id", "exon_number", "cov"])+'\n')
		for line in gtf_file:
			if line[0] != "#":
				entry = line.strip("\n").split("\t")
				if entry[0] in allowed_chromosomes:
					if entry[2] == "transcript":
						detail_info = entry[8]
						gene_id = detail_info.split("gene_id \"")[1].split("\";")[0]
						transcript_id = detail_info.split("transcript_id \"")[1].split("\";")[0]
						cov = detail_info.split("cov \"")[1].split("\";")[0]
						FPKM = detail_info.split("FPKM \"")[1].split("\";")[0]
						TPM = detail_info.split("TPM \"")[1].split("\";")[0]
						transcript_file.write("\t".join([entry[0], entry[3], entry[4], entry[6], gene_id, transcript_id, cov, FPKM, TPM])+'\n')
					if entry[2] == "exon":
						detail_info = entry[8]
						gene_id = detail_info.split("gene_id \"")[1].split("\";")[0]
						transcript_id = detail_info.split("transcript_id \"")[1].split("\";")[0]
						exon_number = detail_info.split("exon_number \"")[1].split("\";")[0]
						cov = detail_info.split("cov \"")[1].split("\";")[0]
						exon_file.write("\t".join([entry[0], entry[3], entry[4], entry[6], gene_id, transcript_id, exon_number, cov])+'\n')


def classfy_transcripts(input_gtf_exon_file, input_gtf_transcript_files, flags):
	""" overlap each exon of each stringtie transcript with repeatmasker """
	""" distinguish TE-promoter, TE-exon and mono-exonic TE transcripts """


	exon_information = pd.read_table(input_gtf_exon_file)
	exon_intersect_with_TE_information = get_column_TE_intersect_parallel(exon_information, get_column_TE_intersect, 5)
	exon_intersect_with_TE_information_gene_information = get_column_gene_intersect_parallel(exon_intersect_with_TE_information, gene_annotation, get_column_gene_intersect, 5)
	del exon_information, exon_intersect_with_TE_information

	print(exon_intersect_with_TE_information_gene_information.head())

	#print("all_exon:",len(exon_information.index))
	#print("TE_exon:", len(exon_intersect_with_TE_information.index))


# https://towardsdatascience.com/make-your-own-super-pandas-using-multiproc-1c04f41944a1
def get_column_TE_intersect_parallel(df, func, n_cores=5):
    df_split = np.array_split(df, n_cores)
    pool = Pool(n_cores)
    df = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    return df

def get_column_TE_intersect(input_dataframe):
	input_dataframe['TE_intersect'] = input_dataframe.apply(intersect_with_TE, axis = 1)
	return(input_dataframe)

def intersect_with_TE(row):

	repeatmasker_file = os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/repeatmasker_sorted.txt.gz"
	repeatmasker = tabix.open(repeatmasker_file)

	results = list(repeatmasker.query(str(row['chr']), int(row['start']), int(row['start'])+1))
	if len(results) == 0:
		return("no")
	elif len(results) == 1:
		return(results[0])		
	elif len(results) > 1:
		## when there're two TE copies overlap with each other
		TE_copy_sizes = []
		for result in results:
			TE_copy_sizes.append(int(result[2])-int(result[1]))
		correct_result_index = results[TE_copy_sizes.index(min(TE_copy_sizes))]
		return(results[correct_result_index])


def get_column_gene_intersect_parallel(df, func, n_cores=5):
    df_split = np.array_split(df, n_cores)
    pool = Pool(n_cores)
    df = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    return df

def get_column_gene_intersect(input_dataframe):
	input_dataframe['gene_intersect'] = input_dataframe.apply(intersect_with_gene, axis = 1)
	return(input_dataframe)

def intersect_with_gene(row):

	gene_annotation_file = os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/gene_annotation_sorted.txt.gz"
	gene_annotation = tabix.open(gene_annotation_file)

	results = list(gene_annotation.query(str(row['chr']), int(row['start']), int(row['start'])+1))
	if len(results) == 0:
		return("no")
	elif len(results) >= 1:
		same_strand_results = [result for result in results if (result[5] == row['strand'])]
		if len(same_strand_results) == 0:
			return("no")
		if len(same_strand_results) == 1:
			return(same_strand_results[0])
		if len(same_strand_results) > 1:
			## coding > pseudo > nonCoding > problem > other
			## if this exon is associated with serveral coding genes, i don't know which one to pick from. Just randomly pick one.
			transcript_types = [result[9] for result in same_strand_results]
			if "coding" in transcript_types:
				for result in same_strand_results:
					if result[9] == "coding":
						return(result)
			elif "pseudo" in transcript_types:
				for result in same_strand_results:
					if result[9] == "pseudo":
						return(result)
			elif "problem" in transcript_types:
				for result in same_strand_results:
					if result[9] == "problem":
						return(result)
			else:
				return("other")






def filter_te_associated_transcripts_with_read_information():
	""" take TE-assocated transcripts, filter them with read information based on three different types of TE-associated transcrits """
	""" for TE-promoter transcripts: chimeric reads and chimeric mate support for the splice junction between exon 1 and exon 2; long read structure support """
	""" for TE-exon transcripts: chimeric reads and chimeric mate support for the splice junction between upstream and downstream exons; long read structure support """
	""" for mono-exonic TE transcripts: read support; long read structure support """






