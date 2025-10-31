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

def get_args():
	""" Fetches the arguments for the program """
	program_desc = """TEProf3 takes aligned files (bam files) and/or assembled data (gtf files)
	and provides you a list of TE-derived transcripts with expression across samples."""
	parser = argparse.ArgumentParser(description=program_desc)

	if len(sys.argv)==1:
	    parser.print_help(sys.stderr)
	    sys.exit(1)

	## for debugging
	parser.add_argument("--test", dest='test',
						help="testing mode", action='store_true')

	## bonus functions
	parser.add_argument("--teprof2", dest='teprof2',
						help="directly use teprof2 output for translation", type=str, default="no")#default="/scratch/yliang/HNSCC/analysis/HNSCC_TSTEA_candidates_v1/assembled/annotatedcufftranscripts.tsv"
	parser.add_argument("--gdcjson", dest='gdcjson',
						help="json file downloaded from GDC for splice junction tsv files to generate sample_manifest.txt file", type=str, default="no")
	parser.add_argument("--gtexjson", dest='gtexjson',
						help="json file downloaded from GTEx for splice junction tsv files to generate sample_manifest.txt file", type=str, default="no")
	
	## prepare reference files
	parser.add_argument("--repeatmasker", dest='repeatmasker',
						help="repeatmasker file, only provide it if you are making a new reference", type=str, default="no")
	parser.add_argument("--geneannotation", dest='geneannotation',
						help="gene annotation file, only provide it if you are making a new reference. GENCODE gtf is recommended", type=str, default="no")
	parser.add_argument("--geneannotationfortranslation", dest='geneannotationfortranslation',
						help="gene annotation file for translation, only provide it if you are making a new reference. GENCODE gtf is recommended", action="store_true")
	parser.add_argument("--hervannotation", dest='hervannotation',
						help="herv annotation file", type=str, default="no")
	parser.add_argument("--geneprotein", dest='geneprotein',
						help="generate protein sequence fasta file from gene annotation", action="store_true")

	## general parameters
	parser.add_argument("-f", "--manifest", dest="manifest",
						help="Dataset manifest file: sample name, short/long read, file name, library strandness (fr or rf or none or leave_empty)(--assemblestrand flag will override this). (tab-delimited)", type=str)
	parser.add_argument("-ki", "--keepintermediate", dest="keepintermediate",
						help="Keep intermediate files (ie, output files from bedtools intersect)", action="store_true")
	parser.add_argument("-s", "--samplenumber", dest="samplenumber",
						help="number of samples to be processed together at the same time (default is 10)", type = int, default=10)
	parser.add_argument("-g", "--guided", dest="guided",
						help="run teprof3 in guided mode. In guided mode, teprof3 will skip process-assemble and TACO step, and pretend to have this provided gtf file as output from TACO. teprof3 will process and filter out non-TE transcripts from this provided gtf file, and use the filtered output for quantification", type=str, default="no")
	parser.add_argument("-rs", "--reset", dest="reset",
						help="Keep intermediate files (ie, output files from bedtools intersect)", action="store_true")
	parser.add_argument("-v", "--version", dest="version",
						help="Print version of teprof3", action="store_true")

	## assemble parameter
	parser.add_argument("-am", "--assemblemode", dest='assemblemode',
						help="how to run transcript de novo assemble (0 for no, 1 for short read, 2 for long read (stringtie), 3 for hybrid (stringtie)) (default is 0)", type = str, default="0")
	parser.add_argument("-at", "--assemblethread", dest='assemblethread',
						help="number of threads used for assemble (default is 4)", type = str, default="4")
	parser.add_argument("-al", "--assemblelength", dest='assemblelength',
						help="minimum transcript length for stringtie assembly (default is 200)", type = str, default="200") # http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
	parser.add_argument("-as", "--assemblesamplenumber", dest='assemblesamplenumber',
						help="number of samples to be processed together at the same time (default is 10)", type = int, default=10)
	parser.add_argument("-ast", "--assemblestrand", dest='assemblestrand',
						help="-1: strandness information is provided in the sample manifest file; 0: the library is not stranded; 1: the library is stranded and firststrand (eg. rf flag in stringtie); 2: the library is stranded and secondstrand (eg. fr flag in stringtie) (default is -1) (this flag overrides information provided in the sample manifest file)", type = str, default="-1") ## 2/fr for Wang lab stranded RNA-seq library
	parser.add_argument("-aj", "--assemblejunctionread", dest='assemblejunctionread',
						help="There should be at least this many spliced reads that align across a junction (i.e. junction coverage). This number can be fractional, since some reads align in more than one place. A read that aligns in n places will contribute 1/n to the junction coverage. Default: 1", type = str, default="1")

	## process assemble parameter
	parser.add_argument("-pt", "--processtpm", dest='processtpm',
						help="tpm cutoff to filter te-derived transcripts right after transcript assembly and before meta assembly (TACO)", type = float, default=0.5)
	parser.add_argument("-ps", "--processsamplenumber", dest='processsamplenumber',
						help="number of samples to be processed together at the same time (default is 10)", type = int, default=10)
	parser.add_argument("-ptn", "--processtranscriptnumber", dest='processtranscriptnumber',
						help="only include samples with >x TE-derived transcripts identified in the sample. Typically, there're >1000 TE-derived transcripts identified in primary tumor samples and cancer cell lines. You can fine-tune this cutoff by plotting the distribution of # of TE-derived transcripts per sample in your dataset.", type = int, default=100)
	#parser.add_argument("-pto", "--processtolerance", dest='processtolerance',
	#					help="when correct the edge of exons, how many bases can be tolerated", type = int, default=3)

	## transcript filtering parameter
	parser.add_argument("-fm", "--filtermode", dest='filtermode',
						help="how to filter TE-derived transcripts (1 for only using short read, 2 for using both short and long read (default is 1)", type = str, default="1")
	#parser.add_argument("-frp", "--filterrepeat", dest='filterrepeat',
	#					help="remove transcripts derived from simple repeats and only focus on transposable elements (default is True)", action='store_false')
	parser.add_argument("-fs", "--filtersamplenumber", dest='filtersamplenumber',
						help="number of samples to be processed together at the same time (default is 10)", type = int, default=10)
	parser.add_argument("-ft", "--filterthread", dest='filterthread',
						help="number of threads used for filtering transcripts (only used in the step where extract reads from bam files)", type = int, default=10)
	#parser.add_argument("-fb", "--filtertranscriptbufferlength", dest='filtertranscriptbufferlength',
	#					help="when collecting reads from bam files, collect reads +- transcriptbufferlength of transcript_start and transcript_stop", type = int, default=3000)
	parser.add_argument("-fa", "--filterannotated", dest='filterannotated',
						help="add this flag to keep annotated transcripts (annotated means exon_1 of TE-derived transcript overlaps with exon_1 of annotated transcript) (when it's annotated, it's usually very noisy, high false positive). Default is to only filter out annotated TE-derived transcripts", action="store_false")
	#parser.add_argument("-fch", "--filterchimeric", dest='filterchimeric',
	#					help="add this flag to only keep TE-derived transcripts that are chimeric (chimeric means exon_1 of TE-derived transcript doesn't overlap with exon of annotated transcript) (when it's not chimeric, it's usually very noisy, high false positive). Default is to keep everything and use chimeric-mate filter to keep \"correct\" non-chimeric transcripts", action="store_true")
	parser.add_argument("-fi", "--filterintronretention", dest='filterintronretention',
						help="number of exons of annotated transcript that one exon of TE-derived transcript overlaps with, this is a cutoff to remove transcripts with intron retention", type = int, default=3)
	parser.add_argument("-fmo", "--filtermonoexon", dest='filtermonoexon',
						help="add this flag to exclude mono-exonic transcripts. (Default is false and keep mono-exonic transcripts)", action='store_true')
	parser.add_argument("-fmot", "--filtermonoexontpm", dest='filtermonoexontpm',
						help="tpm cutoff for mono-exonic transcripts", type=float, default=1)
	parser.add_argument("-fdm", "--filterdownstreammate", dest='filterdownstreammate',
						help="number of reads capturing splicing junction events going to downstream exon. Default is 2", type=int, default=2)
	parser.add_argument("-fr", "--filterratio", dest='filterratio',
						help="ratio of mates that from upstream and splicing downstream. By default, (SJ_uniqlymapped_read_upstream+SJ_multimapped_read_upstream)<=perfect_SJ_uniqlymapped_read_downstream*0.5, 0.5 is this ratio", type=float, default=0.5)
	parser.add_argument("-fncf", "--filternochimericfilter", dest='filternochimericfilter',
						help="disable chimeric mate filter", action='store_true')
	parser.add_argument("-fljt", "--longreadsjtolerance", dest='longreadsjtolerance',
						help="number of bases to tolerate when compare splicing junction from stringtie and LR data", type=int, default=3)
	parser.add_argument("-flmt", "--longreadmtolerance", dest='longreadmtolerance',
						help="number of bases to tolerate when compare mono-exonic from stringtie and LR data", type=int, default=50)

	## mega assembly parameter
	parser.add_argument("-tt", "--tacothread", dest='tacothread',
						help="number of thread for taco", type = str, default="10")
	parser.add_argument("-tto", "--tacotolerance", dest='tacotolerance',
						help="when correct the edge of exons, how many bases can be tolerated", type = int, default=3)

	## quantification parameter
	parser.add_argument("-qm", "--quanmode", dest='quanmode',
						help="how to quantify TE-derived transcripts (1: only using short read, 2: only using SJ (default is 1)", type = str, default="1")
	parser.add_argument("-qs", "--quansamplenumber", dest='quansamplenumber',
						help="number of samples to be processed together at the same time for running stringtie quantification (default is 10)", type = int, default=10)
	parser.add_argument("-qsc", "--quansamplenumbercon", dest='quansamplenumbercon',
						help="number of samples to be processed together at the same time for concatenating stringtie quantification output (default is 50)", type = int, default=50)
	parser.add_argument("-qnp", "--quannoprepde", dest='quannoprepde',
						help="whether to run prepDE.py to generate raw count for deseq2 (default is True)", action = 'store_true')
	parser.add_argument("-qng", "--quannogene", dest='quannogene',
						help="whether to collect tpm information for canonical genes (default is True)", action = 'store_true')
	parser.add_argument("-ql", "--quanreadlength", dest='quanreadlength',
						help="read length of the provided RNA-seq libraries (default is 75)", type = str, default="75")

	args = parser.parse_args()
	return args

def print_version(version):
	""" Save version of each tool used in the pipeline """
	subprocess.run("echo \"TEProf3 version:\" > debug/teprof3_version_info.txt", shell=True, executable='/bin/bash')
	subprocess.run("echo \""+version+"\" >> debug/teprof3_version_info.txt", shell=True, executable='/bin/bash')
	subprocess.run("echo \"Stringtie\" >> debug/teprof3_version_info.txt", shell=True, executable='/bin/bash')
	subprocess.run("stringtie --version >> debug/teprof3_version_info.txt", shell=True, executable='/bin/bash')

def parse_manifest(input_file):
	""" parse the input manifest """
	if os.path.isfile(input_file) is False:
		psys.exit("Can't find " + input_file + ". Please check.")
	input_dict = {}
	with open(input_file, 'r') as file:
		for line in file:
			entry = line.strip('\n').split('\t')
			if len(entry) == 3 or len(entry) == 4:
				if entry[0] not in input_dict:
					## parse strandness information
					if len(entry) == 3:
						input_dict[entry[0]] = {"strandness":""}
					elif len(entry) == 4:
						if entry[3] == "none":
							input_dict[entry[0]] = {"strandness":""}
						elif entry[3] == "fr" or entry[3] == "rf":
							input_dict[entry[0]] = {"strandness":"--"+entry[3]}
						else:
							print_time("Please put one of the following options in the fourth column (strandness) of the sample manifest: fr or rf or none or leave_empty")
							exit()
				if entry[1] == "long":
					input_dict[entry[0]]["bam_LR"] = entry[2].strip()
				elif entry[1] == "short":
					input_dict[entry[0]]["bam_SR"] = entry[2].strip()
				elif entry[1] == "gtf":
					input_dict[entry[0]]["gtf"] = entry[2].strip()
				elif entry[1] == "SJ":
					input_dict[entry[0]]["SJ"] = entry[2].strip()
				else:
					print_time("Please put one of the following options in the second column (data_type) of the sample manifest: long/short/gtf/SJ")
					exit()
			elif entry == ['']:
				continue
			else:
				print_time("please put 3 or 4 arguments in each row in the manifest file")
				exit()
	return(input_dict)

def print_dataset(input_dataset):
	""" print input dataset for confirmation and record """
	for sample in input_dataset:
		for file in input_dataset[sample]:
			print('\t'.join([sample, file, input_dataset[sample][file]]))

def move_file(input_dataset, flags):
	""" Move bam files to bam folder for a cleaner folder space """
	print_time("Move files to folders")
	for sample in input_dataset:
		try:
			subprocess.run("mv "+input_dataset[sample]['bam_LR']+".bai bam/"+input_dataset[sample]['bam_LR']+".bai", shell=True, executable='/bin/bash')
		except:
			None
			#print_time("This is not an error, just a reminder: no index file for long read bam files is provided for "+sample+". Please make sure to make softlinks to .bai files and keep the suffix as .bai")
		try:
			subprocess.run("mv "+input_dataset[sample]['bam_LR']+" bam/"+input_dataset[sample]['bam_LR'], shell=True, executable='/bin/bash')
			input_dataset[sample]['bam_LR'] = "bam/"+input_dataset[sample]['bam_LR']
		except:
			None
			#print_time("This is not an error, just a reminder: no bam file for long read data is provided for "+sample+". Please make sure your input bam files are ended with .bam")
		try:
			subprocess.run("mv "+input_dataset[sample]['bam_SR']+".bai bam/"+input_dataset[sample]['bam_SR']+".bai", shell=True, executable='/bin/bash')
		except:
			if flags.quanmode == "1":
				print_time("no index file for short read bam files is provided for "+sample+". Please make sure to make softlinks to .bai files and keep the suffix as .bai")
				return(1)
			elif flags.quanmode == "2":
				None
		try:
			subprocess.run("mv "+input_dataset[sample]['bam_SR']+" bam/"+input_dataset[sample]['bam_SR'], shell=True, executable='/bin/bash')
			input_dataset[sample]['bam_SR'] = "bam/"+input_dataset[sample]['bam_SR']
		except:
			if flags.quanmode == "1":
				print_time("no bam file for short read data is provided for "+sample+". Please make sure your input bam files are ended with .bam")
				return(1)
			if flags.quanmode == "2":
				None
		try:
			subprocess.run("mv "+input_dataset[sample]['SJ']+" bam/"+input_dataset[sample]['SJ'], shell=True, executable='/bin/bash')
			input_dataset[sample]['SJ'] = "bam/"+input_dataset[sample]['SJ']
		except:
			print_time("no SJ file is provided for "+sample+".")
			return(1)
	return(input_dataset)


def print_time(input_string):
	""" Print message with time stamp """
	program_start_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
	print("[ "+program_start_time+" ] "+input_string )

def run_bash_command(command):
	subprocess.run(command, shell=True, executable='/bin/bash')

def gtf_to_refbed(input_gtf_file):
	""" convert stringtie gtf file to refbed file for washu genome browser visualization """
	# chr, transcript_start, transcript_stop, translation_start, translation_stop, 
	# strand, gene_name, transcript_id, type, exon(including UTR bases) starts, 
	# exon(including UTR bases) stops, and additional gene info (optional)
	output_refbed_file = input_gtf_file.replace("gtf", "refbed")
	input_data = {}
	with open(input_gtf_file, 'r') as input_file:
		for line in input_file:
			if line[0]!="#":
				entry = line.strip('\n').split('\t')
				detail_info = entry[8]
				transcript_id = detail_info.split("transcript_id \"")[1].split("\";")[0]
				try:
					gene_name = detail_info.split("gene_name \"")[1].split("\";")[0]
				except:
					gene_name = detail_info.split("transcript_id \"")[1].split("\";")[0]
				transcript_type = "coding"
				if entry[2] == "transcript":
					input_data[transcript_id] = [[entry[0], str(int(entry[3])-1), entry[4], str(int(entry[3])-1), entry[4], entry[6], gene_name, transcript_id, transcript_type],[],[],[detail_info]]
				elif entry[2] == "exon":
					input_data[transcript_id][1].append(int(str(int(entry[3])-1)))
					input_data[transcript_id][2].append(int(entry[4]))
	with open(output_refbed_file, 'w') as output_file:
		for transcript_id in input_data:
			output_file.write('\t'.join(input_data[transcript_id][0] + [','.join([str(x) for x in sorted(input_data[transcript_id][1])])] + [','.join([str(x) for x in sorted(input_data[transcript_id][2])])] + input_data[transcript_id][3])+'\n')
	run_bash_command("sort -k1,1 -k2,2n "+ str(output_refbed_file) +" | bgzip > "+ str(output_refbed_file) +".sorted.gz")
	run_bash_command("tabix -p bed "+ str(output_refbed_file) +".sorted.gz")
	run_bash_command("rm "+str(output_refbed_file))


def all_gtf_to_refbed(input_dataset, flags):
	gtf_files = []
	if flags.guided=="no":
		for sample in input_dataset:
			gtf_files.append(input_dataset[sample]['gtf'])
		gtf_files += glob.glob("assembled/*.TE.gtf")
	elif flags.guided !="no":
		gtf_files.append("assembled/TE_transcript_consensus.gtf")
	gtf_files.append("assembled/TE_transcript_consensus.filtered.gtf")
	with mp.Pool(20) as pool:
		uesless = pool.map(gtf_to_refbed, gtf_files)


def process_SJ_gdc_json_to_prep_sample_manifest(input_json_file):
	with open("sample_manifest.txt","w") as output_file:
		input_data = json.load(open(input_json_file,"r")) # input_data is a list, with information dictionary for each SJ file
		for sample in tqdm(input_data, desc="Extracting information:"): # entity_id means aliquot_id
			if sample["data_format"] == "TSV":
				output_file.write('\t'.join([sample["associated_entities"][0]["case_id"]+"_"+sample["associated_entities"][0]["entity_id"], "short", sample["analysis"]["input_files"][0]["file_name"]])+"\n")
				output_file.write('\t'.join([sample["associated_entities"][0]["case_id"]+"_"+sample["associated_entities"][0]["entity_id"], "gtf", sample["analysis"]["input_files"][0]["file_name"].replace("bam","gtf")])+"\n")
				output_file.write('\t'.join([sample["associated_entities"][0]["case_id"]+"_"+sample["associated_entities"][0]["entity_id"], "SJ", sample["file_name"].replace(".gz","")])+"\n")

def process_SJ_gtex_json_to_prep_sample_manifest(input_json_file):
	with open("sample_manifest.txt","w") as output_file:
		input_data = json.load(open(input_json_file,"r")) # input_data is a list, with information dictionary for each SJ file
		for sample in tqdm(input_data, desc="Extracting information:"):
			output_file.write('\t'.join([sample["object_id"].split("/")[1], "SJ", sample["file_name"]])+"\n")


def reset_folder(flags):
	input_file = flags.manifest
	if os.path.isfile(input_file) is False:
		print_time("Can't find " + input_file + ". Please check.")
		exit()
	files_in_assembled = []
	files_in_bam = []
	with open(input_file, 'r') as file:
		for line in file:
			entry = line.strip('\n').split('\t')
			if len(entry) == 3 or len(entry) == 4:
				if entry[1] == "long":
					files_in_bam.append(entry[2])
					files_in_bam.append(entry[2]+".bai")
				elif entry[1] == "short":
					files_in_bam.append(entry[2])
					files_in_bam.append(entry[2]+".bai")
				elif entry[1] == "gtf":
					files_in_assembled.append(entry[2])
				elif entry[1] == "SJ":
					files_in_bam.append(entry[2])
			elif entry == ['']:
				continue
			else:
				print_time("please provide a manifest file following the instructions on github.")
				exit()

	print_time("reset files from a failed run")
	for file in files_in_assembled:
		run_bash_command("mv ./assembled/"+file+" ./"+file)
	for file in files_in_bam:
		run_bash_command("mv ./bam/"+file+" ./"+file)

	run_bash_command("rm -rf bam assembled debug")

	return(None)


def cleanup_folder():
	run_bash_command("mkdir assembled/intermediate_files;mv assembled/transcript_count_matrix.csv assembled/intermediate_files/transcript_count_matrix.csv;mv assembled/taco_command.txt assembled/intermediate_files/taco_command.txt;mv assembled/sample_list_for_prepDE.txt assembled/intermediate_files/sample_list_for_prepDE.txt")
	run_bash_command("mv assembled/quantification_command.txt assembled/intermediate_files/quantification_command.txt;mv assembled/TE_transcript_consensus.gtf assembled/intermediate_files/TE_transcript_consensus.gtf;mv assembled/TACO_output assembled/intermediate_files/TACO_output")
	run_bash_command("mv assembled/gtf_to_merge.txt assembled/intermediate_files/gtf_to_merge.txt;mv assembled/gene_count_matrix.csv assembled/intermediate_files/gene_count_matrix.csv; mv assembled/*quantification.gtf assembled/intermediate_files/; mv assembled/*stats assembled/intermediate_files/; mv assembled/*TE.gtf assembled/intermediate_files/")
	run_bash_command("mv assembled/TE_transcript_consensus.filtered.gtf assembled/teprof3_output_TE_transcript_consensus.gtf; mv assembled/TE_transcript_consensus.filtered.refbed.sorted.gz assembled/teprof3_output_TE_transcript_consensus.refbed.sorted.gz; mv assembled/TE_transcript_consensus.filtered.refbed.sorted.gz.tbi assembled/teprof3_output_TE_transcript_consensus.refbed.sorted.gz.tbi")
	run_bash_command("gzip assembled/teprof3_output_quantification.tsv")
	run_bash_command("gzip assembled/teprof3_output_quantification.TE.tsv")


def summarize_stas(input_dataset):
	## clean up, summarize the statistic files
	gtf_to_sample_dict = {}
	input_gtf_files=[]
	for sample in input_dataset:
		gtf_to_sample_dict[input_dataset[sample]["gtf"]] = sample
		input_gtf_files.append(input_dataset[sample]['gtf'])
	stats_dict = {}
	for input_gtf_file in input_gtf_files:
		stats_dict[gtf_to_sample_dict[input_gtf_file]] = []
		header = []
		with open(input_gtf_file+".stats.txt", "r") as stats_file:
			for line in stats_file:
				entry = line.strip("\n").split("\t")
				stats_dict[gtf_to_sample_dict[input_gtf_file]].append(entry[2])
				header.append(entry[1])
	stats_pandas = pd.DataFrame.from_dict(stats_dict, orient='index')
	stats_pandas.columns = header
	stats_pandas = stats_pandas.transpose()
	stats_pandas.to_csv("assembled/teprof3_output_transcript_statistic.tsv", sep="\t")
	run_bash_command("rm assembled/*.gtf.stats.txt")

















