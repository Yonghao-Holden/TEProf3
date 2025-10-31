# TEProf3: TE-derived Promoter Finder 3
# Author: Holden Liang
# de novo assembly script
# -----------------------------------------------------------------------------
# Functions related to processing the bam files and perform de novo assembly
import argparse
import time
import subprocess # https://geekflare.com/python-run-bash/
import os
os.environ['NUMEXPR_MAX_THREADS'] = '4'
os.environ['NUMEXPR_NUM_THREADS'] = '2'
import numexpr as ne 
import multiprocess as mp
from itertools import repeat

import misc

def move_gtf(input_dataset):
	""" move gtf files to assemble folder, if there's no gtf provided """
	for sample in input_dataset:
		try:
			subprocess.run("mv "+input_dataset[sample]['gtf']+" assembled/"+input_dataset[sample]['gtf'], shell=True, executable='/bin/bash')
			input_dataset[sample]['gtf'] = "assembled/"+input_dataset[sample]['gtf']
		except:
			misc.print_time("no gtf file is provided. (1) please change --assemblemode to run de novo assembly with teprof3 (2) please check if you put gtf files in your sample_manifest.txt file")
			exit()
	return(input_dataset)

def start_assemble(input_dataset, flags):
	""" run different assemble commands """
	assemble_mode = flags.assemblemode # 0

	misc.print_time("Run de novo transcript assembly")

	if assemble_mode == "0":
		misc.print_time("You chose NOT to run de novo assembly in TEProf3 and provided gtf files, will move all the gtf files to the assembled folder and proceed to next step")
		input_dataset = move_gtf(input_dataset)
	elif assemble_mode == "1":
		misc.print_time("You chose to run stringtie on your short read mRNA-seq data, will start now")      
		input_dataset = run_stringtie_assemble(input_dataset, flags, "short")
	elif assemble_mode == "2":
		misc.print_time("You chose to run stringtie on your long read mRNA-seq data, will start now")
		input_dataset = run_stringtie_assemble(input_dataset, flags, "long")
	elif assemble_mode == "3":
		misc.print_time("You chose to run stringtie on your short read and long read mRNA-seq data, will start now")
		input_dataset = run_stringtie_assemble(input_dataset, flags, "hybrid")
	else:
		misc.print_time("Please put one of the following options for the --assemblemode flag: 0, 1, 2, 3")
		exit()
	misc.print_time("All files are ready")
	return(input_dataset)

def run_stringtie_assemble(input_dataset, flags, data_type):
	assemble_thread = flags.assemblethread # 4
	assemble_length = flags.assemblelength # 200
	assemble_junctionread = flags.assemblejunctionread # 1
	assemble_samplenumber = flags.assemblesamplenumber # 10

	commands = []
	with open('./assembled/assemble_commands.txt', 'w') as f:
		for sample in input_dataset:
			## process strandness information
			if flags.assemblestrand == '-1':
				assemble_strand = input_dataset[sample]["strandness"]
			elif flags.assemblestrand == '0':
				assemble_strand = ''
			elif flags.assemblestrand == '1':
				assemble_strand = '--rf'
			elif flags.assemblestrand == '2':
				assemble_strand = '--fr'
			else:
				misc.print_time("Unrecognized strand option for flag --assemblestrand")
				exit()
			## make commands
			if data_type == "short":
				input_SR_bam = input_dataset[sample]["bam_SR"]
				command = ' '.join(["samtools view -q 255 -h", input_SR_bam, "| stringtie", "- -o", "assembled/"+input_SR_bam.split("/")[1].replace("bam","gtf"), "-j", assemble_junctionread, "-p", assemble_thread, "-m", assemble_length, assemble_strand])
				input_dataset[sample]["gtf"] = input_SR_bam.replace("bam/","assembled/").replace("bam","gtf")
			elif data_type == "long":
				input_LR_bam = input_dataset[sample]["bam_LR"]
				command = ' '.join(["stringtie -L -o assembled/"+input_LR_bam.split("/")[1].replace("bam","gtf"), "-j", assemble_junctionread, "-p", assemble_thread, "-m", assemble_length, assemble_strand, input_LR_bam])
				input_dataset[sample]["gtf"] = input_LR_bam.replace("bam/","assembled/").replace("bam","gtf")
			elif data_type == "hybrid":
				input_SR_bam = input_dataset[sample]["bam_SR"]
				input_LR_bam = input_dataset[sample]["bam_LR"]
				command = ' '.join(["stringtie --mix -o assembled/"+input_SR_bam.split("/")[1].replace("bam","gtf"), "-j", assemble_junctionread, "-p", assemble_thread, "-m", assemble_length, assemble_strand, input_SR_bam, input_LR_bam])
				input_dataset[sample]["gtf"] = input_SR_bam.replace("bam/","assembled/").replace("bam","gtf")
			f.write(command+"\n")
			commands.append(command)
	with mp.Pool(assemble_samplenumber) as pool:
		uesless = pool.map(misc.run_bash_command, commands)
	return(input_dataset)

def check_strand_info(input_dataset, flags):
	misc.print_time("check if there's strand information in the gtf file")
	input_gtf_files = []
	for sample in input_dataset:
		input_gtf_files.append(input_dataset[sample]['gtf'])

	no_strand_flags = []
	for input_gtf_file in input_gtf_files:
		with open(input_gtf_file, "r") as input_file:
			no_strand_flag = 1
			for line in input_file:
				if line[0] != "#":
					entry = line.strip("\n").split("\t")
					if "GL" not in entry[0] and "KI" not in entry[0] and "chrM" not in entry[0]:
						if entry[2] == "transcript":
							strand = entry[6]
							if strand == "+" or strand == "-":
								no_strand_flag=0
								break
			if no_strand_flag == 1:
				misc.print_time("This gtf file doesn't have strandness information: "+input_gtf_file)
			no_strand_flags.append(no_strand_flag)
	return(max(no_strand_flags))














