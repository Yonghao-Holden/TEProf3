#!/usr/bin/env python3
##########################################################################################################
# TEProf3: TE-derived Promoter Finder 3
# Author: Holden Liang
# main script
# This program takes bam files as input, perform de novo assembly, mega assembly and clean up to output a list of TE-derived transcripts
# teprof3 --test -f sample_manifest.txt
# teprof3 --repeatmasker hg38.fa.out
# teprof3 --geneannotation gencode.v36.annotation.gtf => takes about 10 minutes
# teprof3 -f sample_manifest.txt -am 1 -at 5
# teprof3 --intersect -f sample_manifest_test.txt => will intersect gtf files with TE annotation and output filtered gtf files with only transcripts that starts from TE
# teprof3 --split protein_sequence.fa
# teprof3 --blastp protein_sequence.fa
# teprof3 --gtf2refbed protein_sequence.fa
##########################################################################################################
## delete later:
## difference between ballgown and deseq2: https://support.bioconductor.org/p/107011/
## https://github.com/mortazavilab/TALON/blob/master/src/talon/talon.py used this as template
##########################################################################################################
## import external packages
import argparse
import time
import os
os.environ['NUMEXPR_MAX_THREADS'] = '4'
os.environ['NUMEXPR_NUM_THREADS'] = '2'
import numexpr as ne 
import sys
import subprocess

## import internal packages
import misc
import prepare_ref
import assemble
import process_assemble
import filter_transcripts
import process_meta_assemble
import quantification

def main():
	## step 1. mark program start time for reference
	misc.print_time("Start TEProf3 run")
	version = "v3.2.0"

	## step 2. obtain input arguments
	flags = misc.get_args()

	## reset folder
	if flags.reset==True:
		misc.reset_folder(flags)
		exit()

	if flags.version==True:
		misc.print_time("TEProf3 version is "+version)
		exit()

	if flags.test==True:
		input_dataset = quantification.quantification(input_dataset, flags)
		exit()

	## run teprof3 on gdc files
	if flags.gdcjson != "no":
		misc.process_SJ_gdc_json_to_prep_sample_manifest(flags.gdcjson)
		exit()

	## run teprof3 on gtex files
	if flags.gtexjson != "no":
		misc.process_SJ_gtex_json_to_prep_sample_manifest(flags.gtexjson)
		exit()



	## step 3. prepare reference files
	if flags.repeatmasker != "no":
		prepare_ref.prepare_repeat(flags.repeatmasker)
		exit()
	if flags.geneannotation != "no":
		prepare_ref.prepare_gene_annotation(flags.geneannotation)
		exit()
	if flags.geneannotationfortranslation:
		prepare_ref.prepare_gene_annotation_for_translation()
		exit()
	if flags.hervannotation != "no":
		prepare_ref.prepare_herv_annotation(flags.hervannotation)
		exit()
	if flags.geneprotein:
		prepare_ref.get_fasta_from_gene_annotation(flags)
		exit()
	
	## step 4. check for required files
	if os.path.isfile(os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/repeatmasker_sorted.txt") is False:
		misc.print_time("Please generate a repeatmasker.txt file in the reference folder following instructions on github")
		exit()

	if os.path.isfile(os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/gene_transcript_annotation_sorted.txt") is False:
		misc.print_time("Please generate files for gene annotation in the reference folder following instructions on github")
		exit()

	## step 5. make folders to save output
	if os.path.isdir("./bam") is False:
		os.mkdir("./bam")
	else:
		print("please remove existed folders from last run (teprof3 -f sample_manifest.txt --reset) and re-run teprof3")
		exit()
	if os.path.isdir("./assembled") is False:
		os.mkdir("./assembled")
	else:
		print("please remove existed folders from last run (teprof3 -f sample_manifest.txt --reset) and re-run teprof3")
		exit()
	if os.path.isdir("./debug") is False:
		os.mkdir("./debug")
	else:
		print("please remove existed folders from last run (teprof3 -f sample_manifest.txt --reset) and re-run teprof3")
		exit()

	## step 6. parse input data manifest
	misc.print_version(version=version)

	if flags.samplenumber != 10:
		if flags.assemblesamplenumber == 10:
			flags.assemblesamplenumber = flags.samplenumber
		if flags.processsamplenumber == 10:
			flags.processsamplenumber = flags.samplenumber
		if flags.filtersamplenumber == 10:
			flags.filtersamplenumber = flags.samplenumber
		if flags.quansamplenumber == 10:
			flags.quansamplenumber = flags.samplenumber

	misc.print_time("Read manifest file")
	manifest = flags.manifest # manifest="sample_manifest.txt"
	input_dataset = misc.parse_manifest(manifest)
	if input_dataset == 1:
		exit()
	misc.print_dataset(input_dataset)

	with open("debug/command_used_for_teprof3.txt", "w") as command_file:
		command_file.write(' '.join(sys.argv)+"\n")
	
	
	## step 7. move bam files to bam folder for a cleaner folder
	input_dataset = misc.move_file(input_dataset, flags)
	if input_dataset == 1:
		exit()
	misc.print_dataset(input_dataset)

	if flags.guided=="no":
		## step 8. run assemble	
		input_dataset = assemble.start_assemble(input_dataset, flags)
		misc.print_dataset(input_dataset)
		no_strand_flag = assemble.check_strand_info(input_dataset, flags)
		if no_strand_flag == 1:
			misc.print_time("There's no strandness information in the provided gtf files. Please include this info in the gtf files. (if there's no strandness, we don't know which end is tss)")
			exit()
	
		## step 9. start processing gtf files
		input_dataset = process_assemble.classify_transcripts(input_dataset, flags)

		## step 10. filter transcripts based on read information
		input_dataset = filter_transcripts.filter_transcripts_with_read_info(input_dataset, flags)
		misc.summarize_stas(input_dataset)

		## step 11. merge transcript assembly using TACO
		input_dataset = process_meta_assemble.run_taco(input_dataset, flags)
		
	## step 12. process the TACO output, get TE-derived transcripts, correct tss position and exon position.
	input_dataset = process_meta_assemble.process_taco(input_dataset, flags)
	
	## step 13. quantification across samples
	input_dataset = quantification.quantification(input_dataset, flags)

	## step 14. clean up folder, convert gtf files to refbed files, say it's done
	misc.print_time("Convert all gtf files to refbed files")
	useless = misc.all_gtf_to_refbed(input_dataset, flags)
	misc.print_time("Clean up files")
	useless = misc.cleanup_folder()
	misc.print_time("TEProf3 is done. Upload refbed files to WashU Epigenome Browser (http://epigenomegateway.wustl.edu/browser/) and enjoy.")



if __name__ == '__main__':
	main()


