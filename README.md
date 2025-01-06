<div align="left">
<img src="https://wangftp.wustl.edu/~yliang/TEProf3/teprof3-logo-zip-file/png/logo-no-background_new.png" width="300px"/>
</div>

# TEProf3
TEProf3 (TE-derived Promoter Finder 3) is a user-friendly, one-command tool designed to identify TE-derived promoters and transcripts using transcriptomic data from multiple sources, including short-read RNA-seq data, long-read RNA-seq data and single cell RNA-seq data.

## Table of contents

* [Chapter1 Required Softwares](#Chapter1-Required-softwares)
* [Chapter2 Required Python Packages](#Chapter2-Required-Python-packages)
* [Chapter3 Installation](#Chapter3-installation)
* [Chapter4 Get Started](#Chapter4-get-started)
* [Chapter5 Steps Expained](#Chapter5-Steps-Expained)
* [Chapter6 Output Expained](#Chapter6-Output-Expained)
* [Chapter7 Identify Neoantigen Candidates](#Chapter7-Identify-Neoantigen-Candidates)

## Chapter1-Required Softwares

* Add the following softwares to PATH

stringtie >= 2.2.1 (add both stringtie and prepDE.py3 to your PATH)

bedtools >= v2.17.0

samtools >= 1.2

TACO >= v0.7.3

* After installation of each tool, please run the following commands to make sure they are installed and in the PATH

which stringtie

which prepDE.py3

which bedtools

which samtools

which taco_run

## Chapter2-Required Python Packages

* the following packages can be installed by creating a conda environment with the provided teprof3.yml file

python >= 3.0.0

multiprocess >= 0.70.12.2

pandas >= 2.0.1

tqdm >= 4.65.0

genomepy >= 0.15.0

natsort >= 8.2.0

pysam >= 0.15.4

bgzip >= 1.2.1

tabix >= 1.2.1

## Chapter3-Installation

### 1. Download TEProf3 from github:

```
git clone https://github.com/Yonghao-Holden/TEProf3
```

### 2. Add TEProf3 to your PATH:

```
export PATH="<path_to_TEProf3>/TEProf3/bin/:\$PATH"
```

Or, add the above command to your `~/.bashrc` file so you don't need to export the path everytime you run TEProf3:

```
echo "PATH=\$PATH:<path_to_TEProf3>/TEProf3/bin/" >> ~/.bashrc
```


### 3. Generate reference files in the reference folder.

(1) Make a `reference` folder

```
cd TEProf3
mkdir reference
```

(2) Download repeatmasker (TE annotation) from https://www.repeatmasker.org. (For example: https://www.repeatmasker.org/genomes/hg38/RepeatMasker-rm405-db20140131/hg38.fa.out.gz)  

```
cd reference
wget https://www.repeatmasker.org/genomes/hg38/RepeatMasker-rm405-db20140131/hg38.fa.out.gz; gunzip hg38.fa.out.gz
teprof3 --repeatmasker hg38.fa.out
```

If you want to switch to another genome, remove the existed repeatmasker.txt file, download repeatmasker annotation for your interested genome and re-run `teprof3 --repeatmasker <repeatmasker_annotation_of_your_interested_genome>`.

(3) Download gene annotation from Gencode (https://www.gencodegenes.org/). (For example: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz)

TEProf3 is compatible with Gencode and Ensembl annoation. Please leave a requrest in the Issues tab if you need to use gene annotation from other source.

```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz; gunzip gencode.v36.annotation.gtf.gz
teprof3 --geneannotation gencode.v36.annotation.gtf
```

(4) Download genome reference reference (for TACO)

```
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz
gunzip hg38.fa.gz
ln -s hg38.fa ref_genome.fa
```

(5) Download HERV database from Telescope (only required if you are interested in whether the transcript is overlapped with HERV)

Telescope: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006453#sec019
```
# download S1 File. Annotation of HERV elements from 60 subfamilies in reference genome hg38
wget https://doi.org/10.1371/journal.pcbi.1006453.s006
ln -s journal.pcbi.1006453.s006 HERV.gtf
teprof3 --hervannotation HERV.gtf
```

(6) Download full length LINE1 database from L1Base 2 (only required if you are interested in whether the transcript is overlapped with LINE1)

L1Base 2: http://l1base.charite.de/l1base.php
```
# download S1 File. Annotation of HERV elements from 60 subfamilies in reference genome hg38
wget http://l1base.charite.de/BED/hsflnil1_8438_rm.bed
ln -s hsflnil1_8438_rm.bed LINE1.bed
```

(7) Set up genomepy reference (only required for translation)

```
## in python
genomepy.install_genome("hg38", annotation=True, provider="UCSC", genomes_dir="/bar/yliang/genomes/private/genomepy")

## in TEProf3/reference (name of the softlink needs to be genomepy)
ln -s /bar/yliang/genomes/private/genomepy genomepy
```
(8) Set up Blastp database (only required for classfying proteins) (not neccessary in most cases)

```
## https://www.ncbi.nlm.nih.gov/books/NBK131777/
## https://www.ncbi.nlm.nih.gov/books/NBK52640/
## https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/
## download standalone BLAST local program
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/ncbi-blast-2.13.0+-x64-linux.tar.gz
```

```
## to set up BLAST nr database
## https://danielbruzzese.wordpress.com/2018/03/26/installing-and-querying-a-local-ncbi-nucleotide-database-nt/
mkdir blast_database_nr; cd blast_database_nr
wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.??.tar.gz"
find . -maxdepth 1 -name "*.gz"  | while read file ; do echo "tar -xvf ${file}" >> unzip.txt; done ;
parallel_GNU -j 20 < unzip.txt
```

```
## to set up BLAST database of all ORF from gene annotation
mkdir blast_database_gene_annotation; cd blast_database_gene_annotation
teprof3 --geneprotein --translationlength 20 --translationgenome hg38 2> geneprotein.err &> geneprotein.log
makeblastdb -in gene_protein_sequence.fa \
 -parse_seqids \
 -title "gene annotation derived protein sequences" \
 -dbtype prot \
 -out gene_protein_sequence 2> makeblastdb.err &> makeblastdb.log
```

```
## to set up BLAST database of Swiss-Prot/Uniprot
mkdir blast_database_uniprot; cd blast_database_uniprot
ln -s uniprotkb_AND_reviewed_true_AND_model_o_2024_05_12.fasta uniprot.fa
makeblastdb -in uniprot.fa \
 -parse_seqids \
 -title "Uniprot database" \
 -dbtype prot \
 -out uniprot 2> makeblastdb.err &> makeblastdb.log
```

(9) At the end, you should the following folder structure:

```
--TEProf3
	--bin
	--README.md
	--reference
		--gene_annotation.gtf
		--gene_exon_annotation_sorted.txt
		--gene_exon_coordinates_annotation.txt
		--gene_exon_intron_annotation_sorted.txt
		--gene_start_codon_annotation_sorted.txt
		--gene_transcript_annotation_sorted.txt
		--repeatmasker_sorted.txt
		--ref_genome.fa
		--herv_annotation.txt
		--LINE1.bed
	--teprof3.yml
	--UpdateLog.txt
```

## Chapter4-Get Started

### 1. Create a folder where you want to run TEProf3 in.


### 2. Make softlinks to all the input files (.bam files, .bam.bai files, .SJ files from STAR, and optionally gtf files). Use absolute path when generating softlinks.

(1) Note for input RNA-sea data:

```
(a) Total RNA-seq data is supported and has comparable performance with mRNA-seq data. 
(b) Long read RNA-seq data is supported. TEProf3 will quantify number of reads from LR data that support each novel splice junction, but will NOT use long read information for filtering.
(c) TEProf3 is also compatible with scRNA-seq data. Take pseudobulk file of each cluster as one sample and run TEProf3.
```

(2) Note for bam files:
```
(a) For unstranded short read RNA-seq library, XS tag must be included in alignment step to be compatible with de novo assembly. TopHat and HISAT2 already include this tag in the output bam file, but for STAR, `--outSAMstrandField intronMotif` or `--outSAMattributes XS` must be included in the alignment command.
(b) For stranded short read RNA-seq library, XS tag is not neccessary as the strandness will be specified in the assemble command.
(c) For long read RNA-seq library, `-ax splice` option must be included if using `minimap2` for alignment
(d) For more information and explaination, please refer to manual of stringtie and STAR.
```

(3) there're a few QC metrices that're important and will greatly impact the specificity and sensitivity of TEProf3. 
```
(a) Gene body coverage: because TEProf3 heavily depends on finding reads that support the novel splice junctions at the 5' end of the transcript. Thus, a library with severe 3' bias (for example, a mRNA-seq library with severe RNA degradation during library preparation) is not optimal for TEProf3.
(b) Library complexity: because TE-derived transcripts, in most cases, have a relatively lower expression level. Thus, a library with low library complexity might underestimate the abundance of TE-derived transcripts.
(c) Sequencing depth: it is recommended to have at least 30-40M reads.
(d) Pair-end vs single-end sequencing: TEProf3 heavily depends on the information of reads spanning over splice junction. Single-end sequencing impairs the propability of detecting splice junction, thus underestimate the abundance of TE-derived transcripts.
(e) Read length: longer read length, for example 2x150bp, will have a higher chance having reads detecting splice junctions, thus, increase the sensitivity of the TE-derived transcript discovery.
```
### 3. Create a `sample_manifest.txt` file with information of input files. For example:

```
sample1	short	mRNA_sample1.bam
sample2	short	mRNA_sample2.bam
sample3	short	mRNA_sample3.bam
sample1	SJ	mRNA_sample1.SJ.out.tab
sample2	SJ	mRNA_sample2.SJ.out.tab
sample3	SJ	mRNA_sample3.SJ.out.tab
sample1	gtf	mRNA_sample1.gtf
sample2	gtf	mRNA_sample2.gtf
sample3	gtf	mRNA_sample3.gtf
```

Column 1 is sample name, column 2 is data type (short, long, SJ, and gtf), column 3 is file name. Index files for bam files are not needed in the `sample_manifest.txt` file.

Be careful don't put SPACE in the file name.

Here are some quick commands to automatically genereate the sample_manifest.txt file for you.

```
find . -maxdepth 1 -name "*.bam"  | while read file ; do xbase=$(basename $file); sample_name_temp=${xbase/mRNA_/}; sample_name=${sample_name_temp/.bam/} ; echo -e "${sample_name}\tshort\t${xbase}" >> sample_manifest.txt; done ;
find . -maxdepth 1 -name "*.tab"  | while read file ; do xbase=$(basename $file); sample_name_temp=${xbase/mRNA_/}; sample_name=${sample_name_temp/.SJ.out.tab/} ; echo -e "${sample_name}\tSJ\t${xbase}" >> sample_manifest.txt; done ;
find . -maxdepth 1 -name "*.gtf"  | while read file ; do xbase=$(basename $file); sample_name_temp=${xbase/mRNA_/}; sample_name=${sample_name_temp/.gtf/} ; echo -e "${sample_name}\tgtf\t${xbase}" >> sample_manifest.txt; done ;
```

If you are processing data downloaded from GDC, you can run `teprof3 --gdcjson <downloaded json file>` to generate the sample manifest file.

1. on GDC, add all the `<*.rna_seq.star_splice_junctions.tsv.gz>` files of intersted samples in your cart

2. go into your cart

3. click "Metadata" bottom to download the json file

4. use this downloaded json file to generate manifest file

5. use the following command to filter the sample_manifest.txt file to remove files that weren't downloaded successfully

### 4. Run `teprof3 -f sample_manifest.txt 2> teprof3.err &> teprof3.log`

### 5. If a TEProf3 run errors out or if you want to re-run TEProf3, run `teprof3 -f sample_manifest.txt --reset` to reset the folder

## Chapter5-Steps Explained

<div align="left">
<img src="https://wangftp.wustl.edu/~yliang/TEProf3/teprof3-logo-zip-file/workflow_figure.png" width="800px"/>
</div>


### 1. Data preparation

Short read RNA-seq data give us the highest confidence (both total RNA-seq and mRNA-seq, single-end and paired-end are supported by TEProf3). So, short read bam file is required for TEProf3, along with the index file (.bam.bai) and SJ file from STAR output (.tab). [i.e. short and SJ are required in the `sample_manifest.txt` file]

gtf files are optional and can be generated using TEProf3, which is recommended. [i.e. gtf is optional in the `sample_manifest.txt` file]

Long read RNA-seq data is also supported but given the low sensitivity of long read RNA-seq, information from long read will only be profiled but not used for any filtering. [i.e. long is optional in the `sample_manifest.txt` file]

All steps are designed in a multiprocessing manner with 10 samples processed in parallel by default. Each sample requires one core and ~3Gb memory.

| Flag | Description |
| --- | --- |
| `-f` or `--manifest` | Provide a manifest file contains the following information of each file: sample name, data type (**short**. **long**, **gtf**, **SJ**), file name.|
| `-ki` or `--keepintermediate` | Include this flag to keep intermediate files generated, default is **False**. |
| `-g` or `--guided` | Provide a gtf file and run TEProf3 in guided mode. In guided mode, TEProf3 will parse the provided gtf file, identify TE-derived transcripts, and use them for quantification. |
| `-s` or `--samplenumber` | Number of samples to be processed in parallel for all steps, default is **10**. |
| `-rs` or `--reset` | Include this flag to reset the folder to the original setting. |


### 2. Assemble

TEProf3 uses Stringtie for transcript *de novo* assembly. You can select different assemble mode to decide what data source to use for assembly, e.g. short read RNA-seq, long read RNA-seq or hybrid. 

| Flag | Description |
| --- | --- |
| `-am` or `--assemblemode` | How to run transcript *de novo* assembly: 0 for no, 1 for short read, 2 for long read, 3 for hybrid, default is **0**.|
| `-as` or `--assemblesamplenumber` | Number of samples to be processed together at the same time, default is **10** |
| `-at` or `--assemblethread` | Number of threads used for each sample when running Stringtie, default is **4**. |
| `-al` or `--assemblelength` | Minimum transcript length used for Stringtie, default is **200**. |
| `-ast` or `--assemblestrand` | 0: the library is not stranded; 1: the library is stranded and firststrand (eg. rf flag in stringtie); 2: the library is stranded and secondstrand (eg. fr flag in stringtie), default is **0**. |
| `-aj` or `--assemblejunctionread` | Minimal number of spliced reads that align across a junction (i.e. junction coverage). This number can be fractional, since some reads align in more than one place. A read that aligns in n places will contribute 1/n to the junction coverage, default is **2**. |


### 3. Process assemble

In this step, TEProf3 will collect stringtie output and obtain TE-derived transcripts for each transcripts. 

(1) It will first remove lowly expressed transcripts (default is 0.5 TPM, you can change it with `--processtpm` flag). 

(2) It will classify transcripts into 4 types, "TE transcript", "TE coding gene transcript", "TE noncoding gene transcript" and "TE no gene transcript" based on the target it splices into.

| Flag | Description |
| --- | --- |
| `-ps` or `--processsamplenumber` | Number of samples to be processed in parallel, default is **10**.|
| `-pt` or `--processtpm` | TPM cutoff to filter TE-derived transcripts in each sample before TACO, default is **0.5**. |
| `-ptn` or `--processtranscriptnumber` | Only include samples with >x TE-derived transcripts identified in the sample. Typically, there're >1000 TE-derived transcripts identified in primary tumor samples and cancer cell lines. Extreme low number of TE-derived transcripts usually indicates low library quality. Normal tissue samples and some cancer samples have lower number of TE-derived transcripts. You can fine-tune this cutoff by plotting the distribution of # of TE-derived transcripts per sample in your dataset, default is **100**. |

<div align="left">
<img src="https://wangftp.wustl.edu/~yliang/TEProf3/teprof3-logo-zip-file/different_class_figure_new.png" width="800px"/>
</div>

### 4. Filter transcripts

In this step, TEProf3 will filter transcripts of each sample using the following information and get a list of filtered TE-derived transcript for each sample for TACO mega assembly:

(1) Intron retention: when one exon of a identified TE-derived transcripts overlaps with three exons of an annotated transcript, it usually means this TE-derived transcript is product of intron retention and thus will be removed. You can change this cutoff with the `--filterintronretention` flag.

(2) Novel transcripts: if exon1 of a TE-derived transcript overlaps with exon1 of an annotated transcript, and the annotated TSS is not within the same TE copy the TE-derived transcript is in, it will be removed. 

(3) Chimeric mate filter: a. more than 2 reads contain the splice junction events towards downstream exon. (can change using `--filterdownstreammate` flag) b. reads contain splicing junction events from upstream exon has less than 50% of reads contains splice junction events to downstream exon. (can change using `-filterratio` flag)

(4) Long read filter: TEProf3 will show number of reads that support each novel splice junction. No filtering will be performed based on this information.

It will output a gtf file of TE-derived transcripts and a tsv file for each sample.

| Flag | Description |
| --- | --- |
| `-fs` or `--filtersamplenumber` | Number of samples to be processed in parallel, default is **10**.|
| `-fm` or `--filtermode` | How to filter TE-derived transcripts (0: for no filtering, 1: only using short read, 2: using both short and long read), default is **1**. Switching to filter mode 2 will add information from LR sequencing (for example, number of LR read that support each TE-derived transcript)|
| `-fa` or `--filterannotated` | To include TE-derived transcripts that are NOT novel (i.e. exon_1 of TE-derived transcript overlaps with exon_1 of annotated transcript and the annotated TSS is not within the same TE copy) (very noisy, high false positive, usually comes from 5' noise from stringtie), default is **False**. |
| `-fi` or `--filterintronretention` | Number of exons of annotated transcript that one exon of TE-derived transcript overlaps with, this is a cutoff to remove transcripts with intron retention, default is **3**.|
| `-fmo` or `--filtermonoexon` | Include this flag to exclude mono-exonic transcripts. Default is **False** and keep mono-exonic transcripts.|
| `-fmot` or `--filtermonoexontpm` | TPM cutoff for mono-exonic transcripts (it needs to be >= --processtpm), default is **1**.|
| `-fdm` or `--filterdownstreammate` | Number of reads capturing splicing junction events going to downstream exon, default is **2**.|
| `-fr` or `--filterratio` | Ratio of number of reads splicing from upstream and splicing into downstream, i.e. (SJ_uniqlymapped_read_upstream+SJ_multimapped_read_upstream)<=perfect_SJ_uniqlymapped_read_downstream*0.5, 0.5 is this ratio, default is **0.5**.|
| `-fljt` or `--longreadsjtolerance` | Number of bases to tolerate when comparing splice junction from stringtie and LR data, default is **3**.|
| `-flmt` or `--longreadmtolerance` | Number of bases to tolerate when compare mono-exonic from stringtie and LR data, default is **50**.|


<div align="left">
<img src="https://wangftp.wustl.edu/~yliang/TEProf3/teprof3-logo-zip-file/explain_filtering_figure_new.png" width="800px"/>
</div>

### 5. Mega assembly

In this step, TEProf3 will run TACO to merge transcript assembly information from each sample into one concensus transcript assembly annotation. And collect TE-derived transcripts from it. And apply filter (1) (2) from filter_transcripts again.

Possible flags:

| Flag | Description |
| --- | --- |
| `-tt` or `--tacothread` | Number of thread used for TACO, default is **10**.|

### 6. Quantification

In this step, TEProf3 will use Stringtie to quantify the expression of TE-derived transcripts using TACO output. 

It will output `teprof3_output_quantification_guided_taco_merged.tsv`

TEProf3 will stop at this step as this satisfy most applications.

| Flag | Description |
| --- | --- |
| `-qs` or `--quansamplenumber` | Number of samples to be processed in parallel for running Stringtie quantificaion, default is **10**.|
| `-qsc` or `--quansamplenumbercon` | Number of samples to be processed in parallel for concatenating Stringtie quantification output, default is **100**.|
| `-qm` or `--quantificationmode` | How to quantify TE-derived transcripts (1: use short read and SJ so you will have TPM and SJ support, 2: only use SJ, default is **1**.|
| `-qnp` or `--quannoprepde` | Whether to run prepDE.py to generate raw count for DESeq2, default is **False**. Add this flag to disable prepDE.py command.|
| `-ql` or `--quanreadlength` | Read length of the provided RNA-seq libraries, default is **75**.|


## Chapter6-Output Explained

### 1. `teprof3_output_transcript_statistic.tsv`

Table about number of transcripts that are identified and filtered at each step.

### 2, `teprof3_output_filter_transcript*`

Table with detailed information of each TE-derived transcript in each sample.

There will be duplicated transcript_id in this table as one trancript might be overlapped with multiple full length HERV or LINE1.

### 3. `teprof3_output_filter_transcript_TE_transcript_consensus.tsv`

Table with detailed information of the TE-derived transcripts consensus.

### 4. `teprof3_output_TE_transcript_consensus.gtf`

GTF file of the TE-derived transcripts consensus. You can use it for guided mode when you want to quantify the expression TE-derived transcripts in a different dataset.

### 5. `teprof3_output_quantification.TE.tsv.gz`

This table saves the expression information of each TE-derived transcript in each sample, including splice junction support of each transcript in each sample.

To determine whether a TE-derived transcript is present or not in one sample, we recommned to use the following three criteria:
(1) > 1TPM
(2) >= 1 perfect_SJ_uniqlymapped_read_downstream
(3) SJ_uniqlymapped_read_upstream <= 0.5 * perfect_SJ_uniqlymapped_read_downstream

### 6. `teprof3_output_quantification.tsv.gz`

This table saves the expression information of all transcripts in each sample.

### 7. `*refbed.sorted.gz`

You can directly upload these refbed files to WashU epigenome browser (https://epigenomegateway.wustl.edu/browser/) for visualization.

### 8. `intermediate_files/*.TE.gtf`

These are gtf files of TE-derived transcripts that pass filtering and used as input for mega-assembly. 


## Chapter7-Note for scRNA-seq

TEProf3 on scRNA-seq. Cluster first, and then take each cluster as one sample, take the pseudo-bulk bam file from each cluster as input and run TEPorf3. Then, we should be able to get a consensus transcript assembly across all cells. Then we can use this transcript assembly back to each individual cell for transcript-level quantification and TE quantification. Another way to do this is to treat each cell as one sample, then running stringtie on each cell is too noisy. And running TACO on thousands of samples is messy too. The third way to do this is running stringtie on pseudo-bulk bam file of all cells. Then it will have less sensitivity. Running stringtie on each cluster gives us enough coverage for each cluster, and good sensitivity because those cluster-specific transcripts won't be diluted out. The way Wanqing did it is to take bulk RNA-seq data of the same cell line. But for tumors, we won't have the matched bulk RNA-seq data most of the time. When with high cell number, to quantify the expression of each TE-derived transcript in each cell, use a chimeric gene annotation (e.g. GENCODE + TEProf3 output) to make the genome index and re-run STAR-solo to get the transcript-cell matrix.



## Chapter8 Identify Neoantigen Candidates

This chapter explains how to use teprof3 to identify neoantigen candidates derived from TE-derived transcripts.

Step 8, 9 and 10 are annotating database and characterize protein products translated from TE-derived transcripts. It will output a fasta file that can be used as a database for MassSpec search.

Step 11-15 are to identify neoantigens from TE-derived proteins detected in a MassSpec dataset.

### 7. prepare gene annotation

Go back to the reference folder of TEProf3, run this command to make a gene annotation of only non-TE-derived transcripts:

`teprof3 --geneannotationfortranslation`

### 8. ***in silico*** translation

In this step, TEProf3 will perform ***in silico*** translation on TE-derived transcripts.

This is a standalone step. Please run this after TEProf3 is done:

`teprof3 --translation teprof3_output_TE_transcript_consensus.gtf 2> translation.err &> translation.log`

Possible flags:

| Flag | Description |
| --- | --- |
| `-ti` or `--translation` | for translationmode 1: provide gtf file from teprof3 output for translation (e.g. teprof3 --translation teprof3_output_TE_transcript_consensus.gtf --translationmode 1)\nfor translationmode 2: provide table from teprof3 output for translation (e.g. teprof3 --translation teprof3_output_filter_transcript_TE_transcript_consensus.tsv --translationmode 2).|
| `-tm` or `--translationmode` | mode for translation (1: remove all annotated TE-derived isoforms from the gene annotation, and re-annotate identified TE-derived transcripts using the gene annotation with only non-TE-derived isoforms, this mode can help us differentiate the TE-derived isoform and non-TE-derived isoform of the same gene for translation; 2: this will use the provided gene annotation with no modification to annotate identified TE-derived transcripts for translation), default is **1**.|
| `-tl` or `--translationlength` | Minimum length of translated peptide, default is **20**.|
| `-tg` or `--translationgenome` | Reference genome name (eg. hg38, mm10) (name you used for installing genomepy), default is **hg38**.|

Outputs:

```
teprof3_output_protein_sequence_blastp.fa
teprof3_output_protein_sequence_blastpshort.fa
teprof3_output_protein_sequence.fa
teprof3_output_protein_sequence_gene_blastp.fa
teprof3_output_protein_sequence_gene_blastpshort.fa
teprof3_output_protein_sequence_TE_blastp.fa
teprof3_output_protein_sequence_TE_blastpshort.fa
teprof3_output_protein_information.tsv
teprof3_output_protein_information_light.tsv
command_used_for_teprof3_translation.txt
```

### 9. run blastp to prepare for protein classification

A few limitations:
blastp must be run on a SLURM system. blast must be installed. BLASTDB must be set up.

Run the following commands:

```
mkdir annotate_database; cd annotate_database
teprof3 --split protein_sequence_gene_blastp.fa
teprof3 --split protein_sequence_gene_blastpshort.fa
teprof3 --split protein_sequence_TE_blastp.fa
teprof3 --split protein_sequence_TE_blastpshort.fa

teprof3 --blastp protein_sequence_gene_blastp.fa
teprof3 --blastp protein_sequence_TE_blastp.fa
teprof3 --blastpshort protein_sequence_gene_blastpshort.fa
teprof3 --blastpshort protein_sequence_TE_blastpshort.fa

teprof3 --blastp protein_sequence_gene_blastp.fa --blastpdatabase gene
teprof3 --blastpshort protein_sequence_gene_blastpshort.fa --blastpdatabase gene
```

After running the above blastp commands, you will have a folder structure as followed:

```
--annotate_database
	*.fa
	*blastp_commands.txt
	*split_commands.txt
	--subsampled_fasta
```

### 10. classify protein products

In this step, TEProf3 will classify chimeric protein products into 

```
## run it in the annotate_database folder
teprof3 --classify teprof3_protein_information.tsv
#teprof3 --classify /scratch/yliang/HNSCC/data/CPTAC/2_test2_actual_run_TEProf2/test2_database/teprof3_protein_information.tsv
```

After running the above blastp commands, you will have a folder structure as followed:

```
--annotate_database
	*.fa
	*blastp_commands.txt
	*split_commands.txt
	teprof3_database_blast_to_the_same_protein.tsv
	teprof3_database_classification.tsv
	teprof3_percentage_of_database_protein_type_piechart.pdf
	--subsampled_fasta
```

### 11. run blastp on detected peptides

```
mkdir analysis_use_teprof3; cd analysis_use_teprof3
## run blastp
ln -s <path_to_a_fasta_file_with_detected_peptide> detected_peptide.fa ## ln -s ../detected_peptides.fa .
ln -s <teprof3_database_classification.tsv> teprof3_database_classification.tsv ## ln -s /scratch/yliang/HNSCC/analysis/CPTAC_TSTEA_candidates/2_test2_actual_run_TEProf2/annotate_database_use_teprof3/teprof3_database_classification.tsv
ln -s <teprof3_database_blast_to_the_same_protein.tsv> teprof3_database_blast_to_the_same_protein.tsv ## ln -s /scratch/yliang/HNSCC/analysis/CPTAC_TSTEA_candidates/2_test2_actual_run_TEProf2/annotate_database_use_teprof3/teprof3_database_blast_to_the_same_protein.tsv .
mkdir detected_peptides_blastp_result; cd detected_peptides_blastp_result
ln -s ../detected_peptides.fa .
teprof3 --split detected_peptides.fa
teprof3 --blastpshort detected_peptides.fa
## run blat
cd ..; mkdir detected_peptides_blat_result; cd detected_peptides_blat_result
ln -s <path_to_a_fasta_file_with_detected_peptide> detected_peptide.fa ## ln -s ../detected_peptides.fa .
blat -t=dnax -q=prot -minScore=0 -tileSize=5 -stepSize=1 -minIdentity=5 -repMatch=10000000 /bar/yliang/genomes/private/hg38_autosexchromosom.fa detected_peptides.fa detected_peptides.blat.psl
#pblat -threads=24 -t=dnax -q=prot -minScore=0 -tileSize=5 -stepSize=1 -minIdentity=5 -repMatch=10000000 /bar/yliang/genomes/private/hg38_autosexchromosom.fa detected_peptides.fa detected_peptides.blat.psl
```

right now, the folder structure will be:
```
--analysis_use_teprof3
	detected_peptides.fa
	teprof3_database_classification.tsv
	teprof3_database_blast_to_the_same_protein.tsv
	--detected_peptides_blastp_result
		detected_peptides.fa
		--subsampled_fasta
	--detected_peptides_blat_result
		detected_peptides.fa
		detected_peptides.blat.psl
```

### 12. obtain list of detected TE-derived proteins

```
## run it in the analysis_use_teprof3 folder
teprof3 --detectedprotein
```

right now, the folder structure will be:
```
--analysis_use_teprof3
	teprof3_database_blast_to_the_same_protein.tsv
	teprof3_detected_proteins.fa
	teprof3_detected_proteins.tsv
	teprof3_neoantigen_identification_statistic.tsv
	--detected_peptides_blastp_result
		detected_peptides.fa
		--subsampled_fasta
	--detected_peptides_blat_result
		detected_peptides.fa
		detected_peptides.blat.psl
		detected_peptides.blat.tsv
		teprof3_blat_result.tsv
```


### 13. run netMHCpan to get possible neoantigens

```
netMHCpan -a HLA-A02:01 -f teprof3_detected_proteins.fa -s -xls -xlsfile teprof3_detected_proteins.fa_HLA-A02-01.xls -BA > teprof3_detected_proteins.fa_HLA-A02-01.txt
netMHCpan -a HLA-A02:02 -f teprof3_detected_proteins.fa -s -xls -xlsfile teprof3_detected_proteins.fa_HLA-A02-01.xls -BA > teprof3_detected_proteins.fa_HLA-A02-02.txt
cat teprof3_detected_proteins.fa_HLA* > teprof3_detected_proteins.fa_HLA.txt
teprof3 --parsenetmhcpan <combined_netMHCpan_output> # teprof3 --parsenetmhcpan teprof3_detected_proteins.fa_HLA-A02-01.txt
```

right now, the folder structure will be:
```
--analysis_use_teprof3
	teprof3_detected_proteins.fa
	teprof3_detected_proteins.tsv
	teprof3_detected_proteins.fa_HLA-A02-01.txt
	teprof3_detected_proteins.fa_HLA-A02-01.xls
	teprof3_detected_proteins.fa_HLA-A02-01.txt.SB.txt
	teprof3_detected_proteins.fa_HLA-A02-01.txt.WB.txt
	teprof3_neoantigen_candidates_before_blast.fa
	teprof3_neoantigen_candidates_before_blast.tsv
	teprof3_neoantigen_identification_statistic.tsv
	--detected_peptides_blastp_result
		detected_peptides.fa
		--subsampled_fasta
	--detected_peptides_blat_result
		detected_peptides.fa
		detected_peptides.blat.psl
		detected_peptides.blat.tsv
		teprof3_blat_result.tsv
```

### 14. run blastp to remove neoantigens that could be from other proteins

```
## run blastp
mkdir neoantigen_blastp_result; cd neoantigen_blastp_result
ln -s ../teprof3_neoantigen_candidates_before_blast.fa
teprof3 --split teprof3_neoantigen_candidates_before_blast.fa
teprof3 --blastpshort teprof3_neoantigen_candidates_before_blast.fa
## run blat
cd ..; mkdir neoantigen_blat_result; cd neoantigen_blat_result
blat -t=dnax -q=prot -minScore=0 -tileSize=5 -stepSize=1 -minIdentity=5 -repMatch=10000000 /bar/yliang/genomes/private/hg38_autosexchromosom.fa ../teprof3_neoantigen_candidates_before_blast.fa teprof3_neoantigen_candidates_before_blast.blat.psl
```

right now, the folder structure will be:
```
--analysis_use_teprof3
	teprof3_detected_proteins.fa
	teprof3_detected_proteins.tsv
	teprof3_detected_proteins.fa_HLA-A02-01.txt
	teprof3_detected_proteins.fa_HLA-A02-01.xls
	teprof3_detected_proteins.fa_HLA-A02-01.txt.SB.txt
	teprof3_detected_proteins.fa_HLA-A02-01.txt.WB.txt
	teprof3_neoantigen_candidates_before_blast.fa
	teprof3_neoantigen_candidates_before_blast.tsv
	teprof3_neoantigen_identification_statistic.tsv
	--detected_peptides_blastp_result
		detected_peptides.fa
		--subsampled_fasta
	--detected_peptides_blat_result
		detected_peptides.fa
		detected_peptides.blat.psl
		detected_peptides.blat.tsv
		teprof3_blat_result.tsv
	--neoantigen_blastp_result
		teprof3_neoantigen_candidates_before_blast.fa
		--subsampled_fasta
	--neoantigen_blat_result
		teprof3_neoantigen_candidates_before_blast.blat.psl
```


### 15. get final list of neoantigen candidates

```
## run in the analysis_use_teprof3 folder
teprof3 --neoantigen
```

right now, the folder structure will be:
```
--analysis_use_teprof3
	teprof3_detected_proteins.fa
	teprof3_detected_proteins.tsv
	teprof3_detected_proteins.fa_HLA-A02-01.txt
	teprof3_detected_proteins.fa_HLA-A02-01.xls
	teprof3_detected_proteins.fa_HLA-A02-01.txt.SB.txt
	teprof3_detected_proteins.fa_HLA-A02-01.txt.WB.txt
	teprof3_neoantigen_candidates_before_blast.fa
	teprof3_neoantigen_candidates_before_blast.tsv
	teprof3_neoantigen_candidates_after_blast.tsv
	teprof3_neoantigen_identification_statistic.tsv
	--detected_peptides_blastp_result
		detected_peptides.fa
		--subsampled_fasta
	--detected_peptides_blat_result
		detected_peptides.fa
		detected_peptides.blat.psl
		detected_peptides.blat.tsv
		teprof3_blat_result.tsv
	--neoantigen_blastp_result
		teprof3_neoantigen_candidates_before_blast.fa
		--subsampled_fasta
	--neoantigen_blat_result
		teprof3_neoantigen_candidates_before_blast.blat.psl
	--neoantigen_characteristics
```





















