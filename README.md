<div align="left">
<img src="https://wangcluster.wustl.edu/~yliang/TEProf3/teprof3-logo-zip-file/png/logo-no-background_new.png" width="300px"/>
</div>

#

**TEProf3 (TE-derived Promoter Finder 3)** is a user-friendly, one-command tool for identifying transposable element **(TE)-derived promoters and their resulting transcripts**. It supports transcriptomic data from different modalities, including **short-read RNA-seq**, **long-read RNA-seq**, and **single-cell RNA-seq** datasets.

Thank you for using TEProf3! If you use this tool in your research, please cite the following publication:

* Qu, X., Liang, Y., McCornack, C., Xing, X., Schmidt, H., Tomlinson, C., ... & Wang, T. (2025). Charting the regulatory landscape of TP53 on transposable elements in cancer. Genome Research, 35(6), 1456-1471. [PMID: [40360186](https://pubmed.ncbi.nlm.nih.gov/40360186/)]

## Table of contents

* [Chapter 1 - Required Softwares](#Chapter-1---Required-softwares)
* [Chapter 2 - Required Python Packages](#Chapter-2---Required-Python-packages)
* [Chapter 3 - Installation](#Chapter-3---installation)
* [Chapter 4 - Get Started](#Chapter-4---get-started)
* [Chapter 5 - Steps Explained](#Chapter-5---Steps-Explained)
* [Chapter 6 - Output Explained](#Chapter-6---Output-Explained)
* [Chapter 7 - Note for scRNA-seq](#Chapter-7---Note-for-scRNA-seq)

## Chapter 1 - Required Softwares

* stringtie >= 2.2.1 (both stringtie and prepDE.py3 need to be in the PATH variable)

* bedtools >= v2.17.0

* samtools >= 1.2

* TACO >= v0.7.3

## Chapter 2 - Required Python Packages

* the following packages can be installed by creating a conda environment with the provided teprof3.yml file

* python >= 3.0.0

* multiprocess >= 0.70.12.2

* pandas >= 2.0.1

* tqdm >= 4.65.0

* genomepy >= 0.15.0

* natsort >= 8.2.0

* pysam >= 0.15.4

* bgzip >= 1.2.1

* tabix >= 1.2.1

## Chapter 3 - Installation

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

(2) Download repeatmasker (TE annotation) from the [Repeatmasker](https://www.repeatmasker.org) website. (For example: [hg38](https://www.repeatmasker.org/genomes/hg38/rmsk4.0.5_rb20140131/hg38.fa.out.gz))

```
cd reference
wget https://www.repeatmasker.org/genomes/hg38/RepeatMasker-rm405-db20140131/hg38.fa.out.gz; gunzip hg38.fa.out.gz
teprof3 --repeatmasker hg38.fa.out
```

If you want to switch to another genome, remove the existed repeatmasker.txt file, download the repeatmasker annotation of your interested genome and re-run `teprof3 --repeatmasker <repeatmasker_annotation_of_your_interested_genome>`.

(3) Download gene annotation from [GENCODE](https://www.gencodegenes.org/). (For example: [hg38_GENCODE_v36](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz))

TEProf3 is compatible with GENCODE and Ensembl annoation. Please leave a requrest in the Issues tab if you need to use gene annotation from other source.

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

(5) Download HERV database from [Telescope](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006453#sec019)) (only required if you are interested in whether the transcript is overlapped with HERV)

```
# download S1 File. Annotation of HERV elements from 60 subfamilies in reference genome hg38
wget https://doi.org/10.1371/journal.pcbi.1006453.s006
ln -s journal.pcbi.1006453.s006 HERV.gtf
teprof3 --hervannotation HERV.gtf
```

(6) Download full length LINE1 database from [L1Base 2](http://l1base.charite.de/l1base.php) (only required if you are interested in whether the transcript is overlapped with LINE1)

```
# download S1 File. Annotation of HERV elements from 60 subfamilies in reference genome hg38
wget http://l1base.charite.de/BED/hsflnil1_8438_rm.bed
ln -s hsflnil1_8438_rm.bed LINE1.bed
```

(7) At the end, you should the following folder structure:

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

## Chapter 4 - Get Started

### 1. Create a folder where you want to run TEProf3 in.


### 2. Make softlinks to all the input files (.bam files, .bam.bai files, .SJ files from STAR, and optionally gtf files). Use absolute path when generating softlinks.

(1) Note for input RNA-sea data:

```
(a) Total RNA-seq data is supported and has comparable performance with mRNA-seq data. 
(b) Long read RNA-seq data is supported. TEProf3 will quantify number of reads from LR data that support each novel splice junction, but will NOT use long read information for filtering.
(c) TEProf3 is also compatible with scRNA-seq data. Take pseudobulk file of each cluster as individual sample and run TEProf3.
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
(a) **Gene body coverage**: because TEProf3 heavily depends on finding reads that support the novel splice junctions at the 5' end of the transcript. Thus, a library with severe 3' bias (for example, a mRNA-seq library with severe RNA degradation during library preparation) is not optimal for TEProf3.
(b) **Library complexity**: because TE-derived transcripts, in most cases, have a relatively lower expression level. Thus, a library with low library complexity might underestimate the abundance of TE-derived transcripts.
(c) **Sequencing depth**: it is recommended to have at least 30-40M reads.
(d) **Pair-end vs single-end sequencing**: TEProf3 heavily depends on the information of reads spanning over splice junction. Single-end sequencing impairs the propability of detecting splice junction, thus underestimate the abundance of TE-derived transcripts.
(e) **Read length**: longer read length, for example 2x150bp, will have a higher chance having reads detecting splice junctions, thus, increase the sensitivity of the TE-derived transcript discovery.
```
### 3. Create a `sample_manifest.txt` file with information of input files.

* An example `sample_manifest.txt` file:

```
sample1	short	mRNA_sample1.bam	none
sample2	short	mRNA_sample2.bam	rf
sample3	short	mRNA_sample3.bam	fr
sample1	SJ	mRNA_sample1.SJ.out.tab	none
sample2	SJ	mRNA_sample2.SJ.out.tab	rf
sample3	SJ	mRNA_sample3.SJ.out.tab	fr
sample1	gtf	mRNA_sample1.gtf	none
sample2	gtf	mRNA_sample2.gtf	rf
sample3	gtf	mRNA_sample3.gtf	fr
```

| Column | Description |
| --- | --- |
| 1 | sample name |
| 2 | data type (choose from the following options: short, long, SJ, gtf) |
| 3 | file name |
| 4 | library strandness (choose from the following options: fr, rf, none)(You can leave this cell empty. This cell by default is none)(--assemblestrand flag will override this)(More explanation of fr and rf in [here](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)) |

Here are some quick commands to automatically genereate the sample_manifest.txt file for you.

```
find . -maxdepth 1 -name "*.bam"  | while read file ; do xbase=$(basename $file); sample_name_temp=${xbase/mRNA_/}; sample_name=${sample_name_temp/.bam/} ; echo -e "${sample_name}\tshort\t${xbase}" >> sample_manifest.txt; done ;
find . -maxdepth 1 -name "*.tab"  | while read file ; do xbase=$(basename $file); sample_name_temp=${xbase/mRNA_/}; sample_name=${sample_name_temp/.SJ.out.tab/} ; echo -e "${sample_name}\tSJ\t${xbase}" >> sample_manifest.txt; done ;
find . -maxdepth 1 -name "*.gtf"  | while read file ; do xbase=$(basename $file); sample_name_temp=${xbase/mRNA_/}; sample_name=${sample_name_temp/.gtf/} ; echo -e "${sample_name}\tgtf\t${xbase}" >> sample_manifest.txt; done ;
```

### 4. Run `teprof3 -f sample_manifest.txt 2> teprof3.err &> teprof3.log`

### 5. If a TEProf3 run errors out or if you want to re-run TEProf3, run `teprof3 -f sample_manifest.txt --reset` to reset the folder

## Chapter 5 - Steps Explained

<div align="left">
<img src="https://wangcluster.wustl.edu/~yliang/TEProf3/teprof3-logo-zip-file/teprof3_scheme_2_high_resolution.png" width="800px"/>
</div>


### 1. Data preparation

* TEProf3 supports both total RNA-seq and mRNA-seq, in single-end or paired-end formats.
* A short-read BAM file is required for TEProf3, along with its index file (.bam.bai) and the splice junction (SJ) file generated by STAR (.tab).

> `short` and `SJ` entries are required in the `sample_manifest.txt` file.

* GTF files are optional and can be automatically generated by TEProf3, which is recommended.

> `gtf` is optional in the `sample_manifest.txt` file.

* Long-read RNA-seq data are also supported.

> `long` is optional in the `sample_manifest.txt` file.

All processing steps are designed for multiprocessing, with 10 samples processed in parallel by default. Each sample typically requires one CPU core and approximately 3 GB of memory.

| Flag | Description |
| --- | --- |
| `-f` or `--manifest` | Provide a manifest file containing the following information: sample name, data type (**short**, **long**, **gtf**, **SJ**), file name.|
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
<img src="https://wangcluster.wustl.edu/~yliang/TEProf3/teprof3-logo-zip-file/different_class_figure_v2.png" width="800px"/>
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
<img src="https://wangcluster.wustl.edu/~yliang/TEProf3/teprof3-logo-zip-file/explain_filtering_figure_new.png" width="800px"/>
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


## Chapter 6 - Output Explained

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
> (1) > 1TPM
> (2) >= 1 perfect_SJ_uniqlymapped_read_downstream
> (3) SJ_uniqlymapped_read_upstream <= 0.5 * perfect_SJ_uniqlymapped_read_downstream

### 6. `teprof3_output_quantification.tsv.gz`

This table saves the expression information of all transcripts in each sample.

### 7. `*refbed.sorted.gz`

You can directly upload these refbed files to WashU epigenome browser (https://epigenomegateway.wustl.edu/browser/) for visualization.

### 8. `intermediate_files/*.TE.gtf`

These are gtf files of TE-derived transcripts that pass filtering and used as input for mega-assembly. 


## Chapter 7 - Note for scRNA-seq

TEProf3 on scRNA-seq. Cluster first, and then take each cluster as one sample, take the pseudo-bulk bam file from each cluster as input and run TEPorf3. Then, we should be able to get a consensus transcript assembly across all cells. Then we can use this transcript assembly back to each individual cell for transcript-level quantification and TE quantification. Another way to do this is to treat each cell as one sample, then running stringtie on each cell is too noisy. And running TACO on thousands of samples is messy too. The third way to do this is running stringtie on pseudo-bulk bam file of all cells. Then it will have less sensitivity. Running stringtie on each cluster gives us enough coverage for each cluster, and good sensitivity because those cluster-specific transcripts won't be diluted out. The way Wanqing did it is to take bulk RNA-seq data of the same cell line. But for tumors, we won't have the matched bulk RNA-seq data most of the time. When with high cell number, to quantify the expression of each TE-derived transcript in each cell, use a chimeric gene annotation (e.g. GENCODE + TEProf3 output) to make the genome index and re-run STAR-solo to get the transcript-cell matrix.
















