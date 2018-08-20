# Study gene expression using RNA-seq analysis

### Table of Contents
### 1. Introduction to RNA-seq (Experimental design)
### 2. Set up environment on High Performance Computing and the Genomic Compute Cluster 
### 3. Quality Check/ Alignment/ Estimate transcript abundances of RNA-seq reads
### 4. Expression and Differential Expression analysis
### 5. Visualization of the results
#### * Appendix. Reference-free expression estimation (Kallisto) 

___

### 1. Introduction to RNA-seq (Experimental design)

There are several sequencing analysis methods available based on the **goal of your experiment**:
*	_**RNA-seq**_ : _co-expression networks, differentially expressed genes_ 
*	_**Chip-Seq | ATAC-seq | DNase-seq**_ : _Gene regulation dynamics, chromatin modeling_ 
*	_**Whole-genome Seq**_ : _Rare variant analysis_ 

**RNA-seq** could answer the following experimental questions:
*	_Measure expression variation within or between species_
*	_Transcriptome characterization_ 
*	_Identify splicing sites_
*	_Discover novel transcript_
*	_Differential expression studies of a gene in different conditions, etc._

Before you start, you have to consider the **experiment design**:
*	_What resources do you have already? (reference genome, curated genes, etc)_
*	_Do you need biological replications? (usually yes)_
*	_Do you need technical replications? (mostly not)_
*	_Do you need controls?_
*	_Do you need deep sequencing coverage?_

**Types of RNA-seq reads**

| - | Notes |
| :------------ | :-------------|
| Single-end |* Fast run <br> * less expensive | 
| Paired-end |* More data for each fragment/alignment/assembly <br> * good for isoform-detection <br> * Good for detecting structural variations |
              

* Sanger sequencing - fasta format (1 header followed by any number of sequences lines) 
* NGS sequencing - fastq (Repeated 4 lines) 
    * Note that there are two fastQ files per sample in paired-end sequencing (+ strand, - strand) 
    * forward/reverse reads have almost same headers. 


___
### 2. Set up environment on High Performance Computing and the Genomic Compute Cluster 

If you would like to perform RNA-seq on Quest, you need to first do the followings:

- [ ]	1. Get an account on Quest
- [ ]	2. Request access to the Genomic Cluster (project directory: b1042) on Quest [[link]](https://kb.northwestern.edu/page.php?id=78602) 
- [ ]	3. Log in to Quest [[link]](http://www.it.northwestern.edu/research/user-services/quest/logon.html)
- [ ]	4. In this protocol, we will run an example analysis with chromosome X data of Homo sapiens. (Ref: Nature Protocol 2016) 
    - All necessary data you need are available in the following directory: _**/QuestDownloadPath/**_
      -	`'samples'` directory contains paired-end RNA-seq reads for 6 samples, 3 male and 3 female subjects from YRI (Yoruba from Ibadan, Nigeria) population. 
      -	`‘indexes’` directory contains the indexes for chromosome X for HISAT2. 
      -	`‘genome’` directory contains the sequence of human chromosome X (GrCH38 build 81)
      -	`‘genes’` directory contains human gene annotations for GrCH38 from RefSeq database. 
      -	`‘mergelist.txt’` and ‘geuvadis_phenodata.csv’ are exemplary scripts that you might want to write yourself in a text editor. 
      -	Since it is paired-end reads, each sample has two files: all sequence is in compressed 'fastq' format
        -	(cf) Our analysis only contains the genome of chromosome X, but if someone is interested in the full data sets, these files are ~25 times larger and you can find them: [ftp://ftp.ccb.jhu.edu/pub/RNAseq_protocol](ftp://ftp.ccb.jhu.edu/pub/RNAseq_protocol)
- [ ]	5. Create a working directory for the analysis and copy the data to your directory: 
```bash
	mkdir /YourWorkingDirectory/ 
	cd /YourWorkingDirectory/ 
	cp * /YourWorkingDirectory/ 
```

#### Software setup/installation:
- [ ] Load the necessary modules on Quest: fastqc, trimmomatic, samtools, HISAT2
```bash
	module load fastqc/0.11.5
	module load fastx_toolkit/0.0.14 ### trimmomatic 
	module load hisat2/2.0.4  
	module load samtools/1.6 
	module load stringtie/1.3.4 ### HTseq 
```

- [ ] Install Ballgown package on R 
```bash
	module load R
	R 
```
```R
	library("devtools") 
	source("http://www.bioconductor.org/biocLite.R")
	biocLite(c("alyssafrazee/RSkittleBrewer","ballgown", "genefilter","dplyr","devtools"))

	* Bioconductor version 3.0 or greater and R version 3.1 are required to run this protocol.
```

- [ ] Install StringTie, gffcompare 
    - [ ] To install StringTie, download the latest binary package from http://ccb.jhu.edu/software/stringtie, unpack the StringTie tarfile and cd to the unpacked directory.
    - [ ] To install gffcompare, download the latest binary package from http://github.com/gpertea/gffcompare, and follow the instructions provided in the README.md file.



___
### 3. Quality Check/ Alignment/ Estimate transcript abundances of RNA-seq reads

#### Step1. Analyze raw reads’ quality with FastQC
###### Input: bam/sam/fastQ file	:heavy_minus_sign:	Output: zip/html file contains quality report for each read
-	We have to confirm average quality per read, consistency, GC content (PCR bias), adapter/k-mer content, excessive duplicated reads, etc.
-	You can specify the list of adapter sequences with –adapters (-a) specifier. 
-	HTML ouput generated by FastQC helps visual inspection of overall read quality. Not all yellow and red highlights are problematic, so look through the reports with a grain of salt. 
-	*Other software options: __FastQC__ is for illumina read, __NGSQC__ is for basically every platform*

```bash
* You can find software instruction by typing a commandline: fastqc -h 
	
	mkdir qualitycheck
	fastqc --outdir ./qualitycheck/ *_chrX_*.fastq.gz
```

#### Step2. Filtering raw reads with Trimmomatic
###### Input: fastQ file before filtering	:heavy_minus_sign:	Output: fastQ file after filtering
-	Even if our data looks fine, it is always a good idea to filter out low/poor quality reads. 
-	Appropriate threshold should be determined by each experiment design or organism. 
-	Before running the filtering step, you need to clarify the sequencing method of your reads since the software commands are different. (single-ended or paired-ended?)
-	Our data are paired-ended, so we use ‘PE’ command for Trimmomatic.
-	Trimmomatric removes adapter sequences, low quality reads, too-short reads, etc.
-	Other software options: Trimmomatic, Fastx-toolkit (available on Quest)
-	FastX-toolkit [[Instruction link]](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html) does not accept gzip compressed files, so we would better make pipe and output in compressed format. The following command allows us to throw out any read that fails to meet a threshold of at least 70% of bases with Phred quality score > 20.
	- `gunzip -c ERR188044_chrX_1.fastq.gz | fastq_quality_filter -q 20 -p 70 -i -z -o ERR188044_chrX_1_filtered.fastq` 

```bash
* Paired ended:
	java -jar trimmomatic-0.35.jar PE -phred33 input_forward.fq.gz input_reverse.fq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

* Single ended:
	java -jar trimmomatic-0.35.jar SE -phred33 input.fq.gz output.fq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 

This will perform the following:
•	Remove adapters (ILLUMINACLIP:TruSeq3-PE.fa:2:30:10)
•	Remove leading low quality or N bases (below quality 3) (LEADING:3)
•	Remove trailing low quality or N bases (below quality 3) (TRAILING:3)
•	Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15 (SLIDINGWINDOW:4:15)
•	Drop reads below the 36 bases long (MINLEN:36)
•	Convert quality score to Phred-33 
```

#### Step3. Re-analyzing quality of the filtered reads with FastQC
###### Input: filtered fastQ file	:heavy_minus_sign:	Output: zip/html file contains quality report for each read

-	We have to confirm the read quality after filtering. 
-	You can determine this easily by re-running FastQC on the output fastq files.
	-	How many reads are left?
		-	Count the number of lines in a fastq file, which has 4 lines per entry:
		-	`$ wc -l output.fastq`
			- Divide this number by 4 to get the total number of reads in the fastq file.
	-	What % of raw reads passed the quality filter?
```bash
	fastqc --outdir ./qualitycheck/ *_chrX_*.fastq.gz
```

#### Step4. Alignment of RNA-seq reads to the genome 
###### Input: FastQ reads (2 per sample)	:heavy_minus_sign:	Output: SAM files (1 per sample)
-	HISAT2 (v 2.1.0) maps the reads for each sample to the ref genome: [Help page: `hisat2 –h` ]
-	Note that HISAT2 commands for paired(-1,-2) /unpaired (-U) reads are different. 
-	HISAT2 also provides an option, called --sra-acc, to directly work with NCBI Sequence Read Archive (SRA) data over the internet. This eliminates the need to manually download SRA reads and convert them into fasta/fastq format, without much affecting the run time. 
	-	(e.g.) `--sra-acc SRR353653,SRR353654`
-	Parameter for QC: proportion of mapped read on either genome/transcriptome 
	-	To confirm sequencing accuracy and contaminated DNA 
		- If RNA-seq reads are mapped to _human genome_ – 70~90% + a few multi-mapping reads 
		- If RNA-seq reads are mapped to _transcriptome_ – less mapping % + more multi-mapping reads by sharing same exon among isoforms 
	- If the result screen says that some reads aligned discordantly, it means some occurrences of infusion or translocation. Possibly mismatched/too-far paired-end reads. 
-	_Other software options: **Picard, STAR, PSeQC, Qualimap**_

```bash
	hisat2 -p 2 --dta -x ./indexes/chrX_tran -1 ./samples/ERR188044_chrX_1.fastq.gz -2 ./samples/ERR188044_chrX_2.fastq.gz -S ERR188044_chrX.sam
	hisat2 -p 2 --dta -x ./indexes/chrX_tran -1 ./samples/ERR188104_chrX_1.fastq.gz -2 ./samples/ERR188104_chrX_2.fastq.gz -S ERR188104_chrX.sam
	hisat2 -p 2 --dta -x ./indexes/chrX_tran -1 ./samples/ERR188234_chrX_1.fastq.gz -2 ./samples/ERR188234_chrX_2.fastq.gz -S ERR188234_chrX.sam
	hisat2 -p 2 --dta -x ./indexes/chrX_tran -1 ./samples/ERR188273_chrX_1.fastq.gz -2 ./samples/ERR188273_chrX_2.fastq.gz -S ERR188273_chrX.sam
	hisat2 -p 2 --dta -x ./indexes/chrX_tran -1 ./samples/ERR188454_chrX_1.fastq.gz -2 ./samples/ERR188454_chrX_2.fastq.gz -S ERR188454_chrX.sam
	hisat2 -p 2 --dta -x ./indexes/chrX_tran -1 ./samples/ERR204916_chrX_1.fastq.gz -2 ./samples/ERR204916_chrX_2.fastq.gz -S ERR204916_chrX.sam

* Confirm the QC parameters - % of alignment reads on genome? (e.g. last line in the result screen) 
```

#### Step 5. Sort and convert the SAM file to BAM 
###### Input: SAM file	:heavy_minus_sign:	Output: BAM file 
-	Samtools (v 2.1.0) sorts and converts the SAM file to BAM: [Help page: `samtools –help`]
-	Both SAM/BAM formats represent alignments. BAM is more compressed format. Unmapped reads may also be in the BAM file. Reads that map to multiple location will show up multiple times as well.
-	Exceeding mapping percentage over 100% does not indicate how many reads mapped. They can be inflated by filtering low q reads prior to alignment.

```bash
	samtools sort -@ 2 -o ERR188044_chrX.bam ERR188044_chrX.sam
	samtools sort -@ 2 -o ERR188104_chrX.bam ERR188104_chrX.sam
	samtools sort -@ 2 -o ERR188234_chrX.bam ERR188234_chrX.sam
	samtools sort -@ 2 -o ERR188273_chrX.bam ERR188273_chrX.sam
	samtools sort -@ 2 -o ERR188454_chrX.bam ERR188454_chrX.sam
	samtools sort -@ 2 -o ERR204916_chrX.bam ERR204916_chrX.sam
```

#### Step 6. Assemble and quantify expressed genes and transcripts with StringTie 

- [ ]	(a) Stringtie assembles transcripts for each sample:
###### Input: BAM file + reference GTF file	:heavy_minus_sign:	Output: Assembled GTF file (1 per sample)
```bash
	stringtie -p 2 -G ./genes/chrX.gtf -o ERR188044_chrX.gtf -l ERR188044 ERR188044_chrX.bam
	stringtie -p 2 -G ./genes/chrX.gtf -o ERR188104_chrX.gtf -l ERR188104 ERR188104_chrX.bam
	stringtie -p 2 -G ./genes/chrX.gtf -o ERR188234_chrX.gtf -l ERR188234 ERR188234_chrX.bam
	stringtie -p 2 -G ./genes/chrX.gtf -o ERR188273_chrX.gtf -l ERR188273 ERR188273_chrX.bam
	stringtie -p 2 -G ./genes/chrX.gtf -o ERR188454_chrX.gtf -l ERR188454 ERR188454_chrX.bam
	stringtie -p 2 -G ./genes/chrX.gtf -o ERR204916_chrX.gtf -l ERR204916 ERR204916_chrX.bam
```

- [ ]	(b) Stringtie merges transcripts from all samples:
###### Input: multiple GTF files to be merged + mergelist.txt with filenames	:heavy_minus_sign:	Output: One merged GTF (will be used as a reference for relative comparisons among samples) 
```bash
	stringtie --merge -p 2 -G ./genes/chrX.gtf -o stringtie_merged.gtf ./mergelist.txt
```

- [ ]	(c) Stringtie estimates transcript abundances and create table counts for Ballgown:
###### Input: BAM file of each sample + one merged GTF	:heavy_minus_sign:	Output: several output files for Ballgown-analysis ready 
```bash
	stringtie -e -B -p 2 -G stringtie_merged.gtf -o ./ballgown/ERR188044/ERR188044_chrX.gtf ERR188044_chrX.bam
	stringtie -e -B -p 2 -G stringtie_merged.gtf -o ./ballgown/ERR188104/ERR188104_chrX.gtf ERR188104_chrX.bam
	stringtie -e -B -p 2 -G stringtie_merged.gtf -o ./ballgown/ERR188234/ERR188234_chrX.gtf ERR188234_chrX.bam
	stringtie -e -B -p 2 -G stringtie_merged.gtf -o ./ballgown/ERR188273/ERR188273_chrX.gtf ERR188273_chrX.bam
	stringtie -e -B -p 2 -G stringtie_merged.gtf -o ./ballgown/ERR188454/ERR188454_chrX.gtf ERR188454_chrX.bam
	stringtie -e -B -p 2 -G stringtie_merged.gtf -o ./ballgown/ERR204916/ERR204916_chrX.gtf ERR204916_chrX.bam
```
-	_Other software available: **HTSeq-count, featureCounts**_

___
### 4. Expression and Differential Expression analysis


___
### 5. Visualization of the results
___

### Resources for Further Study 

- [RNA-seq wiki](https://github.com/griffithlab/rnaseq_tutorial/wiki)
- [RNA-seq analysis tutorial with Differential analysis in R DESeq2 package](https://github.com/CandiceChuDVM/RNA-Seq/wiki/RNA-Seq-analysis-tutorial)

##### Online courses
- [Bioconductor for Genomic Data Science (Coursera)](https://www.coursera.org/learn/bioconductor)
- [Command line tools for Genomic Data Science (Coursera)](https://www.coursera.org/learn/genomic-tools)
- [Genomic Data Analysis (edX)](https://www.edx.org/xseries/genomics-data-analysis)

##### Youtube lecture videos
- [Informatics for RNA-seq analysis (5 videos)](https://www.youtube.com/playlist?list=PL3izGL6oi0S849u7OZbX85WTyBxVdcpqx)
- [Introduction to RNA-sequencing (1hr 20mins)](https://youtu.be/Ji9nFCYl7Bk)
- [RPKM, FPKM and TPM](https://youtu.be/TTUrtCY2k-w)
