# Study gene expression using RNA-seq analysis


## Table of Contents
### 1. Introduction to RNA-seq (Experimental design)
### 2. Set up environment on High Performance Computing and the Genomic Compute Cluster 
### 3. Quality Check/ Alignment/ Estimate transcript abundances of RNA-seq reads
### 4. Expression and Differential Expression analysis
### 5. Visualization of the results
* Appendix. Reference-free expression estimation (Kallisto) 

___

### 1. Introduction to RNA-seq (Experimental design)

There are several sequencing analysis methods available based on the **goal of your experiment**:

| Types of Sequencing Analysis | Notes |
| :------------ | :-------------|
|	_**RNA-seq**_ | _co-expression networks, differentially expressed genes_ |
|	_**Chip-Seq \| ATAC-seq \| DNase-seq**_ | _Gene regulation dynamics, chromatin modeling_ |
|	_**Whole-genome Seq**_ | _Rare variant analysis_ |

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
      -	'samples' directory contains paired-end RNA-seq reads for 6 samples, 3 male and 3 female subjects from YRI (Yoruba from Ibadan, Nigeria) population. 
      -	‘indexes’ directory contains the indexes for chromosome X for HISAT2. 
      -	‘genome’ directory contains the sequence of human chromosome X (GrCH38 build 81)
      -	‘genes’ directory contains human gene annotations for GrCH38 from RefSeq database. 
      -	‘mergelist.txt’ and ‘geuvadis_phenodata.csv’ are exemplary scripts that you might want to write yourself in a text editor. 
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
- [ ] _[Optional]_  Install HTSeq, a Python package (http://htseq.readthedocs.io/en/master/install.html#install)


___
### 3. Quality Check/ Alignment/ Estimate transcript abundances of RNA-seq reads
___
### 4. Expression and Differential Expression analysis
___
### 5. Visualization of the results
___

