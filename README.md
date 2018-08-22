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

To quantify expression of transcript/genes among different conditions:
-	Count the number of mapped reads on each transcript  
-	Quantify gene-level expression with GTF (gene transfer format) files. 

###### Input: grouping info of the samples (csv file)	:heavy_minus_sign:	Output: SAM files (1 per sample)

#### Step 7. Run  differential expression analysis with Ballgown

- [ ]	(a) Load relevant R packages (ballgown, RSkittleBrewer, genefilter, dplyr, devtools)
```R
	module load R
	R
	
	# R commands are indicated with '>' below:

	>library("ballgown")
	>library("RSkittleBrewer") # for color setup
	>library("genefilter") # faster calculation of mean/variance
	>library("dplyr") # to sort/arrange results
	>library("devtools")  # reproducibility/installing packages
```

- [ ]	(b) Load phenotype data 
	- In your future experiment, create your own phenotype data specifying the sample conditions you would like to compare. 
	- Each sample information presents on each row of the file 
```R
	> pheno_data = read.csv("phenodata.csv")
	> head(pheno_data)
```

- [ ]	(c) Read in the expression data that were calculated by Stringtie in previous step 6-(c)
	-	IDs of the input files should be matched with the phenotype information. 
	-	Ballgown also supports Cufflinks/RSEM data
```R
	> chrX <- ballgown(dataDir="ballgown", samplePattern="ERR", pData=pheno_data)
	> str(chrX)
```

- [ ]	(e) Filter to remove low-abundance genes 
	-	Genes often have very few or zero counts.
	-	We can apply a variance filter for gene expression analysis. 
```R
	> chrX_filtered <- subset(chrX, "rowVars(texpr(chrX)) >1", genomesubset=TRUE)
	> str(chrX_filtered)
```


#### Step 8. Identify transcripts/genes that show statistically significant differences between groups

- [ ]	(a) Identify **transcripts** that show statistically significant differences between groups 
```R
	* We will look for transcripts that are differentially expressed between sexes, while correcting for any differences in expression due to the condition variable. 

	> results_transcripts <- stattest(chrX_filtered, feature="transcript", covariate="sex", adjustvars=c("condition"), getFC=TRUE, meas="FPKM")
	> head(results_transcripts)

	* Add gene names and gene IDs to the results:

	> results_transcripts <- data.frame(geneNames=ballgown::geneNames(chrX_filtered), geneIDs=ballgown::geneIDs(chrX_filtered), results_transcripts)
	> head(results_transcripts)
```

- [ ] 	(b) Identify **genes** that show statistically significant differences between groups 
```R
	> results_genes <- stattest(chrX_filtered, feature="gene", covariate="sex", adjustvars=c("condition"), getFC=TRUE, meas="FPKM")
	> head(results_genes)
```


-	Several **measurement criteria for RNA-seq quantification**:
	-	**_RPKM (reads per kilobase of exon model per million reads)_** adjusts for feature length and library size by sample normalization 
	-	**_FPKM (fragment per kilobase of exon model per million mapped reads)_** adjusts sample normalization of transcript expression (= similar to RPK) 
		-	RPKM = FPKM (in Single-end sequencing) 
		-	_FPKM can be translated to TPM_
	-	**_TPM (transcript per million)_** is used for measuring RNA-seq gene expression by adjusting transcript differences with overall read # in library: useful in comparing inter-sample comparison with different origins/compositions 
		-	Gene length is not important for inter-sample gene expression comparison, but important in ranking intra-sample gene expression 
		-	All the measures above are not useful for the samples with transcript variances. 

-	In differential expression analysis, we have to eliminate systematic effects that are not due to biological causal differences of interest. _(Normalization)_ We should condition the non-biological differences such as sequencing depth, gene’s length, and variability. Therefore we must calculate the fraction of the reads for each gene compared to the total amount of reads and to the whole RNA library. 
-	Tests for significance must rely on assumptions about the underlying read distributions. The negative binomial distribution is often used for modeling gene expression between biological replicates, because it better accounts for noise than the Poisson distribution, which would otherwise be applicable as an approximation of the binomial distribution (presence of individual reads or not) with a large n (read library) and a small np. Another common assumption is that the majority of the transcriptome is unchanged between the two conditions. If these assumptions are not met by the data, the results will likely be incorrect. This is why it’s important to examine and perform QC on the expression data before running a differential expression analysis!
-	It is highly recommended to have at least two replicates per group. 

-	We can compare the transcripts that are differentially expressed between groups, while correcting for any different expression due to _**‘confounding’**_ variable. 
	-	Ballgown can look at the confounder-adjusted fold change (FC) between the two groups by setting getFC=TRUE parameter in stattest() function. 

-	_**For small sample size (n<4 per group)**_, it is often better to perform regularization than standard linear model-based comparison as Ballgown does. Like “limma-voom” package in Bioconductor, DESeq, edgeR for gene/exon counts are the mostly used ones. (not appropriate for FPKM abundance estimates) 

-	_Other software options: **HTseq-count, DESeq2, edgeR, Kallisto, RSEM** (use expectation maximation to measure TPM value), **NURD** (Transcript expression measure in SE reads with low memory and computing cost), **Cufflinks** (using Tophat mapper for mapping, expectation-maximization algorithm)_

#### Step 9. Explore the results! 
```R
	* Sort the results from the smallest-largest p-value
	
	> results_transcripts <- results_transcripts[order(results_transcripts$pval),]
	> results_genes <- results_genes[order(results_genes$pval),]

	* What are the top transcript/gene expressed differently between sexes? 

	> head(results_transcripts)
	> head(results_genes)

	* (cf) You can also try filtering with q-value (<0.05) with subset() function.
	* Save the analysis results to csv files:

	> write.csv(results_transcripts, file="DifferentialExpressionAnalysis_transcript_results.csv", row.names=FALSE)
	> write.csv(results_genes, file="DifferentialExpressionAnalysis_gene_results.csv", row.names=FALSE)
	> save.image()			# your workspace will be saved as '.RData' in current working directory

```


___
### 5. Visualization of the results

#### Step 10. Choose your environment for Visualization
You can choose to use either __**IGV or UCSC Genome browser**__ for visualizing your overall outcome. Not only for visualizing the expression differences, the step is also essential for checking additional quality control criteria such as PCR duplication caused by variant calling. In our examples, we will use R Ballgown package for RNA-seq analysis specific visualization.

For small and moderately sized interactive analysis: 
-	Go to Rstudio-Quest analytics node on your browser [[https://rstudio.questanalytics.northwestern.edu/auth-sign-in]](https://rstudio.questanalytics.northwestern.edu/auth-sign-in)
-	You might have to re-install the required R packages for differential data analysis described above
```R
	> setwd("/YourWorkingDirectory/")
	> load(".RData")
```

For large sized interactice analysis that might require over 4GB of RAM or more than 4+ cores:
-	Request an interactive session on a compute node [[link]](https://kb.northwestern.edu/69247)
```
	msub -I -l nodes=1:ppn=4 -l walltime=01:00:00 -q genomics -A b1042
```



#### Step 11. 

In our example script, we will explore:
	-	**Boxplot – Distribution of gene abundances across samples**:
		-	Variety of measurements can be compared and visualized other than FPKM values, such as splice junction, exon and gene in the dataset. 
		-	Log transformation is required sometimes to plot some FPKM data = 0. 
	-	**Boxplot – individual expression of a certain transcript between groups**. 
	-	Plot the structure/expression levels in a sample of all transcripts that share the same gene locus.
		-	We can plot their structure and expression levels by passing the gene name and the Ballgown object to the plotTranscripts function.
	-	Plot average expression levels for all transcripts of a gene within different groups

- [ ]	(a). Plot for distribution of gene abundances across samples:
	- In this example, we compare the FPKM measurements for the transcripts colored by 'sex' varaible in phenotype file. 
```R
	> fpkm <- texpr(chrX, meas='FPKM')
	> fpkm <- log2(fpkm +1)
	> boxplot(fpkm, col=as.numeric(pheno_data$sex), las=2,ylab='log2(FPKM+1)')
```

- [ ]	11-2. Plot for individual expression of a certain transcript between groups: 
```R
	* Setup palette with your favorite colors

	> coloring <- c('darkgreen', 'skyblue', 'hotpink', 'orange', 'lightyellow')
	> palette(coloring)

	* Choose your transcript of interest

	*	In this example, by looking head(results_transcripts), I choose to draw the 5th most differientially expressed transcript. (gene name "XIST")
	*	You can also decide the transcript/gene of your interest. If you want to draw 10th transcript in your dataset:
	(ex) > ballgown::transcriptNames(chrX)[10]

	> which(ballgown::geneNames(chrX)=="XIST")	
	
	*	Find the row number of the interested gene in dataset
	# 1492 here 
	
	> ballgown::transcriptNames(bg_chrX)[1492]	# get the transcript name in the gene 
	> plot(fpkm[1492,] ~ pheno_data$sex, border=c(1,2), main=paste(ballgown::geneNames(chrX)[1492], ' : ',ballgown::transcriptNames(chrX)[1492]), pch=19, xlab="sex", ylab='log2(FPKM+1)')
	> points(fpkm[1492,] ~ jitter(as.numeric(pheno_data$sex)), col=as.numeric(pheno_data$sex))
```
-	The output plot shows the name of the transcript (NR_001564) and the name of the gene (XIST) that contains it. 
	- 	Can you tell the exclusive expression of XIST in females? (c.f. In females, the XIST gene is expressed exclusively from the inactive X chromosome, and it is essential for the initiation and spread of X inactivation, which is an early developmental process that transcriptionally silences one of the pair of X chromosomes)

- [ ]	11-3. Plot the structure/expression levels in a sample of all transcripts that share the same gene locus:
```R
	* We choose sample 'ERR204916' for plotting structure/expression level in the same genomic position. 
	
	> plotTranscripts(ballgown::geneIDs(chrX)[1429],chrX, main=c("Gene XIST in sample ERR204916"), sample=c("ERR204916"))
```
	- The output plot shows one transcript per row, colored by its FPKM level. 

- [ ]	11-4. Plot the average expression levels for all transcripts of a gene within different groups:

	- Using plotMeans() function, specify which gene to plot and which variable to group by. 
```R
	> geneIDs(chrX)[1492]
		"MSTRG.501" 
	> plotMeans('MSTRG.501', chrX_filtered, groupvar="sex", legend=FALSE)
	> plotMeans(ballgown::geneIDs(bg_chrX)[1492], chrX, groupvar="sex", legend=FALSE)
```
###### Have fun playing! 


_Other Software Options_:
1.	_Read-level visualization software: **ReadXplorer, UCSC genome browser, integrative Genomics Viewer (IGV), Genome Maps, Savant**_ 
2.	_Gene expression analysis software: **DESeq2, DEXseq **_
3.	_**CummeRbound, Sashimi plot** (junction reads will be more intuitive and aesthetic), **SplicePlot** (can get sashimi, structure, hive plot for sQTL), **TraV** (integrates all the data analysis for visualization but cannot be used for huge genome)_


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
