# dudeML
A python script for the detection of duplications and deletions using machine learning.
By- RIYA KUMARI SINGH

# 1. Requirements
A number of programs are required to install initially.
Most of them are based only on Linux OS, so it is mandatory to run the code in Linux environment
## Python-based:
* Python3
* pandas
* numpy
* scikit-learn
* Biopython
## External
* bedtools
* short-read simulator (e.g. wgsim)
* short-read aligner (e.g. BWA)
* short-read parser (e.g. SAMtools)

ONE TIME SETUPS:
	
	conda create -n dudeml python=3.7 anaconda
	source activate dudeml
	conda install -n dudeml pandas numpy scikit-learn biopython
	conda install -c biopython -n dudeml wgsim bwa samtools bedtools
	conda install -n dudeml -c bioconda samtools=1.11 --force-reinstall
	

# 2. Functions
Subfunctions of the dudeML script, each function provides a specific role and requires differing inputs (described by running the help function of the tool).
1. **winStat**
Finds the average coverage of a window in the chromosome, relative to the average chromosome coverage.
2. **fvecSample**
Reformats the bedfile with coverage information into sets of windows surrounding a focal window.
3. **fvecTrain**
Reformats the bedfile with coverage information into sets of windows surrounding a focal window. Also includes information on if a CNV is present in the window, and the estimated number of copies of that window per chromosome.
4. **simCNV**
Generate coordinates for random CNVs in the fasta file input, after accounting for repetitive content.
5. **recreateTotal**
If already known deletions and duplications are being used in the training data, this function skips the simulation of CNVs and instead generates a file with positions where CNVs should be, for simChrs.
6. **simChr**
Masks known repetitive content in the given chromosome and generates chromosomes with simulated CNVs.
7. **classify**
Given a training file (generated in formatTrain) and a sample file (generated in formatSample), will predicted windows with CNVs based on coverage and standard deviation of coverage.
8. **predict**
Given a training file (generated in formatTrain) and a sample file (generated in formatSample), will predicted windows with CNVs based on coverage and standard deviation of coverage.
9. **simReads**
Simulates read pairs of chosen length to a certain coverage of the chosen chromosome, requires WGsim.
10. **subTrain**
Downsample a training file by a certain percentage or to a certain number of each category. 
11. **summarize**
Combines called CNVs, and if known CNVs are provided, tells you if called CNVs are True-positives or otherwise.
12. **winStatExtra**
Creates summary windows based on known coverage estimates.
13. **covSummary**
Summarizes the coverage of each chromosome in the genomeCoverageBed file.

# 3. Input file formats
## Fasta
The reference sequences for mapping and for generating training files, in the following format:
>\>Chr1

>aagagcctatatca

>\>Chr2

>aagagcctatatca

## BAM files
A binary file containing short reads mapped to a repeat masked reference genome, used as input for genomeCoverageBed within dudeML.


## Duplications
A bed file of known (or simulated) duplications, with the number of copies per chromosome and the frequency of the duplication in the sampled data (e.g. the number of chromosomes with this duplication/ the total number of chromosomes):
>Chr1	1000	1344	dup	3	1.0

>Chr1	2455	6700	dup	2	0.5

>Chr1	34501	36119	dup	2	1.0

>Chr1	45117	48932	dup	4	0.5

## Deletions
A bed file of known (or simulated) deletions, with the number of copies per chromosome and the frequency of the deletion in the sampled data (e.g. the number of chromosomes with this deletion/ the total number of chromosomes):
>Chr1	1000	1344	del	0	1.0

>Chr1	2455	6700	del	0	0.5

>Chr1	34501	36119	del	0	1.0

>Chr1	45117	48932	del	0	0.5

## Fvec files
A type of bed file, containing information on CNVs and copy number if a training file, and containing strain ID if a sample file. Reformats the bedfile from winStat to a feature vector, summarizing the windows around the focal window.
>2L	22240500	22240550	N	1.0	0.64	0.151	0.92	0.071	1.04	0.134	1.04	0.101	1.2	0.075	1.12	0.112	1.44	0.132	1.12	0.168	1.12	0.124	1.6	0.163	1.42	0.145

# 4. A simple walkthrough
## A. Simulate training data
I used the provided DiNV virus genome and repeat locations.
Following that, I simulated CNVs for a homozygous individual, requiring 1 set of chromosomes to be generated. 
It's important that the training data is as similar as possible to the sample being tested.

Cloning the project repository after changing the environment

    source activate dudeml
    git clone https://github.com/riyasingh07/MyDudeML.git
    cd MyDudeML/
    
Masking and Indexing dataset using  Bedtools and BWA
	
    maskFastaFromBed -fi DiNV_CH01M.fa -bed DiNV.sat.bed -fo DiNV_CH01M.fa.masked
    bwa index DiNV_CH01M.fa.masked
 
Creating Train/Test datasets
 
    for i in train test
    do
    mkdir ${i}_sim
    python3 dudeML.py simCNV -fasta DiNV_CH01M.fa -CNV 50 -d ${i}_sim -N 1
    python3 dudeML.py simChr -fasta DiNV_CH01M.fa -cnvBed ${i}_sim/total.1.bed -id ${i} -d ${i}_sim
    done

## B. Estimating coverage in training and test data
  I next simulated reads for the custom chromosomes containing CNVs using WGSIM within dudeML and following mapping, used bedtools to calculate coverage per site, 
  I then used a custom python script to find the mean coverage of each chromosome (to find the relative coverages of each window). 
  In this case a homozygote was simulated.
  
    for i in train test
	do
	python3 dudeML.py simReads -fasta DiNV_CH01M.fa -cov 20 -d ${i}_sim -id ${i} -RL 100
	bwa mem -t 4 DiNV_CH01M.fa.masked ${i}_sim/${i}_20_1.fq ${i}_sim/${i}_20_2.fq | samtools view -Shb - | samtools sort - > ${i}_sim/total.bam
	python3 dudeML.py winStat -fasta DiNV_CH01M.fa -i${i}_sim/total.bam -o ${i}_sim/total_50.bed -w 50 -s 50
	done

## C. Reformatting sample and training datasets.    
  I then removed repetitive regions, reformatted the data (to show the relative coverage of the focal window and the 5 windows on each side). 
  I also prepared the data to filter and extract the regions with known duplications or deletions in training the file. 
  I also labelled CNVs in the test dataset for comparison later. 
  As repetitive regions can be tricky to deal with, I ignored windows with more than 50% of the window masked as an (-c 0.5).

	python3 dudeML.py fvecTrain -i train_sim/total_50.bed -o train_sim/total_50train.bed -w 50 -TE DiNV.sat.bed -dups train_sim/dup.1.bed -dels train_sim/del.1.bed  -windows 5 -c 0.5
	python3 dudeML.py fvecSample -i test_sim/total_50.bed -w 50 -o test_sim/total_50sample.bed -id test_sim -TE DiNV.sat.bed -windows 5 -c 0.5

## D. Predicting CNVs using the generated files.  
Following this, I created a classifier from one of the training features vector files generated and test out predictions of CNVs in the other file.

	python3 dudeML.py classify -i train_sim/total_50train.bed -o train_sim/total_50train.sav
	python3 dudeML.py predict -i test_sim/total_50sample.bed -t train_sim/total_50train.sav -o test_sim/total_50pred.bed
   
   

## E.(OPTIONAL) Training and Testing a real data in dudeML
Real data with known structural variants can be used as a training set using the following pipeline.
	
	python dudeML.py winStat -i knownCNV.bam -o knownCNV_50.bed -w 50
	python3 dudeML.py fvecTrain -i knownCNV_100.bed -o knownCNV_50.bed -w 50 -TE repeats.gff -dups knownDUP.bed -dels knownDEL.bed -c 0.5

Real data can also be used as the test sample to identify unknown CNVs.
	
	python dudeML.py winStat -i unknownCNV.bam -o unknownCNV_50.bed -w 50
	python3 dudeML.py fvecSample -i unknownCNV_50.bed -o unknownCNV_50_sample.bed -w 50 -TE repeats.gff -c 0.5


Then the generated training and test sets can be used to find CNVs.
	
	python3 dudeML.py predict -i unknownCNV_50_sample.bed -t knownCNV_50_train.bed -o unknownCNV_50_pred.bed

