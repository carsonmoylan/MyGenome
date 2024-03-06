# MyGenome
Analyses for CS485G genome assembly

## 1. Analysis of sequence quality
The F1 and R1 sequence datasets using FASTQC: 
```bash
ssh -Y cmo343@cmo343.cs.uky.edu
cd MyGenome
fastqc &
```
Load R1 and R2 datasets into GUI interface.

Screenshot of UFVPY231_1_paired.fastq:
![R1_paired.PNG](/data/R1_paired.PNG)

Screenshot of UFVPY231_2_paired.fastq:
![R2_paired.PNG](/data/R2_paired.PNG)

## 2. Ran trimmomatic
```bash
java -jar ~/sequences/trimmomatic-0.38.jar PE -threads 16 -phred33 -trimlog file.txt UFVPY231_1.fq UFVPY231_2.fq UFVPY231_1_paired.fastq UFVPY231_1_unpaired.fastq UFVPY231_2_paired.fastq UFVPY231_2_unpaired.fastq ILLUMINACLIP:adaptors.fasta:2:30:10 SLIDINGWINDOW:20:20 MINLEN:100
```

## 3. Count the number of forward reads remaining
```bash
grep -c '^@A00261' UFVPY231_1_paired.fastq
```
Output: 9,982,400

## 4. Count the total number of bases in both files
```bash
cat UFVPY231_1_paired.fastq UFVPY231_2_paired.fastq | awk 'NR%4==2 {total += length($0)} END {print total}'
```
The above command tells us that there are 2,989,505,387 total bases
