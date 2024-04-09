# MyGenome (UFVPY231)
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

## 5. Submit first assembly of step size 10 to SLURM
```bash
sbatch velvetoptimiser_noclean.sh UFVPY231 61 131 10
```
Velvet hash value = 101

## 6. Submit second assembly of step size 2 to SLURM
```bash
sbatch velvetoptimiser_noclean.sh UFVPY231 91 111 2
```
Velvet hash value = 99

## 7. Rename sequence headers to a standard format
```bash
perl /project/farman_s24cs485g/SCRIPTs/SimpleFastaHeaders.pl UFVPY231.fasta
```

## 8. Check Genome completeness using BUSCO
```bash
sbatch BuscoSingularity.sh UFVPY231/velvet_UFVPY231_91_111_2_noclean/UFVPY231_nh.fasta
```

## ? other steps??

## 9. Remove contigs less than 200 base pairs long
```bash
perl CullShortContigs.pl UFVPY231_nh.fasta
```

## 10. BLAST the MoMitochondrion.fasta sequence against the final genome assembly
```bash
blastn -query MoMitochondrion.fasta -subject UFVPY231_nh.fasta -evalue 1e-50 -max_target_seqs 20000 -outfmt '6 qseqid sseqid slen length qstart qend sstart send btop' -out MoMitochondrion.UFVPY231.BLAST
```

## 11. Export a list of contigs that mostly comprise mitochondrial sequences
```bash
awk '$4/$3 > 0.9 {print $2 ",mitochondrion"}' MoMitochondrion.UFVPY231.BLAST > UFVPY231_mitochondrion.csv
```

## 12. BLAST the genome assembly against a repeat-masked version of the B71 reference genome
```bash
blastn -query B71v2sh_masked.fasta -subject UFVPY231_Final.fasta -evalue 1e-50 -max_target_seqs 20000 -outfmt '6 qseqid sseqid qstart qend sstart send btop' -out B71v2sh.UFVPY231.BLAST
```

## 13. Identify genetic variants between the B71v2sh genome and the genome assembly
```bash
sbatch CallVariants.sh UFVPY231_BLASTS
```
