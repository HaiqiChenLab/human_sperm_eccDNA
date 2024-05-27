#!/bin/bash

#SBATCH --job-name=job name                                # job name
#SBATCH --partition=xxx                                 # select partion
#SBATCH --nodes=x                                         # number of nodes requested by user
#SBATCH --time=0-999:00:00                                 # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=serialJob.%j.out                         # standard output file name
#SBATCH --error=serialJob.%j.err                         # standard error output file name
#SBATCH --mail-user=your email                            # specify an email address
#SBATCH --mail-type=FAIL,END                              # send email when job status change (start, end, abortion and etc)


# Set the directory containing the genome index files
idx_dir="/your/genome_sequence_file"

# Set the directory for the fastq files
fastq_dir="/your/fastq_files"


module load segemehl/0.3.4
module load bedtools/2.29.2
module load bbmap/38.46 

#########################################

#index the genome
segemehl.x -x "$idx_dir/hg38.idx" -d "$idx_dir/hg38.fa" 

# trim adaptor sequences
bbduk.sh -Xmx1g in1="$fastq_dir/RNA_R1.fastq.gz" in2="$fastq_dir/RNA_R2.fastq.gz" out1="$fastq_dir/RNA_trimmed_R1.fastq.gz" out2="$fastq_dir/RNA_trimmed_R2.fastq.gz" ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo

#alignment
segemehl.x -i "$idx_dir/hg38.idx" -d "$idx_dir/hg38.fa" -q "$fastq_dir/RNA_trimmed_R1.fastq.gz" -p "$fastq_dir/RNA_trimmed_R2.fastq.gz" -D 0 -S > "$fastq_dir/RNA.sam"
  
#Delete reads that are regular, collinear splits
awk '!/A1;R;/' "$fastq_dir/RNA.sngl.bed" > "$fastq_dir/RNA_discord.sngl.bed" 

#Delete the header of the BED file
sed -i '1d' "$fastq_dir/RNA_discord.sngl.bed" 

#Combine and intercept the eccDNA coordinates and the split read summary BED files. Transcript must overlap eccDNA completely.
bedtools intersect -F 1 -wo -a "$fastq_dir/ecDNA_coordinates.bed "-b "$fastq_dir/RNA_discord.sngl.bed" | awk '$2-$5 < 0' | awk '$3-$6 > 0' > "$fastq_dir/intercept.bed"

