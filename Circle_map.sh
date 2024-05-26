#!/bin/bash

#SBATCH --job-name=job name                                # job name
#SBATCH --partition=xxx                                 # select partion
#SBATCH --nodes=x                                         # number of nodes requested by user
#SBATCH --time=0-999:00:00                                 # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=serialJob.%j.out                         # standard output file name
#SBATCH --error=serialJob.%j.err                         # standard error output file name
#SBATCH --mail-user=your email                            # specify an email address
#SBATCH --mail-type=FAIL,END                              # send email when job status change (start, end, abortion and etc)


# Set the directory containing the gzipped fastq files
fastq_dir="/your/directory"


module load BWA/0.7.5
module load samtools/1.6
module load python/3.8.x-anaconda
module load bbmap/38.46

#########################################

#activate the conda environment containing circle-map/1.1.4
conda activate 'env name'

# Loop through all gzipped fastq files in the fastq directory
fastq1 = "$fastq_dir"/*_L2_1.fq.gz

# Extract the base name of the fastq file
base=$(basename "$fastq1" _L2_1.fq.gz)

# Construct the name of the second pair-end gzipped fastq file
fastq2="$fastq_dir/${base}_L2_2.fq.gz"
  
# trim adaptor sequences
bbduk.sh -Xmx1g in1="${fastq1%}" in2="${fastq2%}" out1="$fastq_dir/${base}_trimmed_R1_001.fastq.gz" out2="$fastq_dir/${base}_trimmed_R2_001.fastq.gz" ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
  
#align trimmed fastq files
bwa mem -t 46 hg38.fa "$fastq_dir/${base}_trimmed_R1_001.fastq.gz" "$fastq_dir/${base}_trimmed_R2_001.fastq.gz" > unknown_circle.sam
  
#sort sam files by name and coordinates
samtools sort -n -@ 46 -o qname_unknown_circle.bam unknown_circle.sam
samtools sort -@ 46 -o sorted_unknown_circle.bam unknown_circle.sam
  
#Run circle-map
Circle-Map ReadExtractor -i qname_unknown_circle.bam -o circular_read_candidates.bam
samtools sort -@ 46 -o sort_circular_read_candidates.bam circular_read_candidates.bam
samtools index -@ 46 sort_circular_read_candidates.bam
samtools index -@ 46 sorted_unknown_circle.bam
Circle-Map Realign -i sort_circular_read_candidates.bam -qbam qname_unknown_circle.bam -sbam sorted_unknown_circle.bam -fasta hg38.fa -o my_unknown_circle.bed
  
conda deactivate