#!/bin/bash

#SBATCH --job-name=                # job name
#SBATCH --partition=                                
#SBATCH --nodes=                                        # number of nodes requested by user
#SBATCH --time=                                 # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=                         # standard output file name
#SBATCH --error=                        # standard error output file name
#SBATCH --mail-user=          # specify an email address
#SBATCH --mail-type=                                  # send email when job status change (start, end, abortion and etc)


# Set the directories
WORKDIR=""
fastq_dir=""
GENOME_FASTA=""


module load BWA/0.7.5
module load samtools/1.6
module load python/3.8.x-anaconda
module load bbmap/38.46
module load seqtk/1.2-r94

#########################################

cd "$WORKDIR"

source your_conda_directory/conda.sh

conda activate your_circle_map_virtual_environment

# Loop through all gzipped fastq files in the fastq directory
for fastq1 in "$fastq_dir"/*.fq.gz; do
  # Extract the base name of the fastq file
  base=$(basename "$fastq1" _1.fq.gz)
  # Construct the name of the second pair-end gzipped fastq file
  fastq2="$fastq_dir/${base}_2.fq.gz"
  
  bbduk.sh -Xmx1g in1="${fastq1%}" in2="${fastq2%}" out1="$fastq_dir/${base}_trimmed_R1_001.fastq.gz" out2="$fastq_dir/${base}_trimmed_R2_001.fastq.gz" ref=your_adaptor_list_directory/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
  
  
  bwa mem -t 60 GRCm38.primary_assembly.genome.fa "$fastq_dir/${base}_trimmed_R1_001.fastq.gz" "$fastq_dir/${base}_trimmed_R2_001.fastq.gz" > unknown_circle.sam
 
  echo "Step 1: Creating name-sorted BAM for ReadExtractor..."
  samtools sort -n -@ 60 -o qname_unknown_circle.bam unknown_circle.sam

  echo "Step 2: Extracting circular read candidates with Circle-Map..."
  # We still use Circle-Map's ReadExtractor as it is efficient
  Circle-Map ReadExtractor -i qname_unknown_circle.bam -o circular_read_candidates.bam
  
  echo "Step 3: Sorting candidate reads for our custom script..."
  samtools sort -@ 60 -o sort_circular_read_candidates.bam circular_read_candidates.bam
  samtools index -@ 60 sort_circular_read_candidates.bam

  echo "Step 4: Running custom realignment and eccDNA detection script..."
  # This replaces the `Circle-Map Realign` command.
  ./realign_from_scratch.py \
      -i sort_circular_read_candidates.bam \
      -o custom_eccDNA_calls.bed \
      --min_reads 3 \
      --window 100 \
      --min_mapq 20
done
  
conda deactivate

echo "Custom realignment job finished successfully."