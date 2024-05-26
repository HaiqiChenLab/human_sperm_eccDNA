#!/bin/bash

#SBATCH --job-name=job name                                # job name
#SBATCH --partition=xxx                                 # select partion
#SBATCH --nodes=x                                         # number of nodes requested by user
#SBATCH --time=0-999:00:00                                 # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=serialJob.%j.out                         # standard output file name
#SBATCH --error=serialJob.%j.err                         # standard error output file name
#SBATCH --mail-user=your email                            # specify an email address
#SBATCH --mail-type=FAIL,END                              # send email when job status change (start, end, abortion and etc)



module load homer/4.9

#########################################

cd /your/directory/

# calculate nucleotide frequency at the 20 bp flanking the start and end position of each eccDNA
homerTools freq -format fasta eccDNA_end_junction.fa -gc eccDNA_end_junction_GCdistribution.txt -o sperm_eccDNA_end_junction_positionFrequency.txt

homerTools freq -format fasta eccDNA_start_junction.fa -gc eccDNA_start_junction_GCdistribution.txt -o sperm_eccDNA_start_junction_positionFrequency.txt