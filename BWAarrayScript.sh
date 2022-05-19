#!/bin/bash
#SBATCH --job-name=BWAtestArray 
#SBATCH --ntasks=10
#SBATCH --partition=general          # name of partition to submit job
#SBATCH --time=8:00:00
#SBATCH --mail-type=ALL              # will send email for begin,end,fail
#SBATCH --mail-user=sbw0033@auburn.edu
#SBATCH --output=test_%A_%a.output 		#Changes the output to correspond to each subjob
#SBATCH --error=test_%A_%a.error 		#Changes the error to correspond to each subjob
#SBATCH --array=1,14,141,153,17,18,2,22,23,24,44,47,6,60,64,7,76,9,93

#Load modules using specific version to ensure future compatibility
module load bwa/0.7.17
module load samtools/1.11

########################
#Get into the directory with all the raw reads
cd /scratch/sbw0033/FastqCopies/clean/


bwa mem -M -t 10 \
	/home/sbw0033/TerrapinGenomics/Data/PaintedTurtleResources/BWA_index_chr7/Chr7 \
	/scratch/sbw0033/FastqCopies/clean/"${SLURM_ARRAY_TASK_ID}"/"${SLURM_ARRAY_TASK_ID}"_FW_reads.fq.gz /scratch/sbw0033/FastqCopies/clean/"${SLURM_ARRAY_TASK_ID}"/"${SLURM_ARRAY_TASK_ID}"_RV_reads.fq.gz \
	2> /home/sbw0033/bwaSample_"${SLURM_ARRAY_TASK_ID}".err \
	> /home/sbw0033/TerrapinGenomics/Results/bwaSample_"${SLURM_ARRAY_TASK_ID}".sam
