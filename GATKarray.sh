#!/bin/bash
#SBATCH --job-name=GATKtestArray
#SBATCH --ntasks=20
#SBATCH --partition=general          # name of partition to submit job
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL              # will send email for begin,end,fail
#SBATCH --mail-user=sbw0033@auburn.edu
#SBATCH --output=GATK_test_%A_%a.output 		#Changes the output to correspond to each subjob
#SBATCH --error=GATK_test_%A_%a.error 		#Changes the error to correspond to each subjob
#SBATCH --array=1,14,141,153,17,18,2,22,23,24,44,47,6,60,64,7,76,9,93

referenceChr7="/home/sbw0033/TerrapinGenomics/Data/PaintedTurtleResources/PaintedTurtleChr7.fa.gz"
FinalBamDirectory="/scratch/sbw0033/FinalAlignments"
OutputGvcfDirectory="/scratch/sbw0033/Gvcfs"

#Load modules using specific version to ensure future compatibility
module load gatk/4.1.9.0
module load samtools/1.11

########################
#Get into the directory with all the Cleaned reads
cd /scratch/sbw0033/FinalAlignments

#Run haplotypeCaller on every single sample individually
gatk HaplotypeCaller --input /scratch/sbw0033/FinalAlignments/"${SLURM_ARRAY_TASK_ID}"_0.bam --output /scratch/sbw0033/Gvcfs/"${SLURM_ARRAY_TASK_ID}"_0.g.vcf.gz --reference /home/sbw0033/TerrapinGenomics/Data/PaintedTurtleResources/PaintedTurtleChr7.fa --native-pair-hmm-threads 20

cd /scratch/sbw0033/Gvcfs

#Make a cohort map with samples associated with files
echo "Sample"_"${SLURM_ARRAY_TASK_ID}"'\t'"${SLURM_ARRAY_TASK_ID}"_0.g.vcf.gz >> /home/sbw0033/cohort.sample_map_0
