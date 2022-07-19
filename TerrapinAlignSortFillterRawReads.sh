#!/bin/bash
#SBATCH --job-name=BWAtestArray
#SBATCH --ntasks=14
#SBATCH --partition=general          # name of partition to submit job
#SBATCH --time=96:00:00
#SBATCH --mail-type=ALL              # will send email for begin,end,fail
#SBATCH --mail-user=sbw0033@auburn.edu
#SBATCH --output=test_%A_%a.output 		#Changes the output to correspond to each subjob
#SBATCH --error=test_%A_%a.error 		#Changes the error to correspond to each subjob
#SBATCH --array=1,14,141,153,17,18,2,22,23,24,44,47,6,60,64,7,76,9,93

RawReadDirectory="/scratch/sbw0033/TerrapinGenomics/Data/FastqCopiesFinal/clean"
ReferenceIndexPrefix="/scratch/sbw0033/TerrapinGenomics/Data/TerrapinGenomeBWAIndex/Hap1_index"
AlignmentDirectory="/scratccd h/sbw0033/TerrapinAllAlignments/TerrapinRawAlignments"
SortedBamDirectory="/scratch/sbw0033/TerrapinAllAlignments/TerrapinSortedAlignments"
CleanedBamDirectory="/scratch/sbw0033/TerrapinAllAlignments/TerrapinCleanedAlignments"
FinalAlignmentDirectory="/scratch/sbw0033/TerrapinAllAlignments/TerrapinFinalAlignments"

#Think about having a separate job for forward and reverse reads

#Load modules using specific version to ensure future compatibility
module load bwa/0.7.17
module load samtools/1.11
module load picard/2.23.9
module load java/15.0.1

#index the reference genome
#cd /scratch/sbw0033/TerrapinGenomics/Data/
#bwa index TerrapinReferenceHaplotype1.fasta.gz -p Hap1_index

########################
#Get into the directory with all the raw reads
cd "$RawReadDirectory"

#Run BWA
bwa mem -M -t 14 \
	"${ReferenceIndexPrefix}" \
	"${RawReadDirectory}"/"${SLURM_ARRAY_TASK_ID}"/"${SLURM_ARRAY_TASK_ID}"_FW_reads.fq.gz "${RawReadDirectory}"/"${SLURM_ARRAY_TASK_ID}"/"${SLURM_ARRAY_TASK_ID}"_RV_reads.fq.gz \
	2> /scratch/sbw0033/bwaSample_"${SLURM_ARRAY_TASK_ID}".err \
	> "${AlignmentDirectory}"/"${SLURM_ARRAY_TASK_ID}".sam

#Can Also add read group information this way
#java -Xmx8g -jar picard.jar AddOrReplaceReadGroups I="Sample1_sorted.bam" O="Sample1_IDed.bam" RGPU="RunBarcode1" RGSM="ReadGroupSample1" RGPL="Illumina" RGLB="ReadGroupLibrary"

cd "${AlignmentDirectory}"

#First, convert sam to .bam
samtools view -bS -@ 14 "${SLURM_ARRAY_TASK_ID}".sam > "${SLURM_ARRAY_TASK_ID}".bam

#First, sort all of the bams
samtools sort "${AlignmentDirectory}"/"${SLURM_ARRAY_TASK_ID}".bam --threads 14 -o "${SortedBamDirectory}"/"${SLURM_ARRAY_TASK_ID}"_sorted.bam

cd "${SortedBamDirectory}"

#Now, keep alignments with a map quality score greater than or equal to 30
samtools view --threads 14 -q 30 -b "${SLURM_ARRAY_TASK_ID}"_sorted.bam > "${SortedBamDirectory}"/"${SLURM_ARRAY_TASK_ID}"_sorted.q30.bam

#Remove secondary alignments
samtools view --threads 14 -h -F 0x900 "${SLURM_ARRAY_TASK_ID}"_sorted.q30.bam > "${CleanedBamDirectory}"/"${SLURM_ARRAY_TASK_ID}"_sorted.q30.primary_only.bam

#Get rid of intermediate files
#rm *_sorted.q30.bam
#cd "${AlignmentDirectory}"
#rm *

#Use picard to add read group info
cd "${CleanedBamDirectory}"
java -Xms2g -Xmx16g -jar /tools/picard-2.23.9/libs/picard.jar AddOrReplaceReadGroups \
	I="${SLURM_ARRAY_TASK_ID}"_sorted.q30.primary_only.bam \
	O="${SLURM_ARRAY_TASK_ID}"_IDed.bam \
	RGID=ID_"${SLURM_ARRAY_TASK_ID}" \
	RGLB="lib1" \
	RGPL="ILLUMINA" \
	RGPU=Barcode_"${SLURM_ARRAY_TASK_ID}" \
	RGSM=RGSample_"${SLURM_ARRAY_TASK_ID}"

#Use picard to mark duplicates
java -Xms2g -Xmx16g -jar /tools/picard-2.23.9/libs/picard.jar MarkDuplicates I="${SLURM_ARRAY_TASK_ID}"_IDed.bam O="${FinalAlignmentDirectory}"/"${SLURM_ARRAY_TASK_ID}"_0.bam M="${FinalAlignmentDirectory}"/"${SLURM_ARRAY_TASK_ID}"_marked_dup_metrics.txt CREATE_INDEX=true REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE

cd "${FinalAlignmentDirectory}"

# Index the bam file
samtools index -@ 14 "${SLURM_ARRAY_TASK_ID}"_0.bam

#picard MarkDuplicates I="${SLURM_ARRAY_TASK_ID}"_sorted.q30.primary_only.bam O="${SLURM_ARRAY_TASK_ID}"_marked_duplicates.bam M="${SLURM_ARRAY_TASK_ID}"_marked_dup_metrics.txt
#java -Xms2g -Xmx16g -jar /tools/picard-2.23.9/libs/picard.jar BuildBamIndex I=18_marked_duplicates.bam
