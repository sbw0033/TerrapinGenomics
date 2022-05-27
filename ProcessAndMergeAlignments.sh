#!/bin/bash -l   
#SBATCH --job-name=ProcessAndMergeAlignments           # job name
#SBATCH --mem=9GB            		 # memory
#SBATCH --ntasks=10                  # number of tasks across all nodes
#SBATCH --partition=general          # name of partition to submit job
#SBATCH --time=12:00:00              # Run time (D-HH:MM:SS)
#SBATCH --output=SAMtoolsStuff.out          # Output file. %j is replaced with job ID
#SBATCH --error=SAMtoolsStuff.err         # Error file. %j is replaced with job ID
#SBATCH --mail-type=ALL              # will send email for begin,end,fail
#SBATCH --mail-user=sbw0033@auburn.edu

#ProcessAndMergeAlignments
cd /home/sbw0033/TerrapinGenomics/Results/BWA_out

module load samtools/1.11

#First, rewrite all of them as a bam
for sam in *.sam; do samtools view --threads 10 -S -b "${sam}" > /scratch/sbw0033/AlignmentProcessing/"${sam}".bam; done

cd /scratch/sbw0033/AlignmentProcessing

#Then, sort all of the bams
for bam in *.bam; do samtools sort "${bam}" -o "${bam}"_sorted.bam --threads 10; done

#Do it all in one?
#for sam in *.sam; do samtools view -u "${sam}" --threads 10 | samtools sort --threads 10 -@ 4 -o "${sam}"_sorted.bam; done

#Now we can merge all of our bam files
#samtools merge all.bam *.bam

#First, keep alignments with a map quality score greater than or equal to 30
#samtools view -q 30 -b MergedSamples.sam > Chr7_aligned_reads.q30.bam --threads 10

#Second, remove secondary alignments
#samtools view -h -F 0x900 Chr7_aligned_reads.q30.bam > Chr7_Only_Primary_Alignments.bam --threads 10

#Finally, keep reads where both reads in the read pair mapped.
#samtools view -q 30 -b Chr7_aligned_reads.q30.bam > Chr7_aligned_reads.q30.bam --threads 10

#Now sort the alignment
#samtools sort Chr7_aligned_reads.q30.bam -o Chr7_Sorted_Alignment.bam --threads 10
