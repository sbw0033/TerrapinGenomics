# TerrapinGenomics
Project with Matthew Wolak and Tonia Schwartz at Auburn University

This project starts with pre-processed raw reads. The files containing these reads are re-named and merged with all of the other in the script "TerrapinGenomicsPrepWork"

The major steps in moving from these raw reads to variant sites are as follows
1) Align the raw reads to a reference genome using the program BWA
2) Processing the raw alignments using Picard and Samtools
3) Calling variants using GATK's HaplotypeCaller



###1) RAW READ ALIGNMENT
When using BWA with GATK in mind as the eventual end goal, we need to be sure to include information about the sequencing platform, read group, and sample ID when running BWA.

BWA mem manual (http://bio-bwa.sourceforge.net/bwa.shtml)

###2) RAW READ PROCESSING

Once we have all the alignments we need to do the following things
a) Convert .sam to .bam and sort it 
b) Ditch alignments where quality score is below 30
c) Remove all secondary alignments
d) [Add read group information so GATK doesn't freak out] (https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-)
e) [Mark duplicates- (Also so GATK doesn't freak out)]. (https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-)

###3) CALL VARIANTS

We are going to use GATK's HaplotypeCaller to accomplish this
The BWAarray.sh script


