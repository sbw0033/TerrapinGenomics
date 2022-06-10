# TerrapinGenomics
Project with Matthew Wolak and Tonia Schwartz at Auburn University

This project starts with pre-processed raw reads. The files containing these reads are re-named and merged with all of the other in the script "TerrapinGenomicsPrepWork"

The major steps in moving from these raw reads to variant sites are as follows
1) Align the raw reads to a reference genome using the program BWA
2) Processing the raw alignments using Picard and Samtools
3) Calling variants using GATK's HaplotypeCaller



1) RAW READ ALIGNMENT
When using BWA with GATK in mind as the eventual end goal, we need to be sure to include information about the sequencing platform, read group, and sample ID when running BWA.

BWA mem manual (http://bio-bwa.sourceforge.net/bwa.shtml)


