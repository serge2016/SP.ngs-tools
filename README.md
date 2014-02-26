ngs-tools
========

NGS data diagnose & manipulations tools

## FastQ
* 'fastq_detect.pl': analyze the N first reads from a fastq file (10000) and match with existing FastQ formats. report matching results. 

## SAM
* 'uniq_mappings.pl': filter uniquely mapping reads and multiple mapping reads from a SAM stream to separate BAM files. The Input should be name-sorted (queryname).
