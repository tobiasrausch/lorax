# Lorax: A long-read analysis toolbox for cancer genomics

In cancer genomics, long-read de novo assembly approaches may not be applicable because of tumor heterogeneity, normal cell contamination and aneuploid chromosomes. Likewise, it may not be feasible to generate sufficiently high coverage for each derivative, potentially sub-clonal, chromosome. Lorax is a targeted approach to reveal complex cancer genomic structures such as telomere fusions, templated insertion threads or chromothripsis rearrangements. Lorax is not a long-read SV caller, this functionality is implemented in [delly](https://github.com/dellytools/delly). Lorax requires matched tumor-normal data sequenced using long-reads.

# Telomere fusions

Telomere fusions are a hallmark of cancer and these can be identified with lorax using

`lorax telomere -g hg38.fa -o tumor.bed.gz -m control.bam tumor.bam`

# Templated insertion threads

Templated insertions threads can be identified using

`lorax tithreads -g hg38.fa -o tithreads.bed -m control.bam tumor.bam`




