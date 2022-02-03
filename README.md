# Lorax: A long-read analysis toolbox for cancer genomics

In cancer genomics, long-read de novo assembly approaches may not be applicable because of tumor heterogeneity, normal cell contamination and aneuploid chromosomes. Generating sufficiently high coverage for each derivative, potentially sub-clonal, chromosome is not feasible. Lorax is a targeted approach to reveal complex cancer genomic structures such as telomere fusions, templated insertions or chromothripsis rearrangements. Lorax is NOT a long-read SV caller, this functionality is implemented in [delly](https://github.com/dellytools/delly). Lorax requires matched tumor-normal data sequenced using long-reads.

# Templated insertion threads

Templated insertions threads can be identified using

`lorax tithreads -g hg38.fa -o tithreads.bed -m control.bam tumor.bam`

# Telomere repeats associated with complex rearrangements

Telomere-associated SVs can be identified with lorax using

`lorax telomere -g hg38.fa -o tumor.bed.gz -m control.bam tumor.bam`

# Read selection for target assembly of amplicons

Given a list of amplicon regions and a phased VCF file, lorax can be used to extract amplicon reads for targeted assembly approaches.

`lorax amplicon -g hg38.fa -s sample -v phased.bcf -b amplicons.bed tumor.bam`

# Extracting pairwise matches and FASTA sequences of reads

To get pairwise read to genome matches for a list of reads (`reads.txt`) use

`lorax extract -g hg38.fa -r reads.txt tumor.bam`
