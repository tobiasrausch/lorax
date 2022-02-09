[![C/C++ CI](https://github.com/tobiasrausch/lorax/workflows/C/C++%20CI/badge.svg)](https://github.com/tobiasrausch/lorax/actions)
[![Docker CI](https://github.com/tobiasrausch/lorax/workflows/Docker%20CI/badge.svg)](https://hub.docker.com/r/trausch/lorax/)
[![GitHub license](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/tobiasrausch/lorax/blob/master/LICENSE)
[![GitHub Releases](https://img.shields.io/github/release/tobiasrausch/lorax.svg)](https://github.com/tobiasrausch/lorax/releases)

# Lorax: A long-read analysis toolbox for cancer genomics

In cancer genomics, long-read de novo assembly approaches may not be applicable because of tumor heterogeneity, normal cell contamination and aneuploid chromosomes. Generating sufficiently high coverage for each derivative, potentially sub-clonal, chromosome is not feasible. Lorax is a targeted approach to reveal complex cancer genomic structures such as telomere fusions, templated insertions or chromothripsis rearrangements. Lorax is NOT a long-read SV caller, this functionality is implemented in [delly](https://github.com/dellytools/delly). Lorax requires matched tumor-normal data sequenced using long-reads.

## Templated insertion threads

Templated insertions threads can be identified using

`lorax tithreads -g hg38.fa -o tithreads.bed -m control.bam tumor.bam`

The output file specifies nodes (templated insertion source sequences) and edges (templated insertion adjacencies) of a graph that can be plotted using dot.

`cut -f 4,9 out.bed | sed -e '1s/^/graph {\n/' | sed -e '$a}' > out.dot`

`dot -Tpdf out.dot -o out.pdf`

## Telomere repeats associated with complex rearrangements

Telomere-associated SVs can be identified with lorax using

`lorax telomere -g hg38.fa -o tumor.bed.gz -m control.bam tumor.bam`

The output file clusters reads into distinct telomere junctions and you can trace the alignment matches for every read.

## Read selection for target assembly of amplicons

Given a list of amplicon regions and a phased VCF file, lorax can be used to extract amplicon reads for targeted assembly approaches.

`lorax amplicon -g hg38.fa -s sample -v phased.bcf -b amplicons.bed tumor.bam`

The amplicon subcommand outputs the selected reads and a diagnostic table with amplicon regions and their support by split-reads. Ideally, all amplicon regions are connected and belong to one connected component (one cluster of amplicons). This amplicon graph can be plotted using dot.

`cut -f 4,8 tumor.bed | sed -e '1s/^/graph {\n/' | sed -e '$a}' > out.dot`

`dot -Tpdf out.dot -o out.pdf`


## Extracting pairwise matches and FASTA sequences of reads

To get pairwise read to genome matches for a list of reads (`reads.txt`) use

`lorax extract -g hg38.fa -r reads.txt tumor.bam`
