[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/lorax/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/lorax/badges/downloads.svg)](https://anaconda.org/bioconda/lorax)
[![C/C++ CI](https://github.com/tobiasrausch/lorax/workflows/C/C++%20CI/badge.svg)](https://github.com/tobiasrausch/lorax/actions)
[![Docker CI](https://github.com/tobiasrausch/lorax/workflows/Docker%20CI/badge.svg)](https://hub.docker.com/r/trausch/lorax/)
[![GitHub license](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/tobiasrausch/lorax/blob/master/LICENSE)
[![GitHub Releases](https://img.shields.io/github/release/tobiasrausch/lorax.svg)](https://github.com/tobiasrausch/lorax/releases)

# Lorax: A long-read analysis toolbox for cancer genomics

In cancer genomics, long-read de novo assembly approaches may not be applicable because of tumor heterogeneity, normal cell contamination and aneuploid chromosomes. Generating sufficiently high coverage for each derivative, potentially sub-clonal, chromosome is not feasible. Lorax is a targeted approach to reveal complex cancer genomic structures such as telomere fusions, templated insertions or chromothripsis rearrangements. Lorax is NOT a long-read SV caller, this functionality is implemented in [delly](https://github.com/dellytools/delly). Lorax requires matched tumor-normal data sequenced using long-reads.

## Installing lorax

Lorax is available as a [statically linked binary](https://github.com/tobiasrausch/lorax/releases/), a [singularity container (SIF file)](https://github.com/tobiasrausch/lorax/releases/) or as a [docker container](https://hub.docker.com/r/trausch/lorax/). You can also build lorax from source using a recursive clone and make. Lorax depends on [HTSlib](https://github.com/samtools/htslib) and [Boost](https://www.boost.org/).

`git clone --recursive https://github.com/tobiasrausch/lorax.git`

`cd lorax/`

`make all`

## Templated insertion threads

Templated insertions threads can be identified using

`lorax tithreads -g hg38.fa -o tithreads.bed -m control.bam tumor.bam`

The `out.bed` file specifies nodes (templated insertion source sequences) and edges (templated insertion adjacencies) of a graph that can be plotted using dot.

`cut -f 4,9 out.bed | sed -e '1s/^/graph {\n/' | sed -e '$a}' > out.dot`

`dot -Tpdf out.dot -o out.pdf`

The `out.reads` file lists unique assignments of reads to templated insertion source sequences. To extract the FASTA sequences for all these reads use the `lorax extract` subcommand (below) with the `-a` option.

`tail -n +2 out.reads | cut -f 1 | sort | uniq > reads.lst`

`lorax extract -a -g hg38.fa -r reads.lst tumor.bam`

## Telomere repeats associated with complex rearrangements

Telomere-associated SVs can be identified with lorax using

`lorax telomere -g hg38.fa -o tumor.bed.gz -m control.bam tumor.bam`

The output file clusters reads into distinct telomere junctions and you can trace the alignment matches for every read.

## Read selection for targeted assembly of amplicons

Given a list of amplicon regions and a phased VCF file, lorax can be used to extract amplicon reads for targeted assembly approaches.

`lorax amplicon -g hg38.fa -s sample -v phased.bcf -b amplicons.bed tumor.bam`

The amplicon subcommand outputs the selected reads (as a hash list `out.reads`) and a diagnostic table (`out.bed`) with amplicon regions and their support by split-reads. Ideally, all amplicon regions are connected and belong to one connected component (one cluster of amplicons). This amplicon graph can be plotted using dot.

`cut -f 4,11 out.bed | sed -e '1s/^/graph {\n/' | sed -e '$a}' > out.dot`

`dot -Tpdf out.dot -o out.pdf`

To extract the FASTA sequences for all reads use the `lorax extract` subcommand (below) with the `-a` option.

## Extracting pairwise matches and FASTA sequences of reads

To get FASTA sequences and pairwise read to genome matches for a list of reads (`list.reads`) use

`lorax extract -g hg38.fa -r list.reads tumor.bam`

If the read list contains hashes instead of read names as from the `lorax amplicon` subcommand then please use the `-a` command-line option.

`lorax extract -a -g hg38.fa -r list.reads tumor.bam`

## Citation

Tobias Rausch, Rene Snajder, Adrien Leger, Milena Simovic, Oliver Stegle, Ewan Birney, Marc Jan Bonder, Aurelie Ernst, Jan O. Korbel.          
Long-read sequencing of diagnosis and post-therapy medulloblastoma reveals complex rearrangement patterns and epigenetic signatures.            
[bioRxiv 2022.02.20.480758](https://doi.org/10.1101/2022.02.20.480758)

License
-------
Lorax is distributed under the BSD 3-Clause license. Consult the accompanying [LICENSE](https://github.com/tobiasrausch/lorax/blob/master/LICENSE) file for more details.
