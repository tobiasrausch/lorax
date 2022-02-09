You can build a [lorax](https://github.com/tobiasrausch/lorax) singularity container (SIF file) using

`sudo singularity build lorax.sif lorax.def`

Once you have built the container you can run it using

`singularity exec lorax.sif lorax --help`
