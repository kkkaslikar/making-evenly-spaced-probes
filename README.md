# Introduction

This mini-project attempted to generate PCR probes for [region specific extraction](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-2836-6) of mouse TGFβR1.

Please refer to the  [rendered Github Pages]() for more details.

# File guide

* **oligominer-new-pipeline-explanation.md** contains a step-by-step explanation of the pipeline, the various considerations that went into parameter selection, as well as the troubleshooting steps I took. You may visit the file in this repository directly, or you may look at it as rendered by [GitHub Pages](https://github.com/kkkaslikar/making-evenly-spaced-probes/oligominer-new-pipeline-explanation.html).
* **oligominer-code.md** contains the actual code used the command-line tools such as `OligoMiner` and `bowtie2`. You may visit the file in this repository directly, or you may look at it as rendered by [GitHub Pages](https://github.com/kkkaslikar/making-evenly-spaced-probes/oligominer-code.html).
* **Tgfbr1.fa** contains sequence of the region for which we were attempting to create probes.
* **checking_chromosome_number.R** changes the chromosome nomenclature, updates co-ordinates and spacing of probes, and filters them according to the T~m~. Please refer to this [R notebook](https://kkkaslikar.github.io/making-evenly-spaced-probes/checking_chromosome_number.nb.html) to understand the code output.
* **spacing-probes.Rmd** spaces the probes and enables its visualisation of the probe spacing. Please refer to [this notebook](https://kkkaslikar.github.io/making-evenly-spaced-probes/spacing-probes.nb.html)
* **Tgfbr1.fastq**, **Tgfbr1_probes.sam**,  **Tgfbr1_cleaned_probes.bed**, **Tgfbr1_cleaned_probes.bed**, **Tgfbr1_cleaned_probes_chr.bed**, **Tgfbr1_cleaned_probes_chr_sC.bed** are intermediate output files.
* **Tgfbr1_spaced_probes.bed** is the final output file containing the probes to be used.
* Both the mouse genome sequence (Mus_musculus.GRCm38.dna.primary_assembly.fa; from Ensembl) as well as the bowtie2 index files corresponding to that genome have not been uploaded, due to Github's space constraints.
* Please look to **localtree.txt** to get an idea of the layout of the working directory.
