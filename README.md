# Introduction

This mini-project attempted to generate PCR probes for [region specific extraction](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-2836-6) of mouse TGFβR1.

Please refer to the  [rendered Github Pages]() for more details.

# File guide

* **oligominer-new-pipeline-explanation.md** contains a step-by-step explanation of the pipeline, the various considerations that went into parameter selection, as well as the troubleshooting steps I took. You may visit the file in this repository directly, or you may look at it as rendered by [GitHub Pages](https://github.com/kkkaslikar/making-evenly-spaced-probes/oligominer-new-pipeline-explanation.md).
* **oligominer-code.md** contains the actual code used the command-line tools such as `OligoMiner` and `bowtie2`. You may visit the file in this repository directly, or you may look at it as rendered by [GitHub Pages](https://github.com/kkkaslikar/making-evenly-spaced-probes/oligominer-code.md).
* **Tgfbr1.fa** contains sequence of the region for which we were attempting to create probes.
* **checking_chromosome_number.R** and **spacing-probes.R** are two R scripts used during the amendment, filtering and spacing of probes. Please refer to the GitHub Pages-rendered R notebooks, **checking_chromosome_number.Rmd** and **spacing-probes.Rmd**, corresponding to these scripts, to understand their code.
* **Tgfbr1.fastq**, **Tgfbr1_probes.sam**,  **Tgfbr1_cleaned_probes.bed**, **Tgfbr1_cleaned_probes.bed**, **Tgfbr1_cleaned_probes_chr.bed**, **Tgfbr1_cleaned_probes_chr_sC.bed** are intermediate output files.
* **Tgfbr1_spaced_probes.bed** is the final output file containing the probes to be used.
* Both the mouse genome sequence (Mus_musculus.GRCm38.dna.primary_assembly.fa; from Ensembl) as well as the bowtie2 index files corresponding to that genome have not been uploaded, due to Github's space constraints.
* Please look to **localtree.txt** to get an idea of the layout of the working directory.
