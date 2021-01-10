# Introduction

This mini-project attempted to generate PCR probes for [region specific extraction](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-2836-6) of mouse TGFβR1.

Region-specific extraction is a process which attempts to isolate a large region within a pool of genomic DNA (gDNA). The crux of the strategy is the following:

* Generate multiple short PCR probes specific to the region of interest, which are evenly spaced apart.
* Hybridise the probes to that region during the annealing step.
* Use a PCR DNA polymerase to fill the spaces between the hybridised probes (which amounts to one round of PCR amplification). However, unlike regular PCR which uses regular dNTPs, use a mix of regular and biotinylated dNTPs. Hence, during the extension of probes biotinylated dNTPS are  incorporated into extended product.
* Bind the extended probes (containing biotinylated dNTPs) to streptavidin-coated magnetic beads. Separate  both the extended probes and the hybridised region of interest from the gDNA pool using a magnetic stand.

For PCR probe generation, the authors of the linked paper, Dapprich et al., recommend using the software they created called [Antholigo](https://antholigo.chop.edu/), which generates evenly-spaced specific probes, checks for T_m, and also checks the structure. However, at that point in time, Antholigo did not contain the mouse genome, preventing us from generating probes for mouse TGFβR1.

Hence, I attempted to use [OligoMiner](https://github.com/beliveau-lab/OligoMiner), a software originally designed for generating probes for fluorescence *in situ* hybridisation (FISH), to simulate the output of Antholigo and generate probes for our region of interest.

Since OligoMiner was not designed for this use, I had to resort to several optimisation steps in order to get the desired output.

# File guide

* oligominer-new-pipeline-explanation.md contains a step-by-step explanation of the pipeline, the various considerations that went into parameter selection, as well as the troubleshooting steps I took.
* oligominer-code.md contains the actual code used for using command-line tools such as OligoMiner and bowtie2.
* Tgfbr1.fa contains sequence of the region for which we were attempting to create probes.
* checking_chromosome_number.R and spacing-probes.R are two R scripts used during the generation and filtering of probes. Please refer to the rendered R notebooks, checking_chromosome_number.Rmd and spacing-probes.Rmd, corresponding to these scripts to understand their output.
* Tgfbr1.fastq, Tgfbr1_probes.sam,  Tgfbr1_cleaned_probes.bed, Tgfbr1_cleaned_probes_chr.bed, Tgfbr1_cleaned_probes_chr_sC.bed are intermediate output files.
* Tgfbr1_spaced_probes.bed is the final output file containing the probes to be used.
* Both the mouse genome sequence (Mus_musculus.GRCm38.dna.primary_assembly.fa; from Ensembl) as well as the bowtie2 index files corresponding to that genome have not been uploaded, due to Github's space constraints.
