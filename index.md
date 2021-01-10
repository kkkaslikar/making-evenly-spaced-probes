# Introduction

This mini-project attempted to generate PCR probes for [region specific extraction](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-2836-6) of mouse TGFβR1.

Region-specific extraction is a process which attempts to isolate a large region within a pool of genomic DNA (gDNA). The crux of the strategy is the following:

* Generate multiple short PCR probes specific to the region of interest, which are evenly spaced apart.
* Hybridise the probes to that region during the annealing step.
* Use a PCR DNA polymerase to fill the spaces between the hybridised probes (which amounts to one round of PCR amplification). However, unlike regular PCR which uses regular dNTPs, use a mix of regular and biotinylated dNTPs. Hence, during the extension of probes biotinylated dNTPS are  incorporated into extended product.
* Bind the extended probes (containing biotinylated dNTPs) to streptavidin-coated magnetic beads. Separate  both the extended probes and the hybridised region of interest from the gDNA pool using a magnetic stand.

For PCR probe generation, the authors of the linked paper, Dapprich et al., recommend using the software they created called [Antholigo](https://antholigo.chop.edu/), which generates evenly-spaced specific probes, checks for T~m~, and also checks the structure. However, at that point in time, Antholigo did not contain the mouse genome, preventing us from generating probes for mouse TGFβR1.

Hence, I attempted to use [OligoMiner](https://github.com/beliveau-lab/OligoMiner), a software originally designed for generating probes for fluorescence *in situ* hybridisation (FISH), to simulate the output of Antholigo and generate probes for our region of interest.

Since OligoMiner was not designed for this use, I had to resort to several optimisation steps in order to get the desired output.

For a step-by-step details of the pipeline as well as the rationale behind the selected parameters, look at [this webpage](https://kkkaslikar.github.io/making-evenly-spaced-probes/oligominer-new-pipeline-explanation.html).

For the actual code used for the command-line tools such as `OligoMiner` and `bowtie2`, look at [this link](https://kkkaslikar.github.io/making-evenly-spaced-probes/pipeline.html).

For the R code which changes the chromosome nomenclature from Ensembl to UCSC, modifies the genomic co-ordinates, and filters for T~m~, follow [this link](https://kkkaslikar.github.io/making-evenly-spaced-probes/checking-chromosome-number.nb.html).

For the R code which spaces the probes and visualises the spacing, follow [this link](https://kkkaslikar.github.io/making-evenly-spaced-probes/spacing-probes.nb.html).




