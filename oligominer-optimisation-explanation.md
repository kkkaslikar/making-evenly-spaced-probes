---
title: Obtaining RSE probes for mouse Tgfbr1
---

# General considerations during probe-generation

* Due to the mouse genome not being present on [Antholigo](https://antholigo.chop.edu/) for now, we are attempting to use [OligoMiner](https://github.com/beliveau-lab/OligoMiner) in its place.
* However, since OligoMiner was primarily designed for probe generation for *in situ* hybridisation applications, some of its default parameters and behaviour have to be tweaked for it to emulate Antholigo's output (to the degree to which we can manage).

---

Some of the considerations I am trying to pay attention to during probe generation include:

* Number - There must be sufficient number of probes to cover the region of interest.

* Melting temperature (Tm) - I am attempting to get probes that are not too far apart in Tm, since that may cause problems during PCR.

* Spacing - The probes cannot be spaced too closely together, or too far apart.
  
Getting probes which satisfy these criteria will require a degree of trial and error and optimisation.


# OligoMiner parameter selection


* Since OligoMiner was designed primarily for fluorescence *in situ* hybridisation applications, its default parameters are set for the production of FISH probes.
* For example, OligoMiner assumes an Na^+^ salt concentration of 390 mM and a formamide percentage of 50%, which are the norm for FISH. Both will affect melting temperature (Tm).
* However, as far as I can tell, PCR buffers do not contain formamide and contain only 50 mM of K^+^. K^+^ being a monovalent salt like Na^+^, we can assume that the effects of [Na^+^] and [K^+^] on Tm should be similar.
* In addition, PCR buffers may also contain Mg^2+^, as a cofactor for the polymerase. However, it is not present at a very high concentration, so I do not think it will affect Tm.

# OligoMiner pipeline overview

1. Using the `blockParse.py` to generate probes from the input region (Tgfbr1) to generate a `.fastq` file.
2. Aligning generated probes using `bowtie2` to create a `.sam` file (genome index has to be built first and it takes a long time).

   This will map the probes to the entire mouse genome. Some probes might map to more than just our region of interest. These are not desirable to us.
3. Removing the multi-mapped probes using the `outputClean.py` script. This produces a `BED` file with the genomic co-ordinates (start and end positions) of each of the sequences. 
4. Using the `structureCheck.py` script to check for possible secondary structures within our probes. Probes are filtered further.

# Probe length
  * Dapprich et al. used a probe length of 25 bp in their paper.
  * I correspondingly tried to keep a probe length of between 19 to 25, so as not to be too restrictive, and allow some flexibility, since we want to generate sufficient probes which also account for the spacing and Tm considerations.
   

# Probe spacing

* Dapprich et al. used a probe spacing of 6 to 10 kb in the original paper.
* I used a probe spacing of 6 kb, in order to have a better chance of pulling down our region of interest.
* However, the number of probes that remains after filtering for multimappers, T_m and structure is too low.
* Hence, I resorted to using a different order of operations.

## New probe generation strategy

* Produce probes **first**.
* Filter for multi-mappers
* Filter for those lying within our region of interest
* Filter for secondary structure
* Filter for Tm
* Filter for spacing




 

