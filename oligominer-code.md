
# Building the index

Using bowtie2 to build the index.

```
bowtie2-build -f /mnt/f/OneDrive/GTK-Lab/Region-specific-extraction/Mus_musculus.GRCm38.dna.primary_assembly.fa MM.GRCm38.unmasked

```

## Moving the index files to the directory for indexes

```
mkdir ../Indexes/Mus_musculus_GRCm38_unmasked_bt2_index/

mv ./*.bt2  ../Indexes/Mus_musculus_GRCm38_unmasked_bt2_index/

```

# Generating probes for mouse TGF-beta R1

## Changing directory

```
cd /mnt/f/OneDrive/GTK-Lab/Region-specific-extraction

```

## Generating probes using the `blockParse.py` script.

**Note:** `blockParse.py` has already been added to `PATH` in this case; hence, it can be called directly.

Based on the supplementary material for Dapprich et al. (2016) [The next generation of target capture technologies], the probe length they used was uniformly 25 bp. Hence,a good minimum probe length was deemed to be 20 and the maximum probe length was deemed to be 25. In addition, since a good range of Tm for PCR primers is 55 ℃ to 60 ℃, those were used as the minimum and maximum Tm values respectively.  Since OligoMiner was made specfically for FISH, the default `--salt` concentration of 390 mM does not apply, nor does the default `--formamide` percentage of 50%. Salt concentration in PCR buffers is generally 50 mM of KCl, and as far as I know, no formamide is present in PCR buffers. Hence, the parameters were set accordingly in the call to `blockParse.py`.

Different parameters of spacing were tried, in order to generate probes mapping to the region of interest, with appropriate spacing. Below is an example of the code used to generate probes with 500 bp spacing.

Input:

```
blockParse.py --file Tgfbr1.fa --minLength 19 --maxLength 25 --min_Tm 55 --max_Tm 60 \
 --salt 50 --formamide 0 --Spacing 500
```



Output:

```
0 of 61705
35 candidate probes identified in 59.97 kb yielding 0.58 candidates/kb
```



# Aligning generated probes to create fastq file

Input:

```
bowtie2 -x /mnt/f/Bioinformatics_Data/Indexes/Mus_musculus_GRCm38_unmasked_bt2_index/MM.GRCm38.unmasked -U Tgfbr1.fastq --no-hd -t -k 100 --very-sensitive-local -S Tgfbr1_probes.sam

```


Output:

```
35 reads; of these:
  35 (100.00%) were unpaired; of these:
    0 (0.00%) aligned 0 times
    23 (65.71%) aligned exactly 1 time
    12 (34.29%) aligned >1 times
100.00% overall alignment rate
```

# Removing non-specific probes


Input:

```
outputClean.py --file Tgfbr1_probes.sam --unique --salt 50 --formamide 0 --Report --Meta --output Tgfbr1_cleaned_probes
```

Output:

```
outputClean identified 23 of 35 / 65.7143% candidate probes as unique
```

# Checking probe structure

Input:

```
structureCheck.py --file Tgfbr1_cleaned_probes.bed --formamide 0 --salt 50 --hybTemp 53 --Report
```

Output:

```
structureCheck predicted that 23 of 23 / 100.0000% candidate probes are predicted to have a linear structure with p>0.1000 at 53C in 50 mM Na+ and 0% formamide
```

# Checking probe locations

## Changing probe location nomenclature in `.bed` file

Since the genome I downloaded was from Ensembl, the nomenclature used to annotate the individual chromosome sequences  was in the form of `>1`, `>2`, etc. in the original genome fasta file. Hence, all the files processed by bowtie2, and subsequently OligoMiner, using that genome file used a nomenclature which did not specify the chromosome location in the bed file.

For example: here are the first three lines of `Tgfbr1_cleaned_probes_sC.bed`, which was produced after structure checking:

```
chrom	502	520	CTAAATGCGCCCCGGCCTG	58.31
chrom	2038	2058	AATCATCCAGTGCTGGCGACT	56.20
chrom	2559	2580	TTTAAGCTATGGCCTGCAGGGC	57.28
```

This file cannot be loaded on to the UCSC genome browser (see next heading for the explanation of why we do this), since it does not specify which genome the files came from. What we want to load onto the UCSC browser should look a little like this:


```
chr4	502	520	CTAAATGCGCCCCGGCCTG	58.31
chr4	2038	2058	AATCATCCAGTGCTGGCGACT	56.20
chr4	2559	2580	TTTAAGCTATGGCCTGCAGGGC	57.28
```

The chromosome-mapping information can be retrieved from the original `.sam` file, which is a tab-separated text file that looks like this:



```
chrom:502-520	0	4	47353723	255		19M		*	0	0	CTAAATGCGCCCCGGCCTG	~~~~~~~~~~~~~~~~~~~	AS:i:38	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:19	YT:Z:UU
chrom:2038-2058	0	4	47355259	255	21M	*	0	0	AATCATCCAGTGCTGGCGACT	~~~~~~~~~~~~~~~~~~~~~	AS:i:42	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:21	YT:Z:UU
chrom:2559-2580	0	4	47355780	255	22M	*	0	0	TTTAAGCTATGGCCTGCAGGGC	~~~~~~~~~~~~~~~~~~~~~~	AS:i:44	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:22	YT:Z:UU

```
More information about the fields in the SAM file format can be found here: 

[samtools.github.io/hts-specs/SAMv1.pdf](samtools.github.io/hts-specs/SAMv1.pdf)

As specified by the above pdf file, column 3 in the SAM format is a reference name, or in this case, the name of the chromosome where the probe was aligned to by `bowtie2`. 

To retrieve all the processed probes aligned to chromosome 4 from the original SAM file, the following method was used:

1. Loading all the probes from the original `.sam` file (`Tgfbr1_probes.sam`) (which might contain potential multimappers)
2. Loading  processed probes from the `.bed` file (`Tgfbr1_cleaned_probes_sC.bed`) produced after structure checking
3. Filtering the `.sam`  file probes to limit  them to the ones in the `.bed` file.
4. Verifying which of the `.sam` probes mapped to chromosome 4.
5. Changing the nomenclature of chromosome-4 mapped probes in the `.bed` file from `chrom` to `chr4`, and rewriting that data to a new `.bed` file (`Tgfbr1_cleaned_probes_sC_chr.bed`)

The above steps were performed using the following R script:

```R
library(tidyverse)

sam <- read_tsv("Tgfbr1_probes.sam", col_names = FALSE) # file containing potential multimappers, with chromosome location data
 
bed <- read_tsv("Tgfbr1_cleaned_probes_sC.bed", col_names = FALSE) # file cleaned of all multimappers and checked for structure; contains no chromosome data

sam[which(sam$X10 %in% bed$X4), ] %>% select(X3) %>% distinct() # finding the chromosome number of the probes in bed

bed$X1 <- str_replace(bed$X1, pattern = "chrom", replacement = "chr4") # adding the relevant chromosome name

write_tsv(x = bed, path = "Tgfbr1_cleaned_probes_sC_chr.bed", col_names = FALSE)
```

