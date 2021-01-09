---
title: OligoMiner Pipeline
date: 24/01/2020
---

# Old pipeline

Directory guide

## `/mnt/f/OneDrive/GTK-Lab/Region-specific-extraction`

Contains the mouse TGF-beta R1 sequence fasta file.

## `/mnt/f/Bioinformatics_Data/Ensembl`

Contains the mouse **unmasked and masked genomes**.

## `/mnt/f/Bioinformatics_Data/Indexes/Mus_musculus_GRCm38_unmasked_bt2_index`

Contains the genome index built from the Mus musculus **unmasked genome**.


## Building the index

```
bowtie2-build -f /mnt/f/OneDrive/GTK-Lab/Region-specific-extraction/Mus_musculus.GRCm38.dna.primary_assembly.fa MM.GRCm38.unmasked

```

### Moving the index files to the directory for indexes

```
mv ./*.bt2  ../Indexes/Mus_musculus_GRCm38_unmasked_bt2_index/

```

## Generating probes for mouse TGF-beta R1

### Changing directory

```
cd /mnt/f/OneDrive/GTK-Lab/Region-specific-extraction

```

### Generating probes using the `blockParse.py` script.

**Note:** `blockParse.py` has already been added to `PATH` in this case; hence, it can be called directly.

Based on the supplementary material for Dapprich et al. (2016) [The next generation of target capture technologies], the probe length they used was uniformly 25 bp. Hence,a good minimum probe length was deemed to be 20 and the maximum probe length was deemed to be 25. In addition, since a good range of Tm for PCR primers is 55 ℃ to 60 ℃, those were used as the minimum and maximum Tm values respectively.  Since OligoMiner was made specfically for FISH, the default `--salt` concentration of 390 mM does not apply, nor does the default `--formamide` percentage of 50%. Salt concentration in PCR buffers is generally 50 mM of KCl, and as far as I know, no formamide is present in PCR buffers. Hence, the parameters were set accordingly in the call to `blockParse.py`.

Different parameters of spacing were tried, in order to generate probes mapping to the region of interest, with appropriate spacing. 

```
blockParse.py --file Tgfbr1.fa --minLength 19 --maxLength 25 --min_Tm 55 --max_Tm 60 --salt 50 --formamide 0 --Spacing 500
```



Output:




```
0 of 61705
35 candidate probes identified in 59.97 kb yielding 0.58 candidates/kb
```



## Aligning generated probes to create fastq file


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

## Removing non-specific probes

```
outputClean.py --file Tgfbr1_probes.sam --unique --salt 50 --formamide 0 --Report --Meta --output Tgfbr1_cleaned_probes
```

Output:
```
outputClean identified 23 of 35 / 65.7143% candidate probes as unique
```

## Checking probe structure

```
structureCheck.py --file Tgfbr1_cleaned_probes.bed --formamide 0 --salt 50 --hybTemp 53 --Report
```

Output:

```
structureCheck predicted that 23 of 23 / 100.0000% candidate probes are predicted to have a linear structure with p>0.1000 at 53C in 50 mM Na+ and 0% formamide
```


# New pipeline
## Strategy of new pipeline
	
* [x] Produce probes **first**.
* [x] Filter for multi-mappers
* [x] Filter for those lying within our region of interest
* [x] Filter for Tm
* [x] Filter for secondary structure
* [x] Filter for spacing 

	
## Directory guide

### `/mnt/f/OneDrive/GTK-Lab/Region-specific-extraction`

Contains the mouse TGF-beta R1 sequence fasta file.

### `/mnt/f/Bioinformatics_Data/Ensembl`

Contains the mouse **unmasked and masked genomes**.

### `/mnt/f/Bioinformatics_Data/Indexes/Mus_musculus_GRCm38_unmasked_bt2_index`

Contains the genome index built from the Mus musculus **unmasked genome**.


## Building the index

```
bowtie2-build -f /mnt/f/OneDrive/GTK-Lab/Region-specific-extraction/Mus_musculus.GRCm38.dna.primary_assembly.fa MM.GRCm38.unmasked

```

### Moving the index files to the directory for indexes

```
mv ./*.bt2  ../Indexes/Mus_musculus_GRCm38_unmasked_bt2_index/

```

## Generating probes for mouse TGF-beta R1

### Changing directory

```
cd /mnt/f/OneDrive/GTK-Lab/Region-specific-extraction

```

### Running the `blockParse.py` script.

**Note:** `blockParse.py` has already been added to `PATH` in this case; hence, it can be called directly.

For the new pipeline, it was decided to **not** to space the probes first. The probes were simply generated to have a Tm of between 55 ℃ to 60 ℃, as well as a length of between 19 to 25.

```
blockParse.py --file Tgfbr1.fa --minLength 19 --maxLength 25 --min_Tm 55 --max_Tm 60 --salt 50 --formamide 0 
```

Output:

```
0 of 61705
1261 candidate probes identified in 61.55 kb yielding 20.49 candidates/kb
```
	
## Aligning generated probes to create fastq file


```
bowtie2 -x /mnt/f/Bioinformatics_Data/Indexes/Mus_musculus_GRCm38_unmasked_bt2_index/MM.GRCm38.unmasked -U Tgfbr1.fastq --no-hd -t -k 100 --very-sensitive-local -S Tgfbr1_probes.sam

```


Output:
```
1261 reads; of these:
  1261 (100.00%) were unpaired; of these:
      0 (0.00%) aligned 0 times
	      854 (67.72%) aligned exactly 1 time
		      407 (32.28%) aligned >1 times
			  100.00% overall alignment rate
```

## Removing non-specific probes

```
outputClean.py --file Tgfbr1_probes.sam --unique --salt 50 --formamide 0 --output Tgfbr1_cleaned_probes
```

Output:
```
outputClean identified 854 of 1261 / 67.7240% candidate probes as unique
```

## Checking the chromosome location of uniquely-mapped probes

The following R script was used to check which chromosome did the uniquely mapped probes map to, and for those which were correctly mapped to chromosome 4, the value in the first column was changed to `chr4`, and the resulting data frame was saved to a new bed file to allow convenient visualisation by UCSC genome browser or IGV.


## Checking probe structure

```
structureCheck.py --file Tgfbr1_cleaned_probes.bed --formamide 0 --salt 50 --hybTemp 53 --Report
```

Output:

```
structureCheck predicted that 822 of 831 / 98.9170% candidate probes are predicted to have a linear structure with p>0.1000 at 51C in 50 mM Na+ and 0% formamide
```

## Spacing probes

The following R script was used to space the probes.

# Notes from Greg

1. Run for probe sizes until 25 bp and get a size distribution in R.
2. Check for recommended probe spacing in original paper and emulate that in script.
3. Double-check the probe melting temperatures from the original paper with IDT OligoAnalyzer. Check how far off is OligoAnalyzer compared to Antholigo and OligoMiner.
3. Use UCSC to check the probe location.
