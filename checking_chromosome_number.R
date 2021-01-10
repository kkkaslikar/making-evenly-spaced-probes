library(tidyverse)


# Checking number of the chromosomes that reads aligned to ----------------

## nfields <- system("awk -F'\t' '{print NF}' Tgfbr1_probes.sam", intern = TRUE)

## nfields

## max(nfields)

sam <- read_tsv("Tgfbr1_probes.sam", col_names = FALSE) # file containing potential multimappers, with chromosome location data
 
bed <- read_tsv("Tgfbr1_cleaned_probes.bed", col_names = FALSE) # file cleaned of all multimappers, with no chromosome data

bed

colnames(bed) <- c("chromosome", "start", "end", "probe", "Tm")

sam[which(sam$X10 %in% bed$probe), ] %>% select(X3) %>% distinct() # finding the chromosome ID of the uniquely-mapping probes in bed



# Changing nomenclature of sequences in bed file to the correct chromosome ---------------------------------------

bed$chromosome <- str_replace(bed$chromosome, pattern = "chrom", replacement = "chr4") # adding the relevant chromosome name

bed

# Modifying co-ordinates to genomic co-ordinates for chromosome 4 ---------

tgfbr1_bounds <- c("start" = 47353222, "end" = 47414926) # co-ordinates for Tgfbr1 on chromosome 4

bed_align <- bed   # bed with new co-ordinates

bed_align$start <- (bed_align$start - 1 + tgfbr1_bounds["start"]) # adjusting the start positions

bed_align$end <- (bed_align$end - 1 + tgfbr1_bounds["start"]) # adjusting the end positions

bed_align

# Filtering for consistent Tm ---------------------------------------------

hist(bed_align$Tm)

summary(bed_align$Tm)

bed_align <- bed_align %>% filter(Tm <= 57)

# Saving probes modified for chromosome, co-ordinates and Tm --------

write_tsv(x = bed_align, path = "Tgfbr1_cleaned_probes_chr.bed", col_names = FALSE)

# tgfbr1_range <- GRanges(seqnames = c("chr4"), ranges = IRanges (start = tgfbr1_bounds["start"], end = tgfbr1_bounds["end"]))

# probe_ranges <- IRanges(start = bed$start, end = bed$end)

# x <- bed_align[1, ]

# for(i in 1:nrow(bed_align)){
#   if(bed_align$start[i+1] - bed_align$end[i] > 500){
#     print(bed_align[i, ])
#   }
# }


##########################################################################


# Checking whether co-ordinates from outputClean are w.r.t. Tgfbr1 --------

## library(tidyverse)
## library(GenomicRanges)
## library(Biostrings)

## sequence <- readDNAStringSet("Tgfbr1.fa")

## sequence

## sequence <- sequence[[1]]

## sequence

## probes <- DNAStringSet(x = bed_align$probe)

## al <- lapply(seq_along(probes), FUN = function(i){
  ## matchPattern(pattern = probes[[i]], subject = sequence)
## })

## al

## al <- map(.x = al, .f = ~as(object = .x, Class = "IRanges"))

## al

## al <- IRangesList(al)

## al <- unlist(al)

## al

## al <- IRanges(start = (start(al) + tgfbr1_bounds["start"] - 1), end = (end(al) + tgfbr1_bounds["start"] - 1))

## all(start(al) == bed_align$start)

## all(end(al) == bed_align$end)
