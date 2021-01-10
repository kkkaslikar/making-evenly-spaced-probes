library(tidyverse)
library(Gviz)
library(rafalib)


###### spacing probes

bed_align <- read_tsv(file = "Tgfbr1_cleaned_probes_chr_sC.bed", col_names = FALSE)

colnames(bed_align) <- c("chromosome", "start", "end", "probe", "Tm")

bed_align

x <- vector(mode = "list", length = (nrow(bed_align) - 1)) # initiating a blank list for holding all the probe sets

names(x) <- paste0(as.character(1:(nrow(bed_align) - 1)))

names(x)

min_probe_spacing <- 3000 # minimum space to be kept between end of first probe and start of next

for(i in 1:(nrow(bed_align) - 1)){
  
  j <- i
  
  a <- j+1
  
  while(a <= nrow(bed_align)){
    
    if((bed_align$start[a] - bed_align$end[j]) <= min_probe_spacing){
      a <- a+1
    } else {
      x[[as.character(i)]] <- append(x[[as.character(i)]], c("probe_1" = j, "probe_2" = a))
      j <- a
      a <- j+1
    }
    
  } 
}

# sapply(x, is.null)

# x[(sapply(x, is.null))]

# x <- x[!(sapply(x, is.null))]



length(x)

x$`1`

x$`821`

x <- purrr::compact(x)  # removing null list elements

x <- sapply(x, unique)

x$`1`

names(x) <- paste0("probeset_", names(x))

names(x) %>% head()

# x2 <- lapply(1:(nrow(bed_align)-1), function(i){
#   j <- i
#   a <- j+1
#   while(a <= nrow(bed_align)){
#     if((bed_align$start[a] - bed_align$end[j]) <= 3000){
#       a <- a+1
#     } else {
#       probe_pair <- c("probe_1" = j, "probe_2" = a)
#       j <- a
#       a <- j+1
#       if(!is.null(probe_pair)){
#         return(probe_pair)
#       }
#       }
#     }
#   }
# )

# names(x2) <- as.character( 1:(nrow(bed_align)-1) )


# x2 <-compact(x2)

probeset_lengths <- sapply(x, length)

names(probeset_lengths) <- NULL

probeset_lengths %>% head()

probeset_lengths %>% table()

largest_probesets <- x[which(probeset_lengths == max(probeset_lengths))]

probe_spacing <- vector(mode = "list", length = length(largest_probesets)) # creating null list for holding probe spacing values for all probe sets

names(probe_spacing) <- names(largest_probesets)

mean_probe_spacing <- vector(mode = "numeric", length = length(largest_probesets)) # creating null vector for holding probe spacing values for all probe sets

names(mean_probe_spacing) <- names(largest_probesets)

interquartile <- vector(mode = "numeric", length = length(largest_probesets))  # creating null vector for holding interquartile ranges of all probe sets

names(interquartile) <- names(largest_probesets)

for(i in seq_along(largest_probesets)){  # iterating along all the probe sets
  
  current_set_spacing <- vector(mode = "numeric", length = length(largest_probesets[[i]])) # create a null vector for holding distance between end of one probe and start of next probe (i.e., probe spacing) for current probe set
  
  for(j in 1:(length(current_set_spacing) - 1)){ # iterating along the probes in single probe set
    
    selected_probe_idx <- largest_probesets[[i]][j]  # index of selected probe in probe set
    
    next_probe_idx <- largest_probesets[[i]][j+1] # index of the next probe in probe set
    
    current_set_spacing[j] <- (bed_align$start[next_probe_idx] - bed_align$end[selected_probe_idx]) # start of next probe - end of selected probe
  }
  
  current_set_spacing <- current_set_spacing[-length(current_set_spacing)] # remove final element of current_set_spacing, which is a zero and is redundant
  
  probe_spacing[[i]] <- current_set_spacing # probe spacing for entire probe set
  
  mean_probe_spacing[i] <- mean(current_set_spacing)
  
  interquartile[i] <- IQR(current_set_spacing) # interquartile range of probe spacing of probe set
}

probe_spacing$probeset_1

mean_probe_spacing

interquartile

mypar(1, 1)

boxplot(mean_probe_spacing, main="Mean probe spacing for each set", horizontal=TRUE)

boxplot(interquartile, main = "Interquartile ranges of probe spacing for each set", horizontal=TRUE)

hist(interquartile, main = "Interquartile ranges of probe spacing for each set")

spaced_probes <- bed_align[largest_probesets$probeset_1, ]

write_tsv(x = spaced_probes, path = "Tgfbr1_spaced_probes.bed", col_names = FALSE)


spaced_probe_ranges <- GRanges(seqnames = spaced_probes$chromosome, ranges = IRanges(start = spaced_probes$start, end = spaced_probes$end), probe = spaced_probes$probe, Tm = spaced_probes$Tm)


probe_track <- Gviz::AnnotationTrack(range = spaced_probe_ranges, chromosome = 4, name  = "probes")

gtrack <- GenomeAxisTrack()

tgfbr1_bounds <- c("start" = 47353222, "end" = 47414926)

tgfbr1_range <- GRanges(seqnames = c("chr4"), ranges = IRanges (start = tgfbr1_bounds["start"], end = tgfbr1_bounds["end"]))

tgfbr1_track <- AnnotationTrack(range = tgfbr1_range, name = "Tgfbr1")
  
  
plotTracks(trackList =  list(probe_track, gtrack, tgfbr1_track))

spaced_probes$Tm %>% summary()


######

write_tsv(spaced_probes %>% select(start, probe), path = "oligo_analyzer.txt", col_names = FALSE)

save.image(file = "spacing-probes.RData")
