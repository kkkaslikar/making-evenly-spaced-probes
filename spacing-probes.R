library(tidyverse)
library(Gviz)
library(rafalib)





###### spacing probes


bed_align <- read_tsv(file = "Tgfbr1_cleaned_probes_chr_sC.bed", col_names = FALSE)

colnames(bed_align) <- c("chromosome", "start", "end", "probe", "Tm")

x <- list("blank" = 1) # initiating a blank list


for(i in 1:(nrow(bed_align) - 1)){
  
  j <- i
  
  a <- j+1
  
  while(a <= nrow(bed_align)){
    
    if((bed_align$start[a] - bed_align$end[j]) <= 3000){
      a <- a+1
    } else {
      x[[as.character(i)]] <- append(x[[as.character(i)]], c("probe_1" = j, "probe_2" = a))
      j <- a
      a <- j+1
    }
    
  } 
  
}


names(x)

x <- x[-c(1)]

for(i in names(x)){
  
  x[[i]] <- cbind(x[[i]][names(x[[i]]) == "probe_1"], x[[i]][names(x[[i]]) == "probe_2"])
  
}

for(i in names(x)){
  rownames(x[[i]]) <- NULL
}

sapply(x, nrow) %>% table()

(sapply(x, nrow) %>% summary())



y <- x[which(sapply(x, nrow) >= 10)]

g <- vector(mode = "numeric", length = length(y))

for(i in 1:length(y)){
  
  h <- vector(mode = "numeric", length = nrow(y[[i]]))
  for(j in 1:nrow(y[[i]])){
    h[j] <- (bed_align$start[y[[i]][j, 2]] - bed_align$end[y[[i]][j, 1]])
  }
  
  g[i] <- IQR(h)
}

sapply(y, nrow)[which(sapply(y, nrow) == max(sapply(y, nrow)))]

mypar(1, 2)

boxplot(g, main = "Interquartile ranges of probe spacing")

hist(g, main = "Interquartile ranges of probe spacing")


summary(g)

g[1:27] # IQRs


spaced_probes <- bed_align[y[1] %>% flatten_dbl() %>% unique(), ]

write_tsv(x = spaced_probes, path = "Tgfbr1_spaced_probes.bed", col_names = FALSE)



h <- vector(mode = "numeric", length = nrow(y[[1]]))


for(j in 1:nrow(y[[1]])){
  h[j] <- (bed_align$start[y[[1]][j, 2]] - bed_align$end[y[[1]][j, 1]])
}



summary(h) # probe spacing distribution

boxplot(h, main = "Boxplot of probe spacing distribution")

hist(h, main = "Histogram of probe spacing distribution")

h[which(h == max(h))] # maximum space in selected probe set


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
