options(stringsAsFactors = F)
rm(list = ls())

library(GenomicRanges)

lncrna <- read.table(file = gzfile("I:/Database/GENCODE/gencode.vM2.long_noncoding_RNAs.gtf.gz"), sep = "\t")
head(lncrna)
GRanges(seqnames = "V1", ranges = IRanges(),)