options(stringsAsFactors = F)
rm(list = ls())

data <- list()

dataDir <- "I:/Database/External_Data/SUPeR-seq/"
filenames <- list.files(dataDir, pattern = ".gz$")
for(dataFile in paste0(dataDir,filenames)){
  if(grepl(pattern = "Tang2009", dataFile)){
    data[[dataFile]] <- read.table(con <- gzfile(dataFile, "rt"), header = T,sep = " ")
  }else{
    data[[dataFile]] <- read.table(con <- gzfile(dataFile, "rt"), header = T,sep = " ",quote = "\"")
  }
  close(con)
}

refseqGenes <- data$`I:/Database/External_Data/SUPeR-seq/GSM1290599_SUPeR-seq_ESC1_FPKM_gencode.txt.gz`$gene_name
gencodeGenes <- data$`I:/Database/External_Data/SUPeR-seq/GSM1290602_Tang2009-ESC1_FPKM_gencode.txt.gz`$gene_short_name
gencodeID <- data$`I:/Database/External_Data/SUPeR-seq/GSM1290602_Tang2009-ESC1_FPKM_gencode.txt.gz`$tracking_id
commonGenes <- intersect(refseqGenes,gencodeGenes)

geneBiotype <- read.table(file = "I:/Database/GENCODE/gencode.vM2.long_noncoding_RNAs.geneBiotype.txt", header = F,sep = "\t")
table(geneBiotype[,3])
length(intersect(commonGenes, geneBiotype[,4]))

comMatrix <- matrix(NA,length(commonGenes), 6)
rownames(comMatrix) <- commonGenes
colnames(comMatrix) <- paste0(rep(c("SUPeR-seq_ESC","Tang2009-ESC"),each =3), rep(1:3,times =2))
for(dataFile in names(data)){
  if(grepl(pattern = "Tang2009", dataFile)){
    tmp <- paste0("Tang2009",gsub(pattern = ".*_Tang2009|_FPKM.*",replacement = "",x = dataFile, perl = T))
    comMatrix[,tmp] <- data[[dataFile]]$FPKM[match(commonGenes, data[[dataFile]]$gene_short_name)]
  }else{
    tmp <- paste0("SUPeR-seq",gsub(pattern = ".*_SUPeR-seq|_FPKM.*",replacement = "",x = dataFile, perl = T))
    comMatrix[,tmp] <- data[[dataFile]]$FPKM[match(commonGenes, data[[dataFile]]$gene_name)]
  }
}

meanMatrix <- cbind(indep = rowMeans(comMatrix[,1:3]), polya = rowMeans(comMatrix[,4:6]))
stat <- as.data.frame(meanMatrix > 0)
stat$both <- stat$indep & stat$polya
stat$only_indep <- xor(stat$indep,stat$both)
stat$only_polya <- xor(stat$polya,stat$both)

colSums(stat[intersect(commonGenes, geneBiotype[,4]), ])

colSums(stat)

### polya features of genes or lncrnas
polya <- read.table(gzfile("I:/Database/GENCODE/gencode.vM2.metadata.PolyA_feature.gz"))
polya.lncrna <- intersect(geneBiotype[,2], polya[,1])
length(polya.lncrna)/length(geneBiotype[,2])# polya percentage in lncrna transcripts
geneBiotype$polya <- 0
geneBiotype[match(polya.lncrna, geneBiotype[,2]),"polya"] <- 1
poly.lncrna.genename <- unique(subset(geneBiotype[,c(4,5)], polya == 1)[,1])

colSums(na.omit(stat[poly.lncrna.genename,]))
colSums(na.omit(stat[setdiff(geneBiotype[,4],poly.lncrna.genename),]))
