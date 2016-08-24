options(stringsAsFactors = F)
rm(list = ls())
mrna <- read.table(file = "I:/Database/UCSC/CancerBrowser/TCGA_BRCA_exp_HiSeqV2-2015-02-24/genomicMatrix", header = T,sep = "\t",row.names = 1,quote = "",check.names = F)
mature <- read.table(file = "I:/Database/FireBrowse/gdac.broadinstitute.org_BRCA.miRseq_Mature_Preprocess.Level_3.2016012800.0.0/BRCA.miRseq_mature_RPM.txt",header = T,sep = "\t",quote = "",row.names = 1,check.names = F)
colnames(mature) <- substr(colnames(mature),1,15)
mature[is.na(mature)] <- 0
mature <- log2(mature + 1)

clinical <- read.csv(file = "../sample_TCGA_RNA_Seq_fromWJQ.csv",header = T,row.names = 1)
clinical[clinical[,"N_C_TYPES"] == "C","N_C_TYPES"] <- "Tumor"
clinical[clinical[,"N_C_TYPES"] == "N","N_C_TYPES"] <- "Normal"

samples <- Reduce(intersect,list(colnames(mrna), colnames(mature), rownames(clinical)))

genes <- rownames(mrna)
mirna <- rownames(mature)
brcaData <- cbind(clinical[samples, ], t(mature[,samples]), t(mrna[,samples]))

save(genes, mirna, brcaData, file = "data.Rdata")
