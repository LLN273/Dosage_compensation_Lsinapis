# Differential gene expression using DEseq2

#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
#biocLite("sva")
library("DESeq2")
library(sva)

sessionInfo()

#Clear all states
rm(list=ls(all=TRUE))
dev.off()


#Load gene count matrix and labels ####

setwd("/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/03_Normalized_libs")
countData <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))
head(countData)
nrow(countData)

setwd("/home/luisleal/MYPROJ/3_DosageCompensation_LS/Scripts_NEW")
colData <- read.csv("sample_info.txt", sep="\t", row.names=1)
colData






##### output folders

# standard filtering
OUT_FOLD_STDF <- "/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/10_DESEQ2_sva"








#Pre-filtering I: Only keep genes with >1 samples with non-zero read count ####  
countData <- countData[rowSums(countData > 0) >= 2 , ]
nrow(countData)



#Make separate data sets for each developmental stage ####
Instar_countdata <- subset(countData, select = c("P5052_202_S48", "P5052_210_S51", "P5052_218_S55", 
                                                 "P5052_203_S49", "P5052_211_S52", "P5052_219_S56"))
Instar_coldata <- as.data.frame(colData[c(1, 5, 9, 2, 6, 10), ])

Pupa_countdata <- subset(countData, select = c("P5052_204_S34", "P5052_212_S36", "P5052_220_S57",
                                               "P5052_205_S35", "P5052_213_S53", "P5052_221_S58"))
Pupa_coldata <- as.data.frame(colData[c(3, 7, 11, 4, 8, 12)  , ])

Adult_countdata <- subset(countData, select = c("P5052_233_S60", "P5052_234_S61", "P5052_235_S62",
                                                "P5052_226_S59", "P5052_227_S38", "P7553_325_S10"))
Adult_coldata <- as.data.frame(colData[c(15, 16, 17, 13, 14, 18), ])

Instar_coldata <- Instar_coldata[ , 2:3]
Pupa_coldata <- Pupa_coldata[ , 2:3]
Adult_coldata <- Adult_coldata[ , 2:3]

#Pre-filtering I: Only keep genes with >1 samples with non-zero read count ####    >>> this is now done above, for all samples (as opposed to doing it for each the individual stage)
#Instar_countdata <- Instar_countdata[rowSums(Instar_countdata > 0) >= 2 , ]
#Pupa_countdata <- Pupa_countdata[rowSums(Pupa_countdata > 0) >= 2 , ]
#Adult_countdata <- Adult_countdata[rowSums(Adult_countdata > 0) >= 2 , ]

#Pre-filtering II: for each developmental stage, remove genes with zero counts in all samples ####  
Instar_countdata <- Instar_countdata[rowSums(Instar_countdata > 0) >= 1 , ]
nrow(Instar_countdata)
Pupa_countdata <- Pupa_countdata[rowSums(Pupa_countdata > 0) >= 1 , ]
nrow(Pupa_countdata)
Adult_countdata <- Adult_countdata[rowSums(Adult_countdata > 0) >= 1 , ]
nrow(Adult_countdata)


## Add 1 count to every gene/sample (avoids some false positives)
Instar_countdata <- Instar_countdata + 1
head(Instar_countdata)
Pupa_countdata <- Pupa_countdata + 1
head(Pupa_countdata)
Adult_countdata <- Adult_countdata + 1
head(Adult_countdata)


#Check all sample IDs in colData are also in CountData and match their orders ####
all(rownames(colData) %in% colnames(countData))
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))

all(rownames(Instar_coldata) %in% colnames(Instar_countdata))
Instar_countdata <- Instar_countdata[, rownames(Instar_coldata)]
all(rownames(Instar_coldata) == colnames(Instar_countdata))

all(rownames(Pupa_coldata) %in% colnames(Pupa_countdata))
Pupa_countdata <- Pupa_countdata[, rownames(Pupa_coldata)]
all(rownames(Pupa_coldata) == colnames(Pupa_countdata))

all(rownames(Adult_coldata) %in% colnames(Adult_countdata))
Adult_countdata <- Adult_countdata[, rownames(Adult_coldata)]
all(rownames(Adult_coldata) == colnames(Adult_countdata))




#Create a DESeqDataSet from count matrix and labels (no confounding effects removed at this stage) ####
ddsInstar <- DESeqDataSetFromMatrix(countData = Instar_countdata,
                                    colData = Instar_coldata, design = ~Sex)
ddsInstar

ddsPupa <- DESeqDataSetFromMatrix(countData = Pupa_countdata,
                                  colData = Pupa_coldata, design = ~Sex)
ddsPupa

ddsAdult <- DESeqDataSetFromMatrix(countData = Adult_countdata,
                                   colData = Adult_coldata, design = ~Sex)
ddsAdult

#Pre-filtering IIIs: Only keep rows with baseMean >2 (usually this would be baseMean >1 but we added one count to every gene/sample)
ddsInstar <- ddsInstar[rowMeans(counts(ddsInstar)) > 2, ]
ddsInstar
ddsPupa <- ddsPupa[rowMeans(counts(ddsPupa)) > 2, ]
ddsPupa
ddsAdult <- ddsAdult[rowMeans(counts(ddsAdult)) > 2, ]
ddsAdult





######## MODEL BATCH EFFECTS: Remove hidden batch effects using the sva package

ddsInstar_sva <- ddsInstar
ddsInstar_sva <- estimateSizeFactors(ddsInstar_sva)                   # normalization
dat_Instar  <- counts(ddsInstar_sva, normalized = TRUE)
idx_Instar  <- rowMeans(dat_Instar) > 2                                 # filter
dat_Instar  <- dat_Instar[idx_Instar, ]

ddsPupa_sva <- ddsPupa
ddsPupa_sva <- estimateSizeFactors(ddsPupa_sva)                   
dat_Pupa <- counts(ddsPupa_sva, normalized = TRUE)
idx_Pupa  <- rowMeans(dat_Pupa) > 2                                 
dat_Pupa  <- dat_Pupa[idx_Pupa, ]

ddsAdult_sva <- ddsAdult
ddsAdult_sva <- estimateSizeFactors(ddsAdult_sva)                  
dat_Adult <- counts(ddsAdult_sva, normalized = TRUE)
idx_Adult  <- rowMeans(dat_Adult) > 2                                 
dat_Adult  <- dat_Adult[idx_Adult, ]

# declare model with variable of interest
mod_Instar  <- model.matrix(~ Sex, colData(ddsInstar_sva))  # when looking for DGE across sexes
mod_Pupa  <- model.matrix(~ Sex, colData(ddsPupa_sva))
mod_Adult  <- model.matrix(~ Sex, colData(ddsAdult_sva))

#declare known adjustment variables, or let sva discover batch effects (creates new surrogate variables)
mod0_Instar <- model.matrix(~ 1, colData(ddsInstar_sva))
mod0_Pupa <- model.matrix(~ 1, colData(ddsPupa_sva))
mod0_Adult <- model.matrix(~ 1, colData(ddsAdult_sva))

# estimate number of surrogate variables
#svseq <- svaseq(dat, mod, mod0, n.sv = 3)  ## >>  USE for 03_SDP_LDP (ONLY)
#svseq <- svaseq(dat, mod, mod0, n.sv = 1)  ## >>  USE for 12_LDAF_LDAM (ONLY)
svseq_Instar <- svaseq(dat_Instar, mod_Instar, mod0_Instar)    ##   If  the sva function  is  called  without  the n.sv argument  specified, the  number  
##   of factors will be estimated for you
svseq_Instar$sv

svseq_Pupa <- svaseq(dat_Pupa, mod_Pupa, mod0_Pupa)
svseq_Pupa$sv

svseq_Adult <- svaseq(dat_Adult, mod_Adult, mod0_Adult)
#svseq_Adult <- svaseq(dat_Adult, mod_Adult, mod0_Adult, n.sv = 1)  ## >>  USE for 12_LDAF_LDAM (ONLY)
svseq_Adult$sv

## Use sva to remove any effect of surrogate variables
## (adjust according to number of surrogate variables)
ddssva_Instar_f <- ddsInstar
ddssva_Instar_f$SV1 <- svseq_Instar$sv[,1]
ddssva_Instar_f$SV2 <- svseq_Instar$sv[,2]

design(ddssva_Instar_f) <- ~ SV1 + SV2 + Sex 

ddssva_Pupa_f <- ddsPupa
ddssva_Pupa_f$SV1 <- svseq_Pupa$sv[,1]
ddssva_Pupa_f$SV2 <- svseq_Pupa$sv[,2]

design(ddssva_Pupa_f) <- ~ SV1 + SV2 + Sex 

ddssva_Adult_f <- ddsAdult
ddssva_Adult_f$SV1 <- svseq_Adult$sv[,1]
ddssva_Adult_f$SV2 <- svseq_Adult$sv[,2]

design(ddssva_Adult_f) <- ~ SV1 + SV2 + Sex 


             

               








######### Run the default analysis for DESeq2 ####



# Without batch effect correction

ddsInstar_0 <- DESeq(ddsInstar)
res_Instar_0 <- results(ddsInstar_0, alpha=.05)
res_Instar_0

ddsPupa_0 <- DESeq(ddsPupa)
res_Pupa_0 <- results(ddsPupa_0, alpha=.05)
res_Pupa_0

ddsAdult_0 <- DESeq(ddsAdult)
res_Adult_0 <- results(ddsAdult_0, alpha=.05)
res_Adult_0



#  With sva correction for batch effects
ddsInstar <- DESeq(ddssva_Instar_f)
res_Instar <- results(ddsInstar, alpha=.05)
res_Instar

ddsPupa <- DESeq(ddssva_Pupa_f)
res_Pupa <- results(ddsPupa, alpha=.05)
res_Pupa

ddsAdult <- DESeq(ddssva_Adult_f)
res_Adult <- results(ddsAdult, alpha=.05)
res_Adult







#Results ####


## no batch effects correction

Filtered_P0.05_Instar_0 <- subset(res_Instar_0, padj<0.05)
Filtered_P0.05_Pupa_0 <- subset(res_Pupa_0, padj<0.05)
Filtered_P0.05_Adult_0 <- subset(res_Adult_0, padj<0.05)

summary(Filtered_P0.05_Instar_0)
summary(Filtered_P0.05_Pupa_0)
summary(Filtered_P0.05_Adult_0)


## when correting for batch effects
Filtered_P0.05_Instar <- subset(res_Instar, padj<0.05)
Filtered_P0.05_Pupa <- subset(res_Pupa, padj<0.05)
Filtered_P0.05_Adult <- subset(res_Adult, padj<0.05)

summary(Filtered_P0.05_Instar)
summary(Filtered_P0.05_Pupa)
summary(Filtered_P0.05_Adult)


# Filter logFC > |1.0|) ####
Filtered_P0.05_LFC_1.5_Instar <- subset(Filtered_P0.05_Instar, log2FoldChange >1 | log2FoldChange < -1)
Filtered_P0.05_LFC_1.5_Pupa <- subset(Filtered_P0.05_Pupa, log2FoldChange >1 | log2FoldChange < -1)
Filtered_P0.05_LFC_1.5_Adult <- subset(Filtered_P0.05_Adult, log2FoldChange >1 | log2FoldChange < -1)

summary(Filtered_P0.05_LFC_1.5_Instar)
summary(Filtered_P0.05_LFC_1.5_Pupa)
summary(Filtered_P0.05_LFC_1.5_Adult)


#Filter baseMean > 10 ####
Filtered_P0.05_LFC_1.5_base_Mean_10_Instar <- subset(Filtered_P0.05_LFC_1.5_Instar, baseMean>10)
Filtered_P0.05_LFC_1.5_base_Mean_10_Pupa <- subset(Filtered_P0.05_LFC_1.5_Pupa, baseMean>10)
Filtered_P0.05_LFC_1.5_base_Mean_10_Adult <- subset(Filtered_P0.05_LFC_1.5_Adult, baseMean>10)

summary(Filtered_P0.05_LFC_1.5_base_Mean_10_Instar)
summary(Filtered_P0.05_LFC_1.5_base_Mean_10_Pupa)
summary(Filtered_P0.05_LFC_1.5_base_Mean_10_Adult)




# Select Male/Female biased genes ####
Filtered_P0.05_MBG_Instar <- subset(Filtered_P0.05_LFC_1.5_base_Mean_10_Instar, log2FoldChange >1)
Filtered_P0.05_MBG_Pupa <- subset(Filtered_P0.05_LFC_1.5_base_Mean_10_Pupa , log2FoldChange > 1)
Filtered_P0.05_MBG_Adult <- subset(Filtered_P0.05_LFC_1.5_base_Mean_10_Adult, log2FoldChange > 1)
summary(Filtered_P0.05_MBG_Instar)
summary(Filtered_P0.05_MBG_Pupa)
summary(Filtered_P0.05_MBG_Adult)

Filtered_P0.05_FBG_Instar <- subset(Filtered_P0.05_LFC_1.5_base_Mean_10_Instar, log2FoldChange < -1)
Filtered_P0.05_FBG_Pupa <- subset(Filtered_P0.05_LFC_1.5_base_Mean_10_Pupa, log2FoldChange < -1)
Filtered_P0.05_FBG_Adult <- subset(Filtered_P0.05_LFC_1.5_base_Mean_10_Adult, log2FoldChange < -1)
summary(Filtered_P0.05_FBG_Instar)
summary(Filtered_P0.05_FBG_Pupa)
summary(Filtered_P0.05_FBG_Adult)



# Filter out genes without expression
expressed_instar <- subset(res_Instar, baseMean >0)
expressed_pupa <- subset(res_Pupa, baseMean >0)
expressed_adult <- subset(res_Adult, baseMean >0)



# Non-biased genes
Nonbiased_genes_Instar <- subset(res_Instar, padj>0.05 | baseMean<10 | (log2FoldChange <1 & log2FoldChange > -1))
Nonbiased_genes_Pupa <- subset(res_Pupa, padj>0.05 | baseMean<10 | (log2FoldChange <1 & log2FoldChange > -1))
Nonbiased_genes_Adult <- subset(res_Adult, padj>0.05 | baseMean<10 | (log2FoldChange <1 & log2FoldChange > -1))
summary(Nonbiased_genes_Instar)
summary(Nonbiased_genes_Pupa)
summary(Nonbiased_genes_Adult)

Nonbiased_genes_Instar_f <- subset(res_Instar, padj>0.05 | baseMean<10 | log2FoldChange > -1)
Nonbiased_genes_Pupa_f <- subset(res_Pupa, padj>0.05 | baseMean<10 | log2FoldChange > -1)
Nonbiased_genes_Adult_f <- subset(res_Adult, padj>0.05 | baseMean<10 | log2FoldChange > -1)
summary(Nonbiased_genes_Instar_f)
summary(Nonbiased_genes_Pupa_f)
summary(Nonbiased_genes_Adult_f)

Nonbiased_genes_Instar_m <- subset(res_Instar, padj>0.05 | baseMean<10 | log2FoldChange < 1)
Nonbiased_genes_Pupa_m <- subset(res_Pupa, padj>0.05 | baseMean<10 | log2FoldChange < 1)
Nonbiased_genes_Adult_m <- subset(res_Adult, padj>0.05 | baseMean<10 | log2FoldChange < 1)
summary(Nonbiased_genes_Instar_m)
summary(Nonbiased_genes_Pupa_m)
summary(Nonbiased_genes_Adult_m)







###### Write results to file ####


######### standard filtering (genes with p-adj < 0.05, |logFC| > 1, base_Mean > 10)

setwd(OUT_FOLD_STDF)

write.table(Filtered_P0.05_LFC_1.5_base_Mean_10_Instar,
          file = "Filtered_P0.05_Instar.txt", sep = "\t", quote = FALSE)
write.table(Filtered_P0.05_LFC_1.5_base_Mean_10_Pupa, 
          file = "Filtered_P0.05_Pupa.txt", sep = "\t", quote = FALSE)
write.table(Filtered_P0.05_LFC_1.5_base_Mean_10_Adult, 
          file = "Filtered_P0.05_Adult.txt", sep = "\t", quote = FALSE)

write.table(Filtered_P0.05_MBG_Instar, 
          file = "MBG_Instar.txt", sep = "\t", quote = FALSE)
write.table(Filtered_P0.05_MBG_Pupa, 
          file = "MBG_Pupa.txt", sep = "\t", quote = FALSE)
write.table(Filtered_P0.05_MBG_Adult, 
          file = "MBG_Adult.txt", sep = "\t", quote = FALSE)

write.table(Filtered_P0.05_FBG_Instar, 
          file = "FBG_Instar.txt", sep = "\t", quote = FALSE)
write.table(Filtered_P0.05_FBG_Pupa, 
          file = "FBG_Pupa.txt", sep = "\t", quote = FALSE)
write.table(Filtered_P0.05_FBG_Adult, 
          file = "FBG_Adult.txt", sep = "\t", quote = FALSE)

#All LFC
write.table(expressed_instar, file = "LFC_Instar.txt", sep = "\t", quote = FALSE)
write.table(expressed_pupa, file = "LFC_Pupa.txt", sep = "\t", quote = FALSE)
write.table(expressed_adult, file = "LFC_Adult.txt", sep = "\t", quote = FALSE)


# Non-biased genes
write.table(Nonbiased_genes_Instar, file = "Nonbiased_genes_Instar.txt", sep = "\t", quote = FALSE)
write.table(Nonbiased_genes_Pupa, file = "Nonbiased_genes_Pupa.txt", sep = "\t", quote = FALSE)
write.table(Nonbiased_genes_Adult, file = "Nonbiased_genes_Adult.txt", sep = "\t", quote = FALSE)

write.table(Nonbiased_genes_Instar_f, file = "Nonbiased_genes_Instar_f.txt", sep = "\t", quote = FALSE)
write.table(Nonbiased_genes_Pupa_f, file = "Nonbiased_genes_Pupa_f.txt", sep = "\t", quote = FALSE)
write.table(Nonbiased_genes_Adult_f, file = "Nonbiased_genes_Adult_f.txt", sep = "\t", quote = FALSE)

write.table(Nonbiased_genes_Instar_m, file = "Nonbiased_genes_Instar_m.txt", sep = "\t", quote = FALSE)
write.table(Nonbiased_genes_Pupa_m, file = "Nonbiased_genes_Pupa_m.txt", sep = "\t", quote = FALSE)
write.table(Nonbiased_genes_Adult_m, file = "Nonbiased_genes_Adult_m.txt", sep = "\t", quote = FALSE)






######### 

