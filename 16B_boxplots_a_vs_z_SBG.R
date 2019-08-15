# PART I: Boxplots A vs Z, distributions for sex-biased genes (SBG)
# PART II: Boxplots A vs Z, distributions for upregulated genes (FBG in females; MBG in males)

#Clear all states
rm(list=ls(all=TRUE))
dev.off()



############### Paths and folders

### Folder containind expression data after removing SBG
folder_1 <- "/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/12_Filter_sex_biased_genes"


### Folder containind expression data before removing SBG 
folder_2 <-"/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/03_Normalized_libs"  

### Output folder 
out_dir <-"/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/16B_MWU-test_SBG" 


############## Load data


### First dataset

setwd(folder_1) 

##### Load gene list produced after removing SBG; genes with zero counts removed individually for each sex
instar_V_f_noSBG <- read.delim("nonbiased_genes-instar_V-assigned_A_or_Z_female-filtered.txt", header = TRUE)
instar_V_m_noSBG <- read.delim("nonbiased_genes-instar_V-assigned_A_or_Z_male-filtered.txt", header = TRUE)

pupa_f_noSBG <- read.delim("nonbiased_genes-pupa-assigned_A_or_Z_female-filtered.txt", header = TRUE)
pupa_m_noSBG <- read.delim("nonbiased_genes-pupa-assigned_A_or_Z_male-filtered.txt", header = TRUE)

adult_f_noSBG <- read.delim("nonbiased_genes-adult-assigned_A_or_Z_female-filtered.txt", header = TRUE)
adult_m_noSBG <- read.delim("nonbiased_genes-adult-assigned_A_or_Z_male-filtered.txt", header = TRUE)


####### Load list of female and male biased genes
FBG_Instar <- read.delim("FBG_Instar.txt", header = TRUE)
MBG_Instar <- read.delim("MBG_Instar.txt", header = TRUE)
FBG_Pupa <- read.delim("FBG_Pupa.txt", header = TRUE)
MBG_Pupa <- read.delim("MBG_Pupa.txt", header = TRUE)
FBG_Adult <- read.delim("FBG_Adult.txt", header = TRUE)
MBG_Adult <- read.delim("MBG_Adult.txt", header = TRUE)

FBG_Instar$gene_id <- rownames(FBG_Instar)
MBG_Instar$gene_id <- rownames(MBG_Instar)
FBG_Pupa$gene_id <- rownames(FBG_Pupa)
MBG_Pupa$gene_id <- rownames(MBG_Pupa)
FBG_Adult$gene_id <- rownames(FBG_Adult)
MBG_Adult$gene_id <- rownames(MBG_Adult)



### Second dataset   

setwd(folder_2) 

####### Before removing SBG; genes with zero counts removed individually for each sex
instar_V_f_AllGenes <- read.delim("instar_V-assigned_A_or_Z_female-filtered.txt", header = TRUE)
instar_V_m_AllGenes <- read.delim("instar_V-assigned_A_or_Z_male-filtered.txt", header = TRUE)

pupa_f_AllGenes <- read.delim("pupa-assigned_A_or_Z_female-filtered.txt", header = TRUE)
pupa_m_AllGenes <- read.delim("pupa-assigned_A_or_Z_male-filtered.txt", header = TRUE)

adult_f_AllGenes <- read.delim("adult-assigned_A_or_Z_female-filtered.txt", header = TRUE)
adult_m_AllGenes <- read.delim("adult-assigned_A_or_Z_male-filtered.txt", header = TRUE)



############################################################################ PART I



####### SBG computed as follows (see example below)
####### Example:
####### SBG_instar_v_female = genes in AllGenes_instar_v_female also listed in FBG_instar_v
#######                       + genes in AllGenes_instar_v_female also listed in MBG_instar_v, if their expression in non-zero in females
#######
####### Note: AllGenes refers to gene list before removing SBG

# include all SBG (both regulated and downregulateds)
instar_V_f_allSBG <- instar_V_f_AllGenes[(instar_V_f_AllGenes$gene_id %in% FBG_Instar$gene_id ) | ((instar_V_f_AllGenes$gene_id %in% MBG_Instar$gene_id) & (instar_V_f_AllGenes$FPKM_instar_V_female != "0")), ]
instar_V_m_allSBG <- instar_V_m_AllGenes[(instar_V_m_AllGenes$gene_id %in% MBG_Instar$gene_id) | ((instar_V_m_AllGenes$gene_id %in% FBG_Instar$gene_id) & (instar_V_m_AllGenes$FPKM_instar_V_male != "0")),]

pupa_f_allSBG <- pupa_f_AllGenes[(pupa_f_AllGenes$gene_id %in% FBG_Pupa$gene_id) | ((pupa_f_AllGenes$gene_id %in% MBG_Pupa$gene_id) & (pupa_f_AllGenes$FPKM_pupa_female != "0")), ]
pupa_m_allSBG <- pupa_m_AllGenes[(pupa_m_AllGenes$gene_id %in% MBG_Pupa$gene_id) | ((pupa_m_AllGenes$gene_id %in% FBG_Pupa$gene_id) & (pupa_m_AllGenes$FPKM_pupa_male != "0")), ]

adult_f_allSBG <- adult_f_AllGenes[(adult_f_AllGenes$gene_id %in% FBG_Adult$gene_id) | ((adult_f_AllGenes$gene_id %in% MBG_Adult$gene_id) & (adult_f_AllGenes$FPKM_adult_female != "0")), ]
adult_m_allSBG <- adult_m_AllGenes[(adult_m_AllGenes$gene_id %in% MBG_Adult$gene_id) | ((adult_m_AllGenes$gene_id %in% FBG_Adult$gene_id) & (adult_m_AllGenes$FPKM_adult_male != "0")),]




##################################### PROCESS DATA


############### Dataset I

instar_V_female <- instar_V_f_noSBG[c(4,6)]
instar_V_female$group <- rep("1_instar_V_female", nrow(instar_V_female))
instar_V_female$group <- paste(instar_V_female$group, instar_V_female$chromosome)
instar_V_female <- instar_V_female[c(1,3)]
names(instar_V_female)[1] <-"FPKM"

pupa_female <- pupa_f_noSBG[c(4,6)]
pupa_female$group <- rep("3_pupa_female", nrow(pupa_female))
pupa_female$group <- paste(pupa_female$group, pupa_female$chromosome)
pupa_female <- pupa_female[c(1,3)]
names(pupa_female)[1] <-"FPKM"

adult_female <- adult_f_noSBG[c(4,6)]
adult_female$group <- rep("5_adult_female", nrow(adult_female))
adult_female$group <- paste(adult_female$group, adult_female$chromosome)
adult_female <- adult_female[c(1,3)]
names(adult_female)[1] <-"FPKM"

instar_V_male <- instar_V_m_noSBG[c(5,6)]
instar_V_male$group <- rep("2_instar_V_male", nrow(instar_V_male))
instar_V_male$group <- paste(instar_V_male$group, instar_V_male$chromosome)
instar_V_male <- instar_V_male[c(1,3)]
names(instar_V_male)[1] <-"FPKM"

pupa_male <- pupa_m_noSBG[c(5,6)]
pupa_male$group <- rep("4_pupa_male", nrow(pupa_male))
pupa_male$group <- paste(pupa_male$group, pupa_male$chromosome)
pupa_male <- pupa_male[c(1,3)]
names(pupa_male)[1] <-"FPKM"

adult_male <- adult_m_noSBG[c(5,6)]
adult_male$group <- rep("6_adult_male", nrow(adult_male))
adult_male$group <- paste(adult_male$group, adult_male$chromosome)
adult_male <- adult_male[c(1,3)]
names(adult_male)[1] <-"FPKM"

all_samples <- rbind(instar_V_female, instar_V_male, pupa_female, pupa_male, adult_female, adult_male)


############### Dataset II 

instar_V_female_allSBG <- instar_V_f_allSBG[c(4,6)]
instar_V_female_allSBG$group <- rep("1_instar_V_female_allSBG", nrow(instar_V_female_allSBG))
instar_V_female_allSBG$group <- paste(instar_V_female_allSBG$group, instar_V_female_allSBG$chromosome)
instar_V_female_allSBG <- instar_V_female_allSBG[c(1,3)]
names(instar_V_female_allSBG)[1] <-"FPKM"

pupa_female_allSBG <- pupa_f_allSBG[c(4,6)]
pupa_female_allSBG$group <- rep("3_pupa_female_allSBG", nrow(pupa_female_allSBG))
pupa_female_allSBG$group <- paste(pupa_female_allSBG$group, pupa_female_allSBG$chromosome)
pupa_female_allSBG <- pupa_female_allSBG[c(1,3)]
names(pupa_female_allSBG)[1] <-"FPKM"

adult_female_allSBG <- adult_f_allSBG[c(4,6)]
adult_female_allSBG$group <- rep("5_adult_female_allSBG", nrow(adult_female_allSBG))
adult_female_allSBG$group <- paste(adult_female_allSBG$group, adult_female_allSBG$chromosome)
adult_female_allSBG <- adult_female_allSBG[c(1,3)]
names(adult_female_allSBG)[1] <-"FPKM"

instar_V_male_allSBG <- instar_V_m_allSBG[c(5,6)]
instar_V_male_allSBG$group <- rep("2_instar_V_male_allSBG", nrow(instar_V_male_allSBG))
instar_V_male_allSBG$group <- paste(instar_V_male_allSBG$group, instar_V_male_allSBG$chromosome)
instar_V_male_allSBG <- instar_V_male_allSBG[c(1,3)]
names(instar_V_male_allSBG)[1] <-"FPKM"

pupa_male_allSBG <- pupa_m_allSBG[c(5,6)]
pupa_male_allSBG$group <- rep("4_pupa_male_allSBG", nrow(pupa_male_allSBG))
pupa_male_allSBG$group <- paste(pupa_male_allSBG$group, pupa_male_allSBG$chromosome)
pupa_male_allSBG <- pupa_male_allSBG[c(1,3)]
names(pupa_male_allSBG)[1] <-"FPKM"

adult_male_allSBG <- adult_m_allSBG[c(5,6)]
adult_male_allSBG$group <- rep("6_adult_male_allSBG", nrow(adult_male_allSBG))
adult_male_allSBG$group <- paste(adult_male_allSBG$group, adult_male_allSBG$chromosome)
adult_male_allSBG <- adult_male_allSBG[c(1,3)]
names(adult_male_allSBG)[1] <-"FPKM"

allSBG <- rbind(instar_V_female_allSBG, instar_V_male_allSBG, pupa_female_allSBG, pupa_male_allSBG, adult_female_allSBG, adult_male_allSBG)


######################### BOXPLOTS



boxplot(log2(all_samples$FPKM)~all_samples$group, ylim = c(-6, 18),
        col = c("grey90", "grey90"), notch = TRUE,
        at = c(1,2, 4,5,   8,9, 11,12,  15,16, 18,19),
        outline = FALSE, boxwex = 0.7,  xaxt = "n", frame.plot = FALSE, border="grey65",
        boxlty = 2)
axis(2, labels = "log2 FPKM(>0)", cex.axis = 1.2, at = 4, line = 1.5, tck = 0)

boxplot(log2(allSBG$FPKM)~allSBG$group, ylim = c(-8, 14), add = TRUE, 
        col = c("lightgrey", "darkorange"), notch = FALSE,
        at = c(1.2,2.2, 4.2,5.2,   8.2,9.2, 11.2,12.2,  15.2,16.2, 18.2,19.2),
        outline = FALSE, boxwex = 0.7,  xaxt = "n", frame.plot = FALSE)

text(1.5, -6.5, "\\VE", vfont=c("sans serif","plain"), cex =1.9)
text(4.5, -6.5, "\\MA", vfont=c("sans serif","plain"), cex =1.9)

text(8.5, -6.5, "\\VE", vfont=c("sans serif","plain"), cex =1.9)
text(11.5, -6.5, "\\MA", vfont=c("sans serif","plain"), cex =1.9)

text(15.5, -6.5, "\\VE", vfont=c("sans serif","plain"), cex =1.9)
text(18.5, -6.5, "\\MA", vfont=c("sans serif","plain"), cex =1.9)

axis(1, labels = c("Larva", "Pupa", "Adult"), lwd = 0, cex.axis = 1.5,
     at = c(3, 10, 17), line = 1)
legend(16.5, 18, legend = c("Autosomes","Z"), cex = 1,
       fill = c("lightgrey", "darkorange"), bty = "n")
#legend(0.5, 19, legend = c("Autosomes","Z"), cex = 1,
#       fill = c("lightgrey", "darkorange"), bty = "n")
segments(6.5, -5.5, 6.5, 15, lty = 3, lwd = 1)
segments(13.5, -5.5, 13.5, 15, lty = 3, lwd = 1)
segments(1, -5.5, 19, -5.5, lty = 1, lwd = 1)



#############################################################################

##### SBGs
##### MWU-test of significant difference between Z and A for each sex ####

instar_V_f_mwu <- wilcox.test(instar_V_female_allSBG$FPKM[instar_V_female_allSBG$group=="1_instar_V_female_allSBG A"],
                              instar_V_female_allSBG$FPKM[instar_V_female_allSBG$group=="1_instar_V_female_allSBG Z"], paired = FALSE)

instar_V_m_mwu <- wilcox.test(instar_V_male_allSBG$FPKM[instar_V_male_allSBG$group=="2_instar_V_male_allSBG A"],
                              instar_V_male_allSBG$FPKM[instar_V_male_allSBG$group=="2_instar_V_male_allSBG Z"], paired = FALSE)

pupa_f_mwu <- wilcox.test(pupa_female_allSBG$FPKM[pupa_female_allSBG$group=="3_pupa_female_allSBG A"],
                          pupa_female_allSBG$FPKM[pupa_female_allSBG$group=="3_pupa_female_allSBG Z"], paired = FALSE)

pupa_m_mwu <- wilcox.test(pupa_male_allSBG$FPKM[pupa_male_allSBG$group=="4_pupa_male_allSBG A"],
                          pupa_male_allSBG$FPKM[pupa_male_allSBG$group=="4_pupa_male_allSBG Z"], paired = FALSE)                             

adult_f_mwu <- wilcox.test(adult_female_allSBG$FPKM[adult_female_allSBG$group=="5_adult_female_allSBG A"],
                           adult_female_allSBG$FPKM[adult_female_allSBG$group=="5_adult_female_allSBG Z"], paired = FALSE)

adult_m_mwu <- wilcox.test(adult_male_allSBG$FPKM[adult_male_allSBG$group=="6_adult_male_allSBG A"],
                           adult_male_allSBG$FPKM[adult_male_allSBG$group=="6_adult_male_allSBG Z"], paired = FALSE)


mwu <- matrix(c(instar_V_f_mwu$p.value, pupa_f_mwu$p.value, adult_f_mwu$p.value,
                instar_V_m_mwu$p.value, pupa_m_mwu$p.value, adult_m_mwu$p.value),
              nrow = 3, ncol = 2)
rownames(mwu) <- c("Instar V", "Pupa", "Adult")
colnames(mwu) <- c("Female - A vs Z", "Male - A vs Z")

setwd(out_dir) 
write.table(mwu, file = "mwu_test_p-values_boxplots_SBG.txt", sep = "\t", col.names = NA, row.names = TRUE)



############################################################################ PART II
############################################################################


####### upregulated genes computed as follows (see example below)
####### Example:
####### FBG_instar_v_female = genes in AllGenes_instar_v_female also listed in FBG_instar_v
#######
####### Note: AllGenes refers to gene list before removing SBG

# include only upregulated genes (FMG in females; MBG in males)
instar_V_f_upregulated <- instar_V_f_AllGenes[(instar_V_f_AllGenes$gene_id %in% FBG_Instar$gene_id ), ]
instar_V_m_upregulated <- instar_V_m_AllGenes[(instar_V_m_AllGenes$gene_id %in% MBG_Instar$gene_id), ]

pupa_f_upregulated <- pupa_f_AllGenes[(pupa_f_AllGenes$gene_id %in% FBG_Pupa$gene_id), ]
pupa_m_upregulated <- pupa_m_AllGenes[(pupa_m_AllGenes$gene_id %in% MBG_Pupa$gene_id), ]

adult_f_upregulated <- adult_f_AllGenes[(adult_f_AllGenes$gene_id %in% FBG_Adult$gene_id), ]
adult_m_upregulated <- adult_m_AllGenes[(adult_m_AllGenes$gene_id %in% MBG_Adult$gene_id),]



##################################### PROCESS DATA


instar_V_female_upregulated <- instar_V_f_upregulated[c(4,6)]
instar_V_female_upregulated$group <- rep("1_instar_V_female_upregulated", nrow(instar_V_female_upregulated))
instar_V_female_upregulated$group <- paste(instar_V_female_upregulated$group, instar_V_female_upregulated$chromosome)
instar_V_female_upregulated <- instar_V_female_upregulated[c(1,3)]
names(instar_V_female_upregulated)[1] <-"FPKM"

pupa_female_upregulated <- pupa_f_upregulated[c(4,6)]
pupa_female_upregulated$group <- rep("3_pupa_female_upregulated", nrow(pupa_female_upregulated))
pupa_female_upregulated$group <- paste(pupa_female_upregulated$group, pupa_female_upregulated$chromosome)
pupa_female_upregulated <- pupa_female_upregulated[c(1,3)]
names(pupa_female_upregulated)[1] <-"FPKM"

adult_female_upregulated <- adult_f_upregulated[c(4,6)]
adult_female_upregulated$group <- rep("5_adult_female_upregulated", nrow(adult_female_upregulated))
adult_female_upregulated$group <- paste(adult_female_upregulated$group, adult_female_upregulated$chromosome)
adult_female_upregulated <- adult_female_upregulated[c(1,3)]
names(adult_female_upregulated)[1] <-"FPKM"

instar_V_male_upregulated <- instar_V_m_upregulated[c(5,6)]
instar_V_male_upregulated$group <- rep("2_instar_V_male_upregulated", nrow(instar_V_male_upregulated))
instar_V_male_upregulated$group <- paste(instar_V_male_upregulated$group, instar_V_male_upregulated$chromosome)
instar_V_male_upregulated <- instar_V_male_upregulated[c(1,3)]
names(instar_V_male_upregulated)[1] <-"FPKM"

pupa_male_upregulated <- pupa_m_upregulated[c(5,6)]
pupa_male_upregulated$group <- rep("4_pupa_male_upregulated", nrow(pupa_male_upregulated))
pupa_male_upregulated$group <- paste(pupa_male_upregulated$group, pupa_male_upregulated$chromosome)
pupa_male_upregulated <- pupa_male_upregulated[c(1,3)]
names(pupa_male_upregulated)[1] <-"FPKM"

adult_male_upregulated <- adult_m_upregulated[c(5,6)]
adult_male_upregulated$group <- rep("6_adult_male_upregulated", nrow(adult_male_upregulated))
adult_male_upregulated$group <- paste(adult_male_upregulated$group, adult_male_upregulated$chromosome)
adult_male_upregulated <- adult_male_upregulated[c(1,3)]
names(adult_male_upregulated)[1] <-"FPKM"

upregulated <- rbind(instar_V_female_upregulated, instar_V_male_upregulated, pupa_female_upregulated, pupa_male_upregulated, adult_female_upregulated, adult_male_upregulated)



######################### BOXPLOTS



boxplot(log2(allSBG$FPKM)~allSBG$group, ylim = c(-6, 18),
        col = c("grey90", "grey90"), notch = FALSE,
        at = c(1,2, 4,5,   8,9, 11,12,  15,16, 18,19),
        outline = FALSE, boxwex = 0.7,  xaxt = "n", frame.plot = FALSE, border="grey65",
        boxlty = 2)
axis(2, labels = "log2 FPKM(>0)", cex.axis = 1.2, at = 4, line = 1.5, tck = 0)



boxplot(log2(upregulated$FPKM)~upregulated$group, ylim = c(-6, 18), add = TRUE, 
        col = c("lightgrey", "darkorange"), notch = FALSE,
        at = c(1.2,2.2, 4.2,5.2,   8.2,9.2, 11.2,12.2,  15.2,16.2, 18.2,19.2),
        outline = FALSE, boxwex = 0.7,  xaxt = "n", frame.plot = FALSE)

text(1.5, -6.5, "\\VE", vfont=c("sans serif","plain"), cex =1.9)
text(4.5, -6.5, "\\MA", vfont=c("sans serif","plain"), cex =1.9)

text(8.5, -6.5, "\\VE", vfont=c("sans serif","plain"), cex =1.9)
text(11.5, -6.5, "\\MA", vfont=c("sans serif","plain"), cex =1.9)

text(15.5, -6.5, "\\VE", vfont=c("sans serif","plain"), cex =1.9)
text(18.5, -6.5, "\\MA", vfont=c("sans serif","plain"), cex =1.9)

axis(1, labels = c("Larva", "Pupa", "Adult"), lwd = 0, cex.axis = 1.5,
     at = c(3, 10, 17), line = 1)
legend(16.5, 18, legend = c("Autosomes","Z"), cex = 1,
       fill = c("lightgrey", "darkorange"), bty = "n")
#legend(0.5, 19, legend = c("Autosomes","Z"), cex = 1,
#       fill = c("lightgrey", "darkorange"), bty = "n")
segments(6.5, -5.5, 6.5, 15, lty = 3, lwd = 1)
segments(13.5, -5.5, 13.5, 15, lty = 3, lwd = 1)
segments(1, -5.5, 19, -5.5, lty = 1, lwd = 1)


#############################################################################

##### Upregulated genes
##### MWU-test of significant difference between Z and A for each sex ####

instar_V_f_mwu_up <- wilcox.test(instar_V_female_upregulated$FPKM[instar_V_female_upregulated$group=="1_instar_V_female_upregulated A"],
                                 instar_V_female_upregulated$FPKM[instar_V_female_upregulated$group=="1_instar_V_female_upregulated Z"], paired = FALSE)

instar_V_m_mwu_up <- wilcox.test(instar_V_male_upregulated$FPKM[instar_V_male_upregulated$group=="2_instar_V_male_upregulated A"],
                                 instar_V_male_upregulated$FPKM[instar_V_male_upregulated$group=="2_instar_V_male_upregulated Z"], paired = FALSE)

pupa_f_mwu_up <- wilcox.test(pupa_female_upregulated$FPKM[pupa_female_upregulated$group=="3_pupa_female_upregulated A"],
                             pupa_female_upregulated$FPKM[pupa_female_upregulated$group=="3_pupa_female_upregulated Z"], paired = FALSE)

pupa_m_mwu_up <- wilcox.test(pupa_male_upregulated$FPKM[pupa_male_upregulated$group=="4_pupa_male_upregulated A"],
                             pupa_male_upregulated$FPKM[pupa_male_upregulated$group=="4_pupa_male_upregulated Z"], paired = FALSE)                             

adult_f_mwu_up <- wilcox.test(adult_female_upregulated$FPKM[adult_female_upregulated$group=="5_adult_female_upregulated A"],
                              adult_female_upregulated$FPKM[adult_female_upregulated$group=="5_adult_female_upregulated Z"], paired = FALSE)

adult_m_mwu_up <- wilcox.test(adult_male_upregulated$FPKM[adult_male_upregulated$group=="6_adult_male_upregulated A"],
                              adult_male_upregulated$FPKM[adult_male_upregulated$group=="6_adult_male_upregulated Z"], paired = FALSE)


mwu_up <- matrix(c(instar_V_f_mwu$p.value, pupa_f_mwu$p.value, adult_f_mwu$p.value,
                instar_V_m_mwu$p.value, pupa_m_mwu$p.value, adult_m_mwu$p.value),
              nrow = 3, ncol = 2)
rownames(mwu_up) <- c("Instar V", "Pupa", "Adult")
colnames(mwu_up) <- c("Female - A vs Z", "Male - A vs Z")

setwd(out_dir) 
write.table(mwu_up, file = "mwu_test_p-values_boxplots_upregulated.txt", sep = "\t", col.names = NA, row.names = TRUE)



#############################################################################

############## Additional MWU-tests

# Adult_female_A_SBG vs Adult_female_Z_FBG
Adult_female_A_SBG_vs_Z_FBG <- wilcox.test(adult_female_allSBG$FPKM[adult_female_allSBG$group=="5_adult_female_allSBG A"],
                                 adult_female_upregulated$FPKM[adult_female_upregulated$group=="5_adult_female_upregulated Z"], paired = FALSE)
Adult_female_A_SBG_vs_Z_FBG


# Adult_female_A_nonbiased vs Adult_female_Z_SBG
Adult_female_A_nonbiased_vs_Z_SBG <- wilcox.test(adult_female$FPKM[adult_female$group=="5_adult_female A"],
                                                 adult_female_allSBG$FPKM[adult_female_allSBG$group=="5_adult_female_allSBG Z"], paired = FALSE)
Adult_female_A_nonbiased_vs_Z_SBG



