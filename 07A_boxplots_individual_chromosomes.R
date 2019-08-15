#Boxplots of individual chromosomes

#Clear all states
rm(list=ls(all=TRUE))
dev.off()


##### Indicate whether the analysis is to be performed before (0) or after (1) removing sex-biased genes
FLAG_BG <- 1


############ Paths and folders

##### Input folders
if (FLAG_BG == 0) {
  ##### folder containing datasets before removing biased genes 
  IN_FOLD <- "/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/03_Normalized_libs"
} else if (FLAG_BG == 1) {
  ##### folder containing datasets after removing biased genes (all biased genes, either FBG or MBG)
  IN_FOLD <- "/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/17_Assign_to_chromosomes_m_vs_f_SBG"     
}



############ Load tables of gene expression as mean FPKM-values for each group ####

if (FLAG_BG == 0) {
  
  ### analysis of datasets produced before removing SBG
  setwd(IN_FOLD)
  
  ### genes with zero counts removed individually for each sex
  instar_V_f <- read.delim("instar_V-assigned_to_chromosomes_female-filtered.txt", header = TRUE)
  instar_V_m <- read.delim("instar_V-assigned_to_chromosomes_male-filtered.txt", header = TRUE)
  
  pupa_f <- read.delim("pupa-assigned_to_chromosomes_female-filtered.txt", header = TRUE)
  pupa_m <- read.delim("pupa-assigned_to_chromosomes_male-filtered.txt", header = TRUE)
  
  adult_f <- read.delim("adult-assigned_to_chromosomes_female-filtered.txt", header = TRUE)
  adult_m <- read.delim("adult-assigned_to_chromosomes_male-filtered.txt", header = TRUE)
  
  
} else if (FLAG_BG == 1) {
  
  ### analysis of datasets produced after removing SBG
  setwd(IN_FOLD)
  
  
  ### genes with zero counts removed individually for each sex
  instar_V_f <- read.delim("nonbiased_genes-instar_V-assigned_to_chromosomes_female-filtered.txt", header = TRUE)
  instar_V_m <- read.delim("nonbiased_genes-instar_V-assigned_to_chromosomes_male-filtered.txt", header = TRUE)
  
  pupa_f <- read.delim("nonbiased_genes-pupa-assigned_to_chromosomes_female-filtered.txt", header = TRUE)
  pupa_m <- read.delim("nonbiased_genes-pupa-assigned_to_chromosomes_male-filtered.txt", header = TRUE)
  
  adult_f <- read.delim("nonbiased_genes-adult-assigned_to_chromosomes_female-filtered.txt", header = TRUE)
  adult_m <- read.delim("nonbiased_genes-adult-assigned_to_chromosomes_male-filtered.txt", header = TRUE)
  
}






######################################################################## Analysis

instar_V_female <- instar_V_f[c(4,6)]
instar_V_female_z <- rbind(instar_V_female[instar_V_female$chromosome=="0", ])
instar_V_female_z_median <- median(instar_V_female_z$FPKM_instar_V_female)

instar_V_male <- instar_V_m[c(5,6)]
instar_V_male_z <- rbind(instar_V_male[instar_V_male$chromosome=="0", ])
instar_V_male_z_median <- median(instar_V_male_z$FPKM_instar_V_male)

pupa_female <- pupa_f[c(4,6)]
pupa_female_z <- rbind(pupa_female[pupa_female$chromosome=="0", ])
pupa_female_z_median <- median(pupa_female_z$FPKM_pupa_female)

pupa_male <- pupa_m[c(5,6)]
pupa_male_z <- rbind(pupa_male[pupa_male$chromosome=="0", ])
pupa_male_z_median <- median(pupa_male_z$FPKM_pupa_male)

adult_female <- adult_f[c(4,6)]
adult_female_z <- rbind(adult_female[adult_female$chromosome=="0", ])
adult_female_z_median <- median(adult_female_z$FPKM_adult_female)

adult_male <- adult_m[c(5,6)]
adult_male_z <- rbind(adult_male[adult_male$chromosome=="0", ])
adult_male_z_median <- median(adult_male_z$FPKM_adult_male)



#################################################################### Boxplots

### Female

par(pty="s", mfrow = c(1,3))

boxplot(log2(instar_V_female$FPKM_instar_V_female)~instar_V_female$chromosome, ylim = c(-5, 15), notch = TRUE,
        col = c(rep("darkorange",1),rep("grey",20)), 
        outline = FALSE, xaxt = "n", boxwex = 0.7, frame.plot = FALSE,
        main = "Larva", cex.main = 2)
#axis(3, labels = "A", cex.axis = 3, at = 0, tck = 0)
axis(2, labels = "log2 FPKM(>0)", cex.axis = 1.5, at = 3, line = 1.5, tck = 0)
axis(1, las=2, cex.axis = 1.1, labels = c("Z", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"),
     at = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21))

segments(1, log2(instar_V_female_z_median), 21.5, log2(instar_V_female_z_median), col = "darkorange", lty = 2, lwd = 2)

boxplot(log2(pupa_female$FPKM_pupa_female)~pupa_female$chromosome, ylim = c(-5, 15), notch = TRUE,
        col = c(rep("darkorange",1),rep("grey",20)), 
        outline = FALSE, xaxt = "n", boxwex = 0.7, frame.plot = FALSE,
        main = "Pupa", cex.main = 2, xlab = "Female chromosomes", cex.lab=2)
axis(1, las=2, cex.axis = 1.1, labels = c("Z", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"),
     at = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21))

segments(1, log2(pupa_female_z_median), 21.5, log2(pupa_female_z_median), col = "darkorange", lty = 2, lwd = 2)

#axis(1, labels = c("Female chromosomes"), lwd = 0, at = 11.2, cex.axis = 2.0, tck = 0, line = 4)


boxplot(log2(adult_female$FPKM_adult_female)~adult_female$chromosome, ylim = c(-5, 15), notch = TRUE,
        col = c(rep("darkorange",1),rep("grey",20)), 
        outline = FALSE, xaxt = "n", boxwex = 0.7, frame.plot = FALSE,
        main = "Adult", cex.main = 2)
axis(1, las=2, cex.axis = 1.1, labels = c("Z", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"),
     at = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21))

segments(1, log2(adult_female_z_median), 21.5, log2(adult_female_z_median), col = "darkorange", lty = 2, lwd = 2)



### Male

par(pty="s", mfrow = c(1,3))

boxplot(log2(instar_V_male$FPKM_instar_V_male)~instar_V_male$chromosome, ylim = c(-5, 15), notch = TRUE,
        col = c(rep("darkorange",1),rep("grey",20)), 
        outline = FALSE, xaxt = "n", boxwex = 0.7, frame.plot = FALSE,
        main = "Larva", cex.main = 2)
axis(2, labels = "log2 FPKM(>0)", cex.axis = 1.5, at = 3, line = 1.5, tck = 0)
#axis(3, labels = "B", cex.axis = 3, at = 0, tck = 0)
axis(1, las = 2, cex.axis = 1.1, labels = c("Z", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"),
     at = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21))

segments(1, log2(instar_V_male_z_median), 21.5, log2(instar_V_male_z_median), col = "darkorange", lty = 2, lwd = 2)

boxplot(log2(pupa_male$FPKM_pupa_male)~pupa_male$chromosome, ylim = c(-5, 15), notch = TRUE,
        col = c(rep("darkorange",1),rep("grey",20)), 
        outline = FALSE, xaxt = "n", boxwex = 0.7, frame.plot = FALSE,
        main = "Pupa", cex.main = 2, xlab = "Male chromosomes", cex.lab=2)
axis(1, las=2, cex.axis = 1.1, labels = c("Z", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"),
     at = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21))

segments(1, log2(pupa_male_z_median), 21.5, log2(pupa_male_z_median), col = "darkorange", lty = 2, lwd = 2)

#axis(1, labels = c("Male chromosomes"), lwd = 0, at = 11.2, cex.axis = 2.0, tck = 0, line = 4)

boxplot(log2(adult_male$FPKM_adult_male)~adult_male$chromosome, ylim = c(-5, 15), notch = TRUE,
        col = c(rep("darkorange",1),rep("grey",20)), 
        outline = FALSE, xaxt = "n", boxwex = 0.7, frame.plot = FALSE,
        main = "Adult", cex.main = 2)
axis(1, las=2, cex.axis = 1.1, labels = c("Z", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"),
     at = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21))

segments(1, log2(adult_male_z_median), 21.5, log2(adult_male_z_median), col = "darkorange", lty = 2, lwd = 2)














