# Boxplots A vs Z, before removing sex-biased genes (SBG)

#Clear all states
rm(list=ls(all=TRUE))
dev.off()





##### Input folder
folder_in <- "/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/03_Normalized_libs"           




##################################### Load tables

setwd(folder_in) 

### genes with zero counts removed individually for each sex
instar_V_f <- read.delim("instar_V-assigned_A_or_Z_female-filtered.txt", header = TRUE)
instar_V_m <- read.delim("instar_V-assigned_A_or_Z_male-filtered.txt", header = TRUE)

pupa_f <- read.delim("pupa-assigned_A_or_Z_female-filtered.txt", header = TRUE)
pupa_m <- read.delim("pupa-assigned_A_or_Z_male-filtered.txt", header = TRUE)

adult_f <- read.delim("adult-assigned_A_or_Z_female-filtered.txt", header = TRUE)
adult_m <- read.delim("adult-assigned_A_or_Z_male-filtered.txt", header = TRUE)











##################################### PROCESS DATA

instar_V_female <- instar_V_f[c(4,6)]
instar_V_female$group <- rep("1_instar_V_female", nrow(instar_V_female))
instar_V_female$group <- paste(instar_V_female$group, instar_V_female$chromosome)
instar_V_female <- instar_V_female[c(1,3)]
names(instar_V_female)[1] <-"FPKM"

pupa_female <- pupa_f[c(4,6)]
pupa_female$group <- rep("3_pupa_female", nrow(pupa_female))
pupa_female$group <- paste(pupa_female$group, pupa_female$chromosome)
pupa_female <- pupa_female[c(1,3)]
names(pupa_female)[1] <-"FPKM"

adult_female <- adult_f[c(4,6)]
adult_female$group <- rep("5_adult_female", nrow(adult_female))
adult_female$group <- paste(adult_female$group, adult_female$chromosome)
adult_female <- adult_female[c(1,3)]
names(adult_female)[1] <-"FPKM"

instar_V_male <- instar_V_m[c(5,6)]
instar_V_male$group <- rep("2_instar_V_male", nrow(instar_V_male))
instar_V_male$group <- paste(instar_V_male$group, instar_V_male$chromosome)
instar_V_male <- instar_V_male[c(1,3)]
names(instar_V_male)[1] <-"FPKM"

pupa_male <- pupa_m[c(5,6)]
pupa_male$group <- rep("4_pupa_male", nrow(pupa_male))
pupa_male$group <- paste(pupa_male$group, pupa_male$chromosome)
pupa_male <- pupa_male[c(1,3)]
names(pupa_male)[1] <-"FPKM"

adult_male <- adult_m[c(5,6)]
adult_male$group <- rep("6_adult_male", nrow(adult_male))
adult_male$group <- paste(adult_male$group, adult_male$chromosome)
adult_male <- adult_male[c(1,3)]
names(adult_male)[1] <-"FPKM"

all_samples <- rbind(instar_V_female, instar_V_male, pupa_female, pupa_male, adult_female, adult_male)






######################### BOXPLOTS


boxplot(log2(all_samples$FPKM)~all_samples$group, ylim = c(-6, 18),
        col = c("lightgrey", "darkorange"), notch = TRUE,
        at = c(1,2, 4,5,   8,9, 11,12,  15,16, 18,19),
        outline = FALSE, boxwex = 0.7,  xaxt = "n", frame.plot = FALSE)
axis(2, labels = "log2 FPKM(>0)", cex.axis = 1.2, at = 4, line = 1.5, tck = 0)

text(1.5, -6.5, "\\VE", vfont=c("sans serif","plain"), cex =1.8)
text(4.5, -6.5, "\\MA", vfont=c("sans serif","plain"), cex =1.8)

text(8.5, -6.5, "\\VE", vfont=c("sans serif","plain"), cex =1.8)
text(11.5, -6.5, "\\MA", vfont=c("sans serif","plain"), cex =1.8)

text(15.5, -6.5, "\\VE", vfont=c("sans serif","plain"), cex =1.8)
text(18.5, -6.5, "\\MA", vfont=c("sans serif","plain"), cex =1.8)

axis(1, labels = c("Larva", "Pupa", "Adult"), lwd = 0, cex.axis = 1.5,
     at = c(3, 10, 17), line = 1)
legend(16.5, 18, legend = c("Autosomes","Z"), cex = 1,
       fill = c("lightgrey", "darkorange"), bty = "n")
segments(6.5, -5.5, 6.5, 15, lty = 3, lwd = 1)
segments(13.5, -5.5, 13.5, 15, lty = 3, lwd = 1)
segments(1, -5.5, 19, -5.5, lty = 1, lwd = 1)





