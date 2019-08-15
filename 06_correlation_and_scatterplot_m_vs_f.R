# Correlation between male vs female expression levels

##Install LSD 
#source("http://bioconductor.org/biocLite.R")
#biocLite("LSD")
#biocLite("RColorBrewer")

library(LSD)
library(RColorBrewer)

#Clear all states
rm(list=ls(all=TRUE))

# Input folder
setwd("/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/03_Normalized_libs")


# Load data
instar_V <- read.delim("instar_V-assigned_A_or_Z-filtered.txt")
instar_V_a <- rbind(instar_V[instar_V$chromosome=="A", ])
instar_V_z <- rbind(instar_V[instar_V$chromosome=="Z", ])

pupa <- read.delim("pupa-assigned_A_or_Z-filtered.txt")
pupa_a <- rbind(pupa[pupa$chromosome=="A", ])
pupa_z <- rbind(pupa[pupa$chromosome=="Z", ])

adult <- read.delim("adult-assigned_A_or_Z-filtered.txt", header = TRUE)
adult_a <- rbind(adult[adult$chromosome=="A", ])
adult_z <- rbind(adult[adult$chromosome=="Z", ])

instar_V_female_a <- log2(instar_V_a$FPKM_instar_V_female)
instar_V_female_z <- log2(instar_V_z$FPKM_instar_V_female)
instar_V_male_a <- log2(instar_V_a$FPKM_instar_V_male)
instar_V_male_z <- log2(instar_V_z$FPKM_instar_V_male)

pupa_female_a <- log2(pupa_a$FPKM_pupa_female)
pupa_female_z <- log2(pupa_z$FPKM_pupa_female)
pupa_male_a <- log2(pupa_a$FPKM_pupa_male)
pupa_male_z <- log2(pupa_z$FPKM_pupa_male)

adult_female_a <- log2(adult_a$FPKM_adult_female)
adult_female_z <- log2(adult_z$FPKM_adult_female)
adult_male_a <- log2(adult_a$FPKM_adult_male)
adult_male_z <- log2(adult_z$FPKM_adult_male)

# Analysis
cor.test(instar_V_female_a, instar_V_male_a, alternative = "two.sided", method = "spearman")
cor.test(instar_V_female_z, instar_V_male_z, alternative = "two.sided", method = "spearman")
cor.test(pupa_female_a, pupa_male_a, alternative = "two.sided", method = "spearman")
cor.test(pupa_female_z, pupa_male_z, alternative = "two.sided", method = "spearman")
cor.test(adult_female_a, adult_male_a, alternative = "two.sided", method = "spearman")
cor.test(adult_female_z, adult_male_z, alternative = "two.sided", method = "spearman")



#Scatterplots####
mar.default <- c(5,4,4,2) + 0.1
par(pty="s", mfrow = c(1,3), oma=c(3,3,3,4),mar = mar.default + c(0, 0.5, 0, 0))


# Instar V ####

colscale <- brewer.pal(9, "Greys")

heatscatter(instar_V_male_a, instar_V_female_a, asp = TRUE, pch = 16, cexplot = 0.8,
      nrcol = 5, grid = 1000, colpal = colscale[4:8],
      ylim = c(-6,20), xlim = c(-5,20), main = "", xlab = "",
      ylab = "Female log2 FPKM(>0)", cex.lab = 1.5)
title(main = "Larva", cex.main = 2)
points(instar_V_male_z, instar_V_female_z, pch = 21, bg = "darkorange", col = "black", cex = 0.8)

lines(lowess(instar_V_male_a, instar_V_female_a, f = 1/3, iter = 100), col="dimgrey", lwd = 2) # lowess line (x,y) 
lines(lowess(instar_V_male_z, instar_V_female_z, f = 1/3, iter = 100), col="orange", lwd = 2)
abline(v = c(-10, -5, 0, 5, 10, 15), lty = 3)
abline(h = c(-10, -5, 0, 5, 10, 15), lty = 3)

# Pupa ####

heatscatter(pupa_male_a, pupa_female_a, asp = TRUE, pch = 16, cexplot = 0.8,
     nrcol = 5, grid = 1000, colpal = colscale[4:8],
     ylim = c(-5,20), xlim = c(-5,20), main = "",
     xlab = "", ylab = "", cex.lab = 2)
title(main = "Pupa", cex.main = 2)
points(pupa_male_z, pupa_female_z, col = "black", pch = 21, bg = "darkorange", cex = 0.8)

lines(lowess(pupa_male_a, pupa_female_a, f = 1/3, iter = 100), col="dimgrey", lwd = 2) # lowess line (x,y) 
lines(lowess(pupa_male_z, pupa_female_z, f = 1/3, iter = 100), col="orange", lwd = 2)
abline(v = c(-10, -5, 0, 5, 10, 15), lty = 3)
abline(h = c(-10, -5, 0, 5, 10, 15), lty = 3)

axis(1, labels = "Male log2 FPKM(>0)", cex.axis = 1.5, at = 7.5, line = 4, tck = 0)

# Adult ####

heatscatter(adult_male_a, adult_female_a, asp = TRUE, pch = 16, cex = 0.8,
    nrcol = 5, grid = 1000, colpal = colscale[4:8],
    ylim = c(-6, 20), xlim = c(-5,20), main = "",
    xlab = "", ylab = "", cex.lab = 1.5)
title(main = "Adult", cex.main = 2)
points(adult_male_z, adult_female_z, col = "black", pch = 21, bg = "darkorange", cex = 0.8)

lines(lowess(adult_male_a, adult_female_a, f = 1/3, iter = 100), col="dimgrey", lwd = 2) # lowess line (x,y) 
lines(lowess(adult_male_z, adult_female_z, f = 1/3, iter = 100), col="orange", lwd = 2)
abline(v = c(-10, -5, 0, 5, 10, 15), lty = 3)
abline(h = c(-10, -5, 0, 5, 10, 15), lty = 3)

legend(-6, -14, legend = "Autosomes", cex = 1.5, xpd = TRUE,
       col = "grey", bty = "n", pch = 19)
legend(10, -14, legend = "Z", cex = 1.5, xpd = TRUE,
       col = "darkorange", bty = "n", pch = 19)

