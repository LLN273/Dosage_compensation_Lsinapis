#Density plots m:f ratio


#Clear all states
rm(list=ls(all=TRUE))
dev.off()


## Input folder
setwd("/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/03_Normalized_libs")


## Load data
instar_V <- read.delim("instar_V-assigned_A_or_Z-filtered.txt", header = TRUE)
pupa <- read.delim("pupa-assigned_A_or_Z-filtered.txt", header = TRUE)
adult <- read.delim("adult-assigned_A_or_Z-filtered.txt", header = TRUE)


## Analysis
instar_V$ratio <- instar_V$FPKM_instar_V_male/instar_V$FPKM_instar_V_female
pupa$ratio <- pupa$FPKM_pupa_male/pupa$FPKM_pupa_female
adult$ratio <- adult$FPKM_adult_male/adult$FPKM_adult_female

instar_V_a <- instar_V$ratio[instar_V$chromosome=="A"]
instar_V_z <- instar_V$ratio[instar_V$chromosome=="Z"]
pupa_a <- pupa$ratio[pupa$chromosome=="A"]
pupa_z <- pupa$ratio[pupa$chromosome=="Z"]
adult_a <- adult$ratio[adult$chromosome=="A"]
adult_z <- adult$ratio[adult$chromosome=="Z"]


instar_V_a_median <- median(instar_V_a)
instar_V_z_median <- median(instar_V_z)
mwu_instar_V <- wilcox.test(instar_V_a, instar_V_z, paired = FALSE)
print(mwu_instar_V)

pupa_a_median <- median(pupa_a)
pupa_z_median <- median(pupa_z)
mwu_pupa <- wilcox.test(pupa_a, pupa_z, paired = FALSE)
print(mwu_pupa)

adult_a_median <- median(adult_a)
adult_z_median <- median(adult_z)
mwu_adult <- wilcox.test(adult_a, adult_z, paired = FALSE)
print(mwu_adult)

# Plot ####

par(pty="s", mfrow = c(1,3))

plot_instar_V_a <- density(log2(instar_V_a))
plot(plot_instar_V_a, lwd = 3, ylim = c(0, 0.7), xlim = c(-10, 10), col = "grey", 
     xlab = "", main = "Larva", cex.main = 2, cex.lab = 1.5,
     ylab = "", frame.plot = FALSE)
axis(2, labels = "Density", cex.axis = 1.5, at = 0.35, line = 1.5, tck = 0)
polygon(plot_instar_V_a, border = NA, col = rgb(211, 211, 211, maxColorValue = 255, alpha = 100))
plot_instar_V_z <- density(log2(instar_V_z))
lines(plot_instar_V_z, lwd = 3, col = "darkorange")
polygon(plot_instar_V_z, border = NA, col = rgb(255, 140, 0, maxColorValue = 255, alpha = 25))

segments(log2(instar_V_a_median), 0, log2(instar_V_a_median), 0.68, lty = 2, col = "grey", lwd = 2)
segments(log2(instar_V_z_median), 0, log2(instar_V_z_median), 0.45, lty = 2, col = "darkorange", lwd = 2)

plot_pupa_a <- density(log2(pupa_a))
plot(plot_pupa_a, lwd = 3, ylim = c(0, 0.7), xlim = c(-10, 10), col = "grey", 
     xlab = "", main = "Pupa", cex.main = 2, cex.lab = 1.5,
     ylab = "", frame.plot = FALSE)
polygon(plot_pupa_a,  border = NA, col = rgb(211, 211, 211, maxColorValue = 255, alpha = 100))
plot_pupa_z <- density(log2(pupa_z))
lines(plot_pupa_z, lwd = 3, col = "darkorange")
polygon(plot_pupa_z, border = NA, col = rgb(255, 140, 0, maxColorValue = 255, alpha = 25))

axis(1, labels = "log2 male:female FPKM-ratio", cex.axis = 1.5, at = 0, line = 4, tck = 0)

segments(log2(pupa_a_median), 0, log2(pupa_a_median), 0.45, lty = 2, col = "grey", lwd = 2)
segments(log2(pupa_z_median), 0, log2(pupa_z_median), 0.35, lty = 2, col = "darkorange", lwd = 2)

plot_adult_a <- density(log2(adult_a))
plot(plot_adult_a, lwd = 3, ylim = c(0, 0.7), xlim = c(-10, 10), col = "grey", 
     xlab = "", main = "Adult", cex.main = 2, cex.lab = 1.5,
     ylab = "", frame.plot = FALSE)
polygon(plot_adult_a, border = NA, col = rgb(211, 211, 211, maxColorValue = 255, alpha = 100))
plot_adult_z <- density(log2(adult_z))
lines(plot_adult_z, lwd = 3, col = "darkorange")
polygon(plot_adult_z, border = NA, col = rgb(255, 140, 0, maxColorValue = 255, alpha = 25))

segments(log2(adult_a_median), 0, log2(adult_a_median), 0.29, lty = 2, col = "grey", lwd = 2)
segments(log2(adult_z_median), 0, log2(adult_z_median), 0.32, lty = 2, col = "darkorange", lwd = 2)

text(x=0.370, "*", y=0.35, pos=3, cex=2.5)
text(x=0.370, y = 0.36, '{', srt = -90, cex = 1.5)

legend(-12, -0.2, legend = "Autosomes", cex = 1.5, xpd = TRUE,
       col = "grey", bty = "n", lty = 1, lwd = 3)
legend(3, -0.2, legend = "Z", cex = 1.5, xpd = TRUE,
       col = "darkorange", bty = "n", lty = 1, lwd = 3)

