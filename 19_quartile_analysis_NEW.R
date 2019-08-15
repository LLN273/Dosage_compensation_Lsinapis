#Quartile MWU


#Clear all states
rm(list=ls(all=TRUE))
dev.off()


##### Indicate whether the analysis is to be performed before (0) or after (1) removing sex-biased genes
FLAG_BG <- 1




############ Paths and folders
##### Input folders
if (FLAG_BG == 0) {
  ##### folder containing datasets before removing biased genes 
  IN_FOLDER <- "/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/18A_Quartile_analysis"
} else if (FLAG_BG == 1) {
  ##### folder containing datasets after removing biased genes (all biased genes, either FBG or MBG)
  IN_FOLDER <- "/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/18B_Quartile_analysis_AFTER_REMOVE_SEX_BIASED_GENES"     
}



##### Results files
if (FLAG_BG == 0) {
  table_1_out = "quartile_analysis_p-values.txt"
} else if (FLAG_BG == 1) {
  table_1_out = "quartile_analysis_nonbiased_p-values.txt"
}






############ Load tables of gene expression as mean FPKM-values for each group ####

if (FLAG_BG == 0) {
  
  ### analysis of datasets produced before removing SBG
  setwd(IN_FOLDER)
  
  instar_V <- read.delim("quartiles-sorted-instar_V-assigned_A_or_Z-filtered.txt")
  pupa <- read.delim("quartiles-sorted-pupa-assigned_A_or_Z-filtered.txt")
  adult <- read.delim("quartiles-sorted-adult-assigned_A_or_Z-filtered.txt")
  
  
} else if (FLAG_BG == 1) {
  
  ### analysis of datasets produced after removing SBG
  setwd(IN_FOLDER)
  
  instar_V <- read.delim("quartiles-sorted-nonbiased_genes-instar_V-assigned_A_or_Z-filtered.txt")
  pupa <- read.delim("quartiles-sorted-nonbiased_genes-pupa-assigned_A_or_Z-filtered.txt")
  adult <- read.delim("quartiles-sorted-nonbiased_genes-adult-assigned_A_or_Z-filtered.txt")
  
}





####

#instar_V <- instar_V[grep("[[:digit:]]", instar_V$FPKM_instar_V_female), ]
#pupa <- pupa[grep("[[:digit:]]", pupa$FPKM_pupa_female), ]
#adult <- adult[grep("[[:digit:]]", adult$FPKM_adult_female), ]

#options(digits=10)
#instar_V$FPKM_instar_V_female  <- as.numeric(as.character(instar_V$FPKM_instar_V_female))
#instar_V$FPKM_instar_V_male  <- as.numeric(as.character(instar_V$FPKM_instar_V_male))
#pupa$FPKM_pupa_female  <- as.numeric(as.character(pupa$FPKM_pupa_female))
#pupa$FPKM_pupa_male  <- as.numeric(as.character(pupa$FPKM_pupa_male))
#adult$FPKM_adult_female  <- as.numeric(as.character(adult$FPKM_adult_female))
#adult$FPKM_adult_male  <- as.numeric(as.character(adult$FPKM_adult_male))

####


instar_V_female <- instar_V[c(1,2,4)]
instar_V_female$group <- rep("instar_V_female", nrow(instar_V_female))
instar_V_female$group <- paste(instar_V_female$quartile, instar_V_female$group)
instar_V_female <- instar_V_female[c(1,2,4)]
names(instar_V_female)[2] <- "FPKM"

instar_V_male <- instar_V[c(1,3,4)]
instar_V_male$group <- rep("instar_V_male", nrow(instar_V_male))
instar_V_male$group <- paste(instar_V_male$quartile, instar_V_male$group)
instar_V_male <- instar_V_male[c(1,2,4)]
names(instar_V_male)[2] <- "FPKM"

instar_V_plot <- rbind(instar_V_female, instar_V_male)

pupa_female <- pupa[c(1,2,4)]
pupa_female$group <- rep("pupa_female", nrow(pupa_female))
pupa_female$group <- paste(pupa_female$quartile, pupa_female$group)
pupa_female <- pupa_female[c(1,2,4)]
names(pupa_female)[2] <- "FPKM"

pupa_male <- pupa[c(1,3,4)]
pupa_male$group <- rep("pupa_male", nrow(pupa_male))
pupa_male$group <- paste(pupa_male$quartile, pupa_male$group)
pupa_male <- pupa_male[c(1,2,4)]
names(pupa_male)[2] <- "FPKM"

pupa_plot <- rbind(pupa_female, pupa_male)

adult_female <- adult[c(1,2,4)]
adult_female$group <- rep("adult_female", nrow(adult_female))
adult_female$group <- paste(adult_female$quartile, adult_female$group)
adult_female <- adult_female[c(1,2,4)]
names(adult_female)[2] <- "FPKM"

adult_male <- adult[c(1,3,4)]
adult_male$group <- rep("adult_male", nrow(adult_male))
adult_male$group <- paste(adult_male$quartile, adult_male$group)
adult_male <- adult_male[c(1,2,4)]
names(adult_male)[2] <- "FPKM"

adult_plot <- rbind(adult_female, adult_male)


q1_instar_V_mwu <- wilcox.test(instar_V_female$FPKM[instar_V_female$group=="q_1 instar_V_female"],
                               instar_V_male$FPKM[instar_V_male$group=="q_1 instar_V_male"], paired = FALSE)
q2_instar_V_mwu <- wilcox.test(instar_V_female$FPKM[instar_V_female$group=="q_2 instar_V_female"],
                               instar_V_male$FPKM[instar_V_male$group=="q_2 instar_V_male"], paired = FALSE)
q3_instar_V_mwu <- wilcox.test(instar_V_female$FPKM[instar_V_female$group=="q_3 instar_V_female"],
                               instar_V_male$FPKM[instar_V_male$group=="q_3 instar_V_male"], paired = FALSE)
q4_instar_V_mwu <- wilcox.test(instar_V_female$FPKM[instar_V_female$group=="q_4 instar_V_female"],
                               instar_V_male$FPKM[instar_V_male$group=="q_4 instar_V_male"], paired = FALSE)

q1_pupa_mwu <- wilcox.test(pupa_female$FPKM[pupa_female$group=="q_1 pupa_female"],
                           pupa_male$FPKM[pupa_male$group=="q_1 pupa_male"], paired = FALSE)
q2_pupa_mwu <- wilcox.test(pupa_female$FPKM[pupa_female$group=="q_2 pupa_female"],
                           pupa_male$FPKM[pupa_male$group=="q_2 pupa_male"], paired = FALSE)
q3_pupa_mwu <- wilcox.test(pupa_female$FPKM[pupa_female$group=="q_3 pupa_female"],
                           pupa_male$FPKM[pupa_male$group=="q_3 pupa_male"], paired = FALSE)
q4_pupa_mwu <- wilcox.test(pupa_female$FPKM[pupa_female$group=="q_4 pupa_female"],
                           pupa_male$FPKM[pupa_male$group=="q_4 pupa_male"], paired = FALSE)

q1_adult_mwu <- wilcox.test(adult_female$FPKM[adult_female$group=="q_1 adult_female"],
                            adult_male$FPKM[adult_male$group=="q_1 adult_male"], paired = FALSE)
q2_adult_mwu <- wilcox.test(adult_female$FPKM[adult_female$group=="q_2 adult_female"],
                            adult_male$FPKM[adult_male$group=="q_2 adult_male"], paired = FALSE)
q3_adult_mwu <- wilcox.test(adult_female$FPKM[adult_female$group=="q_3 adult_female"],
                            adult_male$FPKM[adult_male$group=="q_3 adult_male"], paired = FALSE)
q4_adult_mwu <- wilcox.test(adult_female$FPKM[adult_female$group=="q_4 adult_female"],
                            adult_male$FPKM[adult_male$group=="q_4 adult_male"], paired = FALSE)

instar_V_q1_N <- sum(instar_V$quartile=="q_1")
instar_V_q2_N <- sum(instar_V$quartile=="q_2")
instar_V_q3_N <- sum(instar_V$quartile=="q_3")
instar_V_q4_N <- sum(instar_V$quartile=="q_4")
pupa_q1_N <- sum(pupa$quartile=="q_1")
pupa_q2_N <- sum(pupa$quartile=="q_2")
pupa_q3_N <- sum(pupa$quartile=="q_3")
pupa_q4_N <- sum(pupa$quartile=="q_4")
adult_q1_N <- sum(adult$quartile=="q_1")
adult_q2_N <- sum(adult$quartile=="q_2")
adult_q3_N <- sum(adult$quartile=="q_3")
adult_q4_N <- sum(adult$quartile=="q_4")

quartile_table <- matrix(c(instar_V_q1_N, pupa_q1_N, adult_q1_N,
                           q1_instar_V_mwu$p.value, q1_pupa_mwu$p.value, q1_adult_mwu$p.value,
                           instar_V_q2_N, pupa_q2_N, adult_q2_N,
                           q2_instar_V_mwu$p.value, q2_pupa_mwu$p.value, q2_adult_mwu$p.value,
                           instar_V_q3_N, pupa_q3_N, adult_q3_N,
                           q3_instar_V_mwu$p.value, q3_pupa_mwu$p.value, q3_adult_mwu$p.value,
                           instar_V_q4_N, pupa_q4_N, adult_q4_N,
                           q4_instar_V_mwu$p.value, q4_pupa_mwu$p.value, q4_adult_mwu$p.value),
                         nrow = 3, ncol = 8)

row.names(quartile_table) <- c("Instar V", "Pupa", "Adult")
colnames(quartile_table) <- c("N Q1", "p-value", "N Q2", "p-value", "N Q3", "p-value", "N Q4", "p-value")
write.table(quartile_table, file = table_1_out, sep = "\t", col.names = NA, row.names = TRUE)




# Combined plots ####

par(pty="s", mfrow = c(1,3))

boxplot(log2(instar_V_plot$FPKM)~instar_V_plot$group, ylim = c(-3, 15), notch = TRUE,
        at = c(1,2, 4,5, 7,8, 10,11), col = c("red", "blue"), xaxt = "n", outline = FALSE,
        boxwex = 0.7, frame.plot = FALSE, main = "Larva", cex.main = 2)
axis(2, labels = "log2 FPKM(>0)", cex.axis = 1.5, at = 6, line = 1.5, tck = 0)
axis(1, labels = c("Q1", "Q2", "Q3", "Q4"), cex.axis = 1.2, at = c(1.5, 4.5, 7.5, 10.5))
#text(x=4.5, y=6, "**", pos=3, cex=2)
text(x=4.5, y=6, "*", pos=3, cex=2)


boxplot(log2(pupa_plot$FPKM)~pupa_plot$group, ylim = c(-3, 15), notch = TRUE,
        at = c(1,2, 4,5, 7,8, 10,11), col = c("red", "blue"), xaxt = "n", outline = FALSE,
        boxwex = 0.7, frame.plot = FALSE, main = "Pupa", cex.main = 2)
axis(1, labels = c("Q1", "Q2", "Q3", "Q4"), cex.axis = 1.2, at = c(1.5, 4.5, 7.5, 10.5))
axis(1, labels = "Z chromosome quartiles", cex.axis = 1.5, at = 6, line = 4, tck = "0")


boxplot(log2(adult_plot$FPKM)~adult_plot$group, ylim = c(-3, 15), notch = TRUE,
        at = c(1,2, 4,5, 7,8, 10,11), col = c("red", "blue"), xaxt = "n", outline = FALSE,
        boxwex = 0.7, frame.plot = FALSE, main = "Adult", cex.main = 2)
axis(1, labels = c("Q1", "Q2", "Q3", "Q4"), cex.axis = 1.2, at = c(1.5, 4.5, 7.5, 10.5))
#text(x=1.5, y=5, "*", pos=3, cex=2)
#text(x=4.5, y=6.5, "***", pos=3, cex=2)
#text(x=7.5, y=8, "*", pos=3, cex=2)
#text(x=10.5, y=13.5, "**", pos=3, cex=2)
text(x=4.5, y=6.5, "**", pos=3, cex=2)
text(x=7.5, y=8, "**", pos=3, cex=2)

legend(0, -6.5, legend = "Female", cex = 1.3,
       fill = "red", bty = "n", xpd = TRUE)
legend(5, -6.5, legend = "Male", cex = 1.3,
       fill = "blue", bty = "n", xpd = TRUE)

