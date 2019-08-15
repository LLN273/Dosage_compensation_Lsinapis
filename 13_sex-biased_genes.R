# Fisher's exact test

# Distribution of SBG - Z compared to A


# Clear all states
rm(list=ls(all=TRUE))
dev.off()




# Input folder
setwd("/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/12_Filter_sex_biased_genes")



# Results folder
RES_OUT <- "/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/13_Sex-biased_genes"




# Analysis

FBG_instar_V <- read.delim("FisherTest_FBG-instar_V-assigned_A_or_Z.txt", header = TRUE)
FBG_pupa <- read.delim("FisherTest_FBG-pupa-assigned_A_or_Z.txt", header = TRUE)
FBG_adult <- read.delim("FisherTest_FBG-adult-assigned_A_or_Z.txt", header = TRUE)


FBG_instar_V_A <- sum(FBG_instar_V$chromosome=="A")
FBG_instar_V_Z <- sum(FBG_instar_V$chromosome=="Z")
FBG_pupa_A <- sum(FBG_pupa$chromosome=="A")
FBG_pupa_Z <- sum(FBG_pupa$chromosome=="Z")
FBG_adult_A <- sum(FBG_adult$chromosome=="A")
FBG_adult_Z <- sum(FBG_adult$chromosome=="Z")

MBG_instar_V <- read.delim("FisherTest_MBG-instar_V-assigned_A_or_Z.txt", header = TRUE)
MBG_pupa <- read.delim("FisherTest_MBG-pupa-assigned_A_or_Z.txt", header = TRUE)
MBG_adult <- read.delim("FisherTest_MBG-adult-assigned_A_or_Z.txt", header = TRUE)

MBG_instar_V_A <- sum(MBG_instar_V$chromosome=="A")
MBG_instar_V_Z <- sum(MBG_instar_V$chromosome=="Z")
MBG_pupa_A <- sum(MBG_pupa$chromosome=="A")
MBG_pupa_Z <- sum(MBG_pupa$chromosome=="Z")
MBG_adult_A <- sum(MBG_adult$chromosome=="A")
MBG_adult_Z <- sum(MBG_adult$chromosome=="Z")

instar_V <- read.delim("FisherTest-instar_V-assigned_A_or_Z.txt", header = TRUE)
pupa <- read.delim("FisherTest-pupa-assigned_A_or_Z.txt", header = TRUE)
adult <- read.delim("FisherTest-adult-assigned_A_or_Z.txt", header = TRUE)

instar_V_A <- sum(instar_V$chromosome=="A")
instar_V_Z <- sum(instar_V$chromosome=="Z")
pupa_A <- sum(pupa$chromosome=="A")
pupa_Z <- sum(pupa$chromosome=="Z")
adult_A <- sum(adult$chromosome=="A")
adult_Z <- sum(adult$chromosome=="Z")

instar_V_A_non_biased <- instar_V_A-(MBG_instar_V_A+FBG_instar_V_A)
instar_V_Z_non_biased <- instar_V_Z-(MBG_instar_V_Z+FBG_instar_V_Z)
pupa_A_non_biased <- pupa_A-(MBG_pupa_A+FBG_pupa_A)
pupa_Z_non_biased <- pupa_Z-(MBG_pupa_Z+FBG_pupa_Z)
adult_A_non_biased <- adult_A-(MBG_adult_A+FBG_adult_A)
adult_Z_non_biased <- adult_Z-(MBG_adult_Z+FBG_adult_Z)


instar_V_female_test <- c(FBG_instar_V_A, instar_V_A, FBG_instar_V_Z, instar_V_Z)
instar_V_female_test <- matrix(instar_V_female_test, 2, 2)
rownames(instar_V_female_test) <- c("obs","exp")
colnames(instar_V_female_test) <- c("A","Z")
FisherTest_instar_V_female <- fisher.test(instar_V_female_test, alternative = "two.sided")

pupa_female_test <- c(FBG_pupa_A, pupa_A, FBG_pupa_Z, pupa_Z)
pupa_female_test <- matrix(pupa_female_test, 2, 2)
rownames(pupa_female_test) <- c("obs","exp")
colnames(pupa_female_test) <- c("A","Z")
FisherTest_pupa_female <- fisher.test(pupa_female_test, alternative = "two.sided")

adult_female_test <- c(FBG_adult_A, adult_A, FBG_adult_Z, adult_Z)
adult_female_test <- matrix(adult_female_test, 2, 2)
rownames(adult_female_test) <- c("obs","exp")
colnames(adult_female_test) <- c("A","Z")
FisherTest_adult_female <- fisher.test(adult_female_test, alternative = "two.sided")


instar_V_male_test <- c(MBG_instar_V_A, instar_V_A, MBG_instar_V_Z, instar_V_Z)
instar_V_male_test <- matrix(instar_V_male_test, 2, 2)
rownames(instar_V_male_test) <- c("obs","exp")
colnames(instar_V_male_test) <- c("A","Z")
FisherTest_instar_V_male <- fisher.test(instar_V_male_test, alternative = "two.sided")

pupa_male_test <- c(MBG_pupa_A, pupa_A, MBG_pupa_Z, pupa_Z)
pupa_male_test <- matrix(pupa_male_test, 2, 2)
rownames(pupa_male_test) <- c("obs","exp")
colnames(pupa_male_test) <- c("A","Z")
FisherTest_pupa_male <- fisher.test(pupa_male_test, alternative = "two.sided")

adult_male_test <- c(MBG_adult_A, adult_A, MBG_adult_Z, adult_Z)
adult_male_test <- matrix(adult_male_test, 2, 2)
rownames(adult_male_test) <- c("obs","exp")
colnames(adult_male_test) <- c("A","Z")
FisherTest_adult_male <- fisher.test(adult_male_test, alternative = "two.sided")


FisherTest_table <- c(FBG_instar_V_A/(FBG_instar_V_A+FBG_instar_V_Z),
                      FBG_pupa_A/(FBG_pupa_A+FBG_pupa_Z),
                      FBG_adult_A/(FBG_adult_A+FBG_adult_Z),
                      MBG_instar_V_A/(MBG_instar_V_A+MBG_instar_V_Z),
                      MBG_pupa_A/(MBG_pupa_A+MBG_pupa_Z),
                      MBG_adult_A/(MBG_adult_A+MBG_adult_Z),
                      
                      instar_V_A/(instar_V_A+instar_V_Z),
                      pupa_A/(pupa_A+pupa_Z),
                      adult_A/(adult_A+adult_Z),
                      instar_V_A/(instar_V_A+instar_V_Z),
                      pupa_A/(pupa_A+pupa_Z),
                      adult_A/(adult_A+adult_Z),
                      
                      FBG_instar_V_Z/(FBG_instar_V_A+FBG_instar_V_Z),
                      FBG_pupa_Z/(FBG_pupa_A+FBG_pupa_Z),
                      FBG_adult_Z/(FBG_adult_A+FBG_adult_Z),
                      MBG_instar_V_Z/(MBG_instar_V_A+MBG_instar_V_Z),
                      MBG_pupa_Z/(MBG_pupa_A+MBG_pupa_Z),
                      MBG_adult_Z/(MBG_adult_A+MBG_adult_Z),
                      
                      instar_V_Z/(instar_V_A+instar_V_Z),
                      pupa_Z/(pupa_A+pupa_Z),
                      adult_Z/(adult_A+adult_Z),
                      instar_V_Z/(instar_V_A+instar_V_Z),
                      pupa_Z/(pupa_A+pupa_Z),
                      adult_Z/(adult_A+adult_Z),
                      
                      FisherTest_instar_V_female$p.value,
                      FisherTest_pupa_female$p.value,
                      FisherTest_adult_female$p.value,
                      FisherTest_instar_V_male$p.value,
                      FisherTest_pupa_male$p.value,
                      FisherTest_adult_male$p.value
                      )

FisherTest_table <- matrix(FisherTest_table, 6, 5)
row.names(FisherTest_table) <- c("Larva_female", "Pupa_female", "Adult_female",
                                 "larva_male", "Pupa_male", "Adult_male")
colnames(FisherTest_table) <- c("A obs", "A exp", "Z obs", "Z exp", "p-value")



# save results
setwd(RES_OUT)

write.table(FisherTest_table, file = "FisherTest_table.txt", sep = "\t", col.names = NA,
            row.names = TRUE)                               


#Plot
plot_proportion <- matrix(c(NA, NA, NA,
                        FBG_adult_Z/adult_Z, adult_Z_non_biased/adult_Z, MBG_adult_Z/adult_Z,
                        FBG_adult_A/adult_A, adult_A_non_biased/adult_A, MBG_adult_A/adult_A, #adultA
                        NA, NA, NA,
                        
                        FBG_pupa_Z/pupa_Z, pupa_Z_non_biased/pupa_Z, MBG_pupa_Z/pupa_Z, 
                        FBG_pupa_A/pupa_A, pupa_A_non_biased/pupa_A, MBG_pupa_A/pupa_A,
                        NA, NA, NA,
                        
                        FBG_instar_V_Z/instar_V_Z, instar_V_Z_non_biased/instar_V_Z, MBG_instar_V_Z/instar_V_Z,
                        FBG_instar_V_A/instar_V_A, instar_V_A_non_biased/instar_V_A, MBG_instar_V_A/instar_V_A,
                        NA, NA, NA)
                        , 3, 10)

row.names(plot_proportion) <- c("FBG", "UBG", "MBG")
colnames(plot_proportion) <- c(" ",
                           "Z", "A",
                           " ",
                           "Z","A",
                           " ",
                           "Z","A",
                           " ")

par(oma=c(4, 1, 1, 1))
par(mar=c(8, 1, 1, 1))
barplot(plot_proportion, col = c("red", "lightgrey", "blue" ),
        horiz = TRUE, cex.lab = 2, las = 1, ylim = c(1,12))

text(x=0.5, y=11.3, "Larva", cex = 1.5)
text(x=0.5, y=7.6, "Pupa", cex = 1.5)
text(x=0.5, y=4.1, "Adult", cex = 1.5)

legend(0.75, -1, legend = c("Female biased", "Unbiased", "Male biased"), cex = 1.2,
            pch = 22, xpd = TRUE, pt.bg = c("red", "lightgrey", "blue"),
            col = "black", bty = "n", y.intersp = 0.85, text.width = 0.05)
axis(1, labels = "Proportion of genes", cex.axis = 1.5, at = 0.5, line = 3, tck = "0")

