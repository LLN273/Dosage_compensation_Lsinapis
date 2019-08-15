#Kruskal-Wallis test & Dunn's test

#source("https://bioconductor.org/biocLite.R")
#biocLite("dunn.test")
#library(FSA)
library(dunn.test)


#Clear all states
rm(list=ls(all=TRUE))



##### Indicate whether the analysis is to be performed before (0) or after (1) removing sex-biased genes
FLAG_BG <- 1




############ Paths and folders

##### Input folders
if (FLAG_BG == 0) {
  ##### folder containing datasets before removing biased genes 
  IN_FOLD <- "/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/03_Normalized_libs"
} else if (FLAG_BG == 1) {
  ##### folder containing datasets after removing biased genes (all biased genes, either FBG or MBG)
  IN_FOLD <- "/crex1/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/17_Assign_to_chromosomes_m_vs_f_SBG"          
}


##### Output folders
if (FLAG_BG == 0) {
  results_dir <- "/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/08_Kruskal-wallis_dunn"
} else if (FLAG_BG == 1) {
  results_dir <- "/crex1/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/19_Kruskal-wallis_dunn_AFTER_REMOVING_SEX_BIASED_GENES"     
}





########### Load tables of gene expression as mean FPKM-values for each group ####

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

setwd(results_dir)

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

###### Kruskal–Wallis

instar_V_female_kw_test <- kruskal.test(instar_V_female$FPKM_instar_V_female ~ instar_V_female$chromosome)
instar_V_male_kw_test <- kruskal.test(instar_V_male$FPKM_instar_V_male ~ instar_V_male$chromosome)
pupa_female_kw_test <- kruskal.test(pupa_female$FPKM_pupa_female ~ pupa_female$chromosome)
pupa_male_kw_test <- kruskal.test(pupa_male$FPKM_pupa_male ~ pupa_male$chromosome)
adult_female_kw_test <- kruskal.test(adult_female$FPKM_adult_female ~ adult_female$chromosome)
adult_male_kw_test <- kruskal.test(adult_male$FPKM_adult_male ~ adult_male$chromosome)

kw_table <- matrix(c(instar_V_female_kw_test$p.value, instar_V_male_kw_test$p.value,
                     pupa_female_kw_test$p.value, pupa_male_kw_test$p.value,
                     adult_female_kw_test$p.value, adult_male_kw_test$p.value), nrow = 3, ncol = 2)

rownames(kw_table) <- c("Instar_V", "Pupa", "Adult")
colnames(kw_table) <- c("Female", "Male")

write.table(kw_table, file = "kruskal-wallis_test_chromosomes.txt", sep = "\t", col.names = NA, row.names = TRUE)


###### Dunn's test


instar_V_female_dunn_test <- dunn.test(instar_V_female$FPKM_instar_V_female,
                        instar_V_female$chromosome, alpha = 0.05)
instar_V_female_dunn_test_p_values <- data.frame(instar_V_female_dunn_test$comparisons,
                        instar_V_female_dunn_test$P)
colnames(instar_V_female_dunn_test_p_values) <- c("comparison","p_value_instar_V_female")
instar_V_female_dunn_test_p_values <- instar_V_female_dunn_test_p_values[order(instar_V_female_dunn_test_p_values$comparison),]
instar_V_female_significant <- instar_V_female_dunn_test_p_values[instar_V_female_dunn_test_p_values$p_value_instar_V_female<=0.05, ]

instar_V_male_dunn_test <- dunn.test(instar_V_male$FPKM_instar_V_male,
                        instar_V_male$chromosome, alpha = 0.05)
instar_V_male_dunn_test_p_values <- data.frame(instar_V_male_dunn_test$comparisons,
                        instar_V_male_dunn_test$P)
colnames(instar_V_male_dunn_test_p_values) <- c("comparison","p_value_instar_V_male")
instar_V_male_dunn_test_p_values <- instar_V_male_dunn_test_p_values[order(instar_V_male_dunn_test_p_values$comparison),]
instar_V_male_significant <- instar_V_male_dunn_test_p_values[instar_V_male_dunn_test_p_values$p_value_instar_V_male<=0.05, ]


pupa_female_dunn_test <- dunn.test(pupa_female$FPKM_pupa_female,
                        pupa_female$chromosome, alpha = 0.05)
pupa_female_dunn_test_p_values <- data.frame(pupa_female_dunn_test$comparisons,
                        pupa_female_dunn_test$P)
colnames(pupa_female_dunn_test_p_values) <- c("comparison","p_value_pupa_female")
pupa_female_dunn_test_p_values <- pupa_female_dunn_test_p_values[order(pupa_female_dunn_test_p_values$comparison),]
pupa_female_significant <- pupa_female_dunn_test_p_values[pupa_female_dunn_test_p_values$p_value_pupa_female<=0.05, ]

pupa_male_dunn_test <- dunn.test(pupa_male$FPKM_pupa_male,
                        pupa_male$chromosome, alpha = 0.05)
pupa_male_dunn_test_p_values <- data.frame(pupa_male_dunn_test$comparisons,
                        pupa_male_dunn_test$P)
colnames(pupa_male_dunn_test_p_values) <- c("comparison","p_value_pupa_male")
pupa_male_dunn_test_p_values <- pupa_male_dunn_test_p_values[order(pupa_male_dunn_test_p_values$comparison),]
pupa_male_significant <- pupa_male_dunn_test_p_values[pupa_male_dunn_test_p_values$p_value_pupa_male<=0.05, ]


adult_female_dunn_test <- dunn.test(adult_female$FPKM_adult_female,
                        adult_female$chromosome, alpha = 0.05)
adult_female_dunn_test_p_values <- data.frame(adult_female_dunn_test$comparisons,
                        adult_female_dunn_test$P)
colnames(adult_female_dunn_test_p_values) <- c("comparison","p_value_adult_female")
adult_female_dunn_test_p_values <- adult_female_dunn_test_p_values[order(adult_female_dunn_test_p_values$comparison),]
adult_female_significant <- adult_female_dunn_test_p_values[adult_female_dunn_test_p_values$p_value_adult_female<=0.05, ]

adult_male_dunn_test <- dunn.test(adult_male$FPKM_adult_male,
                        adult_male$chromosome, alpha = 0.05)
adult_male_dunn_test_p_values <- data.frame(adult_male_dunn_test$comparisons,
                        adult_male_dunn_test$P)
colnames(adult_male_dunn_test_p_values) <- c("comparison","p_value_adult_male")
adult_male_dunn_test_p_values <- adult_male_dunn_test_p_values[order(adult_male_dunn_test_p_values$comparison),]
adult_male_significant <- adult_male_dunn_test_p_values[adult_male_dunn_test_p_values$p_value_adult_male<=0.05, ]

all_dunn_test <- cbind(instar_V_female_dunn_test_p_values, instar_V_male_dunn_test_p_values,
                        pupa_female_dunn_test_p_values, pupa_male_dunn_test_p_values,
                        adult_female_dunn_test_p_values, adult_male_dunn_test_p_values)

all_dunn_test <- all_dunn_test[c(1, 2, 4, 6, 8, 10, 12)]

write.table(all_dunn_test, file = "dunn_test.txt", sep = "\t", col.names = NA, row.names = TRUE)



### PLOT

split_aux <- data.frame(do.call('rbind', strsplit(as.character(all_dunn_test$comparison),'-',fixed=TRUE)))
all_dunn_test$Chr1 <- as.numeric(as.character(split_aux$X1))
all_dunn_test$Chr2 <- as.numeric(as.character(split_aux$X2))


plot_data <- data.frame(matrix(vector(), 7, 21,
                          dimnames=list(c(), c("Z","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20"))),
                          stringsAsFactors=F)
for (j in 2:7){
   for (i in 0:20){
  
     stage_aux <- all_dunn_test[all_dunn_test[,j] < 0.05, ]
     chr_aux <- stage_aux[stage_aux$Chr1 == i | stage_aux$Chr2 == i, ]
     plot_data[(j-1),(i+1)] <- nrow(chr_aux)
  
  }
}

par(mfrow=c(1,1))

plot_data2 <- plot_data
plot_data2[2,] <- plot_data[3,]
plot_data2[3,] <- plot_data[5,]
plot_data2[4,] <- plot_data[2,]
plot_data2[5,] <- plot_data[4,]

dunn_plot <- c(plot_data2[,1],
               plot_data2[,2],
               plot_data2[,3],
               plot_data2[,4],
               plot_data2[,5],
               plot_data2[,6],
               plot_data2[,7],
               plot_data2[,8],
               plot_data2[,9],
               plot_data2[,10],
               plot_data2[,11],
               plot_data2[,12],
               plot_data2[,13],
               plot_data2[,14],
               plot_data2[,15],
               plot_data2[,16],
               plot_data2[,17],
               plot_data2[,18],
               plot_data2[,19],
               plot_data2[,20],
               plot_data2[,21])




dunn_matrix <- matrix(dunn_plot, nrow = 7, ncol = 21)


barplot(dunn_matrix, beside=TRUE, col = c("red4", "red", "hotpink", "blue3", "dodgerblue2", "cyan2", "white"),
        names.arg = c("Z", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                 "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"))
segments(0, 0, 168, 0, lty = 1, col = "grey70", lwd = 1)
axis(1, labels = "Chromosomes", line = 2, tick = FALSE, cex.axis = 1.5, at=85)

legend(120, 20, legend = c("", "Larva", "Pupa", "Adult"),
       cex = 1.2, xpd = TRUE, pt.cex = 2, pch=c(15, 22, 22, 22), bty = "n", y.intersp = 0.4, x.intersp = 0.4,
       col = c("white", "white", "white", "white"), pt.bg=c("white", "white", "white", "white"))


legend(140, 20, legend = c("", "", "", ""),
       cex = 1.2, xpd = TRUE, pt.cex = 2, pch=c(15, 22, 22, 22), bty = "n", y.intersp = 0.4, x.intersp = 0.4,
       col = c("white", "black", "black", "black"), pt.bg=c("white", "red4", "red", "hotpink"))

legend(134, 20.5, legend = c("Female"),
       cex = 1.2, xpd = TRUE, pt.cex = 2, pch=c(15), bty = "n", y.intersp = 0.4, x.intersp = 0.4,
       col = c("white"), pt.bg=c("white"))


legend(157.5, 20, legend = c("", "", "", ""),
       cex = 1.2, xpd = TRUE, pt.cex = 2, pch=c(15, 22, 22, 22), bty = "n", y.intersp = 0.4, x.intersp = 0,
       col = c("white", "black", "black", "black"), pt.bg=c("white", "blue3", "dodgerblue2", "cyan2"))

legend(155, 20.5, legend = c("Male"),
       cex = 1.2, xpd = TRUE, pt.cex = 2, pch=c(15), bty = "n", y.intersp = 0.4, x.intersp = 0,
       col = c("white"), pt.bg=c("white"))



axis(2, labels="Number of chromosomes with different expression level", at=10, cex.axis=1.2, line = 1.5, tck=0)




####################################################################
