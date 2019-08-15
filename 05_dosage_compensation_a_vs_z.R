#DC analysis - calculate ratios between A and Z median values (before or after removing SBG)


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
  IN_FOLD <- "/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/12_Filter_sex_biased_genes"           
}




##### Output folders
if (FLAG_BG == 0) {
  results_dir <- "/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/05_Dosage_compensation_a_vs_z"
} else if (FLAG_BG == 1) {
  results_dir <- "/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/14_Dosage_compensation_a_vs_z_AFTER_REMOVING_SEX_BIASED_GENES"           
}




##### Results files
if (FLAG_BG == 0) {
  table_1_out = "median_table_FPKM.txt"
  table_2_out = "ratios_of_median_FPKM.txt"
  table_3_out = "mwu_test_p-values.txt"
  table_4_out = "mwu_test_p-values_acrossStages.txt"
} else if (FLAG_BG == 1) {
  table_1_out = "nonbiased_genes_median_table_FPKM.txt"
  table_2_out = "nonbiased_genes_ratios_of_median_FPKM.txt"
  table_3_out = "nonbiased_genes_mwu_test_p-values.txt"
  table_4_out = "nonbiased_genes_mwu_test_p-values_acrossStages.txt"
  table_5_out = "nonbiased_genes_median_table_FPKM-gt-0-in-both-sexes.txt"
}





########### Load tables of gene expression as mean FPKM-values for each group ####

if (FLAG_BG == 0) {
  
  ### analysis of datasets produced before removing SBG
  setwd(IN_FOLD)
  
  ### with FPKM>0 for both sexes
  instar_V <- read.delim("instar_V-assigned_A_or_Z-filtered.txt", header = TRUE)
  pupa <- read.delim("pupa-assigned_A_or_Z-filtered.txt", header = TRUE)
  adult <- read.delim("adult-assigned_A_or_Z-filtered.txt", header = TRUE)
  
  ### genes with zero counts removed individually for each sex
  instar_V_f <- read.delim("instar_V-assigned_A_or_Z_female-filtered.txt", header = TRUE)
  instar_V_m <- read.delim("instar_V-assigned_A_or_Z_male-filtered.txt", header = TRUE)
  
  pupa_f <- read.delim("pupa-assigned_A_or_Z_female-filtered.txt", header = TRUE)
  pupa_m <- read.delim("pupa-assigned_A_or_Z_male-filtered.txt", header = TRUE)
  
  adult_f <- read.delim("adult-assigned_A_or_Z_female-filtered.txt", header = TRUE)
  adult_m <- read.delim("adult-assigned_A_or_Z_male-filtered.txt", header = TRUE)

  
} else if (FLAG_BG == 1) {
  
  ### analysis of datasets produced after removing SBG
  setwd(IN_FOLD)
  

  # with FPKM>0 for both sexes
  instar_V <- read.delim("nonbiased_genes-instar_V-assigned_A_or_Z-filtered.txt", header = TRUE)
  pupa <- read.delim("nonbiased_genes-pupa-assigned_A_or_Z-filtered.txt", header = TRUE)
  adult <- read.delim("nonbiased_genes-adult-assigned_A_or_Z-filtered.txt", header = TRUE)
  
  ### genes with zero counts removed individually for each sex
  instar_V_f <- read.delim("nonbiased_genes-instar_V-assigned_A_or_Z_female-filtered.txt", header = TRUE)
  instar_V_m <- read.delim("nonbiased_genes-instar_V-assigned_A_or_Z_male-filtered.txt", header = TRUE)
  
  pupa_f <- read.delim("nonbiased_genes-pupa-assigned_A_or_Z_female-filtered.txt", header = TRUE)
  pupa_m <- read.delim("nonbiased_genes-pupa-assigned_A_or_Z_male-filtered.txt", header = TRUE)
  
  adult_f <- read.delim("nonbiased_genes-adult-assigned_A_or_Z_female-filtered.txt", header = TRUE)
  adult_m <- read.delim("nonbiased_genes-adult-assigned_A_or_Z_male-filtered.txt", header = TRUE)
  
}









######################################################################## Analysis

setwd(results_dir)

instar_V_A <- rbind(instar_V[instar_V$chromosome == "A", ])
instar_V_Z <- rbind(instar_V[instar_V$chromosome == "Z", ])
pupa_A <- rbind(pupa[pupa$chromosome == "A", ])
pupa_Z <- rbind(pupa[pupa$chromosome == "Z", ])
adult_A <- rbind(adult[adult$chromosome == "A", ])
adult_Z <- rbind(adult[adult$chromosome == "Z", ])

instar_V_A_f <- rbind(instar_V_f[instar_V_f$chromosome == "A", ])
instar_V_A_m <- rbind(instar_V_m[instar_V_m$chromosome == "A", ])
instar_V_Z_f <- rbind(instar_V_f[instar_V_f$chromosome == "Z", ])
instar_V_Z_m <- rbind(instar_V_m[instar_V_m$chromosome == "Z", ])
pupa_A_f <- rbind(pupa_f[pupa_f$chromosome == "A", ])
pupa_A_m <- rbind(pupa_m[pupa_m$chromosome == "A", ])
pupa_Z_f <- rbind(pupa_f[pupa_f$chromosome == "Z", ])
pupa_Z_m <- rbind(pupa_m[pupa_m$chromosome == "Z", ])
adult_A_f <- rbind(adult_f[adult_f$chromosome == "A", ])
adult_A_m <- rbind(adult_m[adult_m$chromosome == "A", ])
adult_Z_f <- rbind(adult_f[adult_f$chromosome == "Z", ])
adult_Z_m <- rbind(adult_m[adult_m$chromosome == "Z", ])





####################### Median and mean values

# Median FPKM (if keeping only genes with FPKM>0 in both sexes) ####
instar_V_male_a_median <- median(instar_V_A$FPKM_instar_V_male)
instar_V_male_z_median <- median(instar_V_Z$FPKM_instar_V_male)
instar_V_female_a_median <- median(instar_V_A$FPKM_instar_V_female)
instar_V_female_z_median <- median(instar_V_Z$FPKM_instar_V_female)

pupa_male_a_median <- median(pupa_A$FPKM_pupa_male)
pupa_male_z_median <- median(pupa_Z$FPKM_pupa_male)
pupa_female_a_median <- median(pupa_A$FPKM_pupa_female)
pupa_female_z_median <- median(pupa_Z$FPKM_pupa_female)

adult_male_a_median <- median(adult_A$FPKM_adult_male)
adult_male_z_median <- median(adult_Z$FPKM_adult_male)
adult_female_a_median <- median(adult_A$FPKM_adult_female)
adult_female_z_median <- median(adult_Z$FPKM_adult_female)

# Median FPKM (if genes with zero counts removed individually for each sex) ####
instar_V_male_a_median_ind <- median(instar_V_A_m$FPKM_instar_V_male)
instar_V_male_z_median_ind <- median(instar_V_Z_m$FPKM_instar_V_male)
instar_V_female_a_median_ind <- median(instar_V_A_f$FPKM_instar_V_female)
instar_V_female_z_median_ind <- median(instar_V_Z_f$FPKM_instar_V_female)

pupa_male_a_median_ind <- median(pupa_A_m$FPKM_pupa_male)
pupa_male_z_median_ind <- median(pupa_Z_m$FPKM_pupa_male)
pupa_female_a_median_ind <- median(pupa_A_f$FPKM_pupa_female)
pupa_female_z_median_ind <- median(pupa_Z_f$FPKM_pupa_female)

adult_male_a_median_ind <- median(adult_A_m$FPKM_adult_male)
adult_male_z_median_ind <- median(adult_Z_m$FPKM_adult_male)
adult_female_a_median_ind <- median(adult_A_f$FPKM_adult_female)
adult_female_z_median_ind <- median(adult_Z_f$FPKM_adult_female)

# Mean of FPKM (if keeping only genes with FPKM>0 in both sexes) ####
if (FLAG_BG == 1) {
   instar_V_male_a_mean_gt0BOTH <- mean(instar_V_A$FPKM_instar_V_male)
   instar_V_male_z_mean_gt0BOTH <- mean(instar_V_Z$FPKM_instar_V_male)
   instar_V_female_a_mean_gt0BOTH <- mean(instar_V_A$FPKM_instar_V_female)
   instar_V_female_z_mean_gt0BOTH <- mean(instar_V_Z$FPKM_instar_V_female)

   pupa_male_a_mean_gt0BOTH <- mean(pupa_A$FPKM_pupa_male)
   pupa_male_z_mean_gt0BOTH <- mean(pupa_Z$FPKM_pupa_male)
   pupa_female_a_mean_gt0BOTH <- mean(pupa_A$FPKM_pupa_female)
   pupa_female_z_mean_gt0BOTH <- mean(pupa_Z$FPKM_pupa_female)

   adult_male_a_mean_gt0BOTH <- mean(adult_A$FPKM_adult_male)
   adult_male_z_mean_gt0BOTH <- mean(adult_Z$FPKM_adult_male)
   adult_female_a_mean_gt0BOTH <- mean(adult_A$FPKM_adult_female)
   adult_female_z_mean_gt0BOTH <- mean(adult_Z$FPKM_adult_female)
}


# Mean FPKM (if genes with zero counts removed individually for each sex) ####
instar_V_male_a_mean <- mean(instar_V_A_m$FPKM_instar_V_male)
instar_V_male_z_mean <- mean(instar_V_Z_m$FPKM_instar_V_male)
instar_V_female_a_mean <- mean(instar_V_A_f$FPKM_instar_V_female)
instar_V_female_z_mean <- mean(instar_V_Z_f$FPKM_instar_V_female)

pupa_male_a_mean <- mean(pupa_A_m$FPKM_pupa_male)
pupa_male_z_mean <- mean(pupa_Z_m$FPKM_pupa_male)
pupa_female_a_mean <- mean(pupa_A_f$FPKM_pupa_female)
pupa_female_z_mean <- mean(pupa_Z_f$FPKM_pupa_female)

adult_male_a_mean <- mean(adult_A_m$FPKM_adult_male)
adult_male_z_mean <- mean(adult_Z_m$FPKM_adult_male)
adult_female_a_mean <- mean(adult_A_f$FPKM_adult_female)
adult_female_z_mean <- mean(adult_Z_f$FPKM_adult_female)

#Table of median and mean FPKM (only for analysis with genes with zero counts removed individually for each sex)

median_and_mean <- matrix(c(instar_V_female_a_median_ind, pupa_female_a_median_ind, adult_female_a_median_ind,
                            instar_V_female_z_median_ind, pupa_female_z_median_ind, adult_female_z_median_ind,
                            instar_V_male_a_median_ind, pupa_male_a_median_ind, adult_male_a_median_ind,
                            instar_V_male_z_median_ind, pupa_male_z_median_ind, adult_male_z_median_ind,
                            instar_V_female_a_mean, pupa_female_a_mean, adult_female_a_mean,
                            instar_V_female_z_mean, pupa_female_z_mean, adult_female_z_mean,
                            instar_V_male_a_mean, pupa_male_a_mean, adult_male_a_mean,
                            instar_V_male_z_mean, pupa_male_z_mean, adult_male_z_mean),
                          nrow = 3, ncol = 8)

rownames(median_and_mean) <- c("Larva","Pupa","Adult")
colnames(median_and_mean) <- c("Female median A", "Female median Z", 
                               "Male median A", "Male median Z",
                               "Female mean A", "Female mean Z",
                               "Male mean A", "Male mean Z")

write.table(median_and_mean, file = table_1_out, sep = "\t", col.names = NA, row.names = TRUE)



#Table of median and mean FPKM (if keeping only genes with FPKM>0 for both sexes)

if (FLAG_BG == 1) {
  median_and_mean_gt0BOTH <- matrix(c(instar_V_female_a_median, pupa_female_a_median, adult_female_a_median,
                                      instar_V_female_z_median, pupa_female_z_median, adult_female_z_median,
                                      instar_V_male_a_median, pupa_male_a_median, adult_male_a_median,
                                      instar_V_male_z_median, pupa_male_z_median, adult_male_z_median,
                                      instar_V_female_a_mean_gt0BOTH, pupa_female_a_mean_gt0BOTH, adult_female_a_mean_gt0BOTH,
                                      instar_V_female_z_mean_gt0BOTH, pupa_female_z_mean_gt0BOTH, adult_female_z_mean_gt0BOTH,
                                      instar_V_male_a_mean_gt0BOTH, pupa_male_a_mean_gt0BOTH, adult_male_a_mean_gt0BOTH,
                                      instar_V_male_z_mean_gt0BOTH, pupa_male_z_mean_gt0BOTH, adult_male_z_mean_gt0BOTH),
                                    nrow = 3, ncol = 8)

   rownames(median_and_mean_gt0BOTH) <- c("Larva","Pupa","Adult")
   colnames(median_and_mean_gt0BOTH) <- c("Female median A", "Female median Z", 
                               "Male median A", "Male median Z",
                               "Female mean A", "Female mean Z",
                               "Male mean A", "Male mean Z")

   write.table(median_and_mean_gt0BOTH, file = table_5_out, sep = "\t", col.names = NA, row.names = TRUE)
   
}





#### Ratios of median expression - male Z:A, female Z:A, A m:f, Z m:f ####
# (only for analysis with genes with zero counts removed individually for each sex)

# Z:A
instar_V_male_ratio <- instar_V_male_z_median_ind/instar_V_male_a_median_ind
instar_V_female_ratio <- instar_V_female_z_median_ind/instar_V_female_a_median_ind
pupa_male_ratio <- pupa_male_z_median_ind/pupa_male_a_median_ind
pupa_female_ratio <- pupa_female_z_median_ind/pupa_female_a_median_ind
adult_male_ratio <- adult_male_z_median_ind/adult_male_a_median_ind
adult_female_ratio <- adult_female_z_median_ind/adult_female_a_median_ind

# m:f 
instar_V_a_ratio <- instar_V_male_a_median/instar_V_female_a_median
instar_V_z_ratio <- instar_V_male_z_median/instar_V_female_z_median
pupa_a_ratio <- pupa_male_a_median/pupa_female_a_median
pupa_z_ratio <- pupa_male_z_median/pupa_female_z_median
adult_a_ratio <- adult_male_a_median/adult_female_a_median
adult_z_ratio <- adult_male_z_median/adult_female_z_median

ratios_of_median <- matrix(c(instar_V_female_ratio, pupa_female_ratio, adult_female_ratio,
                             instar_V_male_ratio, pupa_male_ratio, adult_male_ratio,
                             instar_V_a_ratio, pupa_a_ratio, adult_a_ratio,
                             instar_V_z_ratio, pupa_z_ratio, adult_z_ratio),
                           nrow = 3, ncol = 4)
rownames(ratios_of_median) <- c("Larva", "Pupa", "Adult")
colnames(ratios_of_median) <- c("Female Z:A ratio", "Male Z:A ratio", "A M:F ratio", "Z M:F ratio")


write.table(ratios_of_median, file = table_2_out, sep = "\t", col.names = NA, row.names = TRUE)









##### MWU-test of significant difference between Z and A for each sex, and Z M/F, A M/F ####
instar_V_male_mwu <- wilcox.test(instar_V_A_m$FPKM_instar_V_male, instar_V_Z_m$FPKM_instar_V_male,
                                 paired = FALSE)
instar_V_female_mwu <- wilcox.test(instar_V_A_f$FPKM_instar_V_female, instar_V_Z_f$FPKM_instar_V_female,
                                 paired = FALSE)
pupa_male_mwu <- wilcox.test(pupa_A_m$FPKM_pupa_male, pupa_Z_m$FPKM_pupa_male,
                                 paired = FALSE)
pupa_female_mwu <- wilcox.test(pupa_A_f$FPKM_pupa_female, pupa_Z_f$FPKM_pupa_female,
                             paired = FALSE)
adult_male_mwu <- wilcox.test(adult_A_m$FPKM_adult_male, adult_Z_m$FPKM_adult_male,
                             paired = FALSE)
adult_female_mwu <- wilcox.test(adult_A_f$FPKM_adult_female, adult_Z_f$FPKM_adult_female,
                              paired = FALSE)

instar_V_a_mwu <-wilcox.test(instar_V_A$FPKM_instar_V_male, instar_V_A$FPKM_instar_V_female,
                             paired = FALSE)
instar_V_z_mwu <-wilcox.test(instar_V_Z$FPKM_instar_V_male, instar_V_Z$FPKM_instar_V_female,
                             paired = FALSE)
pupa_a_mwu <-wilcox.test(pupa_A$FPKM_pupa_male, pupa_A$FPKM_pupa_female,
                             paired = FALSE)
pupa_z_mwu <-wilcox.test(pupa_Z$FPKM_pupa_male, pupa_Z$FPKM_pupa_female,
                             paired = FALSE)
adult_a_mwu <-wilcox.test(adult_A$FPKM_adult_male, adult_A$FPKM_adult_female,
                         paired = FALSE)
adult_z_mwu <-wilcox.test(adult_Z$FPKM_adult_male, adult_Z$FPKM_adult_female,
                         paired = FALSE)

mwu <- matrix(c(instar_V_female_mwu$p.value, pupa_female_mwu$p.value, adult_female_mwu$p.value,
                instar_V_male_mwu$p.value, pupa_male_mwu$p.value, adult_male_mwu$p.value,
                instar_V_a_mwu$p.value, pupa_a_mwu$p.value, adult_a_mwu$p.value,
                instar_V_z_mwu$p.value, pupa_z_mwu$p.value, adult_z_mwu$p.value),
              nrow = 3, ncol = 4)
rownames(mwu) <- c("Larva", "Pupa", "Adult")
colnames(mwu) <- c("Female - A vs Z", "Male - A vs Z", "A - F vs M", "Z - F vs M")


write.table(mwu, file = table_3_out, sep = "\t", col.names = NA, row.names = TRUE)





##### MWU-test of significant difference between developmental stages (performed separately for Z and A and for each sex) ####
female_instar_V_vs_pupa_A_mwu <- wilcox.test(instar_V_A_f$FPKM_instar_V_female, pupa_A_f$FPKM_pupa_female,
                                           paired = FALSE)

female_pupa_vs_adult_A_mwu <- wilcox.test(pupa_A_f$FPKM_pupa_female, adult_A_f$FPKM_adult_female,
                                         paired = FALSE)

female_instar_V_vs_adult_A_mwu <- wilcox.test(instar_V_A_f$FPKM_instar_V_female, adult_A_f$FPKM_adult_female,
                                             paired = FALSE)

female_instar_V_vs_pupa_Z_mwu <- wilcox.test(instar_V_Z_f$FPKM_instar_V_female, pupa_Z_f$FPKM_pupa_female,
                                             paired = FALSE)

female_pupa_vs_adult_Z_mwu <- wilcox.test(pupa_Z_f$FPKM_pupa_female, adult_Z_f$FPKM_adult_female,
                                          paired = FALSE)

female_instar_V_vs_adult_Z_mwu <- wilcox.test(instar_V_Z_f$FPKM_instar_V_female, adult_Z_f$FPKM_adult_female,
                                              paired = FALSE)


male_instar_V_vs_pupa_A_mwu <- wilcox.test(instar_V_A_m$FPKM_instar_V_male, pupa_A_m$FPKM_pupa_male,
                                           paired = FALSE)

male_pupa_vs_adult_A_mwu <- wilcox.test(pupa_A_m$FPKM_pupa_male, adult_A_m$FPKM_adult_male,
                                        paired = FALSE)

male_instar_V_vs_adult_A_mwu <- wilcox.test(instar_V_A_m$FPKM_instar_V_male, adult_A_m$FPKM_adult_male,
                                            paired = FALSE)

male_instar_V_vs_pupa_Z_mwu <- wilcox.test(instar_V_Z_m$FPKM_instar_V_male, pupa_Z_m$FPKM_pupa_male,
                                           paired = FALSE)

male_pupa_vs_adult_Z_mwu <- wilcox.test(pupa_Z_m$FPKM_pupa_male, adult_Z_m$FPKM_adult_male,
                                        paired = FALSE)

male_instar_V_vs_adult_Z_mwu <- wilcox.test(instar_V_Z_m$FPKM_instar_V_male, adult_Z_m$FPKM_adult_male,
                                            paired = FALSE)


mwu_stages <- matrix(c(female_instar_V_vs_pupa_A_mwu$p.value, female_pupa_vs_adult_A_mwu$p.value, female_instar_V_vs_adult_A_mwu$p.value,
                       female_instar_V_vs_pupa_Z_mwu$p.value, female_pupa_vs_adult_Z_mwu$p.value, female_instar_V_vs_adult_Z_mwu$p.value,
                       male_instar_V_vs_pupa_A_mwu$p.value, male_pupa_vs_adult_A_mwu$p.value, male_instar_V_vs_adult_A_mwu$p.value,
                       male_instar_V_vs_pupa_Z_mwu$p.value, male_pupa_vs_adult_Z_mwu$p.value, male_instar_V_vs_adult_Z_mwu$p.value),
                     nrow = 3, ncol = 4)
rownames(mwu_stages) <- c("Larva vs Pupa", "Pupa vs Adult", "Larva vs Adult")
colnames(mwu_stages) <- c("Female - A", "Female - Z", "Male - A", "Male - Z")


write.table(mwu_stages, file = table_4_out, sep = "\t", col.names = NA, row.names = TRUE)








# Bootstraping for median ratio confidence intervals ####

library(boot)

# Adult Autosome M:F ratio ####
adult_a_median_ratio_function <- function(data, indices) {
  adult_A <- data[indices,]
  (median(adult_A$FPKM_adult_male))/(median(adult_A$FPKM_adult_female))
}

adult_a_median_boot.res <- boot(data = adult_A, statistic = adult_a_median_ratio_function, R = 10000)
print(adult_a_median_boot.res)
plot(adult_a_median_boot.res)

adult_a_median_CI <- boot.ci(boot.out = adult_a_median_boot.res, type = "basic")
print(adult_a_median_CI)

# Adult Z-chromosome M:F ratio ####

adult_z_median_ratio_function <- function(data, indices) {
  adult_Z <- data[indices,]
  (median(adult_Z$FPKM_adult_male))/(median(adult_Z$FPKM_adult_female))
}

adult_z_median_boot.res <- boot(data = adult_Z, statistic = adult_z_median_ratio_function, R = 10000)
print(adult_z_median_boot.res)
plot(adult_z_median_boot.res)

adult_z_median_CI <- boot.ci(boot.out = adult_z_median_boot.res, type = "basic")
print(adult_z_median_CI)

# Pupa Autosome M:F ratio ####
pupa_a_median_ratio_function <- function(data, indices) {
  pupa_A <- data[indices,]
  (median(pupa_A$FPKM_pupa_male))/(median(pupa_A$FPKM_pupa_female))
}

pupa_a_median_boot.res <- boot(data = pupa_A, statistic = pupa_a_median_ratio_function, R = 10000)
print(pupa_a_median_boot.res)
plot(pupa_a_median_boot.res)

pupa_a_median_CI <- boot.ci(boot.out = pupa_a_median_boot.res, type = "basic")
print(pupa_a_median_CI)

# Pupa Z-chromosome M:F ratio ####
pupa_z_median_ratio_function <- function(data, indices) {
  pupa_Z <- data[indices,]
  (median(pupa_Z$FPKM_pupa_male))/(median(pupa_Z$FPKM_pupa_female))
}

pupa_z_median_boot.res <- boot(data = pupa_Z, statistic = pupa_z_median_ratio_function, R = 10000)
print(pupa_z_median_boot.res)
plot(pupa_z_median_boot.res)

pupa_z_median_CI <- boot.ci(boot.out = pupa_z_median_boot.res, type = "basic")
print(pupa_z_median_CI)


# Instar V Autosome M:F ratio ####
instar_V_a_median_ratio <- function(data, indices) {
  instar_V_A <- data[indices,]
  (median(instar_V_A$FPKM_instar_V_male))/(median(instar_V_A$FPKM_instar_V_female))
}

instar_V_a_median_boot.res <- boot(data = instar_V_A, statistic = instar_V_a_median_ratio, R = 10000)
print(instar_V_a_median_boot.res)
plot(instar_V_a_median_boot.res)

instar_V_a_median_CI <- boot.ci(boot.out = instar_V_a_median_boot.res, type = "basic")
print(instar_V_a_median_CI)

# Instar V Z-chromosome M:F ratio ####
instar_V_z_median_ratio <- function(data, indices) {
  instar_V_Z <- data[indices,]
  (median(instar_V_Z$FPKM_instar_V_male))/(median(instar_V_Z$FPKM_instar_V_female))
}

instar_V_z_median_boot.res <- boot(data = instar_V_Z, statistic = instar_V_z_median_ratio, R = 10000)
print(instar_V_z_median_boot.res)
plot(instar_V_z_median_boot.res)

instar_V_z_median_CI <- boot.ci(boot.out = instar_V_z_median_boot.res, type = "basic")
print(instar_V_z_median_CI)




########################################

# Adult Male Z:A ratio ####
adult_male_median_function <- function(data, indices) {
  adult_m <- data[indices,]
  (median(adult_m$FPKM_adult_male[adult_m$chromosome=="Z"]))/(median(adult_m$FPKM_adult_male[adult_m$chromosome=="A"]))
}

adult_male_median_boot <- boot(data = adult_m, statistic = adult_male_median_function, R = 10000)
print(adult_male_median_boot)
plot(adult_male_median_boot)

adult_male_median_CI <- boot.ci(boot.out = adult_male_median_boot, type = "basic")
print(adult_male_median_CI)

# Adult Female Z:A ratio ####
adult_female_median_function <- function(data, indices) {
  adult_f <- data[indices,]
  (median(adult_f$FPKM_adult_female[adult_f$chromosome=="Z"]))/(median(adult_f$FPKM_adult_female[adult_f$chromosome=="A"]))
}

adult_female_median_boot <- boot(data = adult_f, statistic = adult_female_median_function, R = 10000)
print(adult_female_median_boot)
plot(adult_female_median_boot)

adult_female_median_CI <- boot.ci(boot.out = adult_female_median_boot, type = "basic")
print(adult_female_median_CI)

# Pupa Male Z:A ratio ####
pupa_male_median_function <- function(data, indices) {
  pupa_m <- data[indices,]
  (median(pupa_m$FPKM_pupa_male[pupa_m$chromosome=="Z"]))/(median(pupa_m$FPKM_pupa_male[pupa_m$chromosome=="A"]))
}

pupa_male_median_boot <- boot(data = pupa_m, statistic = pupa_male_median_function, R = 10000)
print(pupa_male_median_boot)
plot(pupa_male_median_boot)

pupa_male_median_CI <- boot.ci(boot.out = pupa_male_median_boot, type = "basic")
print(pupa_male_median_CI)

# Pupa Female Z:A ratio ####
pupa_female_median_function <- function(data, indices) {
  pupa_f <- data[indices,]
  (median(pupa_f$FPKM_pupa_female[pupa_f$chromosome=="Z"]))/(median(pupa_f$FPKM_pupa_female[pupa_f$chromosome=="A"]))
}

pupa_female_median_boot <- boot(data = pupa_f, statistic = pupa_female_median_function, R = 10000)
print(pupa_female_median_boot)
plot(pupa_female_median_boot)

pupa_female_median_CI <- boot.ci(boot.out = pupa_female_median_boot, type = "basic")
print(pupa_female_median_CI)


# Instar V Male Z:A ratio ####
instar_V_male_median_function <- function(data, indices) {
  instar_V_m <- data[indices,]
  (median(instar_V_m$FPKM_instar_V_male[instar_V_m$chromosome=="Z"]))/(median(instar_V_m$FPKM_instar_V_male[instar_V_m$chromosome=="A"]))
}

instar_V_male_median_boot <- boot(data = instar_V_m, statistic = instar_V_male_median_function, R = 10000)
print(instar_V_male_median_boot)
plot(instar_V_male_median_boot)

instar_V_male_median_CI <- boot.ci(boot.out = instar_V_male_median_boot, type = "basic")
print(instar_V_male_median_CI)

# Instar V Female Z:A ratio ####
instar_V_female_median_function <- function(data, indices) {
  instar_V_f <- data[indices,]
  (median(instar_V_f$FPKM_instar_V_female[instar_V_f$chromosome=="Z"]))/(median(instar_V_f$FPKM_instar_V_female[instar_V_f$chromosome=="A"]))
}

instar_V_female_median_boot <- boot(data = instar_V_f, statistic = instar_V_female_median_function, R = 10000)
print(instar_V_female_median_boot)
plot(instar_V_female_median_boot)

instar_V_female_median_CI <- boot.ci(boot.out = instar_V_female_median_boot, type = "basic")
print(instar_V_female_median_CI)

##############################################
