# MWU-test of expression distributions before and after removing SBG, computed separately for A and Z chromosomes


#Clear all states
rm(list=ls(all=TRUE))







############ Paths and folders


##### folder containing datasets before removing biased genes 
IN_FOLD_BEFORE <- "/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/03_Normalized_libs"

##### folder containing datasets after removing biased genes 
IN_FOLD_AFTER <- "/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/12_Filter_sex_biased_genes"           





##### Output folders
results_dir <- "/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/15_MWU-test_before-vs-after_removing-SBG"





##### Results files
table_1_out = "mwu_test_before-vs-after_removing-SBG_p-value.txt"







########### Load tables of gene expression as mean FPKM-values for each group ####


  
######## datasets produced before removing SBG
setwd(IN_FOLD_BEFORE)

### genes with zero counts removed individually for each sex
instar_V_f <- read.delim("instar_V-assigned_A_or_Z_female-filtered.txt", header = TRUE)
instar_V_m <- read.delim("instar_V-assigned_A_or_Z_male-filtered.txt", header = TRUE)
  
pupa_f <- read.delim("pupa-assigned_A_or_Z_female-filtered.txt", header = TRUE)
pupa_m <- read.delim("pupa-assigned_A_or_Z_male-filtered.txt", header = TRUE)
  
adult_f <- read.delim("adult-assigned_A_or_Z_female-filtered.txt", header = TRUE)
adult_m <- read.delim("adult-assigned_A_or_Z_male-filtered.txt", header = TRUE)

  

  
######## datasets produced after removing SBG
setwd(IN_FOLD_AFTER)

### genes with zero counts removed individually for each sex
instar_V_f_AFTER <- read.delim("nonbiased_genes-instar_V-assigned_A_or_Z_female-filtered.txt", header = TRUE)
instar_V_m_AFTER <- read.delim("nonbiased_genes-instar_V-assigned_A_or_Z_male-filtered.txt", header = TRUE)
  
pupa_f_AFTER <- read.delim("nonbiased_genes-pupa-assigned_A_or_Z_female-filtered.txt", header = TRUE)
pupa_m_AFTER <- read.delim("nonbiased_genes-pupa-assigned_A_or_Z_male-filtered.txt", header = TRUE)
  
adult_f_AFTER <- read.delim("nonbiased_genes-adult-assigned_A_or_Z_female-filtered.txt", header = TRUE)
adult_m_AFTER <- read.delim("nonbiased_genes-adult-assigned_A_or_Z_male-filtered.txt", header = TRUE)
  






######################################################################## Analysis

setwd(results_dir)

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

instar_V_A_f_AFTER <- rbind(instar_V_f_AFTER[instar_V_f_AFTER$chromosome == "A", ])
instar_V_A_m_AFTER <- rbind(instar_V_m_AFTER[instar_V_m_AFTER$chromosome == "A", ])
instar_V_Z_f_AFTER <- rbind(instar_V_f_AFTER[instar_V_f_AFTER$chromosome == "Z", ])
instar_V_Z_m_AFTER <- rbind(instar_V_m_AFTER[instar_V_m_AFTER$chromosome == "Z", ])
pupa_A_f_AFTER <- rbind(pupa_f_AFTER[pupa_f_AFTER$chromosome == "A", ])
pupa_A_m_AFTER <- rbind(pupa_m_AFTER[pupa_m_AFTER$chromosome == "A", ])
pupa_Z_f_AFTER <- rbind(pupa_f_AFTER[pupa_f_AFTER$chromosome == "Z", ])
pupa_Z_m_AFTER <- rbind(pupa_m_AFTER[pupa_m_AFTER$chromosome == "Z", ])
adult_A_f_AFTER <- rbind(adult_f_AFTER[adult_f_AFTER$chromosome == "A", ])
adult_A_m_AFTER <- rbind(adult_m_AFTER[adult_m_AFTER$chromosome == "A", ])
adult_Z_f_AFTER <- rbind(adult_f_AFTER[adult_f_AFTER$chromosome == "Z", ])
adult_Z_m_AFTER <- rbind(adult_m_AFTER[adult_m_AFTER$chromosome == "Z", ])





##### MWU-test of significant difference between before and after removing SBG, for each sex and for each chromosomal class ####
instar_V_female_A_mwu <- wilcox.test(instar_V_A_f$FPKM_instar_V_female, instar_V_A_f_AFTER$FPKM_instar_V_female,
                                 paired = FALSE)

instar_V_female_Z_mwu <- wilcox.test(instar_V_Z_f$FPKM_instar_V_female, instar_V_Z_f_AFTER$FPKM_instar_V_female,
                                     paired = FALSE)

pupa_female_A_mwu <- wilcox.test(pupa_A_f$FPKM_pupa_female, pupa_A_f_AFTER$FPKM_pupa_female,
                                     paired = FALSE)

pupa_female_Z_mwu <- wilcox.test(pupa_Z_f$FPKM_pupa_female, pupa_Z_f_AFTER$FPKM_pupa_female,
                                 paired = FALSE)

adult_female_A_mwu <- wilcox.test(adult_A_f$FPKM_adult_female, adult_A_f_AFTER$FPKM_adult_female,
                                 paired = FALSE)

adult_female_Z_mwu <- wilcox.test(adult_Z_f$FPKM_adult_female, adult_Z_f_AFTER$FPKM_adult_female,
                                  paired = FALSE)


instar_V_male_A_mwu <- wilcox.test(instar_V_A_m$FPKM_instar_V_male, instar_V_A_m_AFTER$FPKM_instar_V_male,
                                   paired = FALSE)

instar_V_male_Z_mwu <- wilcox.test(instar_V_Z_m$FPKM_instar_V_male, instar_V_Z_m_AFTER$FPKM_instar_V_male,
                                   paired = FALSE)

pupa_male_A_mwu <- wilcox.test(pupa_A_m$FPKM_pupa_male, pupa_A_m_AFTER$FPKM_pupa_male,
                               paired = FALSE)

pupa_male_Z_mwu <- wilcox.test(pupa_Z_m$FPKM_pupa_male, pupa_Z_m_AFTER$FPKM_pupa_male,
                               paired = FALSE)

adult_male_A_mwu <- wilcox.test(adult_A_m$FPKM_adult_male, adult_A_m_AFTER$FPKM_adult_male,
                                paired = FALSE)

adult_male_Z_mwu <- wilcox.test(adult_Z_m$FPKM_adult_male, adult_Z_m_AFTER$FPKM_adult_male,
                                paired = FALSE)




mwu <- matrix(c(instar_V_female_A_mwu$p.value, pupa_female_A_mwu$p.value, adult_female_A_mwu$p.value,
                instar_V_female_Z_mwu$p.value, pupa_female_Z_mwu$p.value, adult_female_Z_mwu$p.value,
                instar_V_male_A_mwu$p.value, pupa_male_A_mwu$p.value, adult_male_A_mwu$p.value,
                instar_V_male_Z_mwu$p.value, pupa_male_Z_mwu$p.value, adult_male_Z_mwu$p.value),
              nrow = 3, ncol = 4)
rownames(mwu) <- c("Larva", "Pupa", "Adult")
colnames(mwu) <- c("Female A, before vs after",  "Female Z, before vs after", "Male A, before vs after", "Male Z, before vs after")


write.table(mwu, file = table_1_out, sep = "\t", col.names = NA, row.names = TRUE)





##############################################
