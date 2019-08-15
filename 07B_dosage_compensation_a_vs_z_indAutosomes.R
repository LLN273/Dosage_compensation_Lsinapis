#DC analysis - calculate ratios between A and Z median values for each individual autosome

############ Clear all states
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
     IN_FOLD <- "/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/17_Assign_to_chromosomes_m_vs_f_SBG"           
}

##### Output folders
if (FLAG_BG == 0) {
     results_dir <- "/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/05_Dosage_compensation_a_vs_z"
} else if (FLAG_BG == 1) {
     results_dir <- "/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/14_Dosage_compensation_a_vs_z_AFTER_REMOVING_SEX_BIASED_GENES"           
}

##### Results files
if (FLAG_BG == 0) {
    table_1_out = "median_table_FPKM_indAutosomes.txt"
    table_2_out = "ratios_of_median_FPKM_indAutosomes.txt"
    table_3_out = "mwu_test_p-values_indAutosomes.txt"
    table_4_out = "bootstrapping_95CI_indAutosomes.txt"
} else if (FLAG_BG == 1) {
    table_1_out = "nonbiased_genes_median_table_FPKM_indAutosomes.txt"
    table_2_out = "nonbiased_genes_ratios_of_median_FPKM_indAutosomes.txt"
    table_3_out = "nonbiased_genes_mwu_test_p-values_indAutosomes.txt"
    table_4_out = "nonbiased_genes_bootstrapping_95CI_indAutosomes.txt"
}






########### Load tables of gene expression as mean FPKM-values for each group ####

if (FLAG_BG == 0) {
    
     ### analysis of datasets produced before removing SBG
     setwd(IN_FOLD)
  
     ### genes with zero counts removed individually for each sex (Z:A ratios)
     instar_V_f <- read.delim("instar_V-assigned_to_chromosomes_female-filtered.txt", header = TRUE)
     instar_V_m <- read.delim("instar_V-assigned_to_chromosomes_male-filtered.txt", header = TRUE)
     
     pupa_f <- read.delim("pupa-assigned_to_chromosomes_female-filtered.txt", header = TRUE)
     pupa_m <- read.delim("pupa-assigned_to_chromosomes_male-filtered.txt", header = TRUE)
     
     adult_f <- read.delim("adult-assigned_to_chromosomes_female-filtered.txt", header = TRUE)
     adult_m <- read.delim("adult-assigned_to_chromosomes_male-filtered.txt", header = TRUE)

} else if (FLAG_BG == 1) {

     ### analysis of datasets produced after removing SBG
     setwd(IN_FOLD)

     ### genes with zero counts removed individually for each sex (Z:A ratios)
     instar_V_f <- read.delim("nonbiased_genes-instar_V-assigned_to_chromosomes_female-filtered.txt", header = TRUE)
     instar_V_m <- read.delim("nonbiased_genes-instar_V-assigned_to_chromosomes_male-filtered.txt", header = TRUE)
     
     pupa_f <- read.delim("nonbiased_genes-pupa-assigned_to_chromosomes_female-filtered.txt", header = TRUE)
     pupa_m <- read.delim("nonbiased_genes-pupa-assigned_to_chromosomes_male-filtered.txt", header = TRUE)
     
     adult_f <- read.delim("nonbiased_genes-adult-assigned_to_chromosomes_female-filtered.txt", header = TRUE)
     adult_m <- read.delim("nonbiased_genes-adult-assigned_to_chromosomes_male-filtered.txt", header = TRUE)

}






################################################################### Analysis

setwd(results_dir)

##### Collect expression data for each sex, stage, and chromossome
for (i in c(0:20)) {
  assign(paste("instar_V_f_", toString(i), sep=""), rbind(instar_V_f[instar_V_f$chromosome == toString(i), ]))
  assign(paste("instar_V_m_", toString(i), sep=""), rbind(instar_V_m[instar_V_m$chromosome == toString(i), ]))
  assign(paste("pupa_f_", toString(i), sep=""), rbind(pupa_f[pupa_f$chromosome == toString(i), ]))
  assign(paste("pupa_m_", toString(i), sep=""), rbind(pupa_m[pupa_m$chromosome == toString(i), ]))
  assign(paste("adult_f_", toString(i), sep=""), rbind(adult_f[adult_f$chromosome == toString(i), ]))
  assign(paste("adult_m_", toString(i), sep=""), rbind(adult_m[adult_m$chromosome == toString(i), ]))
}

    




###### Median FPKM (genes with zero counts removed individually for each sex) ####

for (i in c(0:20)) {
  assign(paste("instar_V_f_median_", toString(i), sep=""), 
        median(eval(parse(text=paste(paste("instar_V_f_",toString(i), sep = ""),"FPKM_instar_V_female", sep = "$")))))
  
  assign(paste("instar_V_m_median_", toString(i), sep=""), 
         median(eval(parse(text=paste(paste("instar_V_m_",toString(i), sep = ""),"FPKM_instar_V_male", sep = "$")))))
  
  assign(paste("pupa_f_median_", toString(i), sep=""), 
         median(eval(parse(text=paste(paste("pupa_f_",toString(i), sep = ""),"FPKM_pupa_female", sep = "$")))))
  
  assign(paste("pupa_m_median_", toString(i), sep=""), 
         median(eval(parse(text=paste(paste("pupa_m_",toString(i), sep = ""),"FPKM_pupa_male", sep = "$")))))
  
  assign(paste("adult_f_median_", toString(i), sep=""), 
         median(eval(parse(text=paste(paste("adult_f_",toString(i), sep = ""),"FPKM_adult_female", sep = "$")))))
  
  assign(paste("adult_m_median_", toString(i), sep=""), 
         median(eval(parse(text=paste(paste("adult_m_",toString(i), sep = ""),"FPKM_adult_male", sep = "$")))))
}









#Table of median FPKM

median_table <- matrix(c(instar_V_f_median_0, instar_V_f_median_1, instar_V_f_median_2, instar_V_f_median_3, instar_V_f_median_4, instar_V_f_median_5, instar_V_f_median_6, instar_V_f_median_7, instar_V_f_median_8, instar_V_f_median_9, 
                         instar_V_f_median_10, instar_V_f_median_11, instar_V_f_median_12, instar_V_f_median_13, instar_V_f_median_14, instar_V_f_median_15, instar_V_f_median_16, instar_V_f_median_17, instar_V_f_median_18, instar_V_f_median_19, instar_V_f_median_20,
                         instar_V_m_median_0, instar_V_m_median_1, instar_V_m_median_2, instar_V_m_median_3, instar_V_m_median_4, instar_V_m_median_5, instar_V_m_median_6, instar_V_m_median_7, instar_V_m_median_8, instar_V_m_median_9, 
                         instar_V_m_median_10, instar_V_m_median_11, instar_V_m_median_12, instar_V_m_median_13, instar_V_m_median_14, instar_V_m_median_15, instar_V_m_median_16, instar_V_m_median_17, instar_V_m_median_18, instar_V_m_median_19, instar_V_m_median_20,
                         pupa_f_median_0, pupa_f_median_1, pupa_f_median_2, pupa_f_median_3, pupa_f_median_4, pupa_f_median_5, pupa_f_median_6, pupa_f_median_7, pupa_f_median_8, pupa_f_median_9, 
                         pupa_f_median_10, pupa_f_median_11, pupa_f_median_12, pupa_f_median_13, pupa_f_median_14, pupa_f_median_15, pupa_f_median_16, pupa_f_median_17, pupa_f_median_18, pupa_f_median_19, pupa_f_median_20,
                         pupa_m_median_0, pupa_m_median_1, pupa_m_median_2, pupa_m_median_3, pupa_m_median_4, pupa_m_median_5, pupa_m_median_6, pupa_m_median_7, pupa_m_median_8, pupa_m_median_9, 
                         pupa_m_median_10, pupa_m_median_11, pupa_m_median_12, pupa_m_median_13, pupa_m_median_14, pupa_m_median_15, pupa_m_median_16, pupa_m_median_17, pupa_m_median_18, pupa_m_median_19, pupa_m_median_20,
                         adult_f_median_0, adult_f_median_1, adult_f_median_2, adult_f_median_3, adult_f_median_4, adult_f_median_5, adult_f_median_6, adult_f_median_7, adult_f_median_8, adult_f_median_9, 
                         adult_f_median_10, adult_f_median_11, adult_f_median_12, adult_f_median_13, adult_f_median_14, adult_f_median_15, adult_f_median_16, adult_f_median_17, adult_f_median_18, adult_f_median_19, adult_f_median_20,
                         adult_m_median_0, adult_m_median_1, adult_m_median_2, adult_m_median_3, adult_m_median_4, adult_m_median_5, adult_m_median_6, adult_m_median_7, adult_m_median_8, adult_m_median_9, 
                         adult_m_median_10, adult_m_median_11, adult_m_median_12, adult_m_median_13, adult_m_median_14, adult_m_median_15, adult_m_median_16, adult_m_median_17, adult_m_median_18, adult_m_median_19, adult_m_median_20),
                       nrow = 21, ncol = 6)

colnames(median_table) <- c("Instar V Female, median","Instar V Male, median",
                            "Pupa Female, median","Pupa Male, median",
                            "Adult Female, median","Adult Male, median")
rownames(median_table) <- c("Z", "Autosome 1", "Autosome 2", "Autosome 3", "Autosome 4", "Autosome 5", "Autosome 6", "Autosome 7", "Autosome 8", 
                            "Autosome 9", "Autosome 10", "Autosome 11", "Autosome 12", "Autosome 13", "Autosome 14", "Autosome 15", "Autosome 16", 
                            "Autosome 17", "Autosome 18", "Autosome 19", "Autosome 20")

write.table(median_table, file = table_1_out, sep = "\t", col.names = NA, row.names = TRUE)










###### Ratios of median expression - male Z:A, female Z:A ####


# Genes with zero counts removed individually for each sex)



for (i in c(1:20)) {
  
  assign(paste("instar_V_f_ratio_z_to_", toString(i), sep=""), 
         eval(parse(text="instar_V_f_median_0")) / eval(parse(text=paste("instar_V_f_median_",toString(i), sep = ""))))
  
  assign(paste("instar_V_m_ratio_z_to_", toString(i), sep=""), 
         eval(parse(text="instar_V_m_median_0")) / eval(parse(text=paste("instar_V_m_median_",toString(i), sep = ""))))
  
  assign(paste("pupa_f_ratio_z_to_", toString(i), sep=""), 
         eval(parse(text="pupa_f_median_0")) / eval(parse(text=paste("pupa_f_median_",toString(i), sep = ""))))
  
  assign(paste("pupa_m_ratio_z_to_", toString(i), sep=""), 
         eval(parse(text="pupa_m_median_0")) / eval(parse(text=paste("pupa_m_median_",toString(i), sep = ""))))
  
  assign(paste("adult_f_ratio_z_to_", toString(i), sep=""), 
         eval(parse(text="adult_f_median_0")) / eval(parse(text=paste("adult_f_median_",toString(i), sep = ""))))
  
  assign(paste("adult_m_ratio_z_to_", toString(i), sep=""), 
         eval(parse(text="adult_m_median_0")) / eval(parse(text=paste("adult_m_median_",toString(i), sep = ""))))
  
}


## table with ratio of median expressions

ratios_of_median <- matrix(c(instar_V_f_ratio_z_to_1, instar_V_f_ratio_z_to_2, instar_V_f_ratio_z_to_3, instar_V_f_ratio_z_to_4, instar_V_f_ratio_z_to_5, instar_V_f_ratio_z_to_6, instar_V_f_ratio_z_to_7, instar_V_f_ratio_z_to_8, instar_V_f_ratio_z_to_9, 
                                 instar_V_f_ratio_z_to_10, instar_V_f_ratio_z_to_11, instar_V_f_ratio_z_to_12, instar_V_f_ratio_z_to_13, instar_V_f_ratio_z_to_14, instar_V_f_ratio_z_to_15, instar_V_f_ratio_z_to_16, instar_V_f_ratio_z_to_17, instar_V_f_ratio_z_to_18, instar_V_f_ratio_z_to_19, instar_V_f_ratio_z_to_20,
                                 instar_V_m_ratio_z_to_1, instar_V_m_ratio_z_to_2, instar_V_m_ratio_z_to_3, instar_V_m_ratio_z_to_4, instar_V_m_ratio_z_to_5, instar_V_m_ratio_z_to_6, instar_V_m_ratio_z_to_7, instar_V_m_ratio_z_to_8, instar_V_m_ratio_z_to_9, 
                                 instar_V_m_ratio_z_to_10, instar_V_m_ratio_z_to_11, instar_V_m_ratio_z_to_12, instar_V_m_ratio_z_to_13, instar_V_m_ratio_z_to_14, instar_V_m_ratio_z_to_15, instar_V_m_ratio_z_to_16, instar_V_m_ratio_z_to_17, instar_V_m_ratio_z_to_18, instar_V_m_ratio_z_to_19, instar_V_m_ratio_z_to_20,
                                 pupa_f_ratio_z_to_1, pupa_f_ratio_z_to_2, pupa_f_ratio_z_to_3, pupa_f_ratio_z_to_4, pupa_f_ratio_z_to_5, pupa_f_ratio_z_to_6, pupa_f_ratio_z_to_7, pupa_f_ratio_z_to_8, pupa_f_ratio_z_to_9, 
                                 pupa_f_ratio_z_to_10, pupa_f_ratio_z_to_11, pupa_f_ratio_z_to_12, pupa_f_ratio_z_to_13, pupa_f_ratio_z_to_14, pupa_f_ratio_z_to_15, pupa_f_ratio_z_to_16, pupa_f_ratio_z_to_17, pupa_f_ratio_z_to_18, pupa_f_ratio_z_to_19, pupa_f_ratio_z_to_20,
                                 pupa_m_ratio_z_to_1, pupa_m_ratio_z_to_2, pupa_m_ratio_z_to_3, pupa_m_ratio_z_to_4, pupa_m_ratio_z_to_5, pupa_m_ratio_z_to_6, pupa_m_ratio_z_to_7, pupa_m_ratio_z_to_8, pupa_m_ratio_z_to_9, 
                                 pupa_m_ratio_z_to_10, pupa_m_ratio_z_to_11, pupa_m_ratio_z_to_12, pupa_m_ratio_z_to_13, pupa_m_ratio_z_to_14, pupa_m_ratio_z_to_15, pupa_m_ratio_z_to_16, pupa_m_ratio_z_to_17, pupa_m_ratio_z_to_18, pupa_m_ratio_z_to_19, pupa_m_ratio_z_to_20,
                                 adult_f_ratio_z_to_1, adult_f_ratio_z_to_2, adult_f_ratio_z_to_3, adult_f_ratio_z_to_4, adult_f_ratio_z_to_5, adult_f_ratio_z_to_6, adult_f_ratio_z_to_7, adult_f_ratio_z_to_8, adult_f_ratio_z_to_9, 
                                 adult_f_ratio_z_to_10, adult_f_ratio_z_to_11, adult_f_ratio_z_to_12, adult_f_ratio_z_to_13, adult_f_ratio_z_to_14, adult_f_ratio_z_to_15, adult_f_ratio_z_to_16, adult_f_ratio_z_to_17, adult_f_ratio_z_to_18, adult_f_ratio_z_to_19, adult_f_ratio_z_to_20,
                                 adult_m_ratio_z_to_1, adult_m_ratio_z_to_2, adult_m_ratio_z_to_3, adult_m_ratio_z_to_4, adult_m_ratio_z_to_5, adult_m_ratio_z_to_6, adult_m_ratio_z_to_7, adult_m_ratio_z_to_8, adult_m_ratio_z_to_9, 
                                 adult_m_ratio_z_to_10, adult_m_ratio_z_to_11, adult_m_ratio_z_to_12, adult_m_ratio_z_to_13, adult_m_ratio_z_to_14, adult_m_ratio_z_to_15, adult_m_ratio_z_to_16, adult_m_ratio_z_to_17, adult_m_ratio_z_to_18, adult_m_ratio_z_to_19, adult_m_ratio_z_to_20),
                                nrow = 20, ncol = 6)

colnames(ratios_of_median) <- c("Instar V Female, median","Instar V Male, median",
                            "Pupa Female, median","Pupa Male, median",
                            "Adult Female, median","Adult Male, median")
rownames(ratios_of_median) <- c("Z/A1", "Z/A2", "Z/A3", "Z/A4", "Z/A5", "Z/A6", "Z/A7", "Z/A8", 
                                "Z/A9", "Z/A10", "Z/A11", "Z/A12", "Z/A13", "Z/A14", "Z/A15", "Z/A16", 
                                "Z/A17", "Z/A18", "Z/A19", "Z/A20")

write.table(ratios_of_median, file = table_2_out, sep = "\t", col.names = NA, row.names = TRUE)








#################### MWU-test of significant difference between Z and A for each sex ####


for (i in c(1:20)) {
  assign(paste("instar_V_f_mwu_", toString(i), sep=""), 
         wilcox.test(eval(parse(text=paste(paste("instar_V_f_",toString(i), sep = ""),"FPKM_instar_V_female", sep = "$"))),
                     eval(parse(text=paste("instar_V_f_0","FPKM_instar_V_female", sep = "$")))))
  
  assign(paste("instar_V_m_mwu_", toString(i), sep=""), 
         wilcox.test(eval(parse(text=paste(paste("instar_V_m_",toString(i), sep = ""),"FPKM_instar_V_male", sep = "$"))),
                     eval(parse(text=paste("instar_V_m_0","FPKM_instar_V_male", sep = "$")))))
  
  assign(paste("pupa_f_mwu_", toString(i), sep=""), 
         wilcox.test(eval(parse(text=paste(paste("pupa_f_",toString(i), sep = ""),"FPKM_pupa_female", sep = "$"))),
                     eval(parse(text=paste("pupa_f_0","FPKM_pupa_female", sep = "$")))))
  
  assign(paste("pupa_m_mwu_", toString(i), sep=""), 
         wilcox.test(eval(parse(text=paste(paste("pupa_m_",toString(i), sep = ""),"FPKM_pupa_male", sep = "$"))),
                     eval(parse(text=paste("pupa_m_0","FPKM_pupa_male", sep = "$")))))
         
         
  assign(paste("adult_f_mwu_", toString(i), sep=""), 
         wilcox.test(eval(parse(text=paste(paste("adult_f_",toString(i), sep = ""),"FPKM_adult_female", sep = "$"))),
                     eval(parse(text=paste("adult_f_0","FPKM_adult_female", sep = "$")))))
  
  assign(paste("adult_m_mwu_", toString(i), sep=""), 
         wilcox.test(eval(parse(text=paste(paste("adult_m_",toString(i), sep = ""),"FPKM_adult_male", sep = "$"))),
                     eval(parse(text=paste("adult_m_0","FPKM_adult_male", sep = "$")))))
  
  
}






### Table with mwu-test p-values


mwu <- matrix(c(instar_V_f_mwu_1$p.value, instar_V_f_mwu_2$p.value, instar_V_f_mwu_3$p.value, instar_V_f_mwu_4$p.value, instar_V_f_mwu_5$p.value, instar_V_f_mwu_6$p.value, instar_V_f_mwu_7$p.value, instar_V_f_mwu_8$p.value, instar_V_f_mwu_9$p.value, 
                             instar_V_f_mwu_10$p.value, instar_V_f_mwu_11$p.value, instar_V_f_mwu_12$p.value, instar_V_f_mwu_13$p.value, instar_V_f_mwu_14$p.value, instar_V_f_mwu_15$p.value, instar_V_f_mwu_16$p.value, instar_V_f_mwu_17$p.value, instar_V_f_mwu_18$p.value, instar_V_f_mwu_19$p.value, instar_V_f_mwu_20$p.value,
                             instar_V_m_mwu_1$p.value, instar_V_m_mwu_2$p.value, instar_V_m_mwu_3$p.value, instar_V_m_mwu_4$p.value, instar_V_m_mwu_5$p.value, instar_V_m_mwu_6$p.value, instar_V_m_mwu_7$p.value, instar_V_m_mwu_8$p.value, instar_V_m_mwu_9$p.value, 
                             instar_V_m_mwu_10$p.value, instar_V_m_mwu_11$p.value, instar_V_m_mwu_12$p.value, instar_V_m_mwu_13$p.value, instar_V_m_mwu_14$p.value, instar_V_m_mwu_15$p.value, instar_V_m_mwu_16$p.value, instar_V_m_mwu_17$p.value, instar_V_m_mwu_18$p.value, instar_V_m_mwu_19$p.value, instar_V_m_mwu_20$p.value,
                             pupa_f_mwu_1$p.value, pupa_f_mwu_2$p.value, pupa_f_mwu_3$p.value, pupa_f_mwu_4$p.value, pupa_f_mwu_5$p.value, pupa_f_mwu_6$p.value, pupa_f_mwu_7$p.value, pupa_f_mwu_8$p.value, pupa_f_mwu_9$p.value, 
                             pupa_f_mwu_10$p.value, pupa_f_mwu_11$p.value, pupa_f_mwu_12$p.value, pupa_f_mwu_13$p.value, pupa_f_mwu_14$p.value, pupa_f_mwu_15$p.value, pupa_f_mwu_16$p.value, pupa_f_mwu_17$p.value, pupa_f_mwu_18$p.value, pupa_f_mwu_19$p.value, pupa_f_mwu_20$p.value,
                             pupa_m_mwu_1$p.value, pupa_m_mwu_2$p.value, pupa_m_mwu_3$p.value, pupa_m_mwu_4$p.value, pupa_m_mwu_5$p.value, pupa_m_mwu_6$p.value, pupa_m_mwu_7$p.value, pupa_m_mwu_8$p.value, pupa_m_mwu_9$p.value, 
                             pupa_m_mwu_10$p.value, pupa_m_mwu_11$p.value, pupa_m_mwu_12$p.value, pupa_m_mwu_13$p.value, pupa_m_mwu_14$p.value, pupa_m_mwu_15$p.value, pupa_m_mwu_16$p.value, pupa_m_mwu_17$p.value, pupa_m_mwu_18$p.value, pupa_m_mwu_19$p.value, pupa_m_mwu_20$p.value,
                             adult_f_mwu_1$p.value, adult_f_mwu_2$p.value, adult_f_mwu_3$p.value, adult_f_mwu_4$p.value, adult_f_mwu_5$p.value, adult_f_mwu_6$p.value, adult_f_mwu_7$p.value, adult_f_mwu_8$p.value, adult_f_mwu_9$p.value, 
                             adult_f_mwu_10$p.value, adult_f_mwu_11$p.value, adult_f_mwu_12$p.value, adult_f_mwu_13$p.value, adult_f_mwu_14$p.value, adult_f_mwu_15$p.value, adult_f_mwu_16$p.value, adult_f_mwu_17$p.value, adult_f_mwu_18$p.value, adult_f_mwu_19$p.value, adult_f_mwu_20$p.value,
                             adult_m_mwu_1$p.value, adult_m_mwu_2$p.value, adult_m_mwu_3$p.value, adult_m_mwu_4$p.value, adult_m_mwu_5$p.value, adult_m_mwu_6$p.value, adult_m_mwu_7$p.value, adult_m_mwu_8$p.value, adult_m_mwu_9$p.value, 
                             adult_m_mwu_10$p.value, adult_m_mwu_11$p.value, adult_m_mwu_12$p.value, adult_m_mwu_13$p.value, adult_m_mwu_14$p.value, adult_m_mwu_15$p.value, adult_m_mwu_16$p.value, adult_m_mwu_17$p.value, adult_m_mwu_18$p.value, adult_m_mwu_19$p.value, adult_m_mwu_20$p.value),
                           nrow = 20, ncol = 6)

colnames(mwu) <- c("Instar V Female, p.value","Instar V Male, p.value",
                                "Pupa Female, p.value","Pupa Male, p.value",
                                "Adult Female, p.value","Adult Male, p.value")

rownames(mwu) <- c("Z/A1", "Z/A2", "Z/A3", "Z/A4", "Z/A5", "Z/A6", "Z/A7", "Z/A8", 
                                "Z/A9", "Z/A10", "Z/A11", "Z/A12", "Z/A13", "Z/A14", "Z/A15", "Z/A16", 
                                "Z/A17", "Z/A18", "Z/A19", "Z/A20")

write.table(mwu, file = table_3_out, sep = "\t", col.names = NA, row.names = TRUE)











################ Bootstraping for median ratio confidence intervals ####
################ Note: runtime for the next section is about 1h

library(boot)


for (i in c(1:20)) {
  
  # instar_V_f
  assign(paste("instar_V_f_median_function_", toString(i), sep=""), 
         function(data, indices) { instar_V_f <- data[indices,]
         (median(instar_V_f$FPKM_instar_V_female[instar_V_f$chromosome=="0"]))/
         (median(instar_V_f$FPKM_instar_V_female[instar_V_f$chromosome==toString(i)]))
         })  
  
  assign(paste("instar_V_female_median_boot_", toString(i), sep=""), 
         boot(data = instar_V_f, statistic = eval(parse(text=paste("instar_V_f_median_function_", toString(i), sep=""))), R = 10000))

  assign(paste("instar_V_female_median_CI_", toString(i), sep=""), 
         boot.ci(boot.out = eval(parse(text=paste("instar_V_female_median_boot_", toString(i), sep=""))), type = "basic"))
  
  assign(paste("instar_V_female_median_CI_95_lower_", toString(i), sep=""), 
         eval(parse(text=paste(paste(paste("instar_V_female_median_CI_", toString(i), sep=""),
                               "basic", sep = "$"),"[4]", sep=""))))
  
  assign(paste("instar_V_female_median_CI_95_upper_", toString(i), sep=""), 
         eval(parse(text=paste(paste(paste("instar_V_female_median_CI_", toString(i), sep=""),
                                     "basic", sep = "$"),"[5]", sep=""))))
  
  # instar_V_m
  assign(paste("instar_V_m_median_function_", toString(i), sep=""), 
         function(data, indices) { instar_V_m <- data[indices,]
         (median(instar_V_m$FPKM_instar_V_male[instar_V_m$chromosome=="0"]))/
           (median(instar_V_m$FPKM_instar_V_male[instar_V_m$chromosome==toString(i)]))
         })  
  
  assign(paste("instar_V_male_median_boot_", toString(i), sep=""), 
         boot(data = instar_V_m, statistic = eval(parse(text=paste("instar_V_m_median_function_", toString(i), sep=""))), R = 10000))
  
  assign(paste("instar_V_male_median_CI_", toString(i), sep=""), 
         boot.ci(boot.out = eval(parse(text=paste("instar_V_male_median_boot_", toString(i), sep=""))), type = "basic"))
  
  assign(paste("instar_V_male_median_CI_95_lower_", toString(i), sep=""), 
         eval(parse(text=paste(paste(paste("instar_V_male_median_CI_", toString(i), sep=""),
                                     "basic", sep = "$"),"[4]", sep=""))))
  
  assign(paste("instar_V_male_median_CI_95_upper_", toString(i), sep=""), 
         eval(parse(text=paste(paste(paste("instar_V_male_median_CI_", toString(i), sep=""),
                                     "basic", sep = "$"),"[5]", sep=""))))
  
  # pupa_f
  assign(paste("pupa_f_median_function_", toString(i), sep=""), 
         function(data, indices) { pupa_f <- data[indices,]
         (median(pupa_f$FPKM_pupa_female[pupa_f$chromosome=="0"]))/
           (median(pupa_f$FPKM_pupa_female[pupa_f$chromosome==toString(i)]))
         })  
  
  assign(paste("pupa_female_median_boot_", toString(i), sep=""), 
         boot(data = pupa_f, statistic = eval(parse(text=paste("pupa_f_median_function_", toString(i), sep=""))), R = 10000))
  
  assign(paste("pupa_female_median_CI_", toString(i), sep=""), 
         boot.ci(boot.out = eval(parse(text=paste("pupa_female_median_boot_", toString(i), sep=""))), type = "basic"))
  
  assign(paste("pupa_female_median_CI_95_lower_", toString(i), sep=""), 
         eval(parse(text=paste(paste(paste("pupa_female_median_CI_", toString(i), sep=""),
                                     "basic", sep = "$"),"[4]", sep=""))))
  
  assign(paste("pupa_female_median_CI_95_upper_", toString(i), sep=""), 
         eval(parse(text=paste(paste(paste("pupa_female_median_CI_", toString(i), sep=""),
                                     "basic", sep = "$"),"[5]", sep=""))))
  
  
  # pupa_m
  assign(paste("pupa_m_median_function_", toString(i), sep=""), 
         function(data, indices) { pupa_m <- data[indices,]
         (median(pupa_m$FPKM_pupa_male[pupa_m$chromosome=="0"]))/
           (median(pupa_m$FPKM_pupa_male[pupa_m$chromosome==toString(i)]))
         })  
  
  assign(paste("pupa_male_median_boot_", toString(i), sep=""), 
         boot(data = pupa_m, statistic = eval(parse(text=paste("pupa_m_median_function_", toString(i), sep=""))), R = 10000))
  
  assign(paste("pupa_male_median_CI_", toString(i), sep=""), 
         boot.ci(boot.out = eval(parse(text=paste("pupa_male_median_boot_", toString(i), sep=""))), type = "basic"))
  
  assign(paste("pupa_male_median_CI_95_lower_", toString(i), sep=""), 
         eval(parse(text=paste(paste(paste("pupa_male_median_CI_", toString(i), sep=""),
                                     "basic", sep = "$"),"[4]", sep=""))))
  
  assign(paste("pupa_male_median_CI_95_upper_", toString(i), sep=""), 
         eval(parse(text=paste(paste(paste("pupa_male_median_CI_", toString(i), sep=""),
                                     "basic", sep = "$"),"[5]", sep=""))))
  
  # adult_f
  assign(paste("adult_f_median_function_", toString(i), sep=""), 
         function(data, indices) { adult_f <- data[indices,]
         (median(adult_f$FPKM_adult_female[adult_f$chromosome=="0"]))/
           (median(adult_f$FPKM_adult_female[adult_f$chromosome==toString(i)]))
         })  
  
  assign(paste("adult_female_median_boot_", toString(i), sep=""), 
         boot(data = adult_f, statistic = eval(parse(text=paste("adult_f_median_function_", toString(i), sep=""))), R = 10000))
  
  assign(paste("adult_female_median_CI_", toString(i), sep=""), 
         boot.ci(boot.out = eval(parse(text=paste("adult_female_median_boot_", toString(i), sep=""))), type = "basic"))
  
  assign(paste("adult_female_median_CI_95_lower_", toString(i), sep=""), 
         eval(parse(text=paste(paste(paste("adult_female_median_CI_", toString(i), sep=""),
                                     "basic", sep = "$"),"[4]", sep=""))))
  
  assign(paste("adult_female_median_CI_95_upper_", toString(i), sep=""), 
         eval(parse(text=paste(paste(paste("adult_female_median_CI_", toString(i), sep=""),
                                     "basic", sep = "$"),"[5]", sep=""))))
  
  
  # adult_m
  assign(paste("adult_m_median_function_", toString(i), sep=""), 
         function(data, indices) { adult_m <- data[indices,]
         (median(adult_m$FPKM_adult_male[adult_m$chromosome=="0"]))/
           (median(adult_m$FPKM_adult_male[adult_m$chromosome==toString(i)]))
         })  
  
  assign(paste("adult_male_median_boot_", toString(i), sep=""), 
         boot(data = adult_m, statistic = eval(parse(text=paste("adult_m_median_function_", toString(i), sep=""))), R = 10000))
  
  assign(paste("adult_male_median_CI_", toString(i), sep=""), 
         boot.ci(boot.out = eval(parse(text=paste("adult_male_median_boot_", toString(i), sep=""))), type = "basic"))
  
  assign(paste("adult_male_median_CI_95_lower_", toString(i), sep=""), 
         eval(parse(text=paste(paste(paste("adult_male_median_CI_", toString(i), sep=""),
                                     "basic", sep = "$"),"[4]", sep=""))))
  
  assign(paste("adult_male_median_CI_95_upper_", toString(i), sep=""), 
         eval(parse(text=paste(paste(paste("adult_male_median_CI_", toString(i), sep=""),
                                     "basic", sep = "$"),"[5]", sep=""))))
  
}




for (i in c(1:20)) {
  
  assign(paste("instar_V_female_median_CI_95_lower_", toString(i), sep=""), 
         max(0, eval(parse(text=paste("instar_V_female_median_CI_95_lower_", toString(i), sep="")))))
  
  assign(paste("instar_V_male_median_CI_95_lower_", toString(i), sep=""), 
         max(0, eval(parse(text=paste("instar_V_male_median_CI_95_lower_", toString(i), sep="")))))
  
  assign(paste("pupa_female_median_CI_95_lower_", toString(i), sep=""), 
         max(0, eval(parse(text=paste("pupa_female_median_CI_95_lower_", toString(i), sep="")))))
  
  assign(paste("pupa_male_median_CI_95_lower_", toString(i), sep=""), 
         max(0, eval(parse(text=paste("pupa_male_median_CI_95_lower_", toString(i), sep="")))))
  
  assign(paste("adult_female_median_CI_95_lower_", toString(i), sep=""), 
         max(0, eval(parse(text=paste("adult_female_median_CI_95_lower_", toString(i), sep="")))))
  
  assign(paste("adult_male_median_CI_95_lower_", toString(i), sep=""), 
         max(0, eval(parse(text=paste("adult_male_median_CI_95_lower_", toString(i), sep="")))))
  
}






### Table confidence intervals


bootstrap_CI <- matrix(c(instar_V_female_median_CI_95_lower_1, instar_V_female_median_CI_95_lower_2, instar_V_female_median_CI_95_lower_3, instar_V_female_median_CI_95_lower_4, instar_V_female_median_CI_95_lower_5, instar_V_female_median_CI_95_lower_6, instar_V_female_median_CI_95_lower_7, instar_V_female_median_CI_95_lower_8, instar_V_female_median_CI_95_lower_9, 
                         instar_V_female_median_CI_95_lower_10, instar_V_female_median_CI_95_lower_11, instar_V_female_median_CI_95_lower_12, instar_V_female_median_CI_95_lower_13, instar_V_female_median_CI_95_lower_14, instar_V_female_median_CI_95_lower_15, instar_V_female_median_CI_95_lower_16, instar_V_female_median_CI_95_lower_17, instar_V_female_median_CI_95_lower_18, instar_V_female_median_CI_95_lower_19, instar_V_female_median_CI_95_lower_20,
                         instar_V_female_median_CI_95_upper_1, instar_V_female_median_CI_95_upper_2, instar_V_female_median_CI_95_upper_3, instar_V_female_median_CI_95_upper_4, instar_V_female_median_CI_95_upper_5, instar_V_female_median_CI_95_upper_6, instar_V_female_median_CI_95_upper_7, instar_V_female_median_CI_95_upper_8, instar_V_female_median_CI_95_upper_9, 
                         instar_V_female_median_CI_95_upper_10, instar_V_female_median_CI_95_upper_11, instar_V_female_median_CI_95_upper_12, instar_V_female_median_CI_95_upper_13, instar_V_female_median_CI_95_upper_14, instar_V_female_median_CI_95_upper_15, instar_V_female_median_CI_95_upper_16, instar_V_female_median_CI_95_upper_17, instar_V_female_median_CI_95_upper_18, instar_V_female_median_CI_95_upper_19, instar_V_female_median_CI_95_upper_20,
                         instar_V_male_median_CI_95_lower_1, instar_V_male_median_CI_95_lower_2, instar_V_male_median_CI_95_lower_3, instar_V_male_median_CI_95_lower_4, instar_V_male_median_CI_95_lower_5, instar_V_male_median_CI_95_lower_6, instar_V_male_median_CI_95_lower_7, instar_V_male_median_CI_95_lower_8, instar_V_male_median_CI_95_lower_9, 
                         instar_V_male_median_CI_95_lower_10, instar_V_male_median_CI_95_lower_11, instar_V_male_median_CI_95_lower_12, instar_V_male_median_CI_95_lower_13, instar_V_male_median_CI_95_lower_14, instar_V_male_median_CI_95_lower_15, instar_V_male_median_CI_95_lower_16, instar_V_male_median_CI_95_lower_17, instar_V_male_median_CI_95_lower_18, instar_V_male_median_CI_95_lower_19, instar_V_male_median_CI_95_lower_20,
                         instar_V_male_median_CI_95_upper_1, instar_V_male_median_CI_95_upper_2, instar_V_male_median_CI_95_upper_3, instar_V_male_median_CI_95_upper_4, instar_V_male_median_CI_95_upper_5, instar_V_male_median_CI_95_upper_6, instar_V_male_median_CI_95_upper_7, instar_V_male_median_CI_95_upper_8, instar_V_male_median_CI_95_upper_9, 
                         instar_V_male_median_CI_95_upper_10, instar_V_male_median_CI_95_upper_11, instar_V_male_median_CI_95_upper_12, instar_V_male_median_CI_95_upper_13, instar_V_male_median_CI_95_upper_14, instar_V_male_median_CI_95_upper_15, instar_V_male_median_CI_95_upper_16, instar_V_male_median_CI_95_upper_17, instar_V_male_median_CI_95_upper_18, instar_V_male_median_CI_95_upper_19, instar_V_male_median_CI_95_upper_20,
                         pupa_female_median_CI_95_lower_1, pupa_female_median_CI_95_lower_2, pupa_female_median_CI_95_lower_3, pupa_female_median_CI_95_lower_4, pupa_female_median_CI_95_lower_5, pupa_female_median_CI_95_lower_6, pupa_female_median_CI_95_lower_7, pupa_female_median_CI_95_lower_8, pupa_female_median_CI_95_lower_9, 
                         pupa_female_median_CI_95_lower_10, pupa_female_median_CI_95_lower_11, pupa_female_median_CI_95_lower_12, pupa_female_median_CI_95_lower_13, pupa_female_median_CI_95_lower_14, pupa_female_median_CI_95_lower_15, pupa_female_median_CI_95_lower_16, pupa_female_median_CI_95_lower_17, pupa_female_median_CI_95_lower_18, pupa_female_median_CI_95_lower_19, pupa_female_median_CI_95_lower_20,
                         pupa_female_median_CI_95_upper_1, pupa_female_median_CI_95_upper_2, pupa_female_median_CI_95_upper_3, pupa_female_median_CI_95_upper_4, pupa_female_median_CI_95_upper_5, pupa_female_median_CI_95_upper_6, pupa_female_median_CI_95_upper_7, pupa_female_median_CI_95_upper_8, pupa_female_median_CI_95_upper_9, 
                         pupa_female_median_CI_95_upper_10, pupa_female_median_CI_95_upper_11, pupa_female_median_CI_95_upper_12, pupa_female_median_CI_95_upper_13, pupa_female_median_CI_95_upper_14, pupa_female_median_CI_95_upper_15, pupa_female_median_CI_95_upper_16, pupa_female_median_CI_95_upper_17, pupa_female_median_CI_95_upper_18, pupa_female_median_CI_95_upper_19, pupa_female_median_CI_95_upper_20,
                         pupa_male_median_CI_95_lower_1, pupa_male_median_CI_95_lower_2, pupa_male_median_CI_95_lower_3, pupa_male_median_CI_95_lower_4, pupa_male_median_CI_95_lower_5, pupa_male_median_CI_95_lower_6, pupa_male_median_CI_95_lower_7, pupa_male_median_CI_95_lower_8, pupa_male_median_CI_95_lower_9, 
                         pupa_male_median_CI_95_lower_10, pupa_male_median_CI_95_lower_11, pupa_male_median_CI_95_lower_12, pupa_male_median_CI_95_lower_13, pupa_male_median_CI_95_lower_14, pupa_male_median_CI_95_lower_15, pupa_male_median_CI_95_lower_16, pupa_male_median_CI_95_lower_17, pupa_male_median_CI_95_lower_18, pupa_male_median_CI_95_lower_19, pupa_male_median_CI_95_lower_20,
                         pupa_male_median_CI_95_upper_1, pupa_male_median_CI_95_upper_2, pupa_male_median_CI_95_upper_3, pupa_male_median_CI_95_upper_4, pupa_male_median_CI_95_upper_5, pupa_male_median_CI_95_upper_6, pupa_male_median_CI_95_upper_7, pupa_male_median_CI_95_upper_8, pupa_male_median_CI_95_upper_9, 
                         pupa_male_median_CI_95_upper_10, pupa_male_median_CI_95_upper_11, pupa_male_median_CI_95_upper_12, pupa_male_median_CI_95_upper_13, pupa_male_median_CI_95_upper_14, pupa_male_median_CI_95_upper_15, pupa_male_median_CI_95_upper_16, pupa_male_median_CI_95_upper_17, pupa_male_median_CI_95_upper_18, pupa_male_median_CI_95_upper_19, pupa_male_median_CI_95_upper_20,
                         adult_female_median_CI_95_lower_1, adult_female_median_CI_95_lower_2, adult_female_median_CI_95_lower_3, adult_female_median_CI_95_lower_4, adult_female_median_CI_95_lower_5, adult_female_median_CI_95_lower_6, adult_female_median_CI_95_lower_7, adult_female_median_CI_95_lower_8, adult_female_median_CI_95_lower_9, 
                         adult_female_median_CI_95_lower_10, adult_female_median_CI_95_lower_11, adult_female_median_CI_95_lower_12, adult_female_median_CI_95_lower_13, adult_female_median_CI_95_lower_14, adult_female_median_CI_95_lower_15, adult_female_median_CI_95_lower_16, adult_female_median_CI_95_lower_17, adult_female_median_CI_95_lower_18, adult_female_median_CI_95_lower_19, adult_female_median_CI_95_lower_20,
                         adult_female_median_CI_95_upper_1, adult_female_median_CI_95_upper_2, adult_female_median_CI_95_upper_3, adult_female_median_CI_95_upper_4, adult_female_median_CI_95_upper_5, adult_female_median_CI_95_upper_6, adult_female_median_CI_95_upper_7, adult_female_median_CI_95_upper_8, adult_female_median_CI_95_upper_9, 
                         adult_female_median_CI_95_upper_10, adult_female_median_CI_95_upper_11, adult_female_median_CI_95_upper_12, adult_female_median_CI_95_upper_13, adult_female_median_CI_95_upper_14, adult_female_median_CI_95_upper_15, adult_female_median_CI_95_upper_16, adult_female_median_CI_95_upper_17, adult_female_median_CI_95_upper_18, adult_female_median_CI_95_upper_19, adult_female_median_CI_95_upper_20,
                         adult_male_median_CI_95_lower_1, adult_male_median_CI_95_lower_2, adult_male_median_CI_95_lower_3, adult_male_median_CI_95_lower_4, adult_male_median_CI_95_lower_5, adult_male_median_CI_95_lower_6, adult_male_median_CI_95_lower_7, adult_male_median_CI_95_lower_8, adult_male_median_CI_95_lower_9, 
                         adult_male_median_CI_95_lower_10, adult_male_median_CI_95_lower_11, adult_male_median_CI_95_lower_12, adult_male_median_CI_95_lower_13, adult_male_median_CI_95_lower_14, adult_male_median_CI_95_lower_15, adult_male_median_CI_95_lower_16, adult_male_median_CI_95_lower_17, adult_male_median_CI_95_lower_18, adult_male_median_CI_95_lower_19, adult_male_median_CI_95_lower_20,
                         adult_male_median_CI_95_upper_1, adult_male_median_CI_95_upper_2, adult_male_median_CI_95_upper_3, adult_male_median_CI_95_upper_4, adult_male_median_CI_95_upper_5, adult_male_median_CI_95_upper_6, adult_male_median_CI_95_upper_7, adult_male_median_CI_95_upper_8, adult_male_median_CI_95_upper_9, 
                         adult_male_median_CI_95_upper_10, adult_male_median_CI_95_upper_11, adult_male_median_CI_95_upper_12, adult_male_median_CI_95_upper_13, adult_male_median_CI_95_upper_14, adult_male_median_CI_95_upper_15, adult_male_median_CI_95_upper_16, adult_male_median_CI_95_upper_17, adult_male_median_CI_95_upper_18, adult_male_median_CI_95_upper_19, adult_male_median_CI_95_upper_20),
                       nrow = 20, ncol = 12)

colnames(bootstrap_CI) <- c("Instar V Female, lower95", "Instar V Female, upper95", 
                            "Instar V Male, lower95", "Instar V Male, upper95",
                            "Pupa Female, lower95", "Pupa Female, upper95",
                            "Pupa Male, lower95", "Pupa Male, upper95",
                            "Adult Female, lower95", "Adult Female, upper95",
                            "Adult Male, lower95", "Adult Male, upper95")

rownames(bootstrap_CI) <- c("Z/A1", "Z/A2", "Z/A3", "Z/A4", "Z/A5", "Z/A6", "Z/A7", "Z/A8", 
                   "Z/A9", "Z/A10", "Z/A11", "Z/A12", "Z/A13", "Z/A14", "Z/A15", "Z/A16", 
                   "Z/A17", "Z/A18", "Z/A19", "Z/A20")

write.table(bootstrap_CI, file = table_4_out, sep = "\t", col.names = NA, row.names = TRUE)





### Plots

# female

chromosome_list <- c(1:20)
instar_V_f_ratio_z_to_A_median <- rep(-1000,20)
instar_V_f_ratio_z_to_A_lower <- rep(-1000,20)
instar_V_f_ratio_z_to_A_upper <- rep(-1000,20)
pupa_f_ratio_z_to_A_median <- rep(-1000,20)
pupa_f_ratio_z_to_A_lower <- rep(-1000,20)
pupa_f_ratio_z_to_A_upper <- rep(-1000,20)
adult_f_ratio_z_to_A_median <- rep(-1000,20)
adult_f_ratio_z_to_A_lower <- rep(-1000,20)
adult_f_ratio_z_to_A_upper <- rep(-1000,20)


for (i in c(1:20)) {
  instar_V_f_ratio_z_to_A_median[i] <- eval(parse(text = paste("instar_V_f_ratio_z_to_", toString(i), sep="")))
  instar_V_f_ratio_z_to_A_lower[i] <- eval(parse(text = paste("instar_V_female_median_CI_95_lower_", toString(i), sep="")))
  instar_V_f_ratio_z_to_A_upper[i] <- eval(parse(text = paste("instar_V_female_median_CI_95_upper_", toString(i), sep="")))
  
  pupa_f_ratio_z_to_A_median[i] <- eval(parse(text = paste("pupa_f_ratio_z_to_", toString(i), sep="")))
  pupa_f_ratio_z_to_A_lower[i] <- eval(parse(text = paste("pupa_female_median_CI_95_lower_", toString(i), sep="")))
  pupa_f_ratio_z_to_A_upper[i] <- eval(parse(text = paste("pupa_female_median_CI_95_upper_", toString(i), sep="")))
  
  adult_f_ratio_z_to_A_median[i] <- eval(parse(text = paste("adult_f_ratio_z_to_", toString(i), sep="")))
  adult_f_ratio_z_to_A_lower[i] <- eval(parse(text = paste("adult_female_median_CI_95_lower_", toString(i), sep="")))
  adult_f_ratio_z_to_A_upper[i] <- eval(parse(text = paste("adult_female_median_CI_95_upper_", toString(i), sep="")))
  
}

par(pty="s", mfrow = c(1,3))

plot(instar_V_f_ratio_z_to_A_median,
     ylim=range(c(0, 1.3)),
     pch=19, ylab="", xlab="",
     main="Larva", cex.main = 2,
     frame.plot = FALSE, outline = FALSE, xaxt = "n"
     )
axis(2, labels = "Z:A median ratio", cex.axis = 1.5, at = 0.6, line = 1.5, tck = 0)
#axis(3, labels = "C", cex.axis = 3, at = 1, tck = 0)
axis(1, las = 2, cex.axis = 1.1, labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
                                            "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"),
     at = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))

arrows(chromosome_list, instar_V_f_ratio_z_to_A_lower, 
       chromosome_list, instar_V_f_ratio_z_to_A_upper, 
       length=0.05, angle=90, code=3)

segments(0, 0.5, 21, 0.5, col = "darkcyan", lty = 2, lwd = 2)
segments(0, 1, 21, 1, col = "darkcyan", lty = 2, lwd = 2)


plot(pupa_f_ratio_z_to_A_median,
     ylim=range(c(0, 1.3)),
     pch=19, ylab="",
     main="Pupa", cex.main = 2,
     frame.plot = FALSE, outline = FALSE, xaxt = "n", 
     xlab = "Female chromosomes", cex.lab=2
)


axis(1, las = 2, cex.axis = 1.1, labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
                                            "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"),
     at = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))

arrows(chromosome_list, pupa_f_ratio_z_to_A_lower, 
       chromosome_list, pupa_f_ratio_z_to_A_upper, 
       length=0.05, angle=90, code=3)

segments(0, 0.5, 21, 0.5, col = "darkcyan", lty = 2, lwd = 2)
segments(0, 1, 21, 1, col = "darkcyan", lty = 2, lwd = 2)

#axis(1, labels = c("Female autosomes"), lwd = 0, at = 11.2, cex.axis = 2.0, tck = 0, line = 4)


plot(adult_f_ratio_z_to_A_median,
     ylim=range(c(0, 1.3)),
     pch=19, ylab="", xlab="",
     main="Adult", cex.main = 2,
     frame.plot = FALSE, outline = FALSE, xaxt = "n"
)


axis(1, las = 2, cex.axis = 1.1, labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
                                            "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"),
     at = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))

arrows(chromosome_list, adult_f_ratio_z_to_A_lower, 
       chromosome_list, adult_f_ratio_z_to_A_upper, 
       length=0.05, angle=90, code=3)

segments(0, 0.5, 21, 0.5, col = "darkcyan", lty = 2, lwd = 2)
segments(0, 1, 21, 1, col = "darkcyan", lty = 2, lwd = 2)


# male

chromosome_list <- c(1:20)
instar_V_m_ratio_z_to_A_median <- rep(-1000,20)
instar_V_m_ratio_z_to_A_lower <- rep(-1000,20)
instar_V_m_ratio_z_to_A_upper <- rep(-1000,20)
pupa_m_ratio_z_to_A_median <- rep(-1000,20)
pupa_m_ratio_z_to_A_lower <- rep(-1000,20)
pupa_m_ratio_z_to_A_upper <- rep(-1000,20)
adult_m_ratio_z_to_A_median <- rep(-1000,20)
adult_m_ratio_z_to_A_lower <- rep(-1000,20)
adult_m_ratio_z_to_A_upper <- rep(-1000,20)


for (i in c(1:20)) {
  instar_V_m_ratio_z_to_A_median[i] <- eval(parse(text = paste("instar_V_m_ratio_z_to_", toString(i), sep="")))
  instar_V_m_ratio_z_to_A_lower[i] <- eval(parse(text = paste("instar_V_male_median_CI_95_lower_", toString(i), sep="")))
  instar_V_m_ratio_z_to_A_upper[i] <- eval(parse(text = paste("instar_V_male_median_CI_95_upper_", toString(i), sep="")))
  
  pupa_m_ratio_z_to_A_median[i] <- eval(parse(text = paste("pupa_m_ratio_z_to_", toString(i), sep="")))
  pupa_m_ratio_z_to_A_lower[i] <- eval(parse(text = paste("pupa_male_median_CI_95_lower_", toString(i), sep="")))
  pupa_m_ratio_z_to_A_upper[i] <- eval(parse(text = paste("pupa_male_median_CI_95_upper_", toString(i), sep="")))
  
  adult_m_ratio_z_to_A_median[i] <- eval(parse(text = paste("adult_m_ratio_z_to_", toString(i), sep="")))
  adult_m_ratio_z_to_A_lower[i] <- eval(parse(text = paste("adult_male_median_CI_95_lower_", toString(i), sep="")))
  adult_m_ratio_z_to_A_upper[i] <- eval(parse(text = paste("adult_male_median_CI_95_upper_", toString(i), sep="")))
  
}

par(pty="s", mfrow = c(1,3))

plot(instar_V_m_ratio_z_to_A_median,
     ylim=range(c(0, 1.3)),
     pch=19, ylab="", xlab="",
     main="Larva", cex.main = 2,
     frame.plot = FALSE, outline = FALSE, xaxt = "n"
)

axis(2, labels = "Z:A median ratio", cex.axis = 1.5, at = 0.6, line = 1.5, tck = 0)
#axis(3, labels = "D", cex.axis = 3, at = 1, tck = 0)
axis(1, las = 2, cex.axis = 1.1, labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
                                            "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"),
     at = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))

arrows(chromosome_list, instar_V_m_ratio_z_to_A_lower, 
       chromosome_list, instar_V_m_ratio_z_to_A_upper, 
       length=0.05, angle=90, code=3)

segments(0, 0.5, 21, 0.5, col = "darkcyan", lty = 2, lwd = 2)
segments(0, 1, 21, 1, col = "darkcyan", lty = 2, lwd = 2)


plot(pupa_m_ratio_z_to_A_median,
     ylim=range(c(0, 1.3)),
     pch=19, ylab="",
     main="Pupa", cex.main = 2,
     frame.plot = FALSE, outline = FALSE, xaxt = "n",
     xlab = "Male chromosomes", cex.lab=2
)


axis(1, las = 2, cex.axis = 1.1, labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
                                            "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"),
     at = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))

arrows(chromosome_list, pupa_m_ratio_z_to_A_lower, 
       chromosome_list, pupa_m_ratio_z_to_A_upper, 
       length=0.05, angle=90, code=3)

segments(0, 0.5, 21, 0.5, col = "darkcyan", lty = 2, lwd = 2)
segments(0, 1, 21, 1, col = "darkcyan", lty = 2, lwd = 2)

#axis(1, labels = c("Male autosomes"), lwd = 0, at = 11.2, cex.axis = 2.0, tck = 0, line = 4)


plot(adult_m_ratio_z_to_A_median,
     ylim=range(c(0, 1.3)),
     pch=19, ylab="", xlab="",
     main="Adult", cex.main = 2,
     frame.plot = FALSE, outline = FALSE, xaxt = "n"
)


axis(1, las = 2, cex.axis = 1.1, labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
                                            "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"),
     at = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))

arrows(chromosome_list, adult_m_ratio_z_to_A_lower, 
       chromosome_list, adult_m_ratio_z_to_A_upper, 
       length=0.05, angle=90, code=3)

segments(0, 0.5, 21, 0.5, col = "darkcyan", lty = 2, lwd = 2)
segments(0, 1, 21, 1, col = "darkcyan", lty = 2, lwd = 2)


##############################################
