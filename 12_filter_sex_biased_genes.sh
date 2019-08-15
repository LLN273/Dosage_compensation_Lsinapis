#!/bin/bash

#Lars Höök 2017

#Sex-biased genes are removed from genes filtered FPKM>0 - for dosage compensation analysis
#For Fisher test (expected distribution), ALL genes filtered in DEseq2 are assign as A or Z (including SBG)
#Male and female biased genes are picked from list of all genes assigned A or Z, but not filtered FPKM>0
#Gene names are added to sex biased genes for volcano plot




### input folders

# DESEQ2 data
AA=/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/10_DESEQ2_sva

# Normamiled libs with genes assigned to chromosomes
BB=/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/03_Normalized_libs	  



### result folder
RR=/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/12_Filter_sex_biased_genes



### Process data

cd $RR

##

STAGE1="instar_V-assigned_A_or_Z.txt"
FILTERED1="instar_V-assigned_A_or_Z-filtered.txt"
MBG1="MBG_Instar.txt"
FBG1="FBG_Instar.txt"
SBG_STAGE1="Filtered_P0.05_Instar.txt"
ALL_DESEQ1="LFC_Instar.txt"
GENE_NAMES1="instar_V-concatenated.txt"

FILTERED1_f="instar_V-assigned_A_or_Z_female-filtered.txt"
FILTERED1_m="instar_V-assigned_A_or_Z_male-filtered.txt"

##

STAGE2="pupa-assigned_A_or_Z.txt"
FILTERED2="pupa-assigned_A_or_Z-filtered.txt"
MBG2="MBG_Pupa.txt"
FBG2="FBG_Pupa.txt"
SBG_STAGE2="Filtered_P0.05_Pupa.txt"
ALL_DESEQ2="LFC_Pupa.txt"
GENE_NAMES2="pupa-concatenated.txt"

FILTERED2_f="pupa-assigned_A_or_Z_female-filtered.txt"
FILTERED2_m="pupa-assigned_A_or_Z_male-filtered.txt"

##

STAGE3="adult-assigned_A_or_Z.txt"
FILTERED3="adult-assigned_A_or_Z-filtered.txt"
MBG3="MBG_Adult.txt"
FBG3="FBG_Adult.txt"
SBG_STAGE3="Filtered_P0.05_Adult.txt"
ALL_DESEQ3="LFC_Adult.txt"
GENE_NAMES3="adult-concatenated.txt"

FILTERED3_f="adult-assigned_A_or_Z_female-filtered.txt"
FILTERED3_m="adult-assigned_A_or_Z_male-filtered.txt"

##

HEADER_MASTER="gene_id\tgene_name\tscaffold\tFPKM_instar_V_female\tFPKM_instar_V_male\tchromosome"
echo -e $HEADER_MASTER > $RR/FisherTest-$STAGE1
echo -e $HEADER_MASTER > $RR/FisherTest_MBG-$STAGE1
echo -e $HEADER_MASTER > $RR/FisherTest_FBG-$STAGE1

HEADER_MASTER="gene_id\tgene_name\tscaffold\tFPKM_pupa_female\tFPKM_pupa_male\tchromosome"
echo -e $HEADER_MASTER > $RR/FisherTest-$STAGE2
echo -e $HEADER_MASTER > $RR/FisherTest_MBG-$STAGE2
echo -e $HEADER_MASTER > $RR/FisherTest_FBG-$STAGE2

HEADER_MASTER="gene_id\tgene_name\tscaffold\tFPKM_adult_female\tFPKM_adult_male\tchromosome"
echo -e $HEADER_MASTER > $RR/FisherTest-$STAGE3
echo -e $HEADER_MASTER > $RR/FisherTest_MBG-$STAGE3
echo -e $HEADER_MASTER > $RR/FisherTest_FBG-$STAGE3

##



sort $AA/$MBG1 -o $RR/$MBG1
sort $AA/$FBG1 -o $RR/$FBG1

grep -vf <(cut -f1 $AA/$SBG_STAGE1 | sort) $BB/$FILTERED1 > $RR/nonbiased_genes-$FILTERED1
grep -vf <(cut -f1 $AA/$SBG_STAGE1 | sort) $BB/$FILTERED1_f > $RR/nonbiased_genes-$FILTERED1_f
grep -vf <(cut -f1 $AA/$SBG_STAGE1 | sort) $BB/$FILTERED1_m > $RR/nonbiased_genes-$FILTERED1_m
grep -vf <(cut -f1 $AA/$FBG1 | sort) $BB/$FILTERED1_f > $RR/non-female-biased_genes-$FILTERED1_f
grep -vf <(cut -f1 $AA/$MBG1 | sort) $BB/$FILTERED1_m > $RR/non-male-biased_genes-$FILTERED1_m


grep -f <(cut -f1 $AA/$ALL_DESEQ1 | sort) $BB/$STAGE1 >> $RR/FisherTest-$STAGE1

### these two lines are coded differently in pupae and adults...

#grep -f <(cut -f1 $RR/$MBG1) $BB/$STAGE1 >> $RR/FisherTest_MBG-$STAGE1
#grep -f <(cut -f1 $RR/$FBG1) $BB/$STAGE1 >> $RR/FisherTest_FBG-$STAGE1
grep -f <(cut -f1 $RR/$MBG1 | sort) $BB/$STAGE1 >> $RR/FisherTest_MBG-$STAGE1
grep -f <(cut -f1 $RR/$FBG1 | sort) $BB/$STAGE1 >> $RR/FisherTest_FBG-$STAGE1

##

paste <(grep -f <(cut -f1 $RR/$MBG1) <(cut -f1,2 $BB/$GENE_NAMES1)) <(cut -f2-8 $RR/$MBG1) > $RR/gene_names-$MBG1
paste <(grep -f <(cut -f1 $RR/$FBG1) <(cut -f1,2 $BB/$GENE_NAMES1)) <(cut -f2-8 $RR/$FBG1) > $RR/gene_names-$FBG1

#####



sort $AA/$MBG2 -o $RR/$MBG2
sort $AA/$FBG2 -o $RR/$FBG2

grep -vf <(cut -f1 $AA/$SBG_STAGE2 | sort) $BB/$FILTERED2 > $RR/nonbiased_genes-$FILTERED2
grep -vf <(cut -f1 $AA/$SBG_STAGE2 | sort) $BB/$FILTERED2_f > $RR/nonbiased_genes-$FILTERED2_f
grep -vf <(cut -f1 $AA/$SBG_STAGE2 | sort) $BB/$FILTERED2_m > $RR/nonbiased_genes-$FILTERED2_m
grep -vf <(cut -f1 $AA/$FBG2 | sort) $BB/$FILTERED2_f > $RR/non-female-biased_genes-$FILTERED2_f
grep -vf <(cut -f1 $AA/$MBG2 | sort) $BB/$FILTERED2_m > $RR/non-male-biased_genes-$FILTERED2_m


grep -f <(cut -f1 $AA/$ALL_DESEQ2 | sort) $BB/$STAGE2 >> $RR/FisherTest-$STAGE2
grep -f <(cut -f1 $RR/$MBG2 | sort) $BB/$STAGE2 >> $RR/FisherTest_MBG-$STAGE2
grep -f <(cut -f1 $RR/$FBG2 | sort) $BB/$STAGE2 >> $RR/FisherTest_FBG-$STAGE2

paste <(grep -f <(cut -f1 $RR/$MBG2) <(cut -f1,2 $BB/$GENE_NAMES2)) <(cut -f2-8 $RR/$MBG2) > $RR/gene_names-$MBG2
paste <(grep -f <(cut -f1 $RR/$FBG2) <(cut -f1,2 $BB/$GENE_NAMES2)) <(cut -f2-8 $RR/$FBG2) > $RR/gene_names-$FBG2

#####



sort $AA/$MBG3 -o $RR/$MBG3
sort $AA/$FBG3 -o $RR/$FBG3

grep -vf <(cut -f1 $AA/$SBG_STAGE3 | sort) $BB/$FILTERED3 > $RR/nonbiased_genes-$FILTERED3
grep -vf <(cut -f1 $AA/$SBG_STAGE3 | sort) $BB/$FILTERED3_f > $RR/nonbiased_genes-$FILTERED3_f
grep -vf <(cut -f1 $AA/$SBG_STAGE3 | sort) $BB/$FILTERED3_m > $RR/nonbiased_genes-$FILTERED3_m
grep -vf <(cut -f1 $AA/$FBG3 | sort) $BB/$FILTERED3_f > $RR/non-female-biased_genes-$FILTERED3_f
grep -vf <(cut -f1 $AA/$MBG3 | sort) $BB/$FILTERED3_m > $RR/non-male-biased_genes-$FILTERED3_m

grep -f <(cut -f1 $AA/$ALL_DESEQ3 | sort) $BB/$STAGE3 >> $RR/FisherTest-$STAGE3
grep -f <(cut -f1 $RR/$MBG3 | sort) $BB/$STAGE3 >> $RR/FisherTest_MBG-$STAGE3
grep -f <(cut -f1 $RR/$FBG3 | sort) $BB/$STAGE3 >> $RR/FisherTest_FBG-$STAGE3

paste <(grep -f <(cut -f1 $RR/$MBG3) <(cut -f1,2 $BB/$GENE_NAMES3)) <(cut -f2-8 $RR/$MBG3) > $RR/gene_names-$MBG3
paste <(grep -f <(cut -f1 $RR/$FBG3) <(cut -f1,2 $BB/$GENE_NAMES3)) <(cut -f2-8 $RR/$FBG3) > $RR/gene_names-$FBG3

