#!/bin/bash

### Script used to assign non-biased genes to individual chromosomes


### NOTE: might need to run dos2unix:
### dos2unix prepare_stringtie_files.sh

##################################################################

################################################

# 1: Path to folder containing FPKM counts after removing sex biased genes:

################################################

#Folder containing FPKM data after removing sex-biased genes
AA=/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/12_Filter_sex_biased_genes


#OUTPUT FOLDER
RR=/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/17_Assign_to_chromosomes_m_vs_f_SBG

#remember current path
SRCDIR=$(pwd)       

                                      
#######################################

# 2: List with assigned scaffolds

#######################################

ASSIGNED_A_OR_Z=/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/01_Filter_scaffolds/mapped_scaffolds_filtered_above_200_assigned_a_or_z_by_0.95.txt
ASSIGNED_CHR=/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/01_Filter_scaffolds/mapped_scaffolds_filtered_above_200_assigned_to_chromosome_by_0.9.txt


#############################

# 3: Name groups

#############################


STAGE1="instar_V"
STAGE2="pupa"
STAGE3="adult"




################################################

# 4: Reformat input files

################################################


cd $RR

date > prepare_stringtie_files.log

rsync -ah $AA/nonbiased_genes-*-assigned_A_or_Z_female-filtered.txt .
rsync -ah $AA/nonbiased_genes-*-assigned_A_or_Z_male-filtered.txt .

cat $RR/nonbiased_genes-${STAGE1}-assigned_A_or_Z_female-filtered.txt | cut --complement -f 6 > $RR/nonbiased_genes-${STAGE1}-assigned_A_or_Z_female-filtered_CLEAN.txt
cat $RR/nonbiased_genes-${STAGE2}-assigned_A_or_Z_female-filtered.txt | cut --complement -f 6 > $RR/nonbiased_genes-${STAGE2}-assigned_A_or_Z_female-filtered_CLEAN.txt
cat $RR/nonbiased_genes-${STAGE3}-assigned_A_or_Z_female-filtered.txt | cut --complement -f 6 > $RR/nonbiased_genes-${STAGE3}-assigned_A_or_Z_female-filtered_CLEAN.txt

cat $RR/nonbiased_genes-${STAGE1}-assigned_A_or_Z_male-filtered.txt | cut --complement -f 6 > $RR/nonbiased_genes-${STAGE1}-assigned_A_or_Z_male-filtered_CLEAN.txt
cat $RR/nonbiased_genes-${STAGE2}-assigned_A_or_Z_male-filtered.txt | cut --complement -f 6 > $RR/nonbiased_genes-${STAGE2}-assigned_A_or_Z_male-filtered_CLEAN.txt
cat $RR/nonbiased_genes-${STAGE3}-assigned_A_or_Z_male-filtered.txt | cut --complement -f 6 > $RR/nonbiased_genes-${STAGE3}-assigned_A_or_Z_male-filtered_CLEAN.txt

############################################

# 5:Run auxiliary scripts

############################################


echo assigning genes to chromosomes...

cd $RR

perl $SRCDIR/3b_assign_genes_to_chromosomes_IND-SEX.pl $ASSIGNED_CHR nonbiased_genes-${STAGE1}-assigned_A_or_Z_female-filtered_CLEAN.txt
perl $SRCDIR/3b_assign_genes_to_chromosomes_IND-SEX.pl $ASSIGNED_CHR nonbiased_genes-${STAGE1}-assigned_A_or_Z_male-filtered_CLEAN.txt

perl $SRCDIR/3b_assign_genes_to_chromosomes_IND-SEX.pl $ASSIGNED_CHR nonbiased_genes-${STAGE2}-assigned_A_or_Z_female-filtered_CLEAN.txt
perl $SRCDIR/3b_assign_genes_to_chromosomes_IND-SEX.pl $ASSIGNED_CHR nonbiased_genes-${STAGE2}-assigned_A_or_Z_male-filtered_CLEAN.txt

perl $SRCDIR/3b_assign_genes_to_chromosomes_IND-SEX.pl $ASSIGNED_CHR nonbiased_genes-${STAGE3}-assigned_A_or_Z_female-filtered_CLEAN.txt
perl $SRCDIR/3b_assign_genes_to_chromosomes_IND-SEX.pl $ASSIGNED_CHR nonbiased_genes-${STAGE3}-assigned_A_or_Z_male-filtered_CLEAN.txt




echo filtering out non-expressed genes...

perl $SRCDIR/5_filter_genes_individualSexes.pl nonbiased_genes-${STAGE1}-assigned_A_or_Z_female-assigned_to_chromosomes.txt assigned_to_chromosomes female nonbiased_genes-${STAGE1}-assigned_to_chromosomes_female-filtered.txt
perl $SRCDIR/5_filter_genes_individualSexes.pl nonbiased_genes-${STAGE2}-assigned_A_or_Z_female-assigned_to_chromosomes.txt assigned_to_chromosomes female nonbiased_genes-${STAGE2}-assigned_to_chromosomes_female-filtered.txt
perl $SRCDIR/5_filter_genes_individualSexes.pl nonbiased_genes-${STAGE3}-assigned_A_or_Z_female-assigned_to_chromosomes.txt assigned_to_chromosomes female nonbiased_genes-${STAGE3}-assigned_to_chromosomes_female-filtered.txt

perl $SRCDIR/5_filter_genes_individualSexes.pl nonbiased_genes-${STAGE1}-assigned_A_or_Z_male-assigned_to_chromosomes.txt assigned_to_chromosomes male nonbiased_genes-${STAGE1}-assigned_to_chromosomes_male-filtered.txt
perl $SRCDIR/5_filter_genes_individualSexes.pl nonbiased_genes-${STAGE2}-assigned_A_or_Z_male-assigned_to_chromosomes.txt assigned_to_chromosomes male nonbiased_genes-${STAGE2}-assigned_to_chromosomes_male-filtered.txt
perl $SRCDIR/5_filter_genes_individualSexes.pl nonbiased_genes-${STAGE3}-assigned_A_or_Z_male-assigned_to_chromosomes.txt assigned_to_chromosomes male nonbiased_genes-${STAGE3}-assigned_to_chromosomes_male-filtered.txt


### clean up
cd $RR
rm -f *_CLEAN.txt












#############################################################################
