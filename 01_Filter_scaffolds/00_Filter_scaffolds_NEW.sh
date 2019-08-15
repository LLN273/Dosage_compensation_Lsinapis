#!/bin/bash





##########  filter_scaffolds_by_mapping_length


SCRIPT_FOLDER_1=$(pwd)
OUT_FOLDER_1=/crex1/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/01_Filter_scaffolds
LS2HM_CHROMOSOME_MAPPING=/home/luisleal/MYPROJ/3_DosageCompensation_LS/Scripts_NEW/00_Chromosome_mapping_LS2HM/chromosomes_all_info.txt


cd $OUT_FOLDER_1

perl $SCRIPT_FOLDER_1/1_filter_scaffolds_by_mapping_length.pl $LS2HM_CHROMOSOME_MAPPING

#OUTPUT: mapped_scaffolds_filtered_above_200.txt




## assign_scaffolds_a_or_z  (to be used in M:F analysis)

perl $SCRIPT_FOLDER_1/2_A_assign_scaffolds_a_or_z.pl $OUT_FOLDER_1/mapped_scaffolds_filtered_above_200.txt

#OUTPUT: mapped_scaffolds_filtered_above_200_assigned_a_or_z_by_0.95.txt



# assign_scaffolds_a_or_z  (to be used in Z:A analysis)

perl $SCRIPT_FOLDER_1/2_B_mapped_scaffolds_cutoff.pl $OUT_FOLDER_1/mapped_scaffolds_filtered_above_200.txt

#OUTPUT: mapped_scaffolds_filtered_above_200_assigned_to_chromosome_by_0.9.txt










