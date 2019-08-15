#!/bin/bash
#

#module load perl


## INPUT DATA FOLDER
#AA=/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/03_Normalized_libs
AA=/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/12_Filter_sex_biased_genes

## INPUT FILES
#IF_L=instar_V-assigned_A_or_Z-filtered.txt
#IF_P=pupa-assigned_A_or_Z-filtered.txt
#IF_A=adult-assigned_A_or_Z-filtered.txt

IF_L=nonbiased_genes-instar_V-assigned_A_or_Z-filtered.txt
IF_P=nonbiased_genes-pupa-assigned_A_or_Z-filtered.txt
IF_A=nonbiased_genes-adult-assigned_A_or_Z-filtered.txt


## OUTPUT FOLDER
#RR=/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/18A_Quartile_analysis
RR=/crex/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/18B_Quartile_analysis_AFTER_REMOVE_SEX_BIASED_GENES




# current path
SRCDIR_INI=$(pwd)  






cd $RR

rm -f sorted*
rm -f quartiles*

cp $AA/$IF_L .
cp $AA/$IF_P .
cp $AA/$IF_A .


perl $SRCDIR_INI/quartile_analysis_1.pl $IF_L 
perl $SRCDIR_INI/quartile_analysis_1.pl $IF_P
perl $SRCDIR_INI/quartile_analysis_1.pl $IF_A


## The quartile values (eg 66 66 66 66) are the number of genes in file sorted-*****-assigned_A_or_Z-filtered.txt divided by four

## Before removing SBG (m-f_75perc)
#perl $SRCDIR_INI/quartile_analysis_2.pl sorted-${IF_L} 66 66 66 66
#perl $SRCDIR_INI/quartile_analysis_2.pl sorted-${IF_P} 62 63 63 63
#perl $SRCDIR_INI/quartile_analysis_2.pl sorted-${IF_A} 61 61 61 61

## After removing SBG (m-f_75perc)
perl $SRCDIR_INI/quartile_analysis_2.pl sorted-${IF_L} 63 63 63 64
perl $SRCDIR_INI/quartile_analysis_2.pl sorted-${IF_P} 59 60 60 60
perl $SRCDIR_INI/quartile_analysis_2.pl sorted-${IF_A} 45 45 45 46




rm -f $RR/$IF_L
rm -f $RR/$IF_P
rm -f $RR/$IF_A