#!/bin/sh
#$ -S /bin/bash
#$ -cwd
# for PCA
# qsub -pe def_slot 10 1_job_gwas.sh motion_cue_exit_intercept_female

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

DATADIR=../../../../data/1_single_strain/gwas/
IND=${DATADIR}/prep/analyzed_inds.txt
PREFIX=$1
PHENO=${DATADIR}/pheno/df_out_${PREFIX}_scaled_GREML.txt
GWAS_RES=${DATADIR}/result/rawdata
GENO=../../../../data/0_genome/dgrp2
# GCTA: https://yanglab.westlake.edu.cn/software/gcta/
GCTA=${SOFTWARE}/gcta_v1.94.0Beta_linux_kernel_4_x86_64/gcta_v1.94.0Beta_linux_kernel_4_x86_64_static

#$GCTA --bfile $GENO --keep $IND --make-grm --out $GENO
#$GCTA --grm $GENO --make-bK-sparse 0.05 --out $GENO\_sp
#$GCTA --grm $GENO --pca 10 --out $GENO\_pc
$GCTA --bfile $GENO --grm-sparse $GENO\_sp --pheno $PHENO --maf 0.01 --geno 0.1 --qcovar $GENO\_pc.eigenvec --fastGWA-mlm --out $GWAS_RES/$PREFIX --thread-num 10

echo ending at
date

