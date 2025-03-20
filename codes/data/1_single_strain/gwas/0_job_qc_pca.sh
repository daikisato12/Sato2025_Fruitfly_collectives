#!/bin/sh
#$ -S /bin/bash
#$ -cwd
# for PCA
# qsub -pe def_slot 10 0_job_pca.sh

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

DATADIR=../../../../data/1_single_strain/gwas/
IND=${DATADIR}/prep/analyzed_inds.txt
SNP=${DATADIR}/prep/analyzed_SNPs.txt
GENO=../../../../data/0_genome/dgrp2

export OPENBLAS_NUM_THREADS=2

# plink: https://www.cog-genomics.org/plink/1.9/
plink --bfile $GENO --keep analyzed_inds.txt --extract analyzed_SNPs.txt --maf 0.01 --geno 0.1 --make-bed --out dgrp2_QC

# plink2: https://www.cog-genomics.org/plink/2.0/
plink2 --bfile $GENO --pca 10 --threads 10 --out $GENO

echo ending at
date

