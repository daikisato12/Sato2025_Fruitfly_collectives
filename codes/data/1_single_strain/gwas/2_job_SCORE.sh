#!/bin/sh
#$ -S /bin/bash
#$ -cwd
# for SCORE
# qsub -l medium -l s_vmem=20G -l mem_req=20G 2_job_SCORE.sh

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

DATADIR=../../../../data/1_single_strain/gwas/
QC=../../../../data/0_genome/dgrp2_QC
SCORE_PREP=${DATADIR}/SCORE_prep/
SCORE_RES=${DATADIR}/result/heritability/SCORE/
# SCORE: https://github.com/sriramlab/SCORE
SCORE=${SOFTWARE}/SCORE/build/SCORE

LIST=(`cat $DATADIR/pheno/2traits_comb.txt`)
for PREFIX in "${LIST[@]}"; do
	PHENO=${SCORE_PREP}/bivar_${PREFIX}.txt
	$SCORE -g $QC -p $PHENO  -o ${SCORE_RES}/cor_$PREFIX -b 10 -mpheno 1,2
done

echo ending at
date

