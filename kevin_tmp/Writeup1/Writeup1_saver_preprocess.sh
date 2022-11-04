#!/bin/bash
#$ -N bm_saver
#$ -j y
#$ -o ../../../../../out/dvisR_analysis/kevin_tmp/Writeup1/qsub/
#$ -l m_mem_free=25G
#$ -pe openmp 4

Rscript --no-save Writeup1_saver_preprocess.R
