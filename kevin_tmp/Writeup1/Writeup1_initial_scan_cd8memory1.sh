#!/bin/bash
#$ -N initial_scan_cd8memory1
#$ -j y
#$ -o ../../../../../out/dvisR_analysis/kevin_tmp/Writeup1/qsub/
#$ -l m_mem_free=25G

Rscript --no-save Writeup1_initial_scan_cd8memory1.R
