#!/bin/bash
#$ -N initial_scan_cd4memory
#$ -j y
#$ -o ../../../../../out/dvisR_analysis/kevin_tmp/Writeup1/qsub/
#$ -l m_mem_free=25G

Rscript --no-save Writeup1_initial_scan_cd4memory.R