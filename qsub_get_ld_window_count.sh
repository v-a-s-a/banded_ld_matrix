#!/bin/bash
#$ -l h_vmem=16G
#$ -cwd
#$ -t 1-22

source /broad/software/scripts/useuse
use Anaconda3
source activate gear

./get_ld_window_count.py --chr ${SGE_TASK_ID}
