#!/bin/bash
#$ -l h_vmem=64G
#$ -cwd
#$ -t 1-21

source /broad/software/scripts/useuse
use Anaconda3
source activate gear

./ska_chunked_ld.py --chr ${SGE_TASK_ID}
