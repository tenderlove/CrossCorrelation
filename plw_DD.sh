#!/bin/bash
#$ -q free64
#$ -m beas
#$ -R y
#$ -pe one-node-mpi 64
#$ -N plw_DD
module load python/2.7.10
python wpar.py -ncores=64 -cat=u -filt=plw -nbins=16 -wtype=DD > log_plw