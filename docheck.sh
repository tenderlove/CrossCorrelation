#!/bin/bash
python wpar_check.py -filt=check -nbins=14 -ncores=23 -wtype=DD
python wpar_check.py -filt=check -nbins=14 -ncores=23 -wtype=RR
python wpar_check.py -filt=check -nbins=14 -ncores=23 -wtype=DR
