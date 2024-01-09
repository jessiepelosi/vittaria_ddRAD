'''
RunDadiOptimization.py

Purpose: run optimization steps for dadi 
Usage: python RunDadiOptimization.py -i [folded_sfs.txt] -r [number of rounds] -m [dadi demographic model]

Model choices: 
 - Neutral 
 - Exponential Growth
 - Three Epoch
 - Three Epoch with Inbreeding
 - Two Epoch
 - Bottlegrowth 

Note that this must be run in the same directory as Optimize_Functions.py
Optimization script from: Portik et al. 2017, see https://github.com/dportik/dadi_pipeline/tree/master

Jessie Pelosi 2024
'''


import sys
import os
import numpy
import dadi
from datetime import datetime
import Optimize_Functions
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, help="path to FOLDED SFS in dadi format from easysfs", dest="input")
parser.add_argument('-r', '--rounds', required=True, help="number of rounds of optimization", dest="rounds")
parser.add_argument('-m', '--model', help="dadi demographic model", choices=['neutral','exponential_growth', 'three_epoch', 'three_epoch_inbreeding','two_epoch', 'bottlegrowth'], dest="model")

args = parser.parse_args()
INFILE=args.input
ROUNDS=args.rounds
MODEL=args.model

fs = dadi.Spectrum.from_file(INFILE)

#print some useful information about the afs or jsfs
print("\n\n============================================================================")
print("\nData for site frequency spectrum:\n")
#print("Projection: {}".format(proj))
print("Sample sizes: {}".format(fs.sample_sizes))
print("Sum of SFS: {}".format(numpy.around(fs.S(), 2)))
print("\n============================================================================\n")

if MODEL=="two_epoch":

	# Define function - this can be changed and must be changed for each different model 
	func = dadi.Demographics1D.two_epoch
	param_names = ("nu", "T") 

	# For two-epoch model: 
	# nu specifies the ratio of the contemporary population size to the ancient population size 
	# T is the time in the past at which the size change happened in units 2*Na generations 

	#create a prefix to label the output files
	prefix = "TwoEpoch"
                             
	# Define grid points based on sample size
	# For smaller data (largest sample size < 100) [ns+20, ns+30, ns+40] is a good starting point 
	# From dadi manual

	ns = fs.sample_sizes
	pts = [max(ns)+20, max(ns)+30, max(ns)+40]

	#Remember the order for mandatory arguments as below
	#Optimize_Routine(fs, pts, outfile, model_name, func, rounds, param_number, fs_folded)
	Optimize_Functions.Optimize_Routine(fs, pts, prefix, "two_epoch", func, ROUNDS, 2, fs_folded=True)

elif MODEL=="neutral":
	# Define function - this can be changed and must be changed for each different model 
	func = dadi.Demographics1D.snm
	param_names = ("N") 

	# N specifies the constant population size 

	#create a prefix to label the output files
	prefix = "Neutral"
                             
	# Define grid points based on sample size
	# For smaller data (largest sample size < 100) [ns+20, ns+30, ns+40] is a good starting point 
	# From dadi manual

	ns = fs.sample_sizes
	pts = [max(ns)+20, max(ns)+30, max(ns)+40]

	#Remember the order for mandatory arguments as below
	#Optimize_Routine(fs, pts, outfile, model_name, func, rounds, param_number, fs_folded)
	Optimize_Functions.Optimize_Routine(fs, pts, prefix, "Neutral", func, ROUNDS, 1, fs_folded=True)

elif MODEL=="exponential_growth":
	# Define function - this can be changed and must be changed for each different model 
	func = dadi.Demographics1D.growth
	param_names = ("nu", "T") 

	# At time T in the past, an equilibrium population begins growing exponentially, reaching size nu at present

	#create a prefix to label the output files
	prefix = "exponential_growth"
                             
	# Define grid points based on sample size
	# For smaller data (largest sample size < 100) [ns+20, ns+30, ns+40] is a good starting point 
	# From dadi manual

	ns = fs.sample_sizes
	pts = [max(ns)+20, max(ns)+30, max(ns)+40]

	#Remember the order for mandatory arguments as below
	#Optimize_Routine(fs, pts, outfile, model_name, func, rounds, param_number, fs_folded)
	Optimize_Functions.Optimize_Routine(fs, pts, prefix, "exponential_growth", func, ROUNDS, 2, fs_folded=True)

elif MODEL=="three_epoch":
	# Define function - this can be changed and must be changed for each different model 
	func = dadi.Demographics1D.three_epoch
	param_names = ("nuB", "nuF", "TB", "TF") 

	# At time T in the past, an equilibrium population begins growing exponentially, reaching size nu at present

	#create a prefix to label the output files
	prefix = "three_epoch"
                             
	# Define grid points based on sample size
	# For smaller data (largest sample size < 100) [ns+20, ns+30, ns+40] is a good starting point 
	# From dadi manual

	ns = fs.sample_sizes
	pts = [max(ns)+20, max(ns)+30, max(ns)+40]

	#Remember the order for mandatory arguments as below
	#Optimize_Routine(fs, pts, outfile, model_name, func, rounds, param_number, fs_folded)
	Optimize_Functions.Optimize_Routine(fs, pts, prefix, "exponential_growth", func, ROUNDS, 4, fs_folded=True)

elif MODEL=="three_epoch_inbreeding":
	# Define function - this can be changed and must be changed for each different model 
	func = dadi.Demographics1D.three_epoch_inbreeding
	param_names = ("nuB", "nuF", "TB", "TF", "F")  

	# At time T in the past, an equilibrium population begins growing exponentially, reaching size nu at present

	#create a prefix to label the output files
	prefix = "three_epoch_inbreeding"
                             
	# Define grid points based on sample size
	# For smaller data (largest sample size < 100) [ns+20, ns+30, ns+40] is a good starting point 
	# From dadi manual

	ns = fs.sample_sizes
	pts = [max(ns)+20, max(ns)+30, max(ns)+40]

	#Remember the order for mandatory arguments as below
	#Optimize_Routine(fs, pts, outfile, model_name, func, rounds, param_number, fs_folded)
	Optimize_Functions.Optimize_Routine(fs, pts, prefix, "three_epoch_inbreeding", func, ROUNDS, 5, fs_folded=True)

elif MODEL=="bottlegrowth":
	# Define function - this can be changed and must be changed for each different model 
	func = dadi.Demographics1D.bottlegrowth
	param_names = ("nuB", "nuF", "T") 

	# At time T in the past, an equilibrium population begins growing exponentially, reaching size nu at present

	#create a prefix to label the output files
	prefix = "bottlegrowth"
                             
	# Define grid points based on sample size
	# For smaller data (largest sample size < 100) [ns+20, ns+30, ns+40] is a good starting point 
	# From dadi manual

	ns = fs.sample_sizes
	pts = [max(ns)+20, max(ns)+30, max(ns)+40]

	#Remember the order for mandatory arguments as below
	#Optimize_Routine(fs, pts, outfile, model_name, func, rounds, param_number, fs_folded)
	Optimize_Functions.Optimize_Routine(fs, pts, prefix, "bottlegrowth", func, ROUNDS, 3, fs_folded=True)

else:
	print("Model not recognized.")
