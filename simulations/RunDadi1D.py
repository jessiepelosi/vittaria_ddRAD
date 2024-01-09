'''
RunDadi1D.py

Purpose: Run dadi given parameters for a two epoch model 
Usage: python RunDadi1D.py -i [folded SFS] -r [run number] -o [output directory]
								-m [mutation rate] -l [number of sites used to generated SFS, mono- + polymorphic]
								-n [optimized nu] -t [optimized t]

First run RunDadiOptimization.py to get optimized starting parameters 

Modified from: https://github.com/tanyaphung/dadi_tutorial/tree/master

Jessie Pelosi 2024
'''

import matplotlib
matplotlib.use('Agg') # so it doesn't pop up graphics on the cluster
import sys
import argparse
import dadi
from dadi import Numerics, PhiManip, Integration, Spectrum, Demographics1D
from numpy import array
import datetime
import os
todaysdate=datetime.datetime.today().strftime('%Y%m%d')

####################
# Parse input args #
####################

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, help="path to FOLDED SFS in dadi format from easysfs", dest="input")
parser.add_argument('-n', '--iteration', required=True, help="number of iterations", dest="rounds")
parser.add_argument('-m', '--model', help="dadi demographic model", choices=['neutral','exponential_growth', 'three_epoch', 'three_epoch_inbreeding','two_epoch', 'bottlegrowth'], dest="model")
parser.add_argument('-r', '--mutrate', help="mutation rate", dest="rate")
parser.add_argument('-l', '--length', help="number of sites used to generate SFS, mono- and polymorphic", dest = "length")
parser.add_argument('-o', '--output', required=True,help="path to output directory", dest="out")

# Optimized values for each model 
# For neutral model 
parser.add_argument('-N', help="Optimized value of N for neutral model", dest="N")
# For two epoch and exponential growth 
parser.add_argument('--nu', help="Optimized value of nu for two epoch and exponential growth models", dest="NU")
parser.add_argument('--T', help="Optimized value of t for two epoch and exponential growth models", dest="T")
# For bottlegrowth
parser.add_argument('--nuB', help="Optimized value of NuB", dest="NUB")
parser.add_argument('--nuF', help="Optimized value of NuF", dest="NUF")
# For three epoch
parser.add_argument('--TB', help="Optimized value of TB", dest="TB")
parser.add_argument('--TF', help="Optimized value of TF", dest="TF")
# For three epoch with inbreeding 
parser.add_argument('--F', help="Optimized value of F", dest="F")

args = parser.parse_args()
MODEL=args.model
sfs=str(args.input)
runNum=str(args.rounds)
mu=float(args.rate)
L=float(args.length)
outdir=str(args.out)
maxiter=100

# 01. Read input folded SFS
print("Input file:", sfs)
fs = dadi.Spectrum.from_file(sfs)

if MODEL=="neutral":
	N=float(args.N)
	# 02. Define model and extrapolate function
	func = dadi.Demographics1D.snm
	model_ex = dadi.Numerics.make_extrap_func(func)

	# 03. Define parameter space
	popt=(N)
	ns = fs.sample_sizes
	pts = [max(ns)+20, max(ns)+30, max(ns)+40]

	# 04. Calculate the best fit model
	model_fs = model_ex(popt, ns, pts)

	# 05. Compute likelihood of the data given the model allele frequency spectrum 
	ll_model = dadi.Inference.ll_multinom(model_fs, fs)

	# 06. Compute likelihood of data to itself (best possible log-likelihood)
	ll_data = dadi.Inference.ll_multinom(fs, fs)

	# 07. Calculate the best fit theta
	theta = dadi.Inference.optimal_sfs_scaling(model_fs, fs)

	# 08. Model-specific scaling of parameters 
	Nanc = theta / (4*mu*L)
	scaled_param_names = ("Nanc_FromTheta_scaled_dip")
	scaled_popt = (Nanc)

elif MODEL=="two_epoch":
	nu=float(args.NU)
	t=float(args.T)
	# 02. Define model and extrapolate function
	func = dadi.Demographics1D.two_epoch
	model_ex = dadi.Numerics.make_extrap_func(func)

	# 03. Define parameter space
	popt=(nu, t)
	ns = fs.sample_sizes
	pts = [max(ns)+20, max(ns)+30, max(ns)+40]

	# 04. Calculate the best fit model
	model_fs = model_ex(popt, ns, pts)

	# 05. Compute likelihood of the data given the model allele frequency spectrum 
	ll_model = dadi.Inference.ll_multinom(model_fs, fs)

	# 06. Compute likelihood of data to itself (best possible log-likelihood)
	ll_data = dadi.Inference.ll_multinom(fs, fs)

	# 07. Calculate the best fit theta
	theta = dadi.Inference.optimal_sfs_scaling(model_fs, fs)

	# 08. Model-specific scaling of parameters 
	Nanc = theta / (4*mu*L)
	nu_scaled_dip = popt[0]*Nanc
	T_scaled_dip = popt[1]*2*Nanc
	scaled_param_names = ("Nanc_FromTheta_scaled_dip","nu_scaled_dip","T_scaled_gen")
	scaled_popt = (Nanc, nu_scaled_dip, T_scaled_dip)


elif MODEL=="exponential_growth":
	nu=float(args.NU)
	t=float(args.T)
	# 02. Define model and extrapolate function
	func = dadi.Demographics1D.growth
	model_ex = dadi.Numerics.make_extrap_func(func)

	# 03. Define parameter space
	popt=(nu, t)
	ns = fs.sample_sizes
	pts = [max(ns)+20, max(ns)+30, max(ns)+40]

	# 04. Calculate the best fit model
	model_fs = model_ex(popt, ns, pts)

	# 05. Compute likelihood of the data given the model allele frequency spectrum 
	ll_model = dadi.Inference.ll_multinom(model_fs, fs)

	# 06. Compute likelihood of data to itself (best possible log-likelihood)
	ll_data = dadi.Inference.ll_multinom(fs, fs)

	# 07. Calculate the best fit theta
	theta = dadi.Inference.optimal_sfs_scaling(model_fs, fs)

	# 08. Model-specific scaling of parameters 
	Nanc = theta / (4*mu*L)
	nu_scaled_dip = popt[0]*Nanc
	T_scaled_dip = popt[1]*2*Nanc
	scaled_param_names = ("Nanc_FromTheta_scaled_dip","nu_scaled_dip","T_scaled_gen")
	scaled_popt = (Nanc, nu_scaled_dip, T_scaled_dip)

elif MODEL=="three_epoch":
	nuB=float(args.NUB)
	nuF=float(args.NUF)
	tB=float(args.TB)
	tF=float(args.TF)

	# 02. Define model and extrapolate function
	func = dadi.Demographics1D.three_epoch
	model_ex = dadi.Numerics.make_extrap_func(func)

	# 03. Define parameter space
	popt=(nuB, nuF, tB, tF)
	ns = fs.sample_sizes
	pts = [max(ns)+20, max(ns)+30, max(ns)+40]

	# 04. Calculate the best fit model
	model_fs = model_ex(popt, ns, pts)

	# 05. Compute likelihood of the data given the model allele frequency spectrum 
	ll_model = dadi.Inference.ll_multinom(model_fs, fs)

	# 06. Compute likelihood of data to itself (best possible log-likelihood)
	ll_data = dadi.Inference.ll_multinom(fs, fs)

	# 07. Calculate the best fit theta
	theta = dadi.Inference.optimal_sfs_scaling(model_fs, fs)

	# 08. Model-specific scaling of parameters 
	Nanc = theta / (4*mu*L)
	Nanc = theta / (4*mu*L)
	nuB_scaled_dip = popt[0]*Nanc
	nuF_scaled_dip = popt[1]*Nanc
	TB_scaled_dip = popt[2]*2*Nanc
	TF_scaled_dip = popt[3]*2*Nanc
	scaled_param_names = ("Nanc_FromTheta_scaled_dip","nuB_scaled_dip","nuF_scaled_dip", "TB_scaled_gen", "TF_scaled_gen")
	scaled_popt = (Nanc, nuB_scaled_dip,nuF_scaled_dip, TB_scaled_dip, TF_scaled_dip)
	
elif MODEL=="three_epoch_inbreeding":
	nuB=float(args.NUB)
	nuF=float(args.NUF)
	tB=float(args.TB)
	tF=float(args.TF)
	F=float(args.F)

	# 02. Define model and extrapolate function
	func = dadi.Demographics1D.exponential_growth
	model_ex = dadi.Numerics.make_extrap_func(func)

	# 03. Define parameter space
	popt=(nuB, nuF, tB, tF, F)
	ns = fs.sample_sizes
	pts = [max(ns)+20, max(ns)+30, max(ns)+40]

	# 04. Calculate the best fit model
	model_fs = model_ex(popt, ns, pts)

	# 05. Compute likelihood of the data given the model allele frequency spectrum 
	ll_model = dadi.Inference.ll_multinom(model_fs, fs)

	# 06. Compute likelihood of data to itself (best possible log-likelihood)
	ll_data = dadi.Inference.ll_multinom(fs, fs)

	# 07. Calculate the best fit theta
	theta = dadi.Inference.optimal_sfs_scaling(model_fs, fs)

	# 08. Model-specific scaling of parameters 
	Nanc = theta / (4*mu*L)
	nuB_scaled_dip = popt[0]*Nanc
	nuF_scaled_dip = popt[1]*Nanc
	TB_scaled_dip = popt[2]*2*Nanc
	TF_scaled_dip = popt[3]*2*Nanc
	F = popt[4]
	scaled_param_names = ("Nanc_FromTheta_scaled_dip","nuB_scaled_dip","nuF_scaled_dip", "TB_scaled_gen", "TF_scaled_gen", "F")
	scaled_popt = (Nanc, nuB_scaled_dip,nuF_scaled_dip, TB_scaled_dip, TF_scaled_dip, F)

elif MODEL=="bottlegrowth":
	nuB=float(args.NUB)
	nuF=float(args.NUF)
	T=float(args.T)
	# 02. Define model and extrapolate function
	func = dadi.Demographics1D.bottlegrowth
	model_ex = dadi.Numerics.make_extrap_func(func)

	# 03. Define parameter space
	popt=(nuB, nuF, T)
	ns = fs.sample_sizes
	pts = [max(ns)+20, max(ns)+30, max(ns)+40]

	# 04. Calculate the best fit model
	model_fs = model_ex(popt, ns, pts)

	# 05. Compute likelihood of the data given the model allele frequency spectrum 
	ll_model = dadi.Inference.ll_multinom(model_fs, fs)

	# 06. Compute likelihood of data to itself (best possible log-likelihood)
	ll_data = dadi.Inference.ll_multinom(fs, fs)

	# 07. Calculate the best fit theta
	theta = dadi.Inference.optimal_sfs_scaling(model_fs, fs)

	# 08. Model-specific scaling of parameters 
	Nanc = theta / (4*mu*L)
	nuB_scaled_dip = popt[0]*Nanc
	nuF_scaled_dip = popt[1]*Nanc
	T_scaled_dip = popt[2]*2*Nanc
	scaled_param_names = ("Nanc_FromTheta_scaled_dip","nuB_scaled_dip","nuF_scaled_dip","T_scaled_gen")
	scaled_popt = (Nanc, nuB_scaled_dip, nuF_scaled_dip, T_scaled_dip)

else:
	print("Model not recognized.")

############### Write out output (same for any model) ########################
print('##### Writing out parameters #####')

if not os.path.exists(outdir):
    os.mkdir(outdir)

outputFile=open(str(outdir)+"/"+"dadi.inference."+MODEL+".runNum."+str(runNum)+".output","w")
# get all param names:
#param_names_str='\t'.join(str(x) for x in param_names)
scaled_param_names_str=str(scaled_param_names)
header=scaled_param_names_str+"\ttheta\tLL\tLL_data\tmodelFunction\tmu\tL\tmaxiter\trunNumber\trundate\tinitialParameters\tupper_bound\tlower_bound" # add additional parameters theta, log-likelihood, model name, run number and rundate
scaled_popt_str=str(scaled_popt)
# join together all the output fields, tab-separated:
output=[scaled_popt_str,theta,ll_model,ll_data,MODEL,mu,L,maxiter,runNum,todaysdate] # put all the output terms together
output='\t'.join(str(x) for x in output) # write out all the output fields
# this should result in a 2 row table that could be input into R / concatenated with other runs
outputFile.write(('{0}\n{1}\n').format(header,output))
outputFile.close()


############### Output SFS ########################
#print('##### Writing out SFS #####')

#outputSFS=str(outdir)+"/"+".dadi.inference."+str(MODEL)+".runNum."+str(runNum)+".expSFS"
#model.to_file(outputSFS)
#outputSFS.close()

############### Output plot ########################
#print('##### Making plots #####')

#import pylab
#import matplotlib.pyplot as plt
#fig=plt.figure(1)
#pylab.ion()
#outputFigure=str(str(outdir)+"/"+".dadi.inference."+str(MODEL)+".runNum."+str(runNum)+".figure.png")
#dadi.Plotting.plot_1d_comp_multinom(MODEL, fs)
#pylab.show()
#plt.savefig(outputFigure)
#pylab.clf()
