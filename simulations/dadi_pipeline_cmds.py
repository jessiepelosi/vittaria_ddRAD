!{sys.executable} -m pip install --upgrade dadi

import sys
import os
import numpy
import dadi
import dadi.TwoLocus
from datetime import datetime
import Optimize_Functions
import json

pip list

##### BASIC 1D MODELS #####

################################
##                            ##
##          Neutral           ##
##                            ##
################################

### OPTIMIZATION ###

def snm_1d(notused, ns, pts):
    """
    Standard neutral model.

    ns = (n1,)

    n1: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)

    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs
snm_1d.__param_names__ = []
snm = snm_1d

fs = dadi.Spectrum.from_file("new_variants_mapped_October2024/clust1-146.sfs")

#print some useful information about the afs or jsfs
print("\n\n============================================================================")
print("\nData for site frequency spectrum:\n")
#print("Projection: {}".format(proj))
print("Sample sizes: {}".format(fs.sample_sizes))
print("Sum of SFS: {}".format(numpy.around(fs.S(), 2)))
print("\n============================================================================\n")

func = dadi.Demographics1D.snm_1d
param_names = ("n1") 

# nu specifies the ratio of the contemporary population size to the ancient population size 
# T is the time in the past at which the size change happened in units 2*Na generations 

#create a prefix to label the output files
prefix = "neutral"
                             
# Define grid points based on sample size
# For smaller data (largest sample size < 100) [ns+20, ns+30, ns+40] is a good starting point
ns = fs.sample_sizes
pts = [max(ns)+20, max(ns)+30, max(ns)+40]

#Remember the order for mandatory arguments as below
#Optimize_Routine(fs, pts, outfile, model_name, func, rounds, param_number, fs_folded)
#From DadiPipeline 
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "neutral", func, 5, 1, fs_folded=True)

################################
##                            ##
##          Neutral           ##
##                            ##
################################

### SIMULATION ###

# Use optimized values from section above 

popt = [0.3179]
model_ex = dadi.Numerics.make_extrap_func(func)

#Calculate the best fit model
model_fs = model_ex(popt, ns, pts)

#Compute likelihood of the data given the model allele frequency spectrum 
ll_model = dadi.Inference.ll_multinom(model_fs, fs)
print("ll model:", ll_model)

#Compute likelihood of data to itself (best possible log-likelihood)
ll_data = dadi.Inference.ll_multinom(fs, fs)

#Calculate the best fit theta
theta = dadi.Inference.optimal_sfs_scaling(model_fs, fs)
print("Theta:", theta)

#Model-specific scaling of parameters 
mu = 1.172e-8
L = 1609347
Nanc = theta / (4*mu*L)
N = popt[0]*Nanc
print(N)
print(N)

import matplotlib.pyplot as plt
fig = plt.figure(219033)
fig.clear()
dadi.Plotting.plot_1d_comp_multinom(model_fs, fs)
fig.savefig("/blue/barbazuk/jessiepelosi/Vittaria_popgen/demography/new_variants_mapped/neutral-plot.png")

AIC = 2*1-2*ll_model
print("AIC of the model is:", AIC)

################################
##                            ##
##     Exponential Growth     ##
##                            ##
################################

### OPTIMIZATION ###

def growth(params, ns, pts):
    nu, T = params
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)

    nu_func = lambda t : numpy.exp(numpy.log(nu) * t/T)
    phi = Integration.one_pop(phi, xx, T, nu_func)

    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs

# At time T in the past, an equilibrium population begins growing exponentially, reaching size nu at present

fs = dadi.Spectrum.from_file("new_variants_mapped_October2024/clust1-146.sfs")

#print some useful information about the afs or jsfs
print("\n\n============================================================================")
print("\nData for site frequency spectrum:\n")
#print("Projection: {}".format(proj))
print("Sample sizes: {}".format(fs.sample_sizes))
print("Sum of SFS: {}".format(numpy.around(fs.S(), 2)))
print("\n============================================================================\n")

func = dadi.Demographics1D.growth
param_names = ("nu", "T") 

# nu specifies the ratio of the contemporary population size to the ancient population size 
# T is the time in the past at which the size change happened in units 2*Na generations 

#create a prefix to label the output files
prefix = "exp_growth"
                             
# Define grid points based on sample size
# For smaller data (largest sample size < 100) [ns+20, ns+30, ns+40] is a good starting point
ns = fs.sample_sizes
pts = [max(ns)+20, max(ns)+30, max(ns)+40]

#Remember the order for mandatory arguments as below
#Optimize_Routine(fs, pts, outfile, model_name, func, rounds, param_number, fs_folded)
#From DadiPipeline 
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "exp_growth", func, 5, 2, fs_folded=True)



################################
##                            ##
##     Exponential Growth     ##
##                            ##
################################

### SIMULATION ###

# Use optimized values from section above 

popt = [5.9993, 23.2048]
model_ex = dadi.Numerics.make_extrap_func(func)

#Calculate the best fit model
model_fs = model_ex(popt, ns, pts)

#Compute likelihood of the data given the model allele frequency spectrum 
ll_model = dadi.Inference.ll_multinom(model_fs, fs)
print(ll_model)

#Compute likelihood of data to itself (best possible log-likelihood)
ll_data = dadi.Inference.ll_multinom(fs, fs)

#Calculate the best fit theta
theta = dadi.Inference.optimal_sfs_scaling(model_fs, fs)

print("Theta:", theta)
print("nu:", popt[0])
print("T:", popt[1])

#Model-specific scaling of parameters 
mu = 1.172e-8
L = 1609347
Nanc = theta / (4*mu*L)
nu_scaled_dip = popt[0]*Nanc
T_scaled_dip = popt[1]*2*Nanc
scaled_param_names = ("Nanc_FromTheta_scaled_dip","nu_scaled_dip","T_scaled_gen")
scaled_popt = (Nanc, nu_scaled_dip, T_scaled_dip)
print(scaled_param_names)
print(scaled_popt)

import matplotlib.pyplot as plt
fig = plt.figure(219033)
fig.clear()
dadi.Plotting.plot_1d_comp_multinom(model_fs, fs)
fig.savefig("/blue/barbazuk/jessiepelosi/Vittaria_popgen/demography/new_variants_mapped/exp-growth.png")

AIC = 2*2-2*ll_model
print("AIC of the model is:", AIC)

################################
##                            ##
##         Two Epoch          ##
##                            ##
################################

### OPTIMIZATION ###

def two_epoch(params, ns, pts):
    """
    Instantaneous size change some time ago.

    params = (nu,T)
    ns = (n1,)

    nu: Ratio of contemporary to ancient population size
    T: Time in the past at which size change happened (in units of 2*Na 
       generations) 
    n1: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu,T = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    
    phi = Integration.one_pop(phi, xx, T, nu)

    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs
two_epoch.__param_names__ = ['nu', 'T']


fs = dadi.Spectrum.from_file("new_variants_mapped_October2024/clust1-146.sfs")

#print some useful information about the afs or jsfs
print("\n\n============================================================================")
print("\nData for site frequency spectrum:\n")
#print("Projection: {}".format(proj))
print("Sample sizes: {}".format(fs.sample_sizes))
print("Sum of SFS: {}".format(numpy.around(fs.S(), 2)))
print("\n============================================================================\n")

func = dadi.Demographics1D.two_epoch
param_names = ("nu", "T") 

# nu specifies the ratio of the contemporary population size to the ancient population size 
# T is the time in the past at which the size change happened in units 2*Na generations 

#create a prefix to label the output files
prefix = "two_epoch"
                             
# Define grid points based on sample size
# For smaller data (largest sample size < 100) [ns+20, ns+30, ns+40] is a good starting point
ns = fs.sample_sizes
pts = [max(ns)+20, max(ns)+30, max(ns)+40]

#Remember the order for mandatory arguments as below
#Optimize_Routine(fs, pts, outfile, model_name, func, rounds, param_number, fs_folded)
#From DadiPipeline 
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "two_epoch", func, 5, 2, fs_folded=True)

################################
##                            ##
##         Two Epoch          ##
##                            ##
################################

### SIMULATION ###

# Use optimized values from section above 

popt = [29.7, 29.1742]
model_ex = dadi.Numerics.make_extrap_func(func)

#Calculate the best fit model
model_fs = model_ex(popt, ns, pts)

#Compute likelihood of the data given the model allele frequency spectrum 
ll_model = dadi.Inference.ll_multinom(model_fs, fs)
print("ll model:", ll_model)

#Compute likelihood of data to itself (best possible log-likelihood)
ll_data = dadi.Inference.ll_multinom(fs, fs)

#Calculate the best fit theta
theta = dadi.Inference.optimal_sfs_scaling(model_fs, fs)

print("nu", popt[0])
print("T", popt[1])

#Model-specific scaling of parameters 
mu = 1.172e-8
L = 1609347
Nanc = theta / (4*mu*L)
nu_scaled_dip = popt[0]*Nanc
T_scaled_dip = popt[1]*Nanc
scaled_param_names = ("Nanc_FromTheta_scaled_dip","nu_scaled_dip", "T_scaled")
scaled_popt = (Nanc, nu_scaled_dip, T_scaled_dip)
print(scaled_param_names)
print(scaled_popt)

import matplotlib.pyplot as plt
fig = plt.figure(219033)
fig.clear()
dadi.Plotting.plot_1d_comp_multinom(model_fs, fs)
fig.savefig("/blue/barbazuk/jessiepelosi/Vittaria_popgen/demography/two_epoch.png")

AIC = 2*2-2*ll_model
print("AIC of the model is:", AIC)

################################
##                            ##
##        Three Epoch         ##
##                            ##
################################

### OPTIMIZATION ###

def three_epoch(params, ns, pts):
    """
    params = (nuB,nuF,TB,TF)
    ns = (n1,)

    nuB: Ratio of bottleneck population size to ancient pop size
    nuF: Ratio of contemporary to ancient pop size
    TB: Length of bottleneck (in units of 2*Na generations) 
    TF: Time since bottleneck recovery (in units of 2*Na generations) 

    n1: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nuB,nuF,TB,TF = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)

    phi = Integration.one_pop(phi, xx, TB, nuB)
    phi = Integration.one_pop(phi, xx, TF, nuF)

    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs

three_epoch.__param_names__ = ['nuB', 'nuF', 'TB', 'TF']

# At time T in the past, an equilibrium population begins growing exponentially, reaching size nu at present

fs = dadi.Spectrum.from_file("new_variants_mapped_October2024/clust1-146.sfs")

#print some useful information about the afs or jsfs
print("\n\n============================================================================")
print("\nData for site frequency spectrum:\n")
#print("Projection: {}".format(proj))
print("Sample sizes: {}".format(fs.sample_sizes))
print("Sum of SFS: {}".format(numpy.around(fs.S(), 2)))
print("\n============================================================================\n")

func = dadi.Demographics1D.three_epoch
param_names = ("nuB", "nuF", "TB", "TF") 

# nu specifies the ratio of the contemporary population size to the ancient population size 
# T is the time in the past at which the size change happened in units 2*Na generations 

#create a prefix to label the output files
prefix = "three_epoch"
                             
# Define grid points based on sample size
# For smaller data (largest sample size < 100) [ns+20, ns+30, ns+40] is a good starting point
ns = fs.sample_sizes
pts = [max(ns)+20, max(ns)+30, max(ns)+40]

#Remember the order for mandatory arguments as below
#Optimize_Routine(fs, pts, outfile, model_name, func, rounds, param_number, fs_folded)
#From DadiPipeline 
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "three_epoch", func, 5, 4, fs_folded=True)

################################
##                            ##
##        Three Epoch         ##
##                            ##
################################

### SIMULATION ###

# Use optimized values from section above 

popt = [0.3764, 4.7993, 9.0071, 4.5044]
model_ex = dadi.Numerics.make_extrap_func(func)

#Calculate the best fit model
model_fs = model_ex(popt, ns, pts)

#Compute likelihood of the data given the model allele frequency spectrum 
ll_model = dadi.Inference.ll_multinom(model_fs, fs)
print("ll model:", ll_model)

#Compute likelihood of data to itself (best possible log-likelihood)
ll_data = dadi.Inference.ll_multinom(fs, fs)

#Calculate the best fit theta
theta = dadi.Inference.optimal_sfs_scaling(model_fs, fs)

print("nuB", popt[0])
print("nuF", popt[1])
print("TB", popt[2])
print("TF", popt[3])

#Model-specific scaling of parameters 
mu = 1.172e-8
L = 1609347
Nanc = theta / (4*mu*L)
nuB_scaled_dip = popt[0]*Nanc
nuF_scaled_dip = popt[1]*Nanc
TB_scaled_dip = popt[2]*2*Nanc
TF_scaled_dip = popt[3]*2*Nanc
scaled_param_names = ("Nanc_FromTheta_scaled_dip","nuB_scaled_dip","nuF_scaled_dip", "TB_scaled_gen", "TF_scaled_gen")
scaled_popt = (Nanc, nuB_scaled_dip,nuF_scaled_dip, TB_scaled_dip, TF_scaled_dip)
print(scaled_param_names)
print(scaled_popt)

import matplotlib.pyplot as plt
fig = plt.figure(219033)
fig.clear()
dadi.Plotting.plot_1d_comp_multinom(model_fs, fs)
fig.savefig("/blue/barbazuk/jessiepelosi/Vittaria_popgen/demography/dadi/demo_plot.png")

AIC = 2*4-2*ll_model
print("AIC of the model is:", AIC)

################################
##                            ##
##        Three Epoch         ##
##       With Inbreeding      ##
##                            ##  
################################

### OPTIMIZATION ###

def three_epoch_inbreeding(params, ns, pts):
    """
    params = (nuB,nuF,TB,TF,F)
    ns = (n1,)

    nuB: Ratio of bottleneck population size to ancient pop size
    nuF: Ratio of contemporary to ancient pop size
    TB: Length of bottleneck (in units of 2*Na generations) 
    TF: Time since bottleneck recovery (in units of 2*Na generations) 
    F: Inbreeding coefficent

    n1: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nuB,nuF,TB,TF,F = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)

    phi = Integration.one_pop(phi, xx, TB, nuB)
    phi = Integration.one_pop(phi, xx, TF, nuF)

    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx,), (F,), (2,))
    return fs

three_epoch_inbreeding.__param_names__ = ['nuB', 'nuF', 'TB', 'TF', 'F']

fs = dadi.Spectrum.from_file("new_variants_mapped_October2024/clust1-146.sfs")

#print some useful information about the afs or jsfs
print("\n\n============================================================================")
print("\nData for site frequency spectrum:\n")
#print("Projection: {}".format(proj))
print("Sample sizes: {}".format(fs.sample_sizes))
print("Sum of SFS: {}".format(numpy.around(fs.S(), 2)))
print("\n============================================================================\n")

func = dadi.Demographics1D.three_epoch_inbreeding
param_names = ("nuB", "nuF", "TB", "TF", "F") 

#create a prefix to label the output files
prefix = "three_epoch_inbreeding"
                             
# Define grid points based on sample size
# For smaller data (largest sample size < 100) [ns+20, ns+30, ns+40] is a good starting point
ns = fs.sample_sizes
pts = [max(ns)+20, max(ns)+30, max(ns)+40]

#Remember the order for mandatory arguments as below
#Optimize_Routine(fs, pts, outfile, model_name, func, rounds, param_number, fs_folded)
#From DadiPipeline 
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "three_epoch_inbreeding", func, 5, 5, fs_folded=True)


################################
##                            ##
##        Three Epoch         ##
##       With Inbreeding      ##    
##                            ##
################################

### SIMULATION ###

# Use optimized values from section above 

popt = [13.4954, 7.5317, 0.0939, 4.4618, 0.7428]
model_ex = dadi.Numerics.make_extrap_func(func)

#Calculate the best fit model
model_fs = model_ex(popt, ns, pts)

#Compute likelihood of the data given the model allele frequency spectrum 
ll_model = dadi.Inference.ll_multinom(model_fs, fs)
print("ll model:", ll_model)

#Compute likelihood of data to itself (best possible log-likelihood)
ll_data = dadi.Inference.ll_multinom(fs, fs)

#Calculate the best fit theta
theta = dadi.Inference.optimal_sfs_scaling(model_fs, fs)

print("nuB", popt[0])
print("nuF", popt[1])
print("TB", popt[2])
print("TF", popt[3])
print("F", popt[4])

#Model-specific scaling of parameters 
mu = 1.172e-8
L = 1609347
Nanc = theta / (4*mu*L)
nuB_scaled_dip = popt[0]*Nanc
nuF_scaled_dip = popt[1]*Nanc
TB_scaled_dip = popt[2]*2*Nanc
TF_scaled_dip = popt[3]*2*Nanc
F = popt[4]
scaled_param_names = ("Nanc_FromTheta_scaled_dip","nuB_scaled_dip","nuF_scaled_dip", "TB_scaled_gen", "TF_scaled_gen", "F")
scaled_popt = (Nanc, nuB_scaled_dip,nuF_scaled_dip, TB_scaled_dip, TF_scaled_dip, F)
print(scaled_param_names)
print(scaled_popt)

import matplotlib.pyplot as plt
fig = plt.figure(219033)
fig.clear()
dadi.Plotting.plot_1d_comp_multinom(model_fs, fs)
fig.savefig("/blue/barbazuk/jessiepelosi/Vittaria_popgen/demography/three-epoch-inbreeding.png")

AIC = 2*5-2*ll_model
print("AIC of the model is:", AIC)

################################
##                            ##
##        Bottlegrowth        ##
##                            ##
################################

### OPTIMIZATION ###

def bottlegrowth_1d(params, ns, pts):
    """
    Instantanous size change followed by exponential growth.

    params = (nuB,nuF,T)
    ns = (n1,)

    nuB: Ratio of population size after instantanous change to ancient
         population size
    nuF: Ratio of contemporary to ancient population size
    T: Time in the past at which instantaneous change happened and growth began
       (in units of 2*Na generations) 
    n1: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nuB,nuF,T = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)

    nu_func = lambda t: nuB*numpy.exp(numpy.log(nuF/nuB) * t/T)
    phi = Integration.one_pop(phi, xx, T, nu_func)

    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs
bottlegrowth_1d.__param_names__ = ['nuB', 'nuF', 'T']

bottlegrowth = bottlegrowth_1d

fs = dadi.Spectrum.from_file("new_variants_mapped_October2024/clust1-146.sfs")

#print some useful information about the afs or jsfs
print("\n\n============================================================================")
print("\nData for site frequency spectrum:\n")
#print("Projection: {}".format(proj))
print("Sample sizes: {}".format(fs.sample_sizes))
print("Sum of SFS: {}".format(numpy.around(fs.S(), 2)))
print("\n============================================================================\n")

func = dadi.Demographics1D.bottlegrowth_1d
param_names = ("nuB", "nuF", "T") 

#create a prefix to label the output files
prefix = "bottlegrowth"
                             
# Define grid points based on sample size
# For smaller data (largest sample size < 100) [ns+20, ns+30, ns+40] is a good starting point
ns = fs.sample_sizes
pts = [max(ns)+20, max(ns)+30, max(ns)+40]

#Remember the order for mandatory arguments as below
#Optimize_Routine(fs, pts, outfile, model_name, func, rounds, param_number, fs_folded)
#From DadiPipeline 
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "bottlegrowth", func, 5, 3, fs_folded=True)

################################
##                            ##
##        Bottlegrowth        ##
##                            ##
################################

### SIMULATION ###

# Use optimized values from section above 

popt = [28.2469, 1.0772, 0.5297]
model_ex = dadi.Numerics.make_extrap_func(func)

#Calculate the best fit model
model_fs = model_ex(popt, ns, pts)

#Compute likelihood of the data given the model allele frequency spectrum 
ll_model = dadi.Inference.ll_multinom(model_fs, fs)
print("ll model:", ll_model)

#Compute likelihood of data to itself (best possible log-likelihood)
ll_data = dadi.Inference.ll_multinom(fs, fs)

#Calculate the best fit theta
theta = dadi.Inference.optimal_sfs_scaling(model_fs, fs)

print("nuB", popt[0])
print("nuF", popt[1])
print("T", popt[2])

#Model-specific scaling of parameters 
mu = 4.97e-9
L = 40398599
Nanc = theta / (4*mu*L)
nuB_scaled_dip = popt[0]*Nanc
nuF_scaled_dip = popt[1]*Nanc
T_scaled_dip = popt[2]*2*Nanc
scaled_param_names = ("Nanc_FromTheta_scaled_dip","nuB_scaled_dip","nuF_scaled_dip", "T_scaled_gen")
scaled_popt = (Nanc, nuB_scaled_dip,nuF_scaled_dip, T_scaled_dip)
print(scaled_param_names)
print(scaled_popt)

import matplotlib.pyplot as plt
fig = plt.figure(219033)
fig.clear()
dadi.Plotting.plot_1d_comp_multinom(model_fs, fs)
fig.savefig("/blue/barbazuk/jessiepelosi/Vittaria_popgen/demography/new_variants_mapped/bottlegrowth.png")

AIC = 2*3-2*ll_model
print("AIC of the model is:", AIC)