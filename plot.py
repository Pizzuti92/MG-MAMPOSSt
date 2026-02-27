#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 20 15:50:12 2025

@author: lorenzopizzuti
"""

import numpy as np
from math import *
import scipy.optimize as opt
import matplotlib.pyplot as plt

import scipy.interpolate as ipt
from getdist import plots, MCSamples
from scipy.ndimage import gaussian_filter
import models as md
import argparse 
import re #For regular expressions

from collections.abc import Iterable #For checking arrays

from sys import exit
import os


VERSION="2.0"
AUTHORS = "2025 Lorenzo Pizzuti"
description = "This script reads the output from MG-MAMPOSSt,\
            it defines the confidence intervals and plot the results."


# =============================================================================
# Swithces: default values

#if getdata is True, read gas and galaxies for the final plot from external files
getdata = False 

#Exclude the BCG inner region and profiles in the plot? Y/N = 1/0
kbdt = 0

#folder where store the mass profiles
dnamefold = 'MassProfiles'


#Input mg-mamposst file
titfile = 'gomamposst_x.inp' #'gomamposst_x.inp' #'mamComa.inp' 

#Default options file
optFile = 'Options.txt'

plot = True       #plotting

numfac = 10000 # thinning factor for avoid overcomputation of the mass profiles
                         #  for very long chains

median_string = True #Plot the median profiles by default

facext = 1.2      #percentage of rupin up to which plot 

LatexForm = True  #print a Latex formatted output for the constraints

# Compute first the virial mass M200 and plot the marginalized distribution with 
# uncertainties
plotm200 = True  


colore = 'darkgreen'    #color of the contour plot

# Plot the profile of the velocity anisotropy
betaplot = True


#For the Chameleon gravity setup, you can be in f(R) subclass. In that case, we will
# plot directly the scalaron field f_R instead of the chameleon field phi
frcase = 0 

# Default name for gas and galaxies files
namegal = "data/massgal.dat"
namegas = "data/massgas.dat"


#Default name for the gas-galaxies chain file
nameGasGal = "GasGal_mass.txt"

#Plot fonts (Customize here, no input arguments)
fontlabel=22
font=24
#===================== Physical constants =====================================
G = 4.302e-9


#******************************************************************************
parser = argparse.ArgumentParser(description=description)

# main arguments
parser.add_argument("-i", "--Input", help = "[string] Input posterior file")
parser.add_argument("-pf", "--paramfile", help = "[string] Input parameter file") 
parser.add_argument("-op", "--options", help = "[string] Input Options file")

 # Adding optional arguments
parser.add_argument("-gb","--getbaryon", help = "[bool] get baryonic data from additional file (True/False)")
parser.add_argument("-mf","--massfolder", help = "[string] name of the folder where storing the mass profiles. Default is MassProfiles")
parser.add_argument("-pt","--plot", help = "[bool] plot the mass/anisotropy profiles? True/False. Default is True")
parser.add_argument("-m2","--plotm200", help = "[bool] plot the distribution of M200? True/False. Default is True")
parser.add_argument("-c","--color", help = "[string] color of the marginalized distribution. Default is darkgreen")

# Read arguments from command line
args = parser.parse_args()

#Beginning of the script
if (args.Input == None):
    
    name_like = np.loadtxt(titfile,dtype=str)[6]
else:
    name_like = args.Input
    
    
if (args.options != None):
    optFile = args.options
    
if (args.paramfile == None):
    read_input = np.loadtxt(titfile,dtype=str)[1]
else:
    read_input = args.paramfile
    
if (args.getbaryon != None and args.getbaryon.lower() == "True"):
    getdata = True

if (args.massfolder != None):
    dnamefold = args.massfolder

workdir = os.getcwd()

namefold = os.path.join(workdir,dnamefold)
#if the directory doesn't exist, create it
if not os.path.exists(namefold):
    os.mkdir(namefold) 
    
if (args.plot != None and args.plot.lower() == "false"):
    plot = False

if (args.plotm200 != None and args.plotm200.lower() == "false"):
    plotm200 = False
    
if (args.color != None):
    colore = args.color

# =============================================================================

# Main 

print("\nThis is Plot.py version %s\n\t \
      Copyright (C) %s\n" % (VERSION, AUTHORS))
print(" ")
print("marginalization uses the getdist python package")
print(" ")
#Load the posterior and parameters chains

data = np.loadtxt(open(name_like).readlines()[:-4], 
                skiprows=5,dtype=float)


#Create array containing the free parameters

r200=data[:,0].copy()  #the virial radius
rnu=data[:,1].copy()   #scale radius of the number density profile
rs=data[:,2].copy()    #scale radius of the mass profile
beta=data[:,3].copy()  #beta infinity for anisotropy (or r_beta for ML profile) 
mod1=data[:,4].copy()  #first parameter in MG models 
mod2=data[:,5].copy()  #second parameter in MG models
mod3=data[:,6].copy()  #exponent in RG (or gamma of gNFW)
beta2=data[:,7].copy() #second (central) anisotropy parameter
rbeta=data[:,8].copy() #free radius of the anisotropy profile

##NOTE: In the case of general chameleon the RESCALED VARIABLES Q_2 and \phi_2
##are plotted (see Pizzuti et al., 2021)


# Now reading the nmcmc option ================================================


# _nmcmc_re = re.compile(r'^\s*nmcmc\s*(?:=\s*|\s+)([+-]?\d+)\b')

# def read_nmcmc(path=optFile, default=0):
#     with open(path, "r", encoding="utf-8") as f:
#         for raw in f:
#             # togli commenti e newline
#             line = raw.split("#", 1)[0].rstrip()
#             m = _nmcmc_re.match(line)
#             if m:
#                 try:
#                     return int(m.group(1))
#                 except ValueError:
#                     print("parameter nmcmc has undefined value. Switching to default")
#                     return default
#     return default


# nmcmc = read_nmcmc("Options.txt", default=0)
# if nmcmc < 1:
#     print(f"The option for nmcmc is {nmcmc}. No plot for the grid case")
#     raise SystemExit

#==============================================================================
# Reading the parameters from the file (if you want to add parameters, you need 
# to change the read_params block in models)
pars = md.read_params_block(read_input, mod1=mod1, mod2=mod2, r200=r200)

nmcmc = pars["nmcmc"]
if nmcmc < 1:
    print(f"The option for nmcmc is {nmcmc}. No plot for the grid case")
    raise SystemExit

za = pars["za"]
H0 = pars["H0"]
Ol = pars["Olam"]
Om = pars["Omegam"]

kbcg = pars["BCG"]
kbary = pars["Baryons"]
xml = pars["xlumbcg"]

rhogas = pars["rho0gas"]
rsgas  = pars["rsgas"]

rhostar = pars["rhostar"]
rstar   = pars["rstar"]

srhost = pars["srhost"]
srstar = pars["srstar"]

srhoga = pars["srhoga"]
srsga  = pars["srsga"]

eim = pars["eim"]

mass_string = pars["M(r)"]
anis_string = pars["Beta(r)"]
freebeta    = pars["FreeBeta"]

nhone  = pars["nA2"]
mod1low = pars["A1low"]
mod1up  = pars["A1up"]
mod2low = pars["A2low"]
mod2up  = pars["A2up"]

rupin = pars["Rup"]

#If available
if (kbary == 1):
    rhogg, rsgg, rhostarg, rstarg = np.loadtxt(nameGasGal, unpack = True,
                                               skiprows = 1)
    #Consistency check for the lenghts
    lenArray = min(len(rhogg),len(r200))
    rhogas, rsgas, rhostar, rstar = rhogg[:lenArray]*1e13, rsgg[:lenArray], \
    rhostarg[:lenArray]*1e13, rstarg[:lenArray]
    
    
    r200=r200[:lenArray]  #the virial radius
    rnu=rnu[:lenArray]  #scale radius of the number density profile
    rs=rs[:lenArray]   #scale radius of the mass profile
    beta=beta[:lenArray] #beta infinity for anisotropy (or r_beta for ML profile) 
    mod1=mod1[:lenArray] #first parameter in MG models 
    mod2=mod2[:lenArray]  #second parameter in MG models
    mod3=mod3[:lenArray] #exponent in RG (or gamma of gNFW)
    beta2=beta2[:lenArray] #second (central) anisotropy parameter
    rbeta=rbeta[:lenArray]
    
    data = data[:lenArray,:]

# =============================================================================
#Check if you have f(R) as a specific subcase of chameleon:
if (nhone<=0 and isclose(mod2[0], 1/np.sqrt(6), abs_tol=1e-2)):
    if (mass_string=="mNFW_GC" or mass_string=='mBur_GC' 
        or mass_string=='mEin_GC' or mass_string=='mIso_GC' 
        or mass_string =="mTotal_GC" ): 
        
        frcase=1

#Current version: NO BCG in some MG scenarios 
if (mass_string=='mNFW_GC' or  mass_string=='mBur_GC'  or \
    mass_string=='mEin_GC'  or mass_string=='mIso_GC' or  \
    mass_string=='mgNFW_GC' or mass_string=='mNFW_LH' or mass_string=='mNFW_BH'):
    
    kbcg=0 #Force the bcg mass to be excluded
    
#MG Chameleon multicomponent. In this case the BCG MUST be included
if(mass_string == 'mTotal_GC'):
    kbcg = 1

# Two cases if BCG is included or not 
if (kbcg==0 and mass_string != "RG"):
    #Scaled likelihood (normalized s.t. the peak is at 1)
    likemin = np.amin(data[:,9])
    
    like = data[:,9] - np.amin(data[:,9]) 
    

if (kbcg==1 or mass_string=='RG'):
#In this case you have also the additional parameters for the BCG
    #Scaled likelihood (normalized s.t. the peak is at 1)
    likemin = np.amin(data[:,11])
    like = data[:,11] - np.amin(data[:,11]) 
    rjaf = data[:,9].copy()  #rjaf radius
    xmstar = data[:,10].copy() #mass-to-light radius of the BCG
    
    #Store the optimized values
    rjafmin = rjaf[np.where(like==0)]
    xmstarmin = xmstar[np.where(like==0)]


#Store the optimized values
r2min = r200[np.where(like==0)]
rsmin = rs[np.where(like==0)]
betamin = beta[np.where(like==0)]
beta2min = beta2[np.where(like==0)]
mod1min = mod1[np.where(like==0)]
mod2min = mod2[np.where(like==0)]
mod3min = mod3[np.where(like==0)]

#==============================================================================
#find which are the free parameters
total=0

if mass_string == 'gNFW_CDE':
    data[:,5] = data[:,5] - 3.5
    
    mod2low = mod2low - 3.5
    mod2up  = mod2up - 3.5

if(kbcg==0 and mass_string!='RG'):
    #in this case only 9 free parameters
    freep=np.zeros(9)
    for i in range(0,9):
        if(np.any(data[:,i]!=data[0,i])):
            freep[i]=1
else:
    #up to 11 free parameters
    freep=np.zeros(11)
    for i in range(0,11):
        if(np.any(data[:,i]!=data[0,i])):
            freep[i]=1
        
            
il,=np.where(freep==1)
total = len(il)

#==============================================================================

# M200 plot. Do it only if it exists at least one best-fit 
if np.any(il == 0) and plotm200 and mass_string != "RG":

    # --- Costanti e M200 ---
   
    Hz2 = H0**2 * (Om * (1.0 + za)**3 + Ol)
    M200 = (100.0 * Hz2 / G) * r200**3

    # --- best-fit (assumo like==0 identifica il best) ---
    best_idx = np.flatnonzero(like == 0)
    if best_idx.size == 0:
        # sicurezza: nel tuo if sopra controlli il==0, ma qui usi like
        print("No best-fit (like==0) found: skip M200 plot.")
        raise SystemExit

    m_best = float(M200[best_idx[0]])

    # --- GetDist / MCSamples: array (N,1) ---
    names  = ["M200"]
    labels = [r"M_{200}"]
    samples = MCSamples(samples=M200[:, None], names=names, labels=labels)

    # calcola stats UNA volta sola
    stats  = samples.getMargeStats().parWithName("M200")
    lim68  = stats.limits[0]  # 68%
    lim95  = stats.limits[1]  # 95%
    mlow68, mup68 = lim68.lower, lim68.upper
    mlow95, mup95 = lim95.lower, lim95.upper

    # --- stampa intervalli attorno al best-fit (non la mean) ---
    print ("=================================================================")
    print("M200 \nconfidence intervals (around best-fit like==0)")
    print("{0:20s}\t best\t -95%\t -68%\t +68%\t +95%".format("param"))
    print("{0:20s}\t {1:6.3e}\t -{2:6.3e}\t -{3:6.3e}\t {4:6.3e}\t {5:6.3e}".format(
        labels[0],
        m_best,
        m_best - mlow95, m_best - mlow68,
        mup68 - m_best,  mup95 - m_best
    ))
    latex_68 = md.format_uncertainty(m_best/1e14, (mup68-m_best)/1e14, 
                                     (m_best-mlow68)/1e14)
    
    print("written in latex form (units of 1e14 Msun)")
    print(latex_68)
    print ("=================================================================")
    print (" ")
    # --- histogram + smoothing ---
    nbins = 80
    counts, edges = np.histogram(M200, bins=nbins)  # range auto: min/max
    centers = 0.5 * (edges[:-1] + edges[1:])

    hsmooth = gaussian_filter(counts.astype(float), sigma=1.0)
    if hsmooth.max() > 0:
        hsmooth /= hsmooth.max()

    # --- plot ---
    plt.figure(figsize = (9,9))
    plt.plot(centers / 1e14, hsmooth, lw=3, c="black")

    # masks for the intervals filling
    m68 = (centers >= mlow68) & (centers <= mup68)
    m95 = (centers >= mlow95) & (centers <= mup95)

    plt.fill_between((centers[m95] / 1e14), hsmooth[m95], color="black", alpha=0.3)
    plt.fill_between((centers[m68] / 1e14), hsmooth[m68], color="black", alpha=0.3)

    plt.ylabel(r"$P\,[M_{200}]$", fontsize=font)
    plt.xlabel(r"$M_{200}\,[10^{14}\,M_\odot]$", fontsize=font)
    plt.axvline(m_best / 1e14, ls ='--', c='black')
    plt.ylim(0,1.1)
    plt.tick_params(axis="x", labelsize=fontlabel)
    plt.tick_params(axis="y", labelsize=fontlabel)

    plt.show()

#==============================================================================
#Define the plot labels =======================================================
    
labelr2 = r"r_{200}\,\,[Mpc]"
Plabelr2 = r"$P(r_{200})$"

#Current implementation: rho_gas is r200 in Refracted gravity
if mass_string == 'RG':
      labelr2 = r"\rho_g"  

labelm3 = r"\mathcal{Q}"
Plabelm3 = r"$P(\mathcal{Q})$"

#In gNFW or mTotal_GC, the third MG parameter is the gamma of the profile
if (mass_string=='mTotal_GC' or mass_string == 'gNFW_CDE' or mass_string == 'gNFW'):
    labelm3 = r"\gamma"
    Plabelm3 = r"$P(\gamma)$"


labelb = r"r_{\beta}\,\, [Mpc]"
Plabelb = r"P(r_{\beta})"
    
#Collect the anisotropy model from data file
nani = data[0,9]
if (kbcg == 1):
    nani = data[0,11]

if(nani!=1): 
    labelb="\mathcal{A}_\infty"
    Plabelb="P(\mathcal{A}_\infty)"
 

labelrj="r_{jaffe} \,\,[Mpc]"
Plabelrj=  "P(r_{jaffe})"

labelms= "X_L"#"M_{*,BCG}/L"
Plabelms= "P(X_L)" #"P(M_{*,BCG}/L)"
  
labelb2="\mathcal{A}_0"
Plabelb2="P(\mathcal{A}_0)"


labelRb = r"r_{\beta}\,\, [Mpc]"
PlabelRb = r"P(r_{\beta})"

# ===================== modified gravity. Default: Chameleon ==================
labelm1="\phi_2"
Plabelm1=r"$P(\phi_2)$"


labelm2="\mathcal{Q}_2"
Plabelm2=r"$P(\mathcal{Q}_2)$"

#In the sub-case of f(R), plot directly the logarithm of the field
if (frcase==1):
    labelm1="Log_{10}|f_R|"
    Plabelm1=r"$P[Log_{10}|f_R|]$"


if (mass_string == 'gNFW_CDE'):
    labelm1="c_s^2"
    Plabelm1=r"$P[c_s^2]$"
    labelm2="log10 f(w,c_s^2) "
    Plabelm2=r"$P(log10 f)$"
    
# DHOST GRAVITY
if (mass_string=='mNFW_BH'):
    labelm1="Y_1"
    labelm2="Y_2"
    Plabelm1="r$P(Y_1)$"
    Plabelm2="r$P(Y_2)$"

# ===================== NO MORE USED ==========================================
#Please, note that the linear horndeski case is deprecated    
elif (mass_string=='mNFW_LH'):
    labelm1="Log(m/Mpc^{-1})"
    labelm2=r"Q"
    Plabelm1=r"P[Log$(m/Mpc^{-1})$]"
    Plabelm2=r"$P(Q)$"
    
    if (np.any(il==5)):
        data[:,5]=np.log10(data[:,5]) #uses log values for Q in the linear
        #Horndeski case to avoid problems in the getdist marginalizations
        labelm2=r"Log(Q)"
        Plabelm2=r"P[Log$(Q)$]"
# ============================================================================

# REFRACTED GRAVITY
elif (mass_string=='RG'):
	labelm1="\epsilon_0"
	labelm2=r"\log_{10}\rho_c"
	Plabelm2 = r"P[log$_{10}(\rho_c)$]"
	Plabelm1 = r"$P(\epsilon_0)$"

# BOSON STAR HALO 
elif (mass_string=='BS'):
    labelm1=r"\alpha"
    labelm2=r"\delta"
    Plabelm1=r"P[Log$(\alpha)$]"
    Plabelm2=r"$P(\delta)$"
    labelm3 = r"\kappa"
    Plabelm3 = r"$P(\kappa)$"

# =============================================================================
       

# MG columns (0-based)
COL_PHI = 4
COL_Q   = 5
COL_XMA = 6

GC_MODELS = {"mNFW_GC", "mBur_GC", "mEin_GC", "mgNFW_GC", "mIso_GC", "mTotal_GC"}

is_gc_model = mass_string in GC_MODELS

#The MG variables need to be transformed in Chameleon
needs_transform = np.any((il == 4) | (il == 5))

# --- Backup: save original data before any transformation ---
if is_gc_model:
    Qcoup  = data[:, COL_Q].copy()
    phicoup = data[:, COL_PHI].copy()

# --- Now transform (il==4 o il==5) ---
if needs_transform:

    if mass_string == "mNFW_LH":
        # log10 su phi (col 4) + limiti in log10
        data[:, COL_PHI] = np.log10(data[:, COL_PHI])
        mod1up  = np.log10(mod1up)
        mod1low = np.log10(mod1low)

    # f(R) in GC model: log10 of "field" 
    if is_gc_model and frcase == 1:
        a = np.sqrt(2.0 / 3.0) / 1e5

        # tranform: f_R = log10( exp(a*phi) - 1 )
       
        data[:, COL_PHI] = np.log10(np.expm1(a * data[:, COL_PHI]))
        
        #Change the bounds
        mod1up  = np.log10(np.expm1(a * mod1up))
        mod1low = np.log10(np.expm1(a * mod1low))

    # General chameleon (non f(R)): rescaling in [0,1]
    if is_gc_model and frcase == 0:
        data[:, COL_PHI] = 1.0 - np.exp(-0.1 * data[:, COL_PHI])
        data[:, COL_Q]   = data[:, COL_Q] / (1.0 + data[:, COL_Q])

        # limiti “hard-coded”
        mod1up  = 0.99
        mod1low = 0.01
        mod2low = 0.0
        mod2up  = 1.1

# --- For the other models ---
if mass_string == "RG":
    rhocRG   = data[:, COL_Q].copy()
    epsilon0 = data[:, COL_PHI].copy()

elif mass_string == "BS":
    bexp = data[:, COL_PHI].copy()
    gexp = data[:, COL_Q].copy()
    xma  = data[:, COL_XMA].copy()


# ====================== CREATES THE LABELS NOW ===============================
#Array of label for the triangle plots
label = np.array([labelr2,r"r_{\nu}\,\, [Mpc]",r"r_{s}\,\, [Mpc]",labelb,labelm1,
                labelm2,labelm3,labelb2,labelRb,labelrj,labelms])
#Array of labels for the single plots
label1d = np.array([labelr2,r"$r_{\nu}$ [Mpc]",r"$r_{s}$ [Mpc]",r"$"+labelb+"$",
                  r"$"+labelm1+"$",r"$"+labelm2+"$",r"$"+labelm3+"$",
                  r"$"+labelb2+"$",r"$"+labelRb+"$", r"$"+labelrj+"$",r"$"+labelms+"$"])

#Array of y-labels for the single 1d distributions plot
Plabel = np.array([Plabelr2,r"$P(r_{\nu})$",r"$P(r_{s})$",Plabelb,
                 Plabelm1,Plabelm2,Plabelm3,Plabelb2,PlabelRb, Plabelrj,Plabelms])

#==============================================================================
print (" ")
print("Running with  the mass model: \t" + mass_string)
print(" ")

# Prepare the samples 

if total == 1:
# one single free parameter
    labels=label1d[il]
    Ylabels=Plabel[il]
    ndim=1
    names = np.array(["x%s"%i for i in range(ndim)])
    first=data[:,il[0]]
# Marginalize directly and gives confidence intervals
    plotf=np.column_stack((first,np.exp(-like)))
    plotf= plotf[np.argsort(plotf[:, 0])]   
    g = plots.get_single_plotter(width_inch=10, ratio=1)
    plt.grid()
    plt.ylabel(Ylabels[0], fontsize=22)
    plt.xlabel(labels[0], fontsize=22)
    plt.tick_params(axis="x", labelsize=fontlabel)
    plt.tick_params(axis="y", labelsize=fontlabel)
    plt.plot(plotf[:,0],plotf[:,1], lw=3,c='black')
    plt.show()
    
    samples = MCSamples(samples=first,names = names, labels = labels)
    samps=first
   
cols = np.asarray(il, dtype=int)

if total < 1:
    raise ValueError(f"total must be >= 1, got {total}")

if total > cols.size:
    raise ValueError(f"total={total} but len(il)={cols.size}")

ndim = int(total)
samps = data[:, cols[:ndim]]

# =============================================================================
# Get the getdist MCSamples objects for the samples, specifying same parameter
# names and labels; if not specified weights are assumed to all be unity
names = np.array(["x%s"%i for i in range(ndim)])
labels =  label[il]
if(np.any(il==4) and mass_string!='gNFW'):
    ai,=np.where(il==4)
    names[ai]='m1'
    if(np.any(il==5)): 
        ai2,=np.where(il==5)
        names[ai2]='m2'
        samples = MCSamples(samples=samps,names = names, labels = labels, 
         ranges={'m1':[mod1low,mod1up], 'm2':[mod2low,mod2up]})       
    else:
        samples = MCSamples(samples=samps,names = names, labels = labels, 
        ranges={'m1':[mod1low,mod1up]})
elif(np.any(il==5) and mass_string!='gNFW'):
    ai,=np.where(il==5)
    names[ai]='m2'
    samples = MCSamples(samples=samps,names = names, labels = labels,
                        ranges={'m2':[mod2low,mod2up]})
else:
    samples = MCSamples(samples=samps,names = names, labels = labels)

g = plots.get_subplot_plotter(width_inch=10.5)
print(samples)

g.settings.axes_fontsize = 20
g.settings.axes_labelsize=20
g.triangle_plot([samples], filled=True, line_args=[{'ls':'-', 'color':colore}], 
     contour_colors=[colore])
plt.show() 
plt.close('all')   # <--- importantissimo
del g              # opzionale, ma aiuta
#=============================================================================
print("\nconfidence intervals")
print("*******************************************************************")
print("{0:20s}\t best\t low 95\t low 68\t up 68\t up 95".format("param"))

# 1) best-fit index (MAP / min like==0)
best_idx = np.flatnonzero(like == 0)
if best_idx.size == 0:
    raise RuntimeError("No best-fit found: no entry with like==0")

# 2) best-fit point in parameter space (must match names order)
best_point = samps[best_idx[0], :]  # shape (ndim,)

# 3) compute stats once
mstats = samples.getMargeStats()

# 4) collect LaTeX strings (choose 95% by default, like your old code)
latex_cells = []

for i, (parname, label) in enumerate(zip(names, labels)):
    p = mstats.parWithName(parname)

    lim68 = p.limits[0]  # 68%
    lim95 = p.limits[1]  # 95%

    # center on best-fit (NOT mean)
    cen = float(best_point[i])

    # errors w.r.t. best-fit
    low68 = cen - lim68.lower
    low95 = cen - lim95.lower
    up68  = lim68.upper - cen
    up95  = lim95.upper - cen

    print("{0:20s}\t {1: .2f}\t -{2: .2f}\t -{3: .2f}\t +{4: .2f}\t +{5: .2f}".format(
        label, cen, low95, low68, up68, up95
    ))

    # LaTeX: 95% interval (change to up68/low68 if you prefer 68%)
    latex_cells.append(md.format_uncertainty(cen, up95, low95, ndigits=2))

print("*******************************************************************\n")

if LatexForm:
    print("Latex-formatted output - 95% confidence interval\n")
    print(" & ".join([f"${lab}$" for lab in labels]) + r" \\")
    print(" & ".join(latex_cells) + r" \\")


#==============================================================================

# MASS MODELS: The Burkert GC requires the computation of the screening radius
# I will define here a function specifically for the computation of the mass profile

#Burkert profile - field derivative 

#All the MG functions are x**2 dphidr*rs*Q*c**2
def dphib(x,rhos,rs,Q,phinf,expcutoff=False):
    clight=3e5
    limit=0.001 #limiting value for the screening radius
    G=4.302e-9
           
    #eponential cutoff:
    dminf=np.exp(-x*rs*np.sqrt(1e-7*Q/phinf))
    if expcutoff==False:
        dminf=1.0
    
    Facda=Q*rs**2*rhos*8*np.pi*G/clight**2 #factor containing the planck mass     
    #In this case, no analytic solution exists. Need to find numerical solution
    
    def fbur(y):
        return ((1 + y)*(4*phinf - Facda*np.pi + 
      Facda*(2*np.arctan(y)-2*np.log(1+y) + np.log(1 + y**2) )))/(4*y)          
    xsc=limit
    x0=limit
    if fbur(limit)<0:    
        for i in range(0,4000):
            xx=x0+0.01*i
            if fbur(x0)*fbur(xx) < 0:
                sol=opt.root(fbur,xx)
                if sol.success:
                    xsc=sol.x[0]
                else:
                    xsc=xx                
                break
            x0=xx
        if (xsc==limit):
            xsc=xx
    Cb= 1/4*(-4*phinf + Facda*np.pi - 2*Facda*np.log(1 + xsc**2)) 
    
    if (xsc>limit):
        if (x<xsc):
            out=0
        else:
            out=Cb+ Facda/(4)*(-2*np.arctan(x) + 
            2*np.log(1 + x)+np.log(1 + x**2)) #/x**2 
    else:
        out=Facda/(4)*(-2*np.arctan(x) + 
            2*np.log(1 + x)+np.log(1 + x**2))#/x**2
    return out*rs*Q*clight**2/G*dminf
dphiB=np.vectorize(dphib,otypes=[float],excluded=[1,2,3,4])


#generalized Chameleon - Burkert 
def M_bur_CS(r,r200,rs,tmass,Q,z,H0,Om,Ol,expcutoff=True):
    
    c200=r200/rs
    x=r/rs
    G=4.302e-9
    phinf=tmass*1e-5
    
    
    Hz2=H0**2*(Om*(1+z)**3+Ol)
    M200=100*Hz2/G*r200**3
     
    fac200=1.0/(np.log(1+c200*c200)+2.*np.log(1.0+c200)-2.*np.arctan(c200))
    mass=M200*fac200*(np.log(1+x*x)+2.*np.log(1.0+x)-2.*np.arctan(x))
     
    rhos=M200*fac200/(np.pi*rs**3) 
    if (Q!=0 and tmass!=0):  
        dphidr=dphiB(x,rhos,rs,Q,phinf,expcutoff)

    else:
        dphidr=0.0
    return mass+dphidr
     
     
M_bur_CS=np.vectorize(M_bur_CS,otypes=[float])#,excluded=[1,2,3,4,5,6,7,8,9])

#******************************************************************************

# ---------- helper: quantiles in one shot ----------
_QS = np.array([0.16, 0.50, 0.84, 0.025, 0.975])
def qstats(x):
    """Return (q16,q50,q84,q025,q975) for 1D array x."""
    return np.quantile(x, _QS)

# ---------- r grid ----------
r = np.logspace(np.log10(0.05), np.log10(facext * rupin), 60)

# If you want to plot the BCG, make a larger sample
if (kbcg == 1 and kbdt == 0):
    r = np.logspace(np.log10(0.001), np.log10(facext * rupin), 100)

nr = r.size

# ---------- output arrays (DM or main profile) ----------
Mlow   = np.empty(nr)   # 16%
Mmedian= np.empty(nr)   # 50%
Mup    = np.empty(nr)   # 84%
Mlow95 = np.empty(nr)   # 2.5%
Mup95  = np.empty(nr)   # 97.5%

# ---------- output arrays (BCG and totals, only if kbcg==1) ----------
if (kbcg == 1):
    Mlowbcg   = np.empty(nr)
    Mupbcg    = np.empty(nr)
    Mlowbcg95 = np.empty(nr)
    Mupbcg95  = np.empty(nr)

    # total = DM + BCG (only when mass_string != mTotal_GC)
    Mlowtot   = np.empty(nr)
    Muptot    = np.empty(nr)
    Mlowtot95 = np.empty(nr)
    Muptot95  = np.empty(nr)

# ============================================================
# 1) PRECOMPUTE SPECIAL STUFF (mTotal_GC only)
# ============================================================
num = max(1,int(len(r200)/numfac))
if (mass_string == "mTotal_GC"):
    
    rhob_chain = xml * xmstar / (4.0 * np.pi * rjaf**3)  # array, same length as chain

    xsc_chain = md.ScreenTot(
        Qcoup[::num], phicoup[::num] * 1e-5, r200[::num], rs[::num],
        rhostar, rstar, rhogas, rsgas,
        rhob_chain[::num], mod3[::num], rjaf[::num],
        z=za, H0=70, Om=0.3, Ol=0.7
    )

# ============================================================
# 2) DISPATCH: CHAIN PROFILE FUNCTION (computed for each ri)
# ============================================================
if mass_string == "NFW":
    prof_chain = lambda ri: md.M_nfw(ri, r200[::num], rs[::num], za, H0, Om, Ol)
elif mass_string == "Her":
    prof_chain = lambda ri: md.M_Her(ri, r200[::num], rs[::num], za, H0, Om, Ol)
elif mass_string == "Bur":
    prof_chain = lambda ri: md.M_Bur(ri, r200[::num], rs[::num], za, H0, Om, Ol)
elif mass_string == "mNFW_BH":
    prof_chain = lambda ri: md.M_nfw_BH(ri, r200[::num], rs[::num], 
                                        mod1[::num], za, H0, Om, Ol)
elif mass_string == "gNFW":
    prof_chain = lambda ri: md.M_gNFW(ri, r200[::num], rs[::num], 
                                      mod3[::num], za, H0, Om, Ol)
elif mass_string == "mNFW_GC":
    prof_chain = lambda ri: md.M_bnfw_CS(ri, r200[::num], rs[::num], 
                                         phicoup[::num], Qcoup[::num], eim, za, H0, Om, Ol)
elif mass_string == "Eis":
    prof_chain = lambda ri: md.M_Ein(ri, r200[::num], rs[::num], eim, za, H0, Om, Ol)
elif mass_string == "mBur_GC":
    prof_chain = lambda ri: M_bur_CS(ri, r200[::num], rs[::num], 
                                     phicoup[::num], Qcoup[::num], za, H0, Om, Ol)
elif mass_string == "mEin_GC":
    prof_chain = lambda ri: md.M_Ein_CS(ri, r200[::num], rs[::num], 
                                        eim, phicoup[::num], Qcoup[::num], za, H0, Om, Ol)
elif mass_string == "mgNFW_GC":
    prof_chain = lambda ri: md.M_gnfw_CS(ri, r200[::num], rs[::num],
                                         eim, phicoup[::num], Qcoup[::num], za, H0, Om, Ol)
elif mass_string == "mIso_GC":
    prof_chain = lambda ri: md.MIso_CS(ri, r200[::num], rs[::num], 
                                       phicoup[::num], Qcoup[::num], za, H0, Om, Ol)
elif mass_string == "RG":
    #Check if the gas parameters are arrays
    isarray = isinstance(rhogas, Iterable) and not isinstance(rhogas, (str, bytes))
   
    if isarray:
        prof_chain = lambda ri: md.M_RG(ri, epsilon0[::num], rhocRG[::num], mod3[::num], 
                                    xml, xmstar[::num], rjaf[::num], 
                                   rhogas[::num], rhostar[::num],
                                    rsgas[::num], rstar[::num])
    else:
        prof_chain = lambda ri: md.M_RG(ri, epsilon0[::num], rhocRG[::num], mod3[::num], 
                                    xml, xmstar[::num], rjaf[::num], rhogas, rhostar,
                                     rsgas, rstar)
elif mass_string == "BS":
    prof_chain = lambda ri: md.M_BS(ri, rs[::num], bexp[::num], gexp[::num], M200[::num], xma[::num])
elif mass_string == "mTotal_GC":
    prof_chain = lambda ri: md.Mtot(
        ri, Qcoup[::num], phicoup[::num]*1e-5, r200[::num], rs[::num],
        rhostar, rstar, rhogas, rsgas,
        rhob_chain[::num], mod3[::num], rjaf[::num],
        xsc=xsc_chain, z=za, H0=70, Om=0.3, Ol=0.7
    )
elif mass_string == "gNFW_CDE":
    prof_chain = lambda ri: md.M_gNFW(ri, r200[::num], rs[::num], 
        mod3[::num], za, H0, Om, Ol) * (1.0 + 3.0 * 10.0**(mod2 - 3.5) * (1.0 + mod1))
else:
    prof_chain = lambda ri: np.zeros_like(r200)

# ============================================================
# 3) LOOP OVER r: quantiles WITHOUT sorting
# ============================================================
for j, ri in enumerate(r):
    profile = prof_chain(ri)

    q16, q50, q84, q025, q975 = qstats(profile)
    Mlow[j], Mmedian[j], Mup[j] = q16, q50, q84
    Mlow95[j], Mup95[j] = q025, q975

    if kbcg == 1:
        # BCG profile with its own parameter uncertainties
        profbcg = md.jaffe_BCG(ri, xml, xmstar[::num], rjaf[::num])

        b16, b50, b84, b025, b975 = qstats(profbcg)
        Mlowbcg[j], Mupbcg[j] = b16, b84
        Mlowbcg95[j], Mupbcg95[j] = b025, b975

        # total = DM + BCG (only in the "simple sum" cases)
        if mass_string != "mTotal_GC":
            proftot = profbcg + profile
            t16, t50, t84, t025, t975 = qstats(proftot)
            Mlowtot[j], Muptot[j] = t16, t84
            Mlowtot95[j], Muptot95[j] = t025, t975

# ============================================================
# 4) BEST-FIT (minimum) CURVE Mnfw(r) + rm2 (and MDE for gNFW_CDE)
#    (Uses your "*min" arrays like r2min[0], rsmin[0], etc.)
# ============================================================
plot = True
MDE = None  # only for gNFW_CDE

if mass_string == "NFW":
    Mnfw = md.M_nfw(r, r2min[0], rsmin[0], za, H0, Om, Ol)
    rm2 = rs
elif mass_string == "Her":
    Mnfw = md.M_Her(r, r2min[0], rsmin[0], za, H0, Om, Ol)
    rm2 = 0.5 * rs
elif mass_string == "Bur":
    Mnfw = md.M_Bur(r, r2min[0], rsmin[0], za, H0, Om, Ol)
    rm2 = 1.521 * rs
elif mass_string == "mNFW_BH":
    Mnfw = md.M_nfw_BH(r, r2min[0], rsmin[0], mod1min[0], za, H0, Om, Ol)
    rm2 = rs
elif mass_string == "gNFW":
    Mnfw = md.M_gNFW(r, r2min[0], rsmin[0], mod3min[0], za, H0, Om, Ol)
    rm2 = (2.0 - mod3) * rs
elif mass_string == "mNFW_GC":
    Mnfw = md.M_bnfw_CS(r, r2min[0], rsmin[0], mod1min[0], mod2min[0], eim, za, H0, Om, Ol)
    rm2 = rs / eim
elif mass_string == "Eis":
    Mnfw = md.M_Ein(r, r2min[0], rsmin[0], eim, za, H0, Om, Ol)
    rm2 = rs
elif mass_string == "mBur_GC":
    Mnfw = M_bur_CS(r, r2min[0], rsmin[0], mod1min[0], mod2min[0], za, H0, Om, Ol)
    rm2 = 1.521 * rs
elif mass_string == "mEin_GC":
    Mnfw = md.M_Ein_CS(r, r2min[0], rsmin[0], eim, mod1min[0], mod2min[0], za, H0, Om, Ol)
    rm2 = rs
elif mass_string == "RG":
    
    Mnfw = md.M_RG(
        r, mod1min[0], mod2min[0], mod3min[0], xml, xmstarmin[0], rjafmin[0],
        pars["rho0gas"], pars["rhostar"], pars["rsgas"], pars["rstar"])  
    
    rm2 = rs
elif mass_string == "BS":
    Mnfw = md.M_BS(r, rsmin, mod1min[0], mod2min[0], mod3min[0], np.amin(M200))
    rm2 = rs
elif mass_string == "mgNFW_GC":
    Mnfw = md.M_gnfw_CS(r, r2min[0], rsmin[0], eim, mod1min[0], mod2min[0], za, H0, Om, Ol)
    rm2 = (2.0 - mod3) * rs
elif mass_string == "mIso_GC":
    Mnfw = md.MIso_CS(r, r2min[0], rsmin[0], mod1min[0], mod2min[0], za, H0, Om, Ol)
    rm2 = rs
elif mass_string == "mTotal_GC":
    rm2 = (2.0 - mod3) * rs

    rhob_best = xml * xmstarmin[0] / (4.0 * np.pi * rjafmin[0]**3)
    xsc_best = md.ScreenTot(
        mod1min[0], mod2min[0] * 1e-5, r2min[0], rsmin[0],
        rhostar, rstar, rhogas, rsgas, rhob_best, mod3min[0], rjafmin[0],
        z=0.0, H0=70, Om=0.3, Ol=0.7
    )
    Mnfw = md.Mtot(
        r, mod1min[0], mod2min[0] * 1e-5, r2min[0], rsmin[0],
        rhostar, rstar, rhogas, rsgas, rhob_best, mod3min[0], rjafmin[0],
        xsc=xsc_best, z=za, H0=H0, Om=Om, Ol=Ol
    )
elif mass_string == "gNFW_CDE":
    rm2 = (2.0 - mod3) * rs
    base = md.M_gNFW(r, r2min[0], rsmin[0], mod3min[0], za, H0, Om, Ol)

    fac = 3.0 * 10.0**(mod2min[0] - 3.5) * (1.0 + mod1min[0])
    Mnfw = base * (1.0 + fac)
    MDE  = base * fac
else:
    plot = False
    Mnfw = np.zeros_like(r)
    rm2 = rs  # fallback

# ============================================================
def load_baryon_data(namegas = namegas, namegal = namegal):
    tgal = np.loadtxt(namegal) 
    tgas = np.loadtxt(namegas)
    
    # returns arrays in Msun (your files are in 1e13 units)
    gal = (tgal[:, 0], tgal[:, 1] * 1e13)
    gas = (tgas[:, 0], tgas[:, 1] * 1e13)
    gas_up  = (tgas[:, 0], (tgas[:, 1] + tgas[:, 2]) * 1e13)
    gas_low = (tgas[:, 0], (tgas[:, 1] - tgas[:, 2]) * 1e13)
    return gal, gas, gas_up, gas_low

# ============================================================
# PLOT
# ============================================================
if plot:

    # ----------------------------
    # 0) Choose main curve to plot
    # ----------------------------
    Mnfw_plot = Mmedian if median_string else Mnfw

    # ----------------------------
    # 1) Figure / axes cosmetics
    # ----------------------------
    fig, ax = plt.subplots(figsize=(9, 9))
    # ax.grid(True)

    ax.tick_params(axis="x", labelsize=fontlabel)
    ax.tick_params(axis="y", labelsize=fontlabel)

    #ax.set_title(mass_string, fontsize=fontlabel)
    ax.set_xlabel("r [Mpc]", fontsize=font)
    ax.set_ylabel(r"M [$M_{\odot}$]", fontsize=font)

    # ----------------------------
    # 2) Base x-limits (set ONCE, then refine)
    # ----------------------------
    xmin = 0.05
    if kbcg == 1 or mass_string == "mTotal_GC":
        xmin = 0.001 if (kbdt == 0 and mass_string != "RG") else 0.05
    xmax = facext * rupin
    ax.set_xlim(xmin, xmax)

    # ----------------------------
    # 3) DM credible bands (always)
    # ----------------------------
    ax.fill_between(r, Mlow95, Mup95, alpha=0.3, color="black")
    ax.fill_between(r, Mlow,   Mup,   alpha=0.3, color="black")

    # ----------------------------
    # 4) Optional vertical reference
    # ----------------------------
    if mass_string != "RG":
        ax.axvline(r2min[0], color="black", ls="--")

    # ----------------------------
    # 5) BCG + baryons + total (only if kbcg==1)
    # ----------------------------
    # flags
    show_bcg_band = (kbcg == 1) and (mass_string != "RG") and (kbdt == 0)
    can_build_total = (kbcg == 1) and (mass_string not in {"RG", "mTotal_GC"})

    if kbcg == 1:
        # --- BCG best-fit curve ---
        bcg_best = md.jaffe_BCG(r, xml, xmstarmin[0], rjafmin[0])

        # --- Gas+galaxies (if you want to overlay) ---
        # li calcoliamo SOLO se ci serve total, oppure se vuoi comunque plottarli sempre.
        need_baryons = True
        Mgal = Mgas = None

        if need_baryons:
            if getdata:
                (rgal_d, Mgal_d), (rgas_d, Mgas_d), (rgas_u, Mgas_u), \
                    (rgas_l, Mgas_l) = load_baryon_data()
                    
                ax.plot(rgal_d, Mgal_d, label="galaxies")
                ax.plot(rgas_d, Mgas_d, label="gas")
                ax.fill_between(rgas_d, Mgas_l, Mgas_u, alpha=0.4, color="green")

                # p Mgal(r) e Mgas(r) in r grid to build mtotal
                if can_build_total:
                    Mgal = ipt.interp1d(rgal_d, Mgal_d, kind="cubic",
                                        bounds_error=False, fill_value="extrapolate")(r)
                    Mgas = ipt.interp1d(rgas_d, Mgas_d, kind="cubic",
                                        bounds_error=False, fill_value="extrapolate")(r)
            else:
                # analytical profiles (best fit)
                Mgal = md.Mstar_nfw_like(r, pars["rstar"], pars["rhostar"])
                Mgas = md.Mgas_beta_like(r,  pars["rsgas"], pars["rho0gas"])
                ax.plot(r, Mgal, label="galaxies")
                ax.plot(r, Mgas, label="gas")

        # --- BCG curve + BCG uncertainty band ---
        if show_bcg_band:
            ax.loglog(r, bcg_best, ls="--", label="BCG")
            ax.fill_between(r, Mlowbcg95, Mupbcg95, alpha=0.2)

        # --- Total curve + total band (DM+BCG) shifted by baryons (deterministic shift) ---
        if can_build_total:
            Mtot_best = bcg_best + np.asarray(Mnfw_plot, dtype=float)

            # baryonic shift on grid r
            if getdata:
                xmbar = (Mgas + Mgal) if (Mgas is not None and Mgal is not None) else 0.0
            else:
                xmbar = (Mgas + Mgal) if (Mgas is not None and Mgal is not None) else 0.0

            Mtotal2 = Mtot_best + xmbar
            ax.plot(r, Mtotal2, label="total")

            # total CI band (DM+BCG) shifted by baryons
            ax.fill_between(r, Mlowtot95 + xmbar, Muptot95 + xmbar, alpha=0.1, color="green")

    # ----------------------------
    # 6) Main curve(s)
    # ----------------------------
    if (kbcg == 1 and mass_string not in {"RG", "mTotal_GC"}):
        if mass_string == "gNFW_CDE":
            ax.loglog(r, Mnfw_plot, label="Dark Matter+Dark Energy")
            if MDE is not None:
                ax.loglog(r, MDE, label="Dark Energy")
        else:
            ax.loglog(r, Mnfw_plot, label="Dark Matter", color="black")
    else:
        ax.loglog(r, Mnfw_plot)

    # ----------------------------
    # 7) Optional y-limits in kbdt!=0 mode (kept from your logic)
    # ----------------------------
    if (kbdt != 0) and ((kbcg == 1) or (mass_string == "mTotal_GC")) and (mass_string != "RG"):
        ax.set_ylim(1e10, 5e15)

    # ----------------------------
    # 8) Legend (once)
    # ----------------------------
    handles, labels = ax.get_legend_handles_labels()
    if labels:
        ax.legend(fontsize=24)

    plt.show()

    # ============================================================
    # Anisotropy profiles (optimized)
    # ============================================================
    
    # thinning (se non esiste num, fai num=1)
    # num = 1
    
    # b = 1 - 1/beta^2  (come nel tuo codice)
    bb  = 1.0 - 1.0 / (beta**2)
    bb0 = 1.0 - 1.0 / (beta2**2)
    
    # ---- choose nbeta ONCE ----
    ANIS_MAP = {
        "C":   1,
        "ML":  2,
        "OM":  3,
        "T":   4,
        "O":   5,
        "gOM": 6,
        "gT":  7,
        "BP":  8,
    }
    nbeta = ANIS_MAP.get(anis_string)
    if nbeta is None:
        raise ValueError(f"Unknown anis_string='{anis_string}'")
    
    # ---- choose rbetac chain ONCE (same logic as your ifs) ----
    rbetac_chain = rm2  # default
    
    # ML or OM use 'beta' as scale radius (your code: if nbeta==2 or 3 -> rbetac=beta)
    if nbeta in (2, 3):
        rbetac_chain = beta
    
    # gT: if freebeta>0 use rbeta
    elif nbeta == 7 and freebeta > 0:
        rbetac_chain = rbeta
    
    # BP always uses rbeta
    elif nbeta == 8:
        rbetac_chain = rbeta
    
    # ---- thin chains once ----
    rbetac_th = rbetac_chain[::num]
    bb_th     = bb[::num]
    bb0_th    = bb0[::num]
    
    # ---- preallocate outputs ----
    nr = len(r)
    bmlow   = np.empty(nr)
    bmedian = np.empty(nr)
    bmup    = np.empty(nr)
    bmlow95 = np.empty(nr)
    bmup95  = np.empty(nr)
    
    # ---- loop over radii, NO sorting ----
    for j, ri in enumerate(r):
        betaprof = md.betaf(ri, rbetac_th, bb_th, nbeta, bb0_th)
    
        q16, q50, q84, q025, q975 = qstats(betaprof)
        bmlow[j], bmedian[j], bmup[j] = q16, q50, q84
        bmlow95[j], bmup95[j] = q025, q975

        
    plt.figure(figsize=(9,9))
    plt.tick_params(axis="x", labelsize=fontlabel)
    plt.tick_params(axis="y", labelsize=fontlabel)
    plt.fill_between(r, bmup95, bmlow95, alpha=.2,color='royalblue')
    plt.fill_between(r, bmup, bmlow, alpha=.4,color='royalblue')
    
    
    il,= np.where(like==np.amin(like))
    
    plt.xlim(0.05,rupin) #rupin)
    plt.ylim(-1.0,1.0)
    
    plt.title('Anisotropy model: ' + anis_string, fontsize=font)
    plt.minorticks_on()
    plt.plot(r,bmedian,color='black') #, label = r"BP, $r_{\beta}$ free" )
    #plt.plot(r,md.betaf(r,rbetac[il],bb[il],nbeta,bb0[il]))
    plt.xlabel('r [Mpc]',size=font)
    plt.ylabel(r'$\beta(r)$',size=font)
    
    plt.show()
    saveBeta = 'beta_'+anis_string+'_MASS_'+mass_string+'.dat'
    outbeta = np.column_stack((r,bmedian,bmlow,bmup,bmlow95,bmup95))
    head = 'r [Mpc]\t beta\t low68\t up68\t low95\t up95\t'
    np.savetxt(saveBeta,outbeta,header = head, fmt ='%3.4f')

