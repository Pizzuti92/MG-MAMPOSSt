import scipy.special as spec
import numpy as np
import scipy.integrate as ints
from math import *
import scipy.optimize as opt
import sys
import matplotlib.pyplot as mpl
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d as plt3d
import matplotlib.lines as ml
import scipy.interpolate as ipt
from getdist import plots, MCSamples
from scipy.ndimage.filters import gaussian_filter


def matrix(file):
     contents = open(file).read()
     return [item.split() for item in contents.split('\n')[:-1]]
     

titfile='gomamposst_x.inp'
a=(matrix(titfile)[6])
titinp=a[0]


data=np.array(matrix(titinp)[:-4]).astype(float)
r200=data[:,0]
rnu=data[:,1]
rs=data[:,2]
beta=data[:,3]
mod1=data[:,4]
mod2=data[:,5]
like=data[:,6]-np.amin(data[:,6]) #scaled likelihood


#add weigths in the marginalization
ww=np.ones(len(data[:,0]))
weig=1/(r200*rs*beta)

bb=(matrix(titfile)[1])
titinput=bb[0]
massmodel=(matrix(titinput)[22])[0]
massmodel=int(massmodel[0])

#NOTE: In the case of general chameleon the RESCALED VARIABLES Q_2 and \phi_2
#are plotted (see Pizzuti et al., 2021)
#fr case
nhone=(matrix(titinput)[5])[0]
frcase=0
try:
    int(nhone[0])
except:
    frcase=1

f = open('Options.txt', 'r')
for line in map(str.strip,f.read().splitlines()):
    line1 = line.split('=')
    line2 = line.split(' ')
    if 'nmcmc' == line1[0] or 'nmcmc'==line2[0]:
        if(line.find('1')!=-1):
            nmcmc=1
        elif(line.find('2')!=-1):
            nmcmc=2
        elif(line.find('3')!=-1):
            nmcmc=3
        elif(line.find('4')!=-1):
            nmcmc=4            
        else:
            nmcmc=0
f.close()
if (nmcmc<1):
	print("No plot for the grid case")
	exit()

#find which are the free parameters
sigma=1.0
fontlabel=20
font=22
total=0
freep=np.zeros(6)

for i in range(0,6):
    if(np.any(data[:,i]!=data[0,i])):
        freep[i]=1
        
            
il,=np.where(freep==1)
total=len(il)
print(total)

labelb=r"r_{\beta}\,\, [Mpc]"
Plabelb=r"P(r_{\beta})"
if(data[7,0]!=1):
    labelb="\mathcal{A}_\infty"
    Plabelb="P(\mathcal{A}_\infty)"

labelm1="\phi_2"
Plabelm1=r"$P(\phi_2)$"
if (frcase==1):
    labelm1="Log10(\phi/c^2)"
labelm2="\mathcal{Q}_2"
Plabelm2=r"$P(\mathcal{Q}_2)$"

if (frcase==1):
    labelm1="Log10(\phi/c^2)"
    Plabelm1=r"$P[Log10(\phi/c^2)]$"


if(massmodel==8):
    labelm1="Y_1"
    labelm2="Y_2"
    Plabelm1="r$P(Y_1)$"
    Plabelm2="r$P(Y_2)$"
elif (massmodel==7):
	labelm1="Log(m/Mpc^{-1})"
	labelm2=r"Q"
	Plabelm1=r"P[Log$(m/Mpc^{-1})$]"
	Plabelm2=r"$P(Q)$"
elif (massmodel==1):
    massm2=(matrix(titinput)[22])[0]
    if(massm2[1]=='0'):
        labelm1=r"\gamma"
        Plabelm1=r"$P(\gamma)$"
if (massmodel==7 and np.any(il==5)):
    data[:,5]=np.log10(data[:,5]) #uses log values for Q in the linear
    #Horndeski case to avoid problems in the getdist marginalizations
    labelm2=r"Log(Q)"
    Plabelm2=r"P[Log$(Q)$]"
if(np.any(il==4)):
    if (massmodel==7):
        data[:,4]=np.log10(data[:,4])
    if (massmodel==9 and  frcase==1):
        data[:,4]=np.log10(data[:,4])
    if (massmodel==9 and frcase==0):
        data[:,4]=1-np.exp(-data[:,4]*0.1)
        data[:,5]=data[:,5]/(1+data[:,5])  
  
label=np.array(["r_{200}\,\,[Mpc]",r"r_{\nu}\,\, [Mpc]",r"r_{s}\,\, [Mpc]",labelb,labelm1,labelm2])
label1d=np.array([r"$r_{200}$ [Mpc]",r"$r_{\nu}$ [Mpc]",r"$r_{s}$ [Mpc]",labelb,labelm1,labelm2])

Plabel=np.array([r"$P(r_{200})$",r"$P(r_{\nu})$",r"$P(r_{s})$",Plabelb,Plabelm1,Plabelm2])


mod1up=float((matrix(titinput)[38])[0])
mod1low=float((matrix(titinput)[37])[0])

mod2up=float((matrix(titinput)[40])[0])
mod2low=float((matrix(titinput)[39])[0])
if (massmodel==7 or (massmodel==9 and frcase==1)):
    mod1up=np.log10(mod1up)
    mod1low=np.log10(mod1low)
if (massmodel==9 and frcase==0):
    mod1up=1-np.exp(-mod1up*0.1)
    mod1low=1-np.exp(-mod1low*0.1)
    mod2up=0.0
    mod2up=1.1

if (massmodel==7 and np.any(il==5)):
    mod2up=np.log10(mod2up)
    mod2low=np.log10(mod2low)


if total==1:
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
    mpl.grid()
    mpl.ylabel(Ylabels[0], fontsize=22)
    mpl.xlabel(labels[0], fontsize=22)
    mpl.tick_params(axis="x", labelsize=fontlabel)
    mpl.tick_params(axis="y", labelsize=fontlabel)
    mpl.plot(plotf[:,0],plotf[:,1], lw=3,c='black')
    mpl.show()
    print(r"Confidence interval at 2$\sigma$")
    first=np.sort(first)
    print(first[int(len(first)*0.025)],first[int(len(first)*0.975)])
    
elif total==2:
    first=data[:,il[0]]
    second=data[:,il[1]]
    samps=np.column_stack((first,second))
    ndim = 2
    names = np.array(["x%s"%i for i in range(ndim)])
    labels =  label[il]
    samples = MCSamples(samples=samps,names = names, labels = labels)
    if(np.any(il==4)):
        ai,=np.where(il==4)
        names[ai]='m1'
        if(np.any(il==5)): 
            ai2,=np.where(il==5)
            names[ai2]='m2'
            samples = MCSamples(samples=samps,names = names, labels = labels,ranges={'m1':[mod1low,mod1up],
              'm2':[mod2low,mod2up]}) #,weights=1/10**first)#,loglikes=like)
        else:
            samples = MCSamples(samples=samps,names = names, labels = labels,ranges={'m1':[mod1low,mod1up]})
    elif(np.any(il==5)):
        ai,=np.where(il==5)
        names[ai]='m2'
        samples = MCSamples(samples=samps,names = names, labels = labels,ranges={'m2':[mod2low,mod2up]})#,loglikes=like)
    g = plots.get_subplot_plotter(width_inch=7.5)
    g.settings.axes_fontsize = 20
    g.settings.axes_labelsize=20
    g.triangle_plot([samples], filled=True, color='red')
    mpl.show()
    print(r"Confidence interval at 2$\sigma$")
    print(samples.getTable().tableTex())
    
elif total==3:
    first=data[:,il[0]]
    second=data[:,il[1]]
    third=data[:,il[2]]
    samps=np.column_stack((first,second,third))
    ndim = 3
    names = np.array(["x%s"%i for i in range(ndim)])
    labels =  label[il]
    samples = MCSamples(samples=samps,names = names, labels = labels)
    if(np.any(il==4)):
        ai,=np.where(il==4)
        names[ai]='m1'
        if(np.any(il==5)): 
            ai2,=np.where(il==5)
            names[ai2]='m2'
            samples = MCSamples(samples=samps,names = names, labels = labels,ranges={'m1':[mod1low,mod1up],
              'm2':[mod2low,mod2up]}) #,weights=1/10**first)#,loglikes=like)
        else:
            samples = MCSamples(samples=samps,names = names, labels = labels,ranges={'m1':[mod1low,mod1up]})
    elif(np.any(il==5)):
        ai,=np.where(il==5)
        names[ai]='m2'
        samples = MCSamples(samples=samps,names = names, labels = labels,ranges={'m2':[mod2low,mod2up]})#,loglikes=like)
    g = plots.get_subplot_plotter(width_inch=7.5)
    g.settings.axes_fontsize = 20
    g.settings.axes_labelsize=20
    g.triangle_plot([samples], filled=True, color='red')
    mpl.show()
    print(r"Confidence interval at 2$\sigma$")
    print(samples.getTable().tableTex())
    

elif total==4:
    first=data[:,il[0]]
    second=data[:,il[1]]
    third=data[:,il[2]]
    four=data[:,il[3]]
    samps=np.column_stack((first,second,third,four))
    ndim = 4
    names = np.array(["x%s"%i for i in range(ndim)])
    labels =  label[il]
    samples = MCSamples(samples=samps,names = names, labels = labels)
    if(np.any(il==4)):
        ai,=np.where(il==4)
        names[ai]='m1'
        if(np.any(il==5)): 
            ai2,=np.where(il==5)
            names[ai2]='m2'
            samples = MCSamples(samples=samps,names = names, labels = labels, ranges={'m1':[mod1low,mod1up],
              'm2':[mod2low,mod2up]}) #,weights=1/10**first)#,loglikes=like)
        else:
            samples = MCSamples(samples=samps,names = names, labels = labels, ranges={'m1':[mod1low,mod1up]})
    elif(np.any(il==5)):
        ai,=np.where(il==5)
        names[ai]='m2'
        samples = MCSamples(samples=samps,names = names, labels = labels,ranges={'m2':[mod2low,mod2up]})#,loglikes=like)
    g = plots.get_subplot_plotter(width_inch=10.5)
    g.settings.axes_fontsize = 20
    g.settings.axes_labelsize=20
    g.triangle_plot([samples], filled=True, color='red')
    mpl.show()
    print(r"Confidence interval at 2$\sigma$")
    print(samples.getTable().tableTex())
    
elif total==5:
    first=data[:,il[0]]
    second=data[:,il[1]]
    third=data[:,il[2]]
    four=data[:,il[3]]
    fifth=data[:,il[4]]
    samps=np.column_stack((first,second,third,four,fifth))
    ndim = 5
    nsamp = 10000
# Get the getdist MCSamples objects for the samples, specifying same parameter
# names and labels; if not specified weights are assumed to all be unity
    names = np.array(["x%s"%i for i in range(ndim)])
    labels =  label[il]

    samples = MCSamples(samples=samps,names = names,ranges={'x3':[mod1low,mod1up]},
    labels = labels)#, weights=weig,loglikes=like)
    if(np.any(il==4)):
        ai,=np.where(il==4)
        names[ai]='m1'
        if(np.any(il==5)): 
            ai2,=np.where(il==5)
            names[ai2]='m2'
            samples = MCSamples(samples=samps,names = names, labels = labels, ranges={'m1':[mod1low,mod1up],
              'm2':[mod2low,mod2up]}) #,weights=1/10**first)#,loglikes=like)
        else:
            samples = MCSamples(samples=samps,names = names, labels = labels, ranges={'m1':[mod1low,mod1up]})
    elif(np.any(il==5)):
        ai,=np.where(il==5)
        names[ai]='m2'
        samples = MCSamples(samples=samps,names = names, labels = labels,ranges={'m2':[mod2low,mod2up]})#,loglikes=like)
    g = plots.get_subplot_plotter(width_inch=10.5)
    g.settings.axes_fontsize = 20
    g.settings.axes_labelsize=20
    g.triangle_plot([samples], filled=True, color='red',param_limits={'m2':[mod2low,mod2up]})
    mpl.show()
    print(r"Confidence interval at 2$\sigma$")
    print(samples.getTable().tableTex())

elif total==6:
        
    first=data[:,il[0]]
    second=data[:,il[1]]
    third=data[:,il[2]]
    four=data[:,il[3]]
    fifth=data[:,il[4]]
    six=data[:,il[5]]
    samps=np.column_stack((first,second,third,four,fifth,six))
    ndim = 6
    nsamp = 10000
# Get the getdist MCSamples objects for the samples, specifying same parameter
# names and labels; if not specified weights are assumed to all be unity
    names = np.array(["x%s"%i for i in range(ndim)])
    labels =  label[il]
    if(np.any(il==4)):
        ai,=np.where(il==4)
        names[ai]='m1'
        if(np.any(il==5)): 
            ai2,=np.where(il==5)
            names[ai2]='m2'
            samples = MCSamples(samples=samps,names = names, labels = labels, ranges={'m1':[mod1low,mod1up],
              'm2':[mod2low,mod2up]}) #,weights=1/10**first)#,loglikes=like)
        else:
            samples = MCSamples(samples=samps,names = names, labels = labels, ranges={'m1':[mod1low,mod1up]})
    elif(np.any(il==5)):
        ai,=np.where(il==5)
        names[ai]='m2'
        samples = MCSamples(samples=samps,names = names, labels = labels,ranges={'m2':[mod2low,mod2up]})#,loglikes=like)
    g = plots.get_subplot_plotter(width_inch=10.5)
    g.settings.axes_fontsize = 20
    g.settings.axes_labelsize=20
    g.triangle_plot([samples], filled=True, color='red')
    mpl.show()
    print(r"Confidence interval at 2$\sigma$")
    print(samples.getTable().tableTex())

mpl.show()
