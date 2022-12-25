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


data=np.array(matrix(titinp)[3:-4]).astype(float)
r200=data[:,0]
rnu=data[:,1]
rs=data[:,2]
beta=data[:,3]
mod1=data[:,4]
mod2=data[:,5]
beta2=data[:,6]
like=data[:,7]-np.amin(data[:,7]) #scaled likelihood
r2min=r200[np.where(like==0)]
rsmin=rs[np.where(like==0)]
betamin=beta[np.where(like==0)]
beta2min=beta2[np.where(like==0)]
mod1min=mod1[np.where(like==0)]
mod2min=mod2[np.where(like==0)]

#add weigths in the marginalization
ww=np.ones(len(data[:,0]))
weig=1/(r200*rs*beta)

bb=(matrix(titfile)[1])
#titinput='Output/log_for_plot'
#massmodel=(matrix(titinput)[6])[0]
#massmodel=int(massmodel[0])
#
##NOTE: In the case of general chameleon the RESCALED VARIABLES Q_2 and \phi_2
##are plotted (see Pizzuti et al., 2021)
##fr case
#nhone=(matrix(titinput)[5])[0]
frcase=0
#try:
#    int(nhone[0])
#except:
#    frcase=1

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

fpar = open(bb[0],'r') #
for line in map(str.strip,fpar.read().splitlines()):
    line1 = line.split('=')
    line2 = line.split(' ')  #perché posso mettere anche uno spazio prima dell'uguale

    if 'za' == line1[0] or 'za' ==line2[0]:
        try:
            za=float(line1[1]) #è sempre dopo l'uguale!
        except ValueError:
            za=0.0
    if 'H0' == line1[0] or 'H0' ==line2[0]:
        try:
            H0=float(line1[1]) #è sempre dopo l'uguale!
        except ValueError:
            H0=70.0
            
    if 'M(r)' == line1[0] or 'M(r)' ==line2[0]:
        try:
            mass_string=(line1[1].split(' '))[-1] #ultimo elemento
        except ValueError:
            mass_string='NFW'
    if 'nA2'== line1[0] or 'nA2'==line2[0]:
        try:
            nhone=int(line1[1])
        except ValueError: 
            nhone=0
    if 'A1low'== line1[0] or 'A1low'==line2[0]:
        try:
            mod1low=float(line1[1])
        except ValueError: 
            mod1low=np.amin(mod1)
    if 'A1up'== line1[0] or 'A1up'==line2[0]:
        try:
            mod1up=float(line1[1])
        except ValueError: 
            mod1up=np.amax(mod1)
    if 'A2low'== line1[0] or 'A2low'==line2[0]:
        try:
            mod2low=float(line1[1])
        except ValueError: 
            mod2low=np.amin(mod2)
    if 'A2up'== line1[0] or 'A2up'==line2[0]:
        try:
            mod2up=float(line1[1])
        except ValueError: 
            mod2up=np.amax(mod2)
    if 'Rup'== line1[0] or 'Rup'==line2[0]:
        try:
            rupin=float(line1[1])
        except ValueError: 
            rupin=np.amax(r200)

            
fpar.close()
if nhone==-1:
    frcase=1
#find which are the free parameters
sigma=1.0
fontlabel=20
font=22
total=0
freep=np.zeros(7)

for i in range(0,7):
    if(np.any(data[:,i]!=data[0,i])):
        freep[i]=1
        
            
il,=np.where(freep==1)
total=len(il)


labelb=r"r_{\beta}\,\, [Mpc]"
Plabelb=r"P(r_{\beta})"

if(data[0,7]!=1):
    labelb="\mathcal{A}_\infty"
    Plabelb="P(\mathcal{A}_\infty)"
    
labelb2="\mathcal{A}_0"
Plabelb2="P(\mathcal{A}_0)"
labelm1="\phi_2"
Plabelm1=r"$P(\phi_2)$"
if (frcase==1):
    labelm1="Log10(\phi/c^2)"
labelm2="\mathcal{Q}_2"
Plabelm2=r"$P(\mathcal{Q}_2)$"

if (frcase==1):
    labelm1="Log10(\phi/c^2)"
    Plabelm1=r"$P[Log10(\phi/c^2)]$"


if (mass_string=='mNFW_BH'):
    labelm1="Y_1"
    labelm2="Y_2"
    Plabelm1="r$P(Y_1)$"
    Plabelm2="r$P(Y_2)$"
elif (mass_string=='mNFW_LH'):
	labelm1="Log(m/Mpc^{-1})"
	labelm2=r"Q"
	Plabelm1=r"P[Log$(m/Mpc^{-1})$]"
	Plabelm2=r"$P(Q)$"
elif (mass_string=='gNFW'):
#    massm2=(matrix(titinput)[6])[0]
#    if(massm2[1]=='0'):
    labelm1=r"\gamma"
    Plabelm1=r"$P(\gamma)$"
if (mass_string=='mNFW_LH' and np.any(il==5)):
    data[:,5]=np.log10(data[:,5]) #uses log values for Q in the linear
    #Horndeski case to avoid problems in the getdist marginalizations
    labelm2=r"Log(Q)"
    Plabelm2=r"P[Log$(Q)$]"
if(np.any(il==4)):
    if (mass_string=='mNFW_LH'):
        data[:,4]=np.log10(data[:,4])
    if (mass_string=='mNFW_GC' and  frcase==1):
        data[:,4]=np.log10(data[:,4])
    if (mass_string=='mNFW_GC' and frcase==0):
        Qcoup=np.copy(data[:,5])
        phicoup=np.copy(data[:,4])
        data[:,4]=1-np.exp(-data[:,4]*0.1)
        data[:,5]=data[:,5]/(1+data[:,5])  
  
label=np.array(["r_{200}\,\,[Mpc]",r"r_{\nu}\,\, [Mpc]",r"r_{s}\,\, [Mpc]",labelb,labelm1,labelm2,labelb2])
label1d=np.array([r"$r_{200}$ [Mpc]",r"$r_{\nu}$ [Mpc]",r"$r_{s}$ [Mpc]",labelb,labelm1,labelm2,labelb2])

Plabel=np.array([r"$P(r_{200})$",r"$P(r_{\nu})$",r"$P(r_{s})$",Plabelb,Plabelm1,Plabelm2,Plabelb2])


#mod1up=float((matrix(titinput)[17])[0])
#mod1low=float((matrix(titinput)[16])[0])
#
#mod2up=float((matrix(titinput)[19])[0])
#mod2low=float((matrix(titinput)[18])[0])
if (mass_string=='mNFW_LH' or (mass_string=='mNFW_GC' and frcase==1)):
    mod1up=np.log10(mod1up)
    mod1low=np.log10(mod1low)
if (mass_string=='mNFW_GC' and frcase==0):
    mod1up=1-np.exp(-mod1up*0.1)
    mod1low=1-np.exp(-mod1low*0.1)
    mod2up=0.0
    mod2up=1.1

if (mass_string=='mNFW_LH' and np.any(il==5)):
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
#    nsamp = 10000
    names = np.array(["x%s"%i for i in range(ndim)])
    labels =  label[il]
    if(np.any(il==4)):
        ai,=np.where(il==4)
        names[ai]='m1'
        if(np.any(il==5)): 
            ai2,=np.where(il==5)
            names[ai2]='m2'

            samples = MCSamples(samples=samps,names = names, labels = labels, ranges={'m1':[mod1low,mod1up],
              'm2':[mod2low,mod2up]})#,loglikes=like)
        else:
            samples = MCSamples(samples=samps,names = names, labels = labels, ranges={'m1':[mod1low,mod1up]})
    elif(np.any(il==5)):
        ai,=np.where(il==5)
        names[ai]='m2'
        samples = MCSamples(samples=samps,names = names, labels = labels,ranges={'m2':[mod2low,mod2up]})#,loglikes=like)
    g = plots.get_subplot_plotter(width_inch=10.5)
    g.settings.axes_fontsize = 20
    g.settings.axes_labelsize=20
    g.triangle_plot([samples], filled=True, color='red')#,param_limits={'m1':[-0.67,8.2]})#,'m2':[mod2low,mod2up]})
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
#plot of mass profile

def M_nfw(r,r200,rs):
    c200=r200/rs
    x=r/rs
    G=4.302e-9
    z=za
    Hz2=H0**2*(0.3*(1+z)**3+0.7)
    M200=100*Hz2/G*r200**3
    fac200=np.log(1+c200)-c200/(1+c200)
    return (np.log(1+x)-x/(1+x))/fac200*M200

def M_nfw_BH(r,r200,rs,Y1):
    c200=r200/rs
    x=r/rs
    G=4.302e-9
    z=za
    Hz2=H0**2*(0.3*(1+z)**3+0.7)
    M200=100*Hz2/G*r200**3
    fac200=np.log(1+c200)-c200/(1+c200)
    return (Y1*r**2*(rs-r)/(rs+r)**3/4 +(np.log(1+x)-x/(1+x)))/fac200*M200    

def M_Her(r,r200,rs):
    G=4.302e-9
    z=za
    Hz2=H0**2*(0.3*(1+z)**3+0.7)
    M200=100*Hz2/G*r200**3
    return M200*r**2/(r+rs)**2*(r200+rs)*(r200+rs)/(r200*r200)

def M_Bur(r,r200,rs):
    G=4.302e-9
    z=za
    Hz2=H0**2*(0.3*(1+z)**3+0.7)
    M200=100*Hz2/G*r200**3
    trs=r/rs
    rvrs=r200/rs
    fac200=1.0/(np.log(1+rvrs*rvrs)+2.*np.log(1.0+rvrs)-2.*np.arctan(rvrs))
    return M200*fac200*(np.log(1+trs*trs)+2.*np.log(1.0+trs)-2.*np.arctan(trs))

def M_gNFW(r,r200,rs,gamma):
    G=4.302e-9
    z=za
    Hz2=H0**2*(0.3*(1+z)**3+0.7)
    M200=100*Hz2/G*r200**3
    h2x=spec.hyp2f1(3-gamma,3-gamma,4-gamma,-r200/rs)
    h2y=spec.hyp2f1(3-gamma,3-gamma,4-gamma,-r/rs)
    return (r/r200)**(3-gamma)*M200*h2y/h2x
    
def M_nfw_CS(r,r200,rs,tmass,Q):
    c200=r200/rs
    x=r/rs
    G=4.302e-9
    clight=3e5
    z=za
    Hz2=H0**2*(0.3*(1+z)**3+0.7)
    M200=100*Hz2/G*r200**3
    fac200=np.log(1+c200)-c200/(1+c200)  
    phinf=tmass*1e-5
    B=Q*200*c200**3/fac200*Hz2*rs**2/clight**2
    xczero=(B/phinf-1)
    Czero=-B*np.log(1+xczero)+phinf*xczero
    if (xczero<=0.001):
        xczero=0.0   
        Czero=0.0             
    if (x<=xczero):
        dphidr=0.0
    else: 
        cdav=rs*Q*clight**2/(100*Hz2*r200**3)
        dphidr=cdav*(Czero-B*(x/(1.0+x)-np.log(1.0+x)))

        #if the cluster is totally screened, the chameleon contribution is
        #seattled to zero
        
    if (xczero*rs>1.2*rupin): 
        dphidr=0
    #DA FINIRE
    return ((np.log(1+x)-x/(1+x))/fac200+dphidr)*M200

M_nfw_CS=np.vectorize(M_nfw_CS)
#dist=(r200-r2min)**2+(rs-rsmin)**2+(beta-betamin)**2+(mod1-mod1min)**2+(mod1-mod1min)**2
#prova=np.column_stack((dist,r200,rs,beta,mod1,mod2))
#prova=prova[prova[:, 0].argsort()]
#
#r268=prova[0:int(0.68*len(prova)),1]
#rs68=prova[0:int(0.68*len(prova)),2]
#r295=prova[0:int(0.95*len(prova)),1]
#rs95=prova[0:int(0.95*len(prova)),2]
#mod168=prova[0:int(0.68*len(prova)),4]
#mod268=prova[0:int(0.68*len(prova)),5]
#mod195=prova[0:int(0.95*len(prova)),4]
#mod295=prova[0:int(0.95*len(prova)),5]
Mup=[]
Mlow=[]
Mup95=[]
Mlow95=[]
Mmedian=[] 
r=np.linspace(0.05,2.5,50)
median_string=True
plot=True
for ri in r:
    profile=[]
    profile95=[]
#for CS we use the median profile
    for i in range(0,len(r200)):
        if (mass_string=='NFW'):
            profile.append(M_nfw(ri,r200[i],rs[i]))
        elif (mass_string=='Her'):
            profile.append(M_Her(ri,r200[i],rs[i]))
        elif (mass_string=='Bur'):
            profile.append(M_Bur(ri,r200[i],rs[i]))            
        elif (mass_string=='mNFW_BH'):
            profile.append(M_nfw_BH(ri,r200[i],rs[i],mod1[i]))
        elif (mass_string=='gNFW'):
            profile.append(M_gNFW(ri,r200[i],rs[i],mod1[i]))
        elif (mass_string=='mNFW_GC'):
            profile.append(M_nfw_CS(ri,r200[i],rs[i],phicoup[i],Qcoup[i]))
        else:
            profile.append(0.0)
    profile=np.sort(np.array(profile).astype(float))
    Mup.append(profile[int(len(profile)*0.86)]) #0.975
    Mlow.append(profile[int(len(profile)*0.16)]) #0.025
    Mup95.append(profile[int(len(profile)*0.975)]) #0.975
    Mlow95.append(profile[int(len(profile)*0.025)]) #0.025    
    Mmedian.append(profile[int(len(profile)*0.5)])
if (mass_string=='NFW'):
    Mnfw=M_nfw(r,r2min[0],rsmin[0]) 
elif (mass_string=='Her'):
    Mnfw=M_Her(r,r2min[0],rsmin[0]) 
elif (mass_string=='Bur'):
    Mnfw=M_Bur(r,r2min[0],rsmin[0])     
elif (mass_string=='mNFW_BH'):
    Mnfw=M_nfw_BH(r,r2min[0],rsmin[0],mod1min[0])
elif (mass_string=='gNFW'):
    Mnfw=M_gNFW(r,r2min[0],rsmin[0],mod1min[0])
elif (mass_string=='mNFW_GC'):
    M_nfw_CS(ri,r2min[0],rsmin[0],mod1min[0],mod2min[0])
    #Mnfw=M_nfw_CS(r,r2min[0],rsmin[0],mod1min[0],mod2min[0])
else:
    plot=False
    Mnfw=np.zeros(len(r))
if(plot):   
    if median_string==True:    
        Mnfw=Mmedian

    mpl.figure(figsize=(12,12))
    mpl.grid()
    mpl.tick_params(axis="x", labelsize=22)
    mpl.tick_params(axis="y", labelsize=22)
    mpl.fill_between(r, Mup95, Mlow95, alpha=.2)
    mpl.fill_between(r, Mup, Mlow, alpha=.2)
    mpl.xlim(0.05,1.2*rupin)
    mpl.ylim(1e12,1.5*np.amax(Mup95))
    mpl.axvline(r2min[0], 0, 3.0,c='black',ls='--')
    mpl.xlabel('r [Mpc]',size=24)
    mpl.ylabel(r'M [$M_{sun}$]',size=24)
#    if (mass_string=='mNFW_GC'):
#        mpl.semilogy(r,Mnfw)
#    else:
#        mpl.ylim(1e12,3e15)
#        mpl.plot(r,Mnfw)
    mpl.plot(r,Mnfw)
    mpl.show()    
    
    
