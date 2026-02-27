# -*- coding: utf-8 -*-
"""
Created on Sun Nov  6 21:38:25 2022

@author: Lorenzo
"""

import numpy as np
import scipy.special as spec
import scipy.integrate as ints
import csv
import scipy.optimize as opt
from scipy.special import gamma, gammaincc, gammainccinv, betainc
import sys
from scipy.special import binom 
import datetime
import scipy.interpolate as ipt
import re
from dataclasses import dataclass
from typing import Callable, Any
import numpy as np


# Parsing functions 
# prende: key = value  (key può contenere parentesi e spazi; value = primo token non-spazio)
# supporta: key=value, key =value, key= value, key = value
# ignora il resto della riga (commenti ! e testo dopo)
_KV_RE = re.compile(r'^\s*(?P<key>[^=!\n]+?)\s*=\s*(?P<val>[^\s!]+)')

def ffloat(s: str) -> float:
    """float robusto anche per notazione Fortran D/d (es: 1d-3)."""
    return float(s.replace("D", "e").replace("d", "e"))

@dataclass(frozen=True)
class Spec:
    cast: Callable[[str], Any]
    default: Any
    post: Callable[[Any], Any] = lambda x: x  # scaling/transform opzionale

def read_params_block(path: str, mod1=None, mod2=None, r200=None) -> dict[str, Any]:
    # --- defaults (replico i tuoi default) ---
    specs = {
        "za":      Spec(ffloat, 0.0),
        "H0":      Spec(ffloat, 70.0),
        "Olam":    Spec(ffloat, 0.7),
        "Omegam":  Spec(ffloat, 0.3),

        "BCG":     Spec(int,    0),
        "Baryons": Spec(int, 0),

        # scaling 
        "xlumbcg": Spec(ffloat, 1e12,  post=lambda x: x * 1e11),
        "rho0gas": Spec(ffloat, 10e13, post=lambda x: x * 1e13),
        "rsgas":   Spec(ffloat, 0.38),

        "rhostar": Spec(ffloat, 4e13,  post=lambda x: x * 1e13),
        "rstar":   Spec(ffloat, 0.368),

        "srhost":  Spec(ffloat, 0.1e13, post=lambda x: x * 1e13),
        "srstar":  Spec(ffloat, 0.1),

        "srhoga":  Spec(ffloat, 0.5e13, post=lambda x: x * 1e13),
        "srsga":   Spec(ffloat, 0.1),

        "eim":     Spec(ffloat, 1.0),

        "M(r)":    Spec(str,    "NFW"),
        "Beta(r)": Spec(str,    "None"),
        "FreeBeta":Spec(int,    0),

        "nA2":     Spec(int,    0),

        "nmcmc":   Spec(int,   1),
        # bounds: default depedndent from mod1/mod2/r200 
        "A1low":   Spec(ffloat, None),
        "A1up":    Spec(ffloat, None),
        "A2low":   Spec(ffloat, None),
        "A2up":    Spec(ffloat, None),
        "Rup":     Spec(ffloat, None),
    }

    # init with defaults
    out = {k: s.default for k, s in specs.items()}

    with open(path, "r", encoding="utf-8") as f:
        for raw in f:
            # commenti stile Fortran: "!"  (prima li taglio)
            line = raw.split("!", 1)[0].strip()

            # ignora righe vuote e separatori tipo *****
            if not line or set(line) <= {"*"}:
                continue

            m = _KV_RE.match(line)
            if not m:
                continue

            key = " ".join(m.group("key").split())  # normalizza spazi nella chiave
            val = m.group("val").strip()

            if key not in specs:
                continue

            spec = specs[key]
            try:
                out[key] = spec.post(spec.cast(val))
            except ValueError:
                print(f"parameter {key} has undefined value. Switching to default")
                out[key] = spec.default

    # --- fallback per i default dipendenti (replica i tuoi np.amin/amax) ---
    if out["A1low"] is None:
        out["A1low"] = float(np.amin(mod1)) if mod1 is not None else None
    if out["A1up"] is None:
        out["A1up"]  = float(np.amax(mod1)) if mod1 is not None else None
    if out["A2low"] is None:
        out["A2low"] = float(np.amin(mod2)) if mod2 is not None else None
    if out["A2up"] is None:
        out["A2up"]  = float(np.amax(mod2)) if mod2 is not None else None
    if out["Rup"] is None:
        out["Rup"]   = float(np.amax(r200)) if r200 is not None else None

    return out




#Mass concentration relation 
#Ragagnin formula


#Galassie MACS 1206
rhostar = 3.547e+13
rstar = 0.369


rhogas = 2.57018005e+14 
rsgas = 0.38


def format_uncertainty(val, plus, minus, ndigits=2):
    fmt = f"{{:.{ndigits}f}}"
    return rf"${fmt.format(val)}^{{+{fmt.format(plus)}}}_{{-{fmt.format(minus)}}}$"


# Function to format a single uncertainty number for LaTeX
#def format_uncertainty(number):
#    return "${:.2f}^".format(number[0])+"{+"+"{:.2f} ".format(number[1])+"}_{-"+"{:.2f}".format(number[2])+"}$"


#Functions for the general screening Treatment
def fscrtot(r,Q,phinf,rhos,rs,rhostar,rstar,rhog,rg,rhob,gam,rj):
     """
    Auxiliary function for the computation of the screening radius assuming 
    four mass components-modeling in Chameleon Gravity. Dark Matter is modeled 
    as a gNFW, galaxies as a NFW, gas as a beta profile (beta = 1), BCG as a Jaffe
    profile.

    Parameters
    ----------
    r : FLOAT
        Radius in Mpc.
    Q : FLOAT
        Coupling constant 
    phinf : FLOAT
        Value of the field at infinity.
    rhos : FLOAT
        dark matter profile critical density.
    rs : FLOAT
        dark matter profile scale radius (in Mpc).
    rhostar : FLOAT
        galaxies profile critical density.
    rstar : FLOAT
        galaxies profile scale radius (in Mpc).
    rhog : FLOAT
        DESCRIPTION.
    rg : FLOAT
        gas profile scale radius (in Mpc).
    rhob : FLOAT
        BCG profile critical density
    gam : FLOAT
        Slope parameter of the dark matter profile.
    rj : FLOAT
       BCG profile scale radius (in Mpc).

    Returns
    -------
    The value of the function which zero is the screening radius.

    """
     G=4.302e-9
     clight=3e5
     Facda = Q*(8*np.pi*G)/clight**2
#     h2x=spec.hyp2f1(1,1,4-gam,-r/rs)
    
    
     #Mpl2 = 1/(8*np.pi*G)
     
     
     return phinf + Facda * (
         ((rj**3 *rhob)/(r + rj) - (rg**3 *rhog)/np.sqrt(r**2 + rg**2) + (
    rs**2*(1 - (r/(r + rs))**(
       2 - gam)) *rhos)/(-2 + gam) - (
    rstar**3 * rhostar)/(r + rstar)) + (rj**2 *rhob * np.log(r/(r + rj)) ) )

           
           
           
def ScreenTot(Q,phinf,r200,rs,rhostar,rstar,rhog,rg,rhob,gam,rj,
             z = 0.0, H0 = 70,  Om = 0.3, Ol = 0.7):
    """
    This function computes the screening radius of a  
    four mass components-modeled cluster in Chameleon Gravity. Dark Matter is modeled 
    as a gNFW, galaxies as a NFW, gas as a beta profile (beta = 1), BCG as a Jaffe
    profile.

    Parameters
    ----------
    Q : FLOAT
        Coupling constant 
    phinf : FLOAT
        Value of the field at infinity.
    rhos : FLOAT
        dark matter profile r200 (in Mpc).
    rs : FLOAT
        dark matter profile scale radius (in Mpc).
    rhostar : FLOAT
        galaxies profile critical density.
    rstar : FLOAT
        galaxies profile scale radius (in Mpc).
    rhog : FLOAT
        DESCRIPTION.
    rg : FLOAT
        gas profile scale radius (in Mpc).
    rhob : FLOAT
        BCG profile critical density
    gam : FLOAT
        Slope parameter of the dark matter profile.
    rj : FLOAT
       BCG profile scale radius (in Mpc).
    z : FLOAT (optional)
        redhsift of the cluster. Default is 0
    H0 : FLOAT (optional)
        Value of the Hubble parameter today. Default is 70
    Om : FLOAT (optional)
        Value of Omega matter. Default is 0.3
    Ol : FLOAT (optional)
        Value of Omega Lambda. Default is 0.7

    Returns
    -------
    The value of the screening radius.

    """
    rhos = rhos_gNFW(r200, rs, gam, z, H0,  Om, Ol)
    def fgen2(y):
        return fscrtot(y,Q,phinf,rhos,rs,rhostar,rstar,rhog,rg,rhob,gam,rj)
    limit=1e-6
    xsc=limit
    x0=limit   
    
    if(fgen2(x0)<0):
        for i in range(0,4000):
            xx=x0+0.02*i
            
            
            if (fgen2(x0)*fgen2(xx) < 0) :
                sol=opt.root(fgen2,xx)
            
                if sol.success:
                    xsc=sol.x[0]
                else:
                    xsc=xx                
                    break
                x0=xx
        if (xsc==limit):
              xsc=xx
    return xsc 

ScreenTot=np.vectorize(ScreenTot)            
           
           
#Function for the cumulative screening           
def dphiTot(r,Q,phinf,rhos,rs,rhostar,rstar,rhog,rg,
            rhob,gam,rj,xsc = 1e-6, limit = 1e-6):
    
    """
    Computes the value of the fifth force in Chameleon screening assuming 
    a multicomponent mass density distribution of gas, galaxies, DM and BCG.
    The coupling is considered the same for all the components

    Parameters
    ----------
    r : FLOAT
        radius in Mpc.
    Q : FLOAT
        Coupling constant 
    phinf : FLOAT
        Value of the field at infinity (in unit of c^2).
    rhos : FLOAT
        dark matter profile critical density.
    rs : FLOAT
        dark matter profile scale radius (in Mpc).
    rhostar : FLOAT
        galaxies profile critical density.
    rstar : FLOAT
        galaxies profile scale radius (in Mpc).
    rhog : FLOAT
        DESCRIPTION.
    rg : FLOAT
        gas profile scale radius (in Mpc).
    rhob : FLOAT
        BCG profile critical density
    gam : FLOAT
        Slope parameter of the dark matter profile.
    rj : FLOAT
       BCG profile scale radius (in Mpc).
    xsc : FLOAT, optional
        The value of the screening radius. The default is 1e-6.
    limit : FLOAT, optional
        Limiting value for the screening to be considered zero. The default is 1e-6.

    Returns
    -------
    Value of the effective mass at r

    """
    G=4.302e-9
    clight=3e5
    Facda = Q*(8*np.pi*G)/clight**2    
    #limit=1e-6 #limiting value for the screening radius
    h2x= spec.hyp2f1(3-gam,3-gam,4-gam,-r/rs)
    
    h2sc= spec.hyp2f1(3-gam,3-gam,4-gam,-xsc/rs)
    
    Ctot = Facda * ((rj**4 *rhob)/(
    xsc + rj) + ( xsc**3 * (xsc/rs)**(-gam) *rhos * h2sc)/(-3 + gam) + 
    rg**3 *rhog *(xsc/np.sqrt(xsc**2 + rg**2) + np.log(-xsc + np.sqrt(xsc**2 + rg**2))) - 
    rstar**3 *rhostar * (rstar/(xsc + rstar) + np.log(xsc + rstar)))
    

    if (xsc>limit):
        if (r<xsc):
            out=0
        else:
            out= Ctot + Facda* (( r**(3 - gam) * (1 + r/rs)**(-gam) * (r + 
    rs)**gam * rhos* h2x)/( (3 - gam)) + (-((rj**4 *rhob)/(r + rj)) + 
    rg**3 *rhog * (-(r/np.sqrt(r**2 + rg**2) ) - np.log(-r + np.sqrt(r**2 + rg**2))) +
    rstar**3 *rhostar * (rstar/(r + rstar) + np.log(r + rstar) )) ) 
                    
            
    else:

        Ctot= Facda*( (2 *rj**3 *rhob + rg**3 *rhog * np.log(rg**2) - 
   2* rstar**3 * rhostar *(1 + np.log(rstar) )))/(2)
        
        out= Facda* (( r**(3 - gam) * (1 + r/rs)**(-gam) * (r + 
rs)**gam * rhos* h2x)/( (3 - gam)) + ( (-((rj**4 *rhob)/(r + rj)) + 
rg**3 *rhog * (-(r/np.sqrt(r**2 + rg**2) ) - np.log(-r + np.sqrt(r**2 + rg**2))) +
rstar**3 *rhostar * (rstar/(r + rstar) + np.log(r + rstar) ))) ) +Ctot
           
        
    return out*Q*clight**2/G #
dphiTot=np.vectorize(dphiTot,otypes=[float],excluded=[1,2,3,4])
     

def rhos_gNFW(r200, rs, gam, z = 0, H0 = 70,  Om = 0.3, Ol = 0.7):
    """
    Computes the critical density rhos of a gNFW profile.

    Parameters
    ----------
    r200 : FlOAT
        Value of the radius r200 in Mpc.
    rs : FLOAT
        Value of the radius rs in Mpc.
    gam : FLOAT
        Slope parameter.
    z : FLOAT (optional)
        redhsift of the cluster. Default is 0
    H0 : FLOAT (optional)
        Value of the Hubble parameter today. Default is 70
    Om : FLOAT (optional)
        Value of Omega matter. Default is 0.3
    Ol : FLOAT (optional)
        Value of Omega Lambda. Default is 0.7
                                        

    Returns
    -------
    The values of rhos in M_sun/Mpc**3

    """
    G=4.302e-9
    #z=za
    Hz2=H0**2*(Om*(1+z)**3+Ol)
    M200=100*Hz2/G*r200**3
    
    h2x= spec.hyp2f1(3-gam,3-gam,4-gam,-r200/rs)

    return (r200/rs)**(gam-3)*M200 / (4* np.pi * rs**3 * h2x) * ( 3 - gam )
    


def rhotot(r,rhos,rs,rhostar,rstar,rhog,rg,rhob,gam,rj):
    rhogass = rhog/( 1 + (r/rg)**2)**(3/2)
    rhogal = rhostar/( r/rstar * (1 + r/rstar)**2)
    rhoDM = rhos/((r/rs)**gam *(1+(r/rs))**(3-gam))
    rhobcg = rhob/((r/rj)**2  * (1+ r/rj)**2)    
    
    return rhogass + rhogal + rhoDM + rhobcg



def M_tot(x,r200, rs, gamma, xml, xmstar, rjaf, H0, z, Om, Ol, rsgas = 0.16, rhogas = 60.80e13):
    
     xmbar= 4*np.pi*rhogas*rsgas**2*(x-rsgas*np.arctan(x/rsgas))
    
     xmdm = M_gNFW(x,r200,rs,gamma,z,H0,Om,Ol)
     return xmstar*xml*x/rjaf/(1.+x/rjaf)+ xmbar + xmdm
    

def Mtot(r,Q,phinf,r200,rs,rhostar,rstar,rhog,rg,rhob,gam,rj,xsc = 1e-6,
        z = 0, H0 = 70,  Om = 0.3, Ol = 0.7):
    """
    
    Mass profile in Chameleon gravity of a  
    four mass components-modeled cluster in Chameleon Gravity. Dark Matter is modeled 
    as a gNFW, galaxies as a NFW, gas as a beta profile (beta = 1), BCG as a Jaffe
    profile.

    Parameters
    ----------
    r : FLOAT
        value of the radius in Mpc.
    Q : FLOAT
        Coupling constant 
    phinf : FLOAT
        Value of the field at infinity (in unit of c**2).
    rhos : FLOAT
        dark matter profile critical density.
    rs : FLOAT
        dark matter profile scale radius (in Mpc).
    rhostar : FLOAT
        galaxies profile critical density.
    rstar : FLOAT
        galaxies profile scale radius (in Mpc).
    rhog : FLOAT
        DESCRIPTION.
    rg : FLOAT
        gas profile scale radius (in Mpc).
    rhob : FLOAT
        BCG profile critical density
    gam : FLOAT
        Slope parameter of the dark matter profile.
    rj : FLOAT
       BCG profile scale radius (in Mpc).
    xsc : FLOAT, optional
        Value of the screening radius. The default is 1e-6.
    z : FLOAT (optional)
        redhsift of the cluster. Default is 0
    H0 : FLOAT (optional)
        Value of the Hubble parameter today. Default is 70
    Om : FLOAT (optional)
        Value of Omega matter. Default is 0.3
    Ol : FLOAT (optional)
        Value of Omega Lambda. Default is 0.7
                                        

    Returns
    -------
    FLOAT
        Value of the total effective mass in Chameleon gravity.

    """
    
    #compute the value of rhos given r200, rs
    rhos = rhos_gNFW(r200, rs, gam, z= z, H0 = H0,  Om = Om, Ol = Ol)
    
    
    xdm = r/rs
    xg = r/rg
    xstar = r/rstar
    xj = r/rj
    h2x= spec.hyp2f1(3-gam,3-gam,4-gam,-r/rs)
    
    mass = 4* np.pi * r**3 * (rhob/(xj**2 + xj**3) - rhog/(
   xg**2 * np.sqrt(1 + xg**2))  - rhostar/(
   xstar**2 * (1 + xstar)) + (rhog *np.arcsinh(xg))/
   xg**3 + ((xdm)**(
    3 - gam) * rhos *  h2x)/(
   xdm**3 * (3 - gam)) + (rhostar * np.log ( 1 + xstar) )/xstar**3)
    if (Q!=0 and phinf!=0):  
         dphidr=  dphiTot(r,Q,phinf,rhos,rs,rhostar,rstar,rhog,rg,
                          rhob,gam,rj,xsc = xsc)
    else:
         dphidr=0.0
    return mass+dphidr    

Mtot = np.vectorize(Mtot)  


def phiTot(r,Q,phinf,rhos,rs,rhostar,rstar,rhog,rg,rhob,gam,rj,xsc = 1e-6):
    """
    Chameleon field profile of a  
    four mass components-modeled cluster in Chameleon Gravity. Dark Matter is modeled 
    as a gNFW, galaxies as a NFW, gas as a beta profile (beta = 1), BCG as a Jaffe
    profile.

    Parameters
    ----------
    r : FLOAT
        value of the radius in Mpc.
    Q : FLOAT
        Coupling constant 
    phinf : FLOAT
        Value of the field at infinity (in unit of c**2).
    rhos : FLOAT
        dark matter profile critical density.
    rs : FLOAT
        dark matter profile scale radius (in Mpc).
    rhostar : FLOAT
        galaxies profile critical density.
    rstar : FLOAT
        galaxies profile scale radius (in Mpc).
    rhog : FLOAT
        DESCRIPTION.
    rg : FLOAT
        gas profile scale radius (in Mpc).
    rhob : FLOAT
        BCG profile critical density
    gam : FLOAT
        Slope parameter of the dark matter profile.
    rj : FLOAT
       BCG profile scale radius (in Mpc).
    xsc : FLOAT, optional
        Value of the screening radius. The default is 1e-6.

    Returns
    -------
    FLOAT
        Value of the field at r.

    """
    G=4.302e-9
    clight=3e5
    Facda = Q*(8*np.pi*G)/clight**2    
    limit=1e-6 #limiting value for the screening radius
    h2x= spec.hyp2f1(3-gam,3-gam,4-gam,-r/rs)
    
    h2sc= spec.hyp2f1(3-gam,3-gam,4-gam,-xsc/rs)
    
    Ctot = Facda * ((rj**4 *rhob)/(
    xsc + rj) + ( xsc**3 * (xsc/rs)**(-gam) *rhos * h2sc)/(-3 + gam) + 
    rg**3 *rhog *(xsc/np.sqrt(xsc**2 + rg**2) + np.log(-xsc + np.sqrt(xsc**2 + rg**2))) - 
    rstar**3 *rhostar * (rstar/(xsc + rstar) + np.log(xsc + rstar)))
    
    
    if (xsc>limit):
        if (r<xsc):
            out=0
        else:
            out = phinf - Ctot/r + Facda* ( rs**2 *rhos)/( (-2 +gam)) + (Facda  / r)  * (rj**3 * rhob - 
    (r**(3 - gam) * rs**2 * (r + rs)**(-2 +gam ) * rhos)/(-2 + gam) - 
    rstar**3 * rhostar - ( (r/rs)**(3 - gam) * rs**3 *rhos * h2x)/(3 - gam) + 
    rg**3 *rhog * np.log(-r + np.sqrt(r**2 + rg**2) ) + 
    r * rj**2 * rhob *np.log(r/(r + rj)) - rstar**3 * rhostar * np.log(r + rstar) )
                    
            
    else:    
        Ctot= Facda*( (2 *rj**3 *rhob + rg**3 *rhog * np.log(rg**2) - 
   2* rstar**3 * rhostar *(1 + np.log(rstar) )))/(2)
        
        out = phinf - Ctot/r + Facda* ( rs**2 *rhos)/( (-2 +gam)) + (Facda  / r)  * (rj**3 * rhob - 
    (r**(3 - gam) * rs**2 * (r + rs)**(-2 +gam ) * rhos)/(-2 + gam) - 
    rstar**3 * rhostar + ( (r/rs)**(3 - gam) * rs**3 *rhos * h2x)/(3 - gam) + 
    rg**3 *rhog * np.log(-r + np.sqrt(r**2 + rg**2) ) + 
    r * rj**2 * rhob *np.log(r/(r + rj)) - rstar**3 * rhostar * np.log(r + rstar) )


    return out

phiTot=np.vectorize(phiTot,otypes=[float],excluded=[1,2,3,4])

def Cm_Rag(M200,h0,zred,om,omb,sig8=0.809):

    #sig=0.384
    zp=1./0.877-1
    mp=17.4e13
    
    a0=np.exp(1.24)
    b0=-0.05
    c0=0.20
    am=0.632 
    ab=-0.246  
    asw=0.561
    ah=-0.026
    bm=-0.118
    bb=0.112 
    bs=0.056
    bh=-0.044
    cm=0.352
    cb=-0.039
    cs=0.767 
    ch=-0.276
    omp=0.272 
    ombp=0.0465 
    sig8p=0.809
    h0p=70.4
    #omb=ombp
    #sig8=sig8p
    
    #omb=0.049
    #h0=h
    
    a=np.exp(np.log(a0)+am*np.log(om/omp)+ab*np.log(omb/ombp)+asw*np.log(sig8/sig8p)+ah*np.log(h0/h0p))
    b=b0+bm*np.log(om/omp)+bb*np.log(omb/ombp)+bs*np.log(sig8/sig8p)+bh*np.log(h0/h0p)
    c=c0+cm*np.log(om/omp)+cb*np.log(omb/ombp)+cs*np.log(sig8/sig8p)+ch*np.log(h0/h0p)
    
    return np.exp(np.log(a)+b*np.log(M200/mp)+c*np.log((1+zp)/(1+zred)))#+sig) the scattering
# is not needed in the computation of the concentration


def C_Merten(M200,H0,Om,z):
    A=np.random.normal(3.66, 0.16)
    B=np.random.normal(-0.14,0.52)
    C=np.random.normal(-0.32,0.18) 
        
    return A/(1.37/(1+z))**B*(M200/(8e14/(0.01*H0)))**C



def FindM200(r200, z, H0=70, Om = 0.3, Ol = 0.7):
    """
    Find M200 given r200 and a cosmology

    Parameters
    ----------
    r200 : Float
        r200 of the halo.
    z : Float
        redshift.
    H0 : Float, optional
        Hubble constant at z=0. The default is 70.
    Om : Float, optional
        Omega matter. The default is 0.3.
    Ol : Float, optional
        Omega Lambda. The default is 0.7.

    Returns
    -------
    M200 : Float
        Value of the mass.

    """
    G=4.302e-9
    #z=za
    Hz2=H0**2*(Om*(1+z)**3+Ol)
    M200=100*Hz2/G*r200**3
    return M200

#Define and vectorize the mass functions
#******************************************************************************


def M_nfw(r,r200,rs,z,H0,Om,Ol):
    c200=r200/rs
    x=r/rs
    G=4.302e-9
    #z=za
    Hz2=H0**2*(Om*(1+z)**3+Ol)
    M200=100*Hz2/G*r200**3
    fac200=np.log(1+c200)-c200/(1+c200)
    return (np.log(1+x)-x/(1+x))/fac200*M200
M_nfw=np.vectorize(M_nfw)


def M_Bur(r,r200,rs,z,H0,Om,Ol):
    G=4.302e-9
    #z=za
    Hz2=H0**2*(Om*(1+z)**3+Ol)
    M200=100*Hz2/G*r200**3
    trs=r/rs
    rvrs=r200/rs
    fac200=1.0/(np.log(1+rvrs*rvrs)+2.*np.log(1.0+rvrs)-2.*np.arctan(rvrs))
    return M200*fac200*(np.log(1+trs*trs)+2.*np.log(1.0+trs)-2.*np.arctan(trs))
M_Bur=np.vectorize(M_Bur)


def M_Her(r,r200,rs,z,H0,Om,Ol):
    G=4.302e-9
    #z=za
    Hz2=H0**2*(Om*(1+z)**3+Ol)
    M200=100*Hz2/G*r200**3
    return M200*r**2/(r+rs)**2*(r200+rs)*(r200+rs)/(r200*r200)
M_Her= np.vectorize(M_Her)


def M_Ein(r,r200,rs,n,z,H0,Om,Ol):
    G=4.302e-9
    #z=za
    Hz2=H0**2*(Om*(1+z)**3+Ol)
    M200=100*Hz2/G*r200**3 
    c200=r200/rs
    x=r/rs
    fac200=(1-gammaincc(3*n,(2*n*(c200)**(1/n))))*gamma(3*n)   
    return M200/fac200*(1-gammaincc(3*n,(2*n*(x)**(1/n))))*gamma(3*n) # 

def M_gNFW(r,r200,rs,gamma,z,H0,Om,Ol):
    G=4.302e-9
    #z=za
    Hz2=H0**2*(Om*(1+z)**3+Ol)
    M200=100*Hz2/G*r200**3
    h2x=spec.hyp2f1(3-gamma,3-gamma,4-gamma,-r200/rs)
    h2y=spec.hyp2f1(3-gamma,3-gamma,4-gamma,-r/rs)
    return (r/r200)**(3-gamma)*M200*h2y/h2x
M_gNFW=np.vectorize(M_gNFW)


#Jaffe profile (for the BCG)
def jaffe_BCG(r,xml,mxs,rj):
    return mxs*xml*r/rj/(1.+r/rj)


jaffe_BCG=np.vectorize(jaffe_BCG)



# ---------- helpers for baryons ----------
def Mstar_nfw_like(r, rs, rhos):
    # your analytic stellar cumulative mass
    return 4.0 * np.pi * rs**3 * rhos * (np.log(1.0 + r/rs) - (r/rs)/(1.0 + r/rs))

def Mgas_beta_like(r, rs, rhos):
    # your analytic gas cumulative mass 
    return rs**2 * 4.0 * np.pi * rhos * (r - rs*np.arctan(r/rs))



def M_nfw_BH(r,r200,rs,Y1,z,H0,Om,Ol):
    c200=r200/rs
    x=r/rs
    G=4.302e-9
    #z=za
    Hz2=H0**2*(Om*(1+z)**3+Ol)
    M200=100*Hz2/G*r200**3
    fac200=np.log(1+c200)-c200/(1+c200)
    return (Y1*r**2*(rs-r)/(rs+r)**3/4 +(np.log(1+x)-x/(1+x)))/fac200*M200
M_nfw_BH=np.vectorize(M_nfw_BH)


def M_Bur_BH(r,r200,rs,Y1,z,H0,Om,Ol):
    G=4.302e-9
    #z=za
    Hz2=H0**2*(Om*(1+z)**3+Ol)
    M200=100*Hz2/G*r200**3
    trs=r/rs
    rvrs=r200/rs
    fac200=1.0/(np.log(1+rvrs*rvrs)+2.*np.log(1.0+rvrs)-2.*np.arctan(rvrs))
    
    
    return fac200*M200*((r**3*(-r**3 + r*rs**2 + 2*rs**3)*Y1)/((r + rs)**2*(r**2 + rs**2)**2) 
                        -2*np.arctan(r/rs) + np.log(1 + r**2/rs**2) + 2*np.log(1+trs))
    
M_Bur_BH=np.vectorize(M_Bur_BH)   


def M_bnfw_BH(r,r200,rs,Y1,b,z,H0,Om,Ol):
     c200=r200/rs
     x=r/rs
     G=4.302e-9
    
     Hz2=H0**2*(Om*(1+z)**3+Ol)
     M200=100*Hz2/G*r200**3             
     if (b!=1 and b!=2):
         f200=(c200+1)**(-b)*((c200+1)**b-b*c200*(c200+1)+c200**2-1)
         
         return (M200*(1 + x)**(-b)*(4*rs**2*((r + rs)/rs)**b - 
                4*(r + rs)*((-1 + b)*r + rs) + r**2*(2 - 2*b - 
           (1 - b)*b*((r + rs)/rs)**(-2 + b))*Y1))/(4*f200*rs**2)
     elif b==2:         
         return M_nfw_BH(r,r200,rs,Y1,z,H0,Om,Ol)         
     else:
         return 0
     
M_bnfw_BH=np.vectorize(M_bnfw_BH)



def M_Ein_BH(r,r200,rs,Y1,n,z,H0,Om,Ol):
    
    G=4.302e-9

    Hz2=H0**2*(Om*(1+z)**3+Ol)
    M200=100*Hz2/G*r200**3 
    c200=r200/rs
    fac200=(1-gammaincc(3*n,(2*n*(c200)**(1/n))))*gamma(3*n) 
    fmod=-(8**n * np.exp(-2*n*(r/rs)**(1/n))*
           (-1 +(r/rs)**(1/n))*(n*(r/rs)**(1/n))**(3*n)*Y1)/n
  
    return M200*(fmod+gamma(3*n)*(2-2*gammaincc(3*n, 2*n*(r/rs)**(1/n)) ) )/(2*fac200)


M_Ein_BH=np.vectorize(M_Ein_BH)
#******************************************************************************



def M_nfw_CS(r,r200,rs,tmass,Q,z,H0,Om,Ol):
    c200=r200/rs
    x=r/rs
    G=4.302e-9
    clight=3e5
    #z=za
    Hz2=H0**2*(Om*(1+z)**3+Ol)
    M200=100*Hz2/G*r200**3
    fac200=np.log(1+c200)-c200/(1+c200)  
    phinf=tmass*1e-5
    B=Q*200*c200**3/fac200*Hz2*rs**2/clight**2
    xczero=(B/phinf-1)

    Czero=-B*np.log(1+xczero)+phinf*xczero
    
    if (xczero<=0.001):
        #for very small screening radius, the integration constant should be different
        xczero=0.0   
        Czero=-B              
    if (x<=xczero):
        dphidr=0.0
    else: 
        cdav=rs*Q*clight**2/(100*Hz2*r200**3)
        dphidr=cdav*(Czero-B*(x/(1.0+x)-np.log(1.0+x)))

        #if the cluster is totally screened, the chameleon contribution is
        #seattled to zero
    rupin=r200 #define a cutoff for the fifth force   
    if (xczero*rs>1.2*rupin): 
        dphidr=0
    return ((np.log(1+x)-x/(1+x))/fac200+dphidr)*M200

M_nfw_CS=np.vectorize(M_nfw_CS)

#=======================================================================
#Generalized NFW profile 

def dphigNFW(x,rs,rhos,gam=1.0,Q=0.5,phinf=1e-4,expcutoff=False):
    """
    Gradient of the chameleon field assuming a gNFW density profile    
    
    Parameters
    ----------
    x : Float
        radius, in unit of r/rs.
    rs : Float
        scale radius of the profile.
    rhos : Float
        Characteristic density of the profile.
    gam : Float
        Inner slope exponent of the density profile.
    Q : Float
        Coupling constant.
    phinf : Float
        Value of the field at infinty.
    expcutoff : BOOL, optional
        Add an exponential cutoff to the effective mass?. The default is True.

    Returns
    -------
    Float
        The value of the gradient of the field at x.

    """
    clight= 2.99792458e5
    limit=1e-6 #limiting value for the screening radius
    G=4.302e-9
    
    
    #eponential cutoff:
    dminf=np.exp(-x*rs*np.sqrt(1e-7*Q/phinf))
    if expcutoff==False:
        dminf=1.0
        
    h2x=spec.hyp2f1(3-gam,3-gam,4-gam,-x)  

    Facda=Q*rs**2*rhos*8*np.pi*G/clight**2 #factor containing the planck mass     
    #In this case, analytic solution exists. No need to find numerical solution
    
    if (1+(phinf/Facda)*(gam - 2) < 0 ):
        xsc = limit
    else:
        xsc = 1/(1 - (1+(phinf/Facda)*(gam - 2) )**(1/(2-gam)) ) -1 
    if(np.isnan(xsc)):
        xsc=limit
    if (xsc<0):
        xsc = limit #1e6 #completely screend
    Clog = 1 #+ np.log(rs)
     #+np.log(rs) 

    if (xsc>limit):
        hxsc = spec.hyp2f1(3-gam,3-gam,4-gam,-xsc) 
        Cgnfw = - Facda  * ( Clog + xsc**(3-gam)/(3 - gam)*hxsc)
        if (x<xsc):
            out=0.0
        else:            
            out= Cgnfw + Facda *( Clog  + x**(3- gam)/ (3-gam)*h2x) #+ np.log(rs)
    else:
        Ccg = - Facda  * (Clog)
        out=Ccg + Facda  * ( Clog + x**(3- gam)/ (3-gam)*h2x)
        
    return out*Q*rs*clight**2/G #*dminf


dphigNFW=np.vectorize(dphigNFW,otypes=[float],excluded=[1,2,3,4,5])




def M_gnfw_CS(r,r200,rs,gamma,tmass,Q,z,H0,Om,Ol,expcutoff = False):
    """
    Returns the gNFW effective mass profile in Chameleon gravity

    Parameters
    ----------
    r : Float
        input radius [Mpc].
    r200 : Float
        Cluster virial radius.
    rs : Float
        scale radius of the mass profile.
    gamma : Float
        exponent of the density profile's slope.
    tmass : Float
        field value at infinity (in unit of 1e-5 c**2).
    Q : Float
        Couoling constant.
    z : Float
        cluster redshift.
    H0 : Float
        Hubble parameter at z = 0.
    Om : Float
        Omega matter.
    Ol : Float
        Omega lambda.
    expcutoff : BOOL, optional
        Add an exponential cutoff to the effective mass?. The default is True.

    Returns
    -------
    Float
        DESCRIPTION.

    """
    x=r/rs
    G=4.302e-9
    phinf=tmass*1e-5
     
    Hz2=H0**2*(Om*(1+z)**3+Ol)
    M200=100*Hz2/G*r200**3
    h2x=spec.hyp2f1(3-gamma,3-gamma,4-gamma,-r200/rs)
     
    rhos=((M200 * r200**(-3 + gamma)* rs**(-gamma)*(3 - gamma))/(4*np.pi*h2x))

    h2x=spec.hyp2f1(3-gamma,3-gamma,4-gamma,-r200/rs)
    h2y=spec.hyp2f1(3-gamma,3-gamma,4-gamma,-r/rs)
    mass = (r/r200)**(3-gamma)*M200*h2y/h2x
     
    if (Q!=0 and tmass!=0):  
         dphidr=dphigNFW(x,rs,rhos,gamma,Q,phinf,expcutoff)
    else:
         dphidr=0.0
    return mass+dphidr     
     

M_gnfw_CS=np.vectorize(M_gnfw_CS)



#modified gravity functions for other mass profiles ***************************
def dphibn(x,rhos,rs,b,Q,phinf,expcutoff=False):
    #This function compute Meff=Q/G *r**2 * dphidr (where phi is Phi/Mpl)
    #rhos is the normalization. It can be expressed in term of M200 and c200
    #for each profile
    clight=3e5
    limit=0.001 #limiting value for the screening radius
    G=4.302e-9
    
    Facda=Q*rs**2*rhos*8*np.pi*G/clight**2 #factor containing the planck mass
    
    #eponential cutoff:
    dminf=np.exp(-x*rs*np.sqrt(1e-7*Q/phinf))
    if expcutoff==False:
        dminf=1.0
        
    if (b!=1 and b!=2):
        cdav=Q**2*rs**3*rhos*8*np.pi/clight**2 #factor in front of everything
    #general case
        xcs=(phinf*(b-1)/Facda)**(1/(1-b))-1
        if xcs>=limit:
            #already put away the constant factor
            Cs=(1+xcs)**(1-b)*(1-xcs+b*xcs)/((b-1)*(b-2)) 
        else: #case of field always in the unscreened regime
            xcs=0            
            Cs=1/((b-1)*(b-2)) 
        if (x>xcs):
            out=((x+1)**(1-b)*(x-b*x-1)/((b-1)*(b-2))+Cs)*cdav*clight**2*dminf
        else:
            out=0
    #NFW:
    elif b==2:
        
        xcs=Facda/phinf-1

        if xcs>=limit:
            Cs=-phinf-Facda*np.log(Facda/phinf)
        else:
            xcs=0
            Cs=-Facda
        if (x>xcs):
            out=(Facda*(1/(x+1.0)+np.log(1+x))+Cs)*rs*Q*clight**2/G*dminf
        else:
            out=0
    else:
        out=0
    return out

dphibn=np.vectorize(dphibn)

#Effective mass
def M_bnfw_CS(r,r200,rs,tmass,Q,b,z,H0,Om,Ol,expcutoff=True):
     c200=r200/rs
     x=r/rs
     G=4.302e-9
     phinf=tmass*1e-5
    
     Hz2=H0**2*(Om*(1+z)**3+Ol)
     M200=100*Hz2/G*r200**3
     
     
     if (b!=1 and b!=2):
         f200=(c200+1)**(-b)*((c200+1)**b-b*c200*(c200+1)+c200**2-1)
         mass=M200*(1 + (1 + x)**(1 - b)*(-1 + x - b*x))/f200
         
         
         rhos=M200*(b-1)*(b-2)/(4*np.pi*rs**3*f200)
         if (Q!=0 and tmass!=0):  
             dphidr=dphibn(x,rhos,rs,b,Q,phinf,expcutoff)
         else:
             dphidr=0.0
         return mass+dphidr
     elif b==2:
         
         fac200=np.log(1+c200)-c200/(1+c200) 
         mass=M200*(np.log(1+x)-x/(1+x))/fac200
         
         rhos=M200/(4*np.pi*rs**3*fac200)
         if (Q!=0 and tmass!=0):  
             dphidr=dphibn(x,rhos,rs,b,Q,phinf,expcutoff)
         else:
             dphidr=0.0
         
         return mass+dphidr
         
     else:
         return 0
     
M_bnfw_CS=np.vectorize(M_bnfw_CS)     

def dphiGen(r,rhos,rhog,rhob,rs,rg,rj,Q,phinf,
            xsc = 0.001, expcutoff = True):
    clight=3e5
    limit=0.001 #limiting value for the screening radius
    G=4.302e-9    
    Facda=Q*8*np.pi*G/clight**2

    #eponential cutoff:
    dminf = np.exp(-r*np.sqrt(1e-7*Q/phinf))
    if expcutoff == False:
        dminf = 1.0
    if(xsc>limit):
        if (r<= xsc):
            return 0
        else:
           Cgen = - Facda*(rg**2*rhog*xsc - rhob*rj**4/(4*np.pi*(rj+xsc))
            + rhos*rs**4/(rs+xsc)-rg**3*rhog*np.arctan(xsc/rg)+rhos*rs**3*np.log(rs+xsc))
           
           dout = Cgen + Facda*(rg**2*rhog*r - rhob*rj**4/(4*np.pi*(rj+r))
            + rhos*rs**4/(rs+r) -
            rg**3*rhog*np.arctan(r/rg)+rhos*rs**3*np.log(rs+r))
    else:
        Cgen = Facda/(rhob*rj**3/(4*np.pi) - rhos*rs**3*(1+np.log(rs)))
        dout = Cgen + Facda*(rg**2*rhog*r - rhob*rj**4/(4*np.pi*(rj+r))
         + rhos*rs**4/(rs+r)-
         rg**3*rhog*np.arctan(r/rg)+rhos*rs**3*np.log(rs+r))
    return dout*Q*clight**2/G*dminf
dphiGen=np.vectorize(dphiGen)

#Effective mass
def M_Gen_CS(r,r200,rs,M,rj,rhog,rg,tmass,Q,z,H0,Om,Ol,
             xsc = 0.001, expcutoff=False):
    """
    
    Parameters
    ----------
    r : FLOAT
        input radius in Mpc.
    r200 : FLOAT
        virial radius of the DM profile.
    rs : FLOAT
        scale radius of the density profile.
    M : FLOAT
        Stellar mass of the BCG.
    rj : FLOAT
        scale radius of the BCG.
    rhog : FOLAT
        inner density of the baryonic profile.
    rg : FLOAT
        scale radius of the baryonic profile.
    tmass : FLOAT
        Field at infinity (in unit of 1e-5c**2).
    Q : FLOAT
        Coupling with matter.
    z : FLOAT
        redshift of the cluster.
    H0 : FLOAT
        Hubble constant at z=0.
    Om : FLOAT
        Omega matter.
    Ol : FLOAT
        Omega Lambda.
    xsc : FLOAT, optional
        Value of the screening radius. The default is 0.001.
    expcutoff : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    FLOAT 
    the value of the dynamical mass at that radius
    

    """    
    c200=r200/rs
    x=r/rs
    G=4.302e-9
    phinf=tmass*1e-5
    
    Hz2=H0**2*(Om*(1+z)**3+Ol)
    M200=100*Hz2/G*r200**3
    
    fac200=np.log(1+c200)-c200/(1+c200)     
    rhos=M200/(4*np.pi*rs**3*fac200)   
    
    rho0BCG=M/(np.pi*4*rj**3)
    
    Mbcg = M*r/rj/(1.+r/rj)
    
    MDM= M200*(np.log(1+x)-x/(1+x) )/fac200
    
    Mgas = 4*np.pi*rhog*rg**2*(r-rg*np.arctan(r/rg))
    
    mass = Mbcg + MDM + Mgas
    if (Q!=0 and tmass!=0):  
        dphidr= dphiGen(r,rhos,rhog,rho0BCG,rs,rg,rj,Q,phinf,
                    xsc = xsc, expcutoff = True)
    else:
        dphidr=0.0
    return dphidr+mass

M_Gen_CS = np.vectorize(M_Gen_CS)



#************** Modified gravity profile: Burkert ***************************
#All the MG functions are x**2 dphidr*rs*Q*c**2
def dphib(x,rhos,rs,Q,phinf,xsc = 0.001, expcutoff = True):
    clight=3e5
    limit=0.001 #limiting value for the screening radius
    G=4.302e-9
    
    #eponential cutoff:
    dminf = np.exp(-x*rs*np.sqrt(1e-7*Q/phinf))
    if expcutoff == False:
        dminf = 1.0
    
    Facda=Q*rs**2*rhos*8*np.pi*G/clight**2 #factor containing the planck mass     
    #In this case, no analytic solution exists. Need to find numerical solution
    
    Cb= 1/4*(-4*phinf + Facda*np.pi - 2*Facda*np.log(1 + xsc**2)) 
    
    if (xsc>limit):
        if (x<xsc):
            out=0
        else:
            out=Cb+ Facda/(4)*(-2*np.arctan(x) + 
            2*np.log(1 + x)+np.log(1 + x**2)) #/x**2 
    else:
        #test!
        #Cb= 1/4*(-4*phinf + Facda*np.pi) NO! Viene negativa la derivata 
        out=Facda/(4)*(-2*np.arctan(x) + 
            2*np.log(1 + x)+np.log(1 + x**2)) #+ Cb #/x**2
        
    return out*rs*Q*clight**2/G*dminf
dphiB=np.vectorize(dphib,otypes=[float],excluded=[1,2,3,4])


#Effective mass
def M_bur_CS(r,r200,rs,tmass,Q,z,H0,Om,Ol,xsc = 0.001, expcutoff=False):
    
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
        dphidr=dphiB(x,rhos,rs,Q,phinf,xsc,expcutoff)

    else:
        dphidr=0.0
    return dphidr+mass
     
     
M_bur_CS=np.vectorize(M_bur_CS,otypes=[float],excluded=[1,2,3,4,5,6,7,8,9])     


# Modified gravity Chameleon Einasto profile ********************************
def dphiE(x,rs,rhos,n,Q,phinf,expcutoff=True):
    clight=3e5
    limit=0.001 #limiting value for the screening radius
    G=4.302e-9
    
    
    #eponential cutoff:
    dminf=np.exp(-x*rs*np.sqrt(1e-7*Q/phinf))
    if expcutoff==False:
        dminf=1.0
    
    Facda=Q*rs**2*rhos*8*np.pi*G/clight**2 #factor containing the planck mass     
    #In this case, analytic solution exists. No need to find numerical solution
    
    xsc=(gammainccinv(2*n,4**n*phinf*n**(2*n-1)/(Facda*np.exp(2*n)*gamma(2*n)))/(2*n))**n
    if(np.isnan(xsc)):
        xsc=limit
    CE=8**(-n)*np.exp(2*n)*Facda*n**(1-3*n)*gammaincc(3*n,(2*n*(xsc)**(1/n)))*gamma(3*n)
    if (xsc>limit):
        if (x<xsc):
            out=0
        else:
            out=CE-(8)**(-n)*np.exp(2*n)*Facda*n**(1-3*n)*gammaincc(3*n, 2*n*x**(1/n))*gamma(3*n) #/x**2
    else:
        Cce=8**(-n)*np.exp(2*n)*n**(1 - 3*n)*gamma(3*n)*Facda
        out=Cce-(8)**(-n)*np.exp(2*n)*Facda*n**(1-3*n)*gammaincc(3*n, 2*n*x**(1/n))*gamma(3*n) #/x**2
    
    return out*rs*Q*clight**2/G*dminf
dphiE=np.vectorize(dphiE,otypes=[float])#,excluded=[1,2,3,4,5])


# Effective mass
def M_Ein_CS(r,r200,rs,n,tmass,Q,z,H0,Om,Ol,expcutoff=False):
    G=4.302e-9
    #z=za
    Hz2=H0**2*(Om*(1+z)**3+Ol)
    M200=100*Hz2/G*r200**3 
    c200=r200/rs
    phinf=tmass*1e-5
    
    x=r/rs
    fac200=(1-gammaincc(3*n,(2*n*(c200)**(1/n))))*gamma(3*n)   
    mass= M200/fac200*(1-gammaincc(3*n,(2*n*(x)**(1/n))))*gamma(3*n)  
    
    rhos=M200*(8**n)*n**(3*n-1)/(4*np.pi*rs**3*fac200*np.exp(2*n))
    if (Q!=0 and tmass!=0):  
        
        dphidr=dphiE(x,rs,rhos,n,Q,phinf,expcutoff)
    else:
        dphidr=0.0
    return dphidr+mass

M_Ein_CS=np.vectorize(M_Ein_CS,otypes=[float])#,excluded=[1,2,3,4,5,6,7,8,9])    


def dphiso(x,rhos,rs,Q,phinf,expcutoff=False):
    #This function compute Meff=Q/G *r**2 * dphidr (where phi is Phi/Mpl)
    #rhos is the normalization. It can be expressed in term of M200 and c200
    #for each profile
    clight=3e5
    limit=0.001 #limiting value for the screening radius
    G=4.302e-9
    
    #eponential cutoff:
    dminf = np.exp(-x*rs*np.sqrt(1e-7*Q/phinf))
    if expcutoff == False:
        dminf = 1.0
    
    Facda=Q*rs**2*rhos*8*np.pi*G/clight**2 #factor containing the planck mass     
    #In this case, no analytic solution exists. Need to find numerical solution
    
    if (Facda/phinf>1):
        
        xsc = np.sqrt((Facda/phinf)**2 - 1) 
        Ciso = Facda*( xsc/np.sqrt(1+xsc**2) - np.log(np.sqrt(1+xsc**2) + xsc )  )
        
    else:
        xsc = limit
        Ciso = 0

    if (x >= xsc):
        
        outpdd = rs*(Ciso - Facda*(x/np.sqrt(1+x**2) - np.log(np.sqrt(1+x**2)  + x ) ))
        #return rs*(Ciso - (Facda* x)/np.sqrt(1 + x**2) - Facda* np.log(1/(x + 
       #np.sqrt(1 + x**2) )))*Q*clight**2/G*dminf
         
    else:
        outpdd=  0.0
    return outpdd*Q*clight**2/G
dphiso = np.vectorize(dphiso)


def phi_iso(x,rhos,rs,Q,phinf):
    #rhos is the normalization. It can be expressed in term of M200 and c200
    #for each profile
    clight=3e5
    limit=0.001 #limiting value for the screening radius
    G=4.302e-9
    
    
    Facda=Q*rs**2*rhos*8*np.pi*G/clight**2 #factor containing the planck mass     
    #In this case, no analytic solution exists. Need to find numerical solution
    
    if (Facda/phinf>1):
        
        xsc = np.sqrt((Facda/phinf)**2 - 1) 
        #Ciso = phinf*xsc - Facda*np.log(np.sqrt(1+xsc**2)+xsc)
        
        Ciso = Facda*( xsc/np.sqrt(1+xsc**2) - np.log(np.sqrt(1+xsc**2) + xsc )  )
        
    else:
        xsc = limit
        Ciso = 0
    
    if (x>=xsc):
        
        outpdd = -Ciso/x - Facda*(np.log(np.sqrt(1+x**2)  + x ) )/x + phinf
        
    else:
        outpdd = 0.0
        
    return outpdd    
phi_iso = np.vectorize(phi_iso)


def MIso_CS(r,r200,rs,tmass,Q,z,H0,Om,Ol,expcutoff=False):
    G=4.302e-9
    #z=za
    Hz2=H0**2*(Om*(1+z)**3+Ol)
    M200=100*Hz2/G*r200**3 
    c200=r200/rs
    phinf=tmass*1e-5
    x = r/rs

    sqc = np.sqrt(1+c200**2)           
    rhoss= - M200*sqc/((c200 + sqc*np.log(np.sqrt(c200**2 +1) - c200) ) )    
    mGR = rhoss*( -x/(np.sqrt( 1 + x**2 ) ) + np.log( np.sqrt( 1 + x**2 ) + x  )     )
    
    rhos = rhoss/(np.pi*4*rs**3)

    if (Q!=0 and tmass!=0):  
        
        dphidr=dphiso(x,rhos,rs,Q,phinf,expcutoff=expcutoff)
    else:
        dphidr=0.0
    return dphidr+mGR

MIso_CS = np.vectorize(MIso_CS)


#===================== SCREENING FUNCTIONS ====================================



# The routines below computes the screening radius as a function of the model
# parameters in Chameleon Gravity, depending on the mass profile model


# Burkert *********************************************************************


def ScreenBur(r200,rs,Q,phinf,H0,z,Om,Ol):
    clight=3e5
    limit=0.0001 #limiting value for the screening radius
    G=4.302e-9
    c200=r200/rs
    Hz2=H0**2*(Om*(1+z)**3+Ol)
    M200=100*Hz2/G*r200**3
     
    fac200=1.0/(np.log(1+c200*c200)+2.*np.log(1.0+c200)-2.*np.arctan(c200))
     
    rhos=M200*fac200/(np.pi*rs**3)
    Facda=Q*rs**2*rhos*8*np.pi*G/clight**2
# No analytic solution for the screening radius. Find a numerical solution to 
# the following equation:     
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
    return xsc
ScreenBur=np.vectorize(ScreenBur)


# NFW *************************************************************************
def ScreenNFW(r200,rs,Q,phinf,H0,z,Om,Ol):
    clight=3e5
    limit=0.001 #limiting value for the screening radius
    G=4.302e-9
    c200=r200/rs
    Hz2=H0**2*(Om*(1+z)**3+Ol)
    M200=100*Hz2/G*r200**3
    fac200=np.log(1+c200)-c200/(1+c200) 
    rhos=M200/(4*np.pi*rs**3*fac200)

    Facda=Q*rs**2*rhos*8*np.pi*G/clight**2
    xcs=Facda/phinf-1
    if xcs <=limit:
        return limit
    return xcs
ScreenNFW=np.vectorize(ScreenNFW) 


def ScreengNFW(r200,rs,gamma,Q,phinf,H0,z,Om,Ol):
    clight=3e5
    limit=1e-6 #limiting value for the screening radius
    G=4.302e-9
    c200=r200/rs

    Hz2=H0**2*(Om*(1+z)**3+Ol)
    M200=100*Hz2/G*r200**3
    h2x=spec.hyp2f1(3-gamma,3-gamma,4-gamma,-r200/rs)
    
    rhos=((M200 * r200**(-3 + gamma)* rs**(-gamma)*(3 - gamma))/(4*np.pi*h2x))

    Facda=Q*rs**2*rhos*8*np.pi*G/clight**2
    
    if ((phinf/Facda - 1)*(gamma - 2) < 0 ):
        return limit
    
    xsc = 1/(1 - ( (phinf/Facda - 1)*(gamma - 2) )**(1/(2-gamma)) ) -1
    
    
    if(np.isnan(xsc)):
        xsc=limit
    
    return xsc
ScreengNFW=np.vectorize(ScreengNFW)     
    

# bNFW (b!=1,2) ***************************************************************
def ScreenbNFW(r200,rs,b,Q,phinf,H0,z,Om,Ol):
     clight=3e5
     limit=0.001 #limiting value for the screening radius
     G=4.302e-9
     c200=r200/rs

     Hz2=H0**2*(Om*(1+z)**3+Ol)
     M200=100*Hz2/G*r200**3
     if (b!=1 and b!=2):
        f200=(c200+1)**(-b)*((c200+1)**b-b*c200*(c200+1)+c200**2-1)    
        rhos=M200*(b-1)*(b-2)/(4*np.pi*rs**3*f200)
        Facda=Q*rs**2*rhos*8*np.pi*G/clight**2
        xcs=(phinf*(b-1)/Facda)**(1/(1-b))-1
        if xcs<=limit:
            xcs=limit
        return xcs
     elif b == 2:
         return ScreenNFW(r200,rs,Q,phinf,H0,z,Om,Ol)
     else:
        return 0
ScreenbNFW=np.vectorize(ScreenbNFW) 

# Isothermal ==================================================================
def ScreenIso(r200,rs,Q,phinf,H0,z,Om,Ol):
    clight=3e5
    limit=0.001 #limiting value for the screening radius
    G=4.302e-9
    c200=r200/rs
    Hz2=H0**2*(Om*(1+z)**3+Ol)
    M200=100*Hz2/G*r200**3
    rhos=M200/(4*np.pi*rs**3*(c200/(np.sqrt(1+c200**2)) + np.log(np.sqrt(c200**2 +1) + c200) ) )

    Facda=Q*rs**2*rhos*8*np.pi*G/clight**2
    if (Facda/phinf>1):       
    
        xcs = np.sqrt((Facda/phinf)**2 - 1)
    else:
        xcs = limit  
    
    if xcs <=limit:
        return limit
    return xcs
ScreenIso=np.vectorize(ScreenIso) 




# Einasto *********************************************************************

def ScreenEin(r200,rs,n,Q,phinf,H0,z,Om,Ol):
    clight=3e5
    limit=0.001 #limiting value for the screening radius
    G=4.302e-9
    c200=r200/rs

    Hz2=H0**2*(Om*(1+z)**3+Ol)
    M200=100*Hz2/G*r200**3
    fac200=(1-gammaincc(3*n,(2*n*(c200)**(1/n))))*gamma(3*n) 
    rhos=M200*(8**n)*n**(3*n-1)/(4*np.pi*rs**3*fac200*np.exp(2*n))
    Facda=Q*rs**2*rhos*8*np.pi*G/clight**2

    xsc=(gammainccinv(2*n,4**n*phinf*n**(2*n-1)/(Facda*np.exp(2*n)*gamma(2*n)))/(2*n))**n
    if(np.isnan(xsc)):
        xsc=limit
    
    return xsc
ScreenEin=np.vectorize(ScreenEin) 




# This function compute the screening radius for Einasto model but using the
#numerical solution of the algebric equation 
def ScreenEinNum(r200,rs,n,Q,phinf,H0,z,Om,Ol):
    clight=3e5
    limit=0.001 #limiting value for the screening radius
    G=4.302e-9
    c200=r200/rs

    Hz2=H0**2*(Om*(1+z)**3+Ol)
    M200=100*Hz2/G*r200**3
    fac200=(1-gammaincc(3*n,(2*n*(c200)**(1/n))))*gamma(3*n) 
    rhos=M200*(8**n)*n**(3*n-1)/(4*np.pi*rs**3*fac200*np.exp(2*n))
    Facda=Q*rs**2*rhos*8*np.pi*G/clight**2


    def fein(y):
          return 4**n*phinf*n**(2*n-1)/(Facda*np.exp(2*n)*
                 gamma(2*n))-gammaincc(2*n,2*n*y**(1/n))          
       
    xsc=limit
    x0=limit   
    for i in range(0,4000):
        xx=x0+0.01*i
        if fein(x0)*fein(xx) < 0:
            sol=opt.root(fein,xx)
            if sol.success:
                xsc=sol.x[0]
            else:
                xsc=xx                
            break
        x0=xx
        if (xsc==limit):
              xsc=xx
    return xsc 



    

    

#BCG VDP functions
#Only for DarkMatter gNFW
def fintegrand(xl,Q,phinf,r200, rs, gamma, xml, xmstar, rjaf, H0, z, Om, Ol, 
               rsgas = 0.16, rhogas = 60.80e13, rstar = 0.18, rhostar = 1e14,
               xsc = 1e-6, rabcg = 20.0, anibcg ='OM'):
      x=np.exp(xl)
      
      rhob = xmstar*xml/(np.pi*4*rjaf**3)
      
      fmass = Mtot(x,Q,phinf,r200,rs,rhostar,rstar,rhogas,
           rsgas,rhob, gamma,rjaf,xsc = xsc,
              z= z, H0= H0,  Om = Om, Ol  = Ol )
      
      fint=frho(x,rjaf)*fmass/x*fbetaint(x,rabcg,anibcg) #integral in dlog!
      return fint
  
    
def fintegrandOLD(xl,r200, rs, gamma, xml, xmstar, rjaf, H0, z, Om, Ol, 
               rsgas = 0.16, rhogas = 60.80e13,
               rabcg = 20.0, anibcg ='OM'):
      x=np.exp(xl)
      
      fmass = M_tot(x,r200, rs, gamma, xml, xmstar, rjaf, H0, z, Om, Ol,
              rsgas =rsgas, rhogas = rhogas)
      
      fint=frho(x,rjaf)*fmass/x*fbetaint(x,rabcg,anibcg) #integral in dlog!
      return fint


def fintegrand2(xl,eps, ldens, Q, xml, xmstar, rjaf, 
               rsgas = 0.16, rhogas = 60.80e13, rabcg = 20.0, anibcg ='OM'):
      x=np.exp(xl)
      fmass = M_RG(x, eps, ldens, Q, xml, xmstar, rjaf, rsgas, rhogas)
      fint=frho(x,rjaf)*fmass/x*fbetaint(x,rabcg,anibcg) #integral in dlog!
      return fint


#     exp(2*int(beta/x dx)) - OM model / constant model
#     !!! Note that the constant part of the integration is in the MAIN pgm

def fbetaint(x,rabcg,anibcg = 'OM'):
    if (anibcg == 'OM'):
        fb=rabcg*rabcg+x*x
    else: # (anibcg == 'C') then 
        bcgbetac=rabcg
        fb=x**(2.*bcgbetac) # beta=bcgbetac=constant at all radii 
    
    return fb



def frho(x,rjaf):
    return rjaf/(4.*np.pi*x*x*(rjaf+x)*(rjaf+x))
    

#     projected density profile (N(R) or Sigma*(R) or I(R))
#     Jaffe or Hernquist

def frhoproj(x, rjaf):

#      parameter (pi=3.1415927d0, Dmj=0.184565d0)
      
   s=x/rjaf

   if (abs(s-1.) < 0.001): 

     s1=0.999
     cs1=np.arccosh(1.0/s1)
     s2=1.001
     cs2=np.arccos(1.0/s2)
     
     ticaz1=(1.0-(2.0-s1*s1)*cs1/np.sqrt(abs(s1*s1-1.0)))
     ticaz1=ticaz1/2.0/(s1*s1-1.0)  
     frhoproj1=np.pi/(4.0*s1)-ticaz1
 
     ticaz2=(1.0-(2.0-s2*s2)*cs2/np.sqrt(abs(s2*s2-1.0)))
     ticaz2=ticaz2/2.0/(s2*s2-1.0)               

     frhoproj2=np.pi/(4.0*s2)-ticaz2

     frhoproj=(frhoproj2-frhoproj1)/(s2-s1)*(s-s1)+frhoproj1

   else:

     if (s < 1.):
         cs=np.arccosh(1.0/s)
     else:
         cs=np.arccos(1.0/s)

     pipp=(1.0-(2.0-s*s)*cs/np.sqrt(abs(s*s-1.0)))/2.0/(s*s-1.0)
     frhoproj=np.pi/(4.0*s)-pipp
  
   return frhoproj/(np.pi*rjaf*rjaf)

frhoproj = np.vectorize(frhoproj)


def fbeta(x,rabcg, anibcg):
    if (anibcg == 'OM'):
        fb=x*x/(rabcg*rabcg+x*x)
    else: # (anibcg == 'C') then 
        fb=rabcg
        
    
    return fb
    




def sigmaLOS_BCG(R, r200, rs, Y1, Y2, Y3, xml, xmstar, rjaf,z, H0,  Om, Ol, 
               rsgas = 0.16, rhogas = 60.80e13, rstar = 0.18, rhostar = 1e14,
               xsc = 1e-6, rabcg = 20.0, anibcg ='OM', rtbcg=20.0, 
               nout = 25, nmass = 'GR', b2 = 0.3):
    """
    
    Parameters
    ----------
    
    r : FLOAT
        value of the radius in Mpc.
    Q : FLOAT
        Coupling constant 
    phinf : FLOAT
        Value of the field at infinity (in unit of c**2).
    rhos : FLOAT
        dark matter profile critical density.
    rs : FLOAT
        dark matter profile scale radius (in Mpc).
    rhostar : FLOAT
        galaxies profile critical density.
    rstar : FLOAT
        galaxies profile scale radius (in Mpc).
    rhog : FLOAT
        DESCRIPTION.
    rg : FLOAT
        gas profile scale radius (in Mpc).
    rhob : FLOAT
        BCG profile critical density
    gam : FLOAT
        Slope parameter of the dark matter profile.
    rj : FLOAT
       BCG profile scale radius (in Mpc).
    xsc : FLOAT, optional
        Value of the screening radius. The default is 1e-6.
    
    
    R : FLOAT
        value of the projected radius in Mpc.
    r200 : FLOAT
        dark matter r200 (in Mpc).
    rs : FLOAT
        dark matter profile scale radius (in Mpc).
    Y1 : FLOAT
        Value of the Chameleon field at infinity (in unit of 1e-5c**2),
        or value of the permittivity in vacuum for RG framework.
    Y2 : FLOAT
        Coupling constant for Chameleon gravity, or value of the logarithm 
        of the critical density for RG framework
    Y3 : FLOAT
        gamma parameter of the gNFW profile for Dark Matter (Chameleon or GR)
        or exponent Q for the RG framework.
    xml : FLOAT
        Luminosity of the BCG in L_sun.
    xmstar : FLOAT
        Stellar Mass to light ratio of the BCG.
    rjaf : FLOAT
        jaffe radius for the stellar profile of the BCG.
    z : FLOAT
        redshift.
    H0 : FLOAT
        Hubble parameter at z = 0.
    Om : FLOAT
        Omega matter.
    Ol : FLOAT
        Omega Lambda.
    rsgas : FLOAT, optional
        Scale radius of the gas profile. The default is 0.16.
    rhogas : FLOAT, optional
        Critical density of the gas profile. The default is 60.80e13.
    rstar : FLOAT, optional
        Scale radius of the galaxy stellar profile. The default is 0.18.
    rhostar : TYPE, optional
        Critical density of the galaxy stellar profile. The default is 1e14.
    xsc : FLOAT, optional
        Value of the screening in Chameleon gravity. The default is 1e-6.
    rabcg : FLOAT, optional
        Scale radius of the anisotropy profile of the BCG. The default is 20.0.
    anibcg : STRING, optional
         Model of the. The default is 'OM'. Other option is 'C'
    rtbcg : FLOAT, optional
        Integration radius of the BCG. The default is 20.0.
    nout : TYPE, optional
        DESCRIPTION. The default is 25.
    nmass : TYPE, optional
        DESCRIPTION. The default is 'GR'.
    b2 : TYPE, optional
        DESCRIPTION. The default is 0.3.

    Returns
    -------
    TYPE
        DESCRIPTION.
    
    EXAMPLE OF USAGE
    
    Y1 = 10
    Y2 = 0.3
    rstar = 0.2
    rhostar = 1.00e+14
    rjaf = 0.038
    z= 0.44
    H0 = 70
    Om = 0.3
    Ol = 0.7
    rhogas = 2.57018005e+14
    rsgas = 0.38
    rs = 0.5
    rhob = xml*xmstar/(4*np.pi*rjaf**3)
    gam = 0.7
    xsc = ScreenTot(Y2,Y1*1e-5,r200,rs,rhostar,rstar,rhog,rg,rhob,gam,rjaf,
                 z = z, H0 = 70,  Om = 0.3, Ol = 0.7)
    testbcgGR2 = sigmaLOS_BCG(Rbcg, r200, rs, Y1, Y2, 0.7, xml, xmstar, 0.038,z, H0,  Om, Ol, 
                   rsgas = rsgas, rhogas = rhogas, rstar = rstar, rhostar = cc.rhostar,
                   xsc = xsc, rabcg = 20.0, anibcg ='OM', rtbcg=20.0, 
                   nout = 40, nmass = 'GR', b2 = 0.3)
    
    
    """
    grav = 4.302e-9
    bx=np.log(rtbcg/10.0) #rtbcg
    
    #a1=np.log(1e-4)

    #nout = 40
    #dlrr = np.linspace(a1,bx,nout)
    #deltalog=(bx-a1)/(nout)

    if nmass == 'RG':
        def fint(x):
            return fintegrand2(x,Y1, Y2, Y3, xml, xmstar, rjaf, 
                           rsgas, rhogas, rabcg, anibcg)
        
    elif nmass == 'GR':
        gam = Y3
        def fint(x):
            return fintegrand(x,Y2,Y1*1e-5, r200, rs, gam, xml, xmstar, 
                              rjaf, H0, z, Om, Ol, 
                           rsgas, rhogas, rstar, rhostar, xsc, rabcg, anibcg)        
    
    #elif nmass == 'OLD':
    #    def fint(x):
    #        gam = Y3
    #        return fintegrandOLD(x,r200, rs, gam, xml, xmstar, rjaf, H0, z, Om, Ol, 
    #                              rsgas = rsgas, rhogas = rhogas,
    #                              rabcg =rabcg, anibcg =anibcg) 
        
    def integral(x,b):
        
        if anibcg == 'C':
            denbeta=np.exp(2.*rabcg*x)
        elif (anibcg == 'OM'):
            denbeta=rabcg*rabcg+np.exp(2.*x)
        return ints.quad(fint,x,b,
                        epsabs=1.49e-8, epsrel=1.49e-8)[0]*grav/denbeta

    
    integral = np.vectorize(integral)
    
    #lsrnu = np.log(integral(dlrr,2*bx))
    #finterp = ipt.interp1d(dlrr,lsrnu, kind = 'cubic')
    
    def fabbeta(xl,ax):
        #s1 = np.exp(finterp(xl))
        s1 = integral(xl,2*bx)
        #print(s1)
        x=np.exp(xl)
        b1=fbeta(x,rabcg, anibcg)

        a=np.exp(ax)

        if (abs(1.-x/a) < 0.001):
          fab=0.0
        else:

          fab=(s1)*x*x/np.sqrt(x*x-a*a)*(1.-b1*a*a/(x*x))
        return fab    
        
    
    def integral2(x,b):
        a=np.exp(x)
        fir=frhoproj(a,rjaf)
        return np.sqrt(2.0*ints.quad(fabbeta,x,b, args = (x),
                        epsabs=1.49e-6, epsrel=1.49e-6)[0]/fir)
    integral2 = np.vectorize(integral2)
    
    
    

    return integral2(np.log(R),b2)
        
sigmaLOS_BCG = np.vectorize(sigmaLOS_BCG)
  



def M_RG(x, eps, ldens, Q, xml, xmstar, rjaf, 
         rhogas, rhostar, rsgas, rstar, kalpha = 0 ):
    """
    Function computing the effective mass in refracted gravity
    baryonic (gas+stellar) component is modeled as a beta profile with rho0gas and 
    rsgas as fixed parameters.

    Parameters
    ----------
    x : FLOAT (1D-Array)
        The value of the radius at which the effective potential is computed (Mpc).
    ldens : FLOAT
        The logarithm of the critical density, in M_sun/Mpc**3.
    eps : FLOAT
        Gravitational permittivity in vacuum.
    Q : FLOAT
        Exponent that defines the sharpness of the transition.
    xml : FLOAT
        BCG luminosity (in L_sun).
    xmstar : FLOAT
        BCG masss-to-light ratio.
    rjaf : FLOAT
        BCG scale radius.

    Returns
    -------
    FLOAT (1D-Array)
       The value of the effective mass at x.

    """
    #r200gas=17.066 
    #rsgas= 0.16 #MACS 1206 0.12   
    
    #Critical density: use the log of screen
    arhoc=10**(ldens)
    
    #gas and galaxies: beta profile
    #rho0gas= 60.80e13   # MACS 1206 89.08e13   
    
    xg = x/rsgas
    if (kalpha==1):
        #New: Gas (this is a beta profile with beta = 1)
        fbb1 = (-1.0/(xg**2*np.sqrt(1+xg**2))+ (np.arcsinh(xg))/xg**3)                  
    else:
        #profile with beta = 2/3 
        fbb1 = rsgas**2 * (x - rsgas*np.arctan(xg)) / x**3
     
    xstar = x/rstar
    fbb2 = (-1.0/(xstar**2*(1+xstar))+(np.log(1+xstar))/xstar**3)
       
    xmbar = 4*np.pi*x**3*(rhostar*fbb2 + rhogas*fbb1) 
    
    
    dmtot=xmstar*xml*x/rjaf/(1.+x/rjaf)+ xmbar
    
    #densità
    rbc= x/rjaf 
    # density of the BCG
    rhobcg=xmstar*xml/(4.0*np.pi*rjaf)/((rbc)**2*(1.0+rbc)**2)
    #density baryons
    #NEW! Density gas and galaxies:
    if (kalpha == 1):
        rhogass = rhogas/(1.0 +(x/rsgas)**2)**(3.0/2)
    else:
        rhogass = rhogas/(1.0 +(x/rsgas)**2)

    rhogals = rhostar/( x/rstar * (1.00 + x/rstar)**2)
    
    rhotot=rhobcg + rhogass + rhogals       

    epsi=eps+(1-eps)*(1.0/2)*(np.tanh(Q*np.log(rhotot/arhoc))+1)

    return (dmtot) /(epsi)
    
M_RG = np.vectorize(M_RG)   


#Boson star approximation

def M_BS(r,rs,b,g,Minf,xM0):
    """
    Function to approximate a "cluster Boson star"

    Parameters
    ----------
    r : FLOAT
        radius (Mpc).
    rs : TYPE
        scale radius of the profile (Mpc).
    b : FOLAT
        First exponent.
    g : FLOAT
        Second exponent.
    Minf : FLOAT
        BS Mass at infinity (M_sun).
    xM0 : FLOAT
        Fraction of the boson star mass in the inner part.

    Returns
    -------
    FLOAT
        Value of the mass.

    """
    if (r/rs <= 2.0): 
        return  Minf*(np.tanh((r/rs)**b)+xM0*(r/rs)**g)
    else:
       return  Minf*(np.tanh((2.0)**b)+xM0*(2.0)**g) 
M_BS = np.vectorize(M_BS)



#***** DEFINE NUMBER DENISTY PROFILE FUNCTIONS ******************************
#******************* N(R) NFW model *****************************************
def NpNFW(R,rnu,r200,Npar):
    
    c=r200/rnu
    Mrs=Npar/(np.log(1+c)-c/(1+c))
    xx=R/rnu
    try:
        len(xx)
        npp=[]
        for xxx in xx:
            if xxx<1:
                npp.append( Mrs*(np.arccosh(1/xxx)/(np.sqrt(1-xxx**2))+np.log(xxx/2)))
            elif xxx==1:
                npp.append( Mrs*(1-np.log(2)))
            elif xxx>1:
                npp.append( Mrs*(np.arccos(1/xxx)/(np.sqrt(xxx**2-1))+np.log(xxx/2)))
        return  np.array(npp)
        
    except TypeError:
        if xx<1:
            return Mrs*(np.arccosh(1/xx)/(np.sqrt(1-xx**2))+np.log(xx/2))
        elif xx==1:
            return Mrs*(1-np.log(2))
        elif xx>1:
            return Mrs*(np.arccos(1/xx)/(np.sqrt(xx**2-1))+np.log(xx/2))


#******  N(R) - projected Hernquist - see Hernquist (1990)

def NpHer(R,rnu,r200,Npar):
      C=Npar*(rnu+r200)**2/r200**2
      s=R/rnu
      xs=1
      if s>1:
          xs=np.arccos(1./s)/np.sqrt(s*s-1.)
          sigmar2=s**2*(xs-1)/(1-s**2)
      elif s<1:
          xs=np.log((1.+np.sqrt(1.-s*s))/s)/np.sqrt(1.-s*s)
          sigmar2=s**2*(xs-1)/(1-s**2)
      else:
          sigmar2=1/3
     
      return sigmar2*C
NpHer=np.vectorize(NpHer)



#3D Number density profile **************************************************
def N3d(r,rs,r2,Npar,nnu):
#3d number profile   
    c=r2/rs 
    if nnu==1: #NFW
        fc=1/(np.log(1+c)-c/(1+c))
        return Npar*fc*(np.log(1+r/rs)-(r/rs)/(1+r/rs))
    elif nnu==2: #Hernquist
        fc=(rs+r2)**2/r2**2
        return Npar*fc*(r**2/(r+rs)**2)



#****************************************************************************
#Function for the anisotropy profile. Note: beta is the real anisotropy parameter
def betaf(r,rbetac,betaiic,anis,beta0=0): 
   if anis==1:
       return betaiic #constant anisotropy case
   elif anis==2:
       return r/(2.0*(r+rbetac)) #M&L case 
   elif anis==3:
       return r**2/(r**2+rbetac**2) #OM case
   elif anis==4:  #T case
       bec=betaiic#1.-1./(betaiic*betaiic)
       #if (bec>=1.):
       #      bec=0.9990  # avoid unphysical values
       return r*bec/((r+rbetac)) 
   elif anis==5:
       bec=betaiic #1.-1./(betaiic*betaiic)
       return (r-rbetac)*bec/((r+rbetac)) #O case
   elif anis==6:
       bec=betaiic #1.-1./(betaiic*betaiic)
       bec0=beta0 #1.-1./(beta0*beta0)
       return bec0+(bec-bec0)*r**2/(r**2+rbetac**2) #gOM
   elif anis==7:
       bec=betaiic #1.-1./(betaiic*betaiic)
       bec0=beta0 #1.-1./(beta0*beta0)
       return bec0+(bec-bec0)*r/(r+rbetac) #gT
   elif anis==8:
       bec=betaiic #1.-1./(betaiic*betaiic)
       bec0=beta0 #1.-1./(beta0*beta0)
       return bec0 + (bec - bec0)*(r/(r + rbetac) + \
                (r**2/rbetac**2)*np.exp(-(r/rbetac)**2)) #BP
   else:
       return 0
   
betaf=np.vectorize(betaf)


# ========================== Lensing functions ================================
#Probably I will add them in a separate python module

#COSMOLOGY

#E(z)^-1 time evolution of the Hubble factor
def Ez(z,Om,Ol): 
    return 1/(Om*(1+z)**3+Ol)**(1/2)

#comoving distance (Mpc)
def Chi(z,H0,Om,Ol):
    cl=2.99792458e5 #speed of light
    prova=ints.quad(Ez,0,z,args=(Om,Ol),
                    epsabs=1.49e-8, epsrel=1.49e-8)
    return  cl/H0*prova[0]


def Age_univ(param):
#given a set of cosmological parameter Om, Or, Ok, Ol, H0
#it returns the age of the Universe    
    Om=param[0]
    Or=param[1]
    Ok=param[2]
    Ol=param[3]
    H0=param[4]
    def Eza(a):
        return 1/(a*(Om*a**(-3)+Or*a**(-4)+Ok*a**(-2)+Ol)**(1/2))
    prova=ints.quad(Eza,0,1,
                    epsabs=1.49e-8, epsrel=1.49e-8)
    agesec=1/(3.24078e-20*H0)*prova[0]
    return agesec/(3600*24*365)



#angular diameter distance between twp objects at redshift z1, z2 (Mpc) 
#z2>z1
def Da(z2,z1=0,H0=70,Om=0.3,Ol=0.7):
    return (Chi(z2,H0,Om,Ol)-Chi(z1,H0,Om,Ol))/(1+z2)  

Da =np.vectorize(Da)
#****************************************************************************
#lensing
#sources redshift distributions
def n(z,ngal):
    #ngal=30  #According to Chen et al., 2020, expetcted number from Euclid
    z0=1/3*(ngal/30)**(1/3)
    return z**2*np.exp(-z/z0)



def n_z(z,ngal=30,zm=10): #normalized
    #normalize the profile
    #zm=10 #maximum z at which the sources are expected to be localised
    normfac=ints.quad(n,0,zm,args=(ngal),epsabs=1.49e-8, epsrel=1.49e-8)[0]
    
    return n(z,ngal)/normfac

n_z=np.vectorize(n_z)


#Beta factors of lensing 
#defined as the ratio between Dls and Ds: to be averaged
#over the redshift distribution  
def beta_lens(zs,zc,H0=70,Om=0.3,Ol=0.7,ngal=30,zm=10):
    return Da(zs,zc,H0,Om,Ol)/Da(zs,0,H0,Om,Ol)*n_z(zs,ngal,zm)


#Squared: needed to compute the lensing strenght
def beta_square(zs,zc,H0=70,Om=0.3,Ol=0.7,ngal=30,zm=10):
    return (Da(zs,zc,H0,Om,Ol)/Da(zs,0,H0,Om,Ol))**2*n_z(zs,ngal,zm)

def lensfactor(zc,H0=70,Om=0.3,Ol=0.7,ngal=30,zm=10):
    
#This fucntion computes the critical surface density and the lensing factor 
#which goes in front of <k> in the expression of the reduced
#tangential shear (see Umetsu 2020)   
    betave=ints.quad(beta_lens,zc,zm,args=(zc,H0,Om,Ol,ngal,zm),
                     epsabs=1.49e-8, epsrel=1.49e-8)
    betave_squared=ints.quad(beta_square,zc,zm,args=(zc,H0,Om,Ol,ngal,zm),
                             epsabs=1.49e-8, epsrel=1.49e-8)
            
    sigcrit_ave=betave[0]*Da(zc,0,H0,Om,Ol) #this is <Dl*Dls/Ds>

#this is the factor which goes in front of <k> in the expression of the reduced
#tangential shear (see Umetsu 2020)
    fl=betave_squared[0]/betave[0]**2
    #print('the fl factor in averaged shear is ', fl)
    
    return sigcrit_ave,fl

#beta_inf=Da2(zc,zm)/Da(zm)  #value of beta infinity to put in the definition
                            #of the averaged quantities (note: it is unnecessary)


#root function for NFW
def fmin(x,eps2):
    return x**2-eps2*(np.log(1+x)-x/(1+x))
fmin=np.vectorize(fmin)


#root function for bNFW
def fminb(x,eps2,b):
    return x**2-eps2*((1+x)**(1-b)*(1+(b-1)*x))/((b-1)*(b-2))
fminb=np.vectorize(fminb)

#root function for gNFW
def fming(x,eps2,b):
    h2y=spec.hyp2f1(3-b,3-b,4-b,-x)
    return x**2-eps2*x**(3-b)*h2y
fminb=np.vectorize(fming)


def Re_sim(M,zc,H0=70,Om=0.3,Ol=0.7,zm=5):
    cl=2.99792458e5
    G=4.301e-9
    return np.sqrt(Da(zm,zc,H0,Om,Ol)/Da(zc,0,H0,Om,Ol)/
                   Da(zm,0,H0,Om,Ol)/cl**2*4*G*M)*180/np.pi*3600

def Einstein_radius(r200,rs,zc,H0=70,Om=0.3,Ol=0.7,zm=5,deg=False,
                    model='NFW',b=3):
    """
    This function compute the Einstein radius for a lens deflector at a given 
    redshfit, assuming a NFW mass profile

    Parameters
    ----------
    r200 : REAL
        virial radius of the halo.
    rs : REAL
        scale radius of the halo.
    zc : REAL
        average redhsift of the halo.
    H0 : REAL, optional
        Hubble parameter at z=0. The default is 70.
    Om : REAL, optional
        Matter density parameter. The default is 0.3.
    Ol : REAL, optional
        Dark Energy density parameter. The default is 0.7.
    zm : REAL, optional
        redshift of the source. The default is 5.
    deg : BOOL, optional
        Do you want the result in arcseconds?. The default is False.
    model : STR, optional
        The mass profile is a NFW, bNFW, gNFW. The default is NFW.

    Returns
    -------
    REAL
        The Einstein radius of the halo in unit of Mpc.

    """
    
    cl=2.99792458e5 #speed of light
    #sigmave2,_=lensfactor(zc,H0,Om,Ol,30,zm)
    sigmave=Da(zm,zc,H0,Om,Ol)*Da(zc,0,H0,Om,Ol)/Da(zm,0,H0,Om,Ol)
    Dd=Da(zc,0,H0,Om,Ol)
    
    Hz2=H0**2*(Om*(1+zc)**3+Ol)    
    M200=100*Hz2*r200**3
    eps=sigmave*4*M200/cl**2
    c200=r200/rs
    if model=='NFW' or (model=='bNFW' and b==2):
        

      fac200=np.log(1+c200)-c200/(1+c200)  
      eps2=eps/fac200/rs**2

      xtest=np.linspace(0.001,100,1000)
      ftest=fmin(xtest,eps2)
      idx2 = np.where(np.sign(ftest[:-1]) != np.sign(ftest[1:]))[0] + 1
      
      if len(idx2)!= 0:
        
          xguess=xtest[idx2]
        
          roots=opt.fsolve(fmin,xguess,args=eps2)*rs  
          rootfinal=roots[0]
      else:
          rootfinal=1e-3
          
       
      if deg==True:
          rootfinal=rootfinal*180/np.pi*3600/Dd
      return rootfinal
    
    elif model=='bNFW':
      f200=(1 + (1 + c200)**(1 - b)*(-1 + c200 - b*c200))
        
      eps2=eps*(b-1)*(b-2)/(f200)/rs**2
      xguess=13.2
      roots=opt.fsolve(fminb,xguess,args=(eps2,b))*rs   
      rootfinal=roots[0]
      
      return rootfinal

  
    elif model=='gNFW':
        h2x=spec.hyp2f1(3-b,3-b,4-b,-r200/rs) 
        fac200=h2x*c200**(3-b)
        eps2=eps/fac200/rs**2
    
        
        xguess=1
        xtest=np.linspace(0.001,100,1000)
        ftest=fming(xtest,eps2,b)
        idx2 = np.where(np.sign(ftest[:-1]) != np.sign(ftest[1:]))[0] + 1  
        if len(idx2)!= 0:
        
            xguess=xtest[idx2]
        
            roots=opt.fsolve(fming,xguess,args=(eps2,b))*rs  
            rootfinal=roots[0]
        else:
            rootfinal=1e-3
        if deg==True:
           rootfinal=rootfinal*180/np.pi*3600/Dd
        return rootfinal       
          
        
    
    else:
        res=eps**0.5
        if deg==True:
            res=res*180/np.pi*3600/Dd
        return res


# SUFRACE DENISTY PROFILES ====================================================
# NFW
#
#this is the function of the surface density profile in GR (NFW)
def f_nfw(x):
    if (x<1):
        return 1/(1-x**2)*(-1+2/(np.sqrt(1-x**2))*np.arctanh(np.sqrt((1-x)/(1+x))))
    elif (x==1):
        return 1/3
    else:
        return 1/(x**2-1)*(1-2/(np.sqrt(x**2-1))*np.arctan(np.sqrt((x-1)/(1+x))))
    
f_nfw=np.vectorize(f_nfw)

#******************************************************************************
#this is the function for the average projected mass density M_pj GR (NFW)
def g_nfw(x):
    if (x<1):
        return 2/(x**2)*(2/(np.sqrt(1-x**2))*np.arctanh(np.sqrt((1-x)/(1+x)))+
                  np.log(x/2))
    elif (x==1):
        return 2*(1+np.log(1/2))
    else:
        return 2/(x**2)*(2/(np.sqrt(x**2-1))*np.arctan(np.sqrt((x-1)/(1+x)))+
                  np.log(x/2))
    
    
g_nfw=np.vectorize(g_nfw)

def gamma_NFW(R,r200,rs,zc,H0=70,Om=0.3,Ol=0.7,ngal=30,zm=10):
    
    sigmave,fl=lensfactor(zc,H0,Om,Ol,ngal,zm)
    
    c200=r200/rs
    x=R/rs
    cl=2.99792458e5 #speed of light
   
    Hz2=H0**2*(Om*(1+zc)**3+Ol)
    fac200=np.log(1+c200)-c200/(1+c200)  
    sigcrit_ave,_=lensfactor(zc,H0,Om,Ol,ngal,zm)
    facfront=200*sigcrit_ave*r200**3*Hz2/(rs**2*cl**2*fac200)
    return facfront*(g_nfw(x)-f_nfw(x))
gamma_NFW=np.vectorize(gamma_NFW) 



def kappa_NFW(R,r200,rs,zc,H0=70,Om=0.3,Ol=0.7,ngal=30,zm=10):
    """"
    This function compute the redshift-averaged convergence map <k> 
    for a NFW profile model of the lens mass distribution
    """
    sigmave,fl=lensfactor(zc,H0,Om,Ol,ngal,zm)    
    ct=r200/rs #concentration
    fac200t=np.log(1+ct)-ct/(1+ct)
    x=R/rs
    cl=2.99792458e5 #speed of light
    Hz2=H0**2*(Om*(1+zc)**3+Ol)

    facfront=200*sigmave*r200**3*Hz2
    facfront=facfront/(rs**2*cl**2*fac200t)
    ak_true=facfront*f_nfw(x)
    
    return ak_true
kappa_NFW=np.vectorize(kappa_NFW)



def gt_aveNFW(R,r200,rs,zc,H0=70,Om=0.3,Ol=0.7,ngal=30,zm=10):
    
    """
    This routine computes the average tangential shear profile,
    eq. (34) in Pizzuti et al., 2021, at the projected radius R for a NFW
     model described by the external parameters `r200, rs`.
    """
    sigmave,fl=lensfactor(zc,H0,Om,Ol,ngal,zm)    
    stro=(1-fl*kappa_NFW(R,r200,rs,zc,H0,Om,Ol,ngal,zm))
    gt_true=gamma_NFW(R,r200,rs,zc,H0,Om,Ol,ngal,zm)/stro
    
    return gt_true
gt_aveNFW=np.vectorize(gt_aveNFW)



         
#******************************************************************************
# bNFW
#Da sistemare come si deve...
#da capire come mai non ci va il fattore 4pi... Bisogna fare il conto da zero
def Integrand_bNFW(y,x,b):
    return (1 + y)**(-b)/np.sqrt(y**2 - x**2)

def Integrand2_bNFW(y,x,b):
    return  (1 + y)**(-b)*(y - np.sqrt(y**2 - x**2))


#Sigma(R/rs)/(rs*rhos)
def f_bNFW(x,b):
    pp=ints.quad(Integrand_bNFW,x,np.inf,args=(x,b),
                             epsabs=1.49e-8, epsrel=1.49e-8)
    return 2*pp[0]

f_bNFW=np.vectorize(f_bNFW)   
    
#Sigmaproj/rs^3rhos la differenza con il g_nfw è in un fattore 2 pi
#g_bNFW/(2pi*x**2)=g_nfw(x)
def g_bNFW(x,b):     
    if (b!=1 and b!=2):
        mass= (1 + (1 + x)**(1 - b)*(-1 + x - b*x))/((b-1)*(b-2))
    elif (b==2):
        mass=(-(x/(1 + x)) + np.log(1 + x))
    else:
        return 0
    pp=ints.quad(Integrand2_bNFW,x,np.inf,args=(x,b),
                             epsabs=1.49e-8, epsrel=1.49e-8)
    return 4*np.pi*(mass+pp[0])
        
g_bNFW=np.vectorize(g_bNFW) 



def gamma_bNFW(R,r200,rs,b,zc,H0=70,Om=0.3,Ol=0.7,ngal=30,zm=10):
    #gamma_t averaged over the redshift distribution with ngal per arcmin 
    #square and zm as a maximum redshift, for a cluster lens at z=zc
    #For b=2 it coincides with the expression above gamma_NFW, which
    #is recomanded as it has analytical expression
    sigmave,fl=lensfactor(zc,H0,Om,Ol,ngal,zm)
    
    c200=r200/rs
    x=R/rs
    cl=2.99792458e5 #speed of light
   
    Hz2=H0**2*(Om*(1+zc)**3+Ol)
    M200=100*Hz2*r200**3
    
    if (b!=1 and b!=2):
        f200=(1 + (1 + c200)**(1 - b)*(-1 + c200 - b*c200))
        
        rhos=M200*(b-1)*(b-2)/(4*np.pi*rs**3*f200)

    elif b==2:       
        fac200=np.log(1+c200)-c200/(1+c200)        
        rhos=M200/(4*np.pi*rs**3*fac200)
    else:
         return 0
    
    sigmar=f_bNFW(x,b)*rhos*rs 
    sigmamed=g_bNFW(x,b)*rhos*rs**3/(R**2*np.pi)
    
    return sigmave*4*np.pi/cl**2*(sigmamed-sigmar)  #*4*np.pi


def kappa_bNFW(R,r200,rs,b,zc,H0=70,Om=0.3,Ol=0.7,ngal=30,zm=10):
    """"
    This function compute the redshift-averaged convergence map <k> 
    for a bNFW profile model of the lens mass distribution
    """
    sigmave,fl=lensfactor(zc,H0,Om,Ol,ngal,zm)    
    c200=r200/rs
    x=R/rs
    cl=2.99792458e5 #speed of light
   
    Hz2=H0**2*(Om*(1+zc)**3+Ol)
    M200=100*Hz2*r200**3
    
    if (b!=1 and b!=2):
        f200=(1 + (1 + c200)**(1 - b)*(-1 + c200 - b*c200))
        
        rhos=M200*(b-1)*(b-2)/(4*np.pi*rs**3*f200)

    elif b==2:       
        fac200=np.log(1+c200)-c200/(1+c200)        
        rhos=M200/(4*np.pi*rs**3*fac200)
    else:
         return 0
    sigmar=f_bNFW(x,b)*rhos*rs 
    ak_true=sigmave*4*np.pi/cl**2*(sigmar) 
    
    return ak_true
kappa_bNFW=np.vectorize(kappa_bNFW)



def gt_bNFW(R,r200,rs,b,zc,H0=70,Om=0.3,Ol=0.7,ngal=30,zm=10):
    
    """
    This routine computes the average tangential shear profile,
    eq. (34) in Pizzuti et al., 2021, at the projected radius R for a bNFW
     model described by the external parameters `r200, rs`.
    """
    sigmave,fl=lensfactor(zc,H0,Om,Ol,ngal,zm)    
    stro=(1-fl*kappa_bNFW(R,r200,rs,b,zc,H0,Om,Ol,ngal,zm))
    gt_true=gamma_bNFW(R,r200,rs,b,zc,H0,Om,Ol,ngal,zm)/stro
    
    return gt_true
gt_ave_bNFW=np.vectorize(gt_bNFW)

#******************************************************************************
#gNFW
# Sono qui. Devo verificare la formula per sigma crit e per gNFW da mettere 


def Integrand_gNFW(y,x,b):
    return y**(1-b)*(1 + y)**(b-3)/(y**2 - x**2)**(0.5)

def Integrand2_gNFW(y,x,b):
    return  y**(1-b)*(1 + y)**(b-3)*(y - np.sqrt(y**2 - x**2))

#RICORDA: Sigma(R)=(rs*rhos)f(x,b) (x=R/rs)
def f_gNFW(x,b):
    pp=ints.quad(Integrand_gNFW,x,np.inf,args=(x,b),
                              epsabs=1.49e-8, epsrel=1.49e-8, limit=100)
    return 2*pp[0]

f_gNFW=np.vectorize(f_gNFW)

#Ricorda: Mproj(R)=rhos*rs**3*g(x,b)
def g_gNFW(x,b):     
    #if (b!=1):
    h2y=x**(3-b)*spec.hyp2f1(3-b,3-b,4-b,-x)/(3-b)
    mass=h2y
    #else:
    #    mass=(-(x/(1 + x)) + np.log(1 + x))

    pp=ints.quad(Integrand2_gNFW,x,np.inf,args=(x,b),
                             epsabs=1.49e-8, epsrel=1.49e-8)
    return 4*np.pi*(mass+pp[0])
        
g_gNFW=np.vectorize(g_gNFW) 



def gamma_gNFW(R,r200,rs,b,zc,H0=70,Om=0.3,Ol=0.7,ngal=30,zm=10):
    #Eq 35 in Pizzuti et al, 2021
    #gamma_t averaged over the redshift distribution with ngal per arcmin 
    #square and zm as a maximum redshift, for a cluster lens at z=zc
    #For b=2 it coincides with the expression above gamma_NFW, which
    #is recomanded as it has analytical expression
    sigmave,fl=lensfactor(zc,H0,Om,Ol,ngal,zm)
    
    c200=r200/rs
    x=R/rs
    cl=2.99792458e5 #speed of light
    
    Hz2=H0**2*(Om*(1+zc)**3+Ol)
    M200=100*Hz2*r200**3
    h2x=spec.hyp2f1(3-b,3-b,4-b,-r200/rs)   
    rhosrscubic=M200*(3-b)*c200**(b-3)/(4*np.pi*h2x) #rhos*rs**3
     
    
    sigmar=f_gNFW(x,b)*rhosrscubic/rs**2 
    sigmamed=g_gNFW(x,b)*rhosrscubic/(R**2*np.pi)
    
    return 4*np.pi*sigmave/cl**2*(sigmamed-sigmar) 

gamma_gNFW=np.vectorize(gamma_gNFW)


def kappa_gNFW(R,r200,rs,b,zc,H0=70,Om=0.3,Ol=0.7,ngal=30,zm=10):
    """"
    This function compute the redshift-averaged convergence map <k> 
    for a NFW profile model of the lens mass distribution
    """
    sigmave,fl=lensfactor(zc,H0,Om,Ol,ngal,zm)    
    c200=r200/rs #concentration
    x=R/rs
    cl=2.99792458e5 #speed of light
    Hz2=H0**2*(Om*(1+zc)**3+Ol)
    M200=100*Hz2*r200**3
    h2x=spec.hyp2f1(3-b,3-b,4-b,-r200/rs)   
    rhosrscubic=M200*(3-b)*c200**(b-3)/(4*np.pi*h2x) #rhos*rs**3


    sigmar=f_gNFW(x,b)*rhosrscubic/rs**2 

    ak_true=4*np.pi*sigmave/cl**2*sigmar
    
    return ak_true

kappa_gNFW=np.vectorize(kappa_gNFW)


def gt_ave_gNFW(R,r200,rs,b,zc,H0=70,Om=0.3,Ol=0.7,ngal=30,zm=10):
    
    """
    This routine computes the average tangential shear profile,
    eq. (34) in Pizzuti et al., 2021, at the projected radius R for a NFW
     model described by the external parameters `r200, rs`.
    """
    sigmave,fl=lensfactor(zc,H0,Om,Ol,ngal,zm)    
    stro=(1-fl*kappa_gNFW(R,r200,rs,b,zc,H0,Om,Ol,ngal,zm))
    gt_true=gamma_gNFW(R,r200,rs,b,zc,H0,Om,Ol,ngal,zm)/stro
    
    return gt_true
gt_ave_gNFW=np.vectorize(gt_ave_gNFW)


#******************************************************************************

class Model:
    def __init__(self,param,r200,rs,beta,rnu=0,Y1=0,Y2=0,Y3=0.0,
                 gamma=1,beta0=0,b=2,expc = False, xml = 1e12, xmstar = 4.6, rjaf = 0.029,
                 rsgas = 0.17, rhogas = 4.58229220e+14, rstar = 0.369, rhostar = 3.547e+13):
        

        #first: number identifying the models of mass, beta and nu(r)
        self.nmass=param[0]
        self.nbeta=param[1]
        self.nnu=param[2]
        #Then: cosmological parameters (Needed to compute M200)
        self.H0=param[3]
        self.z=param[4]
        self.Om=param[5]
        self.Ol=param[6] 
        
        self.r200=r200
        self.rs=rs
        self.beta=beta
        self.beta0=beta0
        if(rnu==0):
            self.rnu=rs
        else:
            self.rnu=rnu
        self.gamma=gamma        
        self.Y1=Y1
        self.Y2=Y2
        self.Y3 = Y3
        self.b=b
        self.expc = expc
        
        if self.nmass==6 and self.b<=1:
            print('value of b not acceptable. Switching to b=2')
            self.b=2
        
        self.xml = xml
        self.xmstar = xmstar
        self.rjaf = rjaf
        self.rsgas = rsgas
        self.rhogas = rhogas
        self.rstar = rstar
        self.rhostar = rhostar
        
        #define useful quantities:
        self.hz = self.H0*np.sqrt(self.Om*(1.+self.z)**3+self.Ol) 
        self.v200 = 10.*self.hz*r200
        
        #prepare the screening radius 
        self.xsc = 0.001
        
        #Get physical value of the anisotropy:
        if not (self.nbeta==2 or self.nbeta==3):

            #avoid unphysical values
            if self.beta>=1.000:
                print('')
                print("WARNING: Model.beta not physical")
                print("for this model beta should be less than 1!")
                print("re-defining maximum value to be 1")
                print('')
                
                self.beta=0.99999
            if self.beta0>=1.000:
                print('')
                print("WARNING: Model.beta0 not physical")
                print("for this model beta0 should be less than 1!")
                print("re-defining maximum value to be 1")
                print('')
                self.beta0=0.99999
            #avoid unphysical values for the case of beta O
            if self.nbeta==5 and self.beta<=-1.000:
                self.beta=-0.99999
        
        self.get_sigmar=np.vectorize(self.__get_sigma)
        self.get_sigmaLOS=np.vectorize(self.__get_sigmaLOS)
        self.get_Np = np.vectorize(self.__numdens)
        
        if (self.nmass == 7):
            self.xsc = ScreenBur(self.r200,self.rs,self.Y2,self.Y1*1e-5,self.H0,self.z,self.Om,self.Ol)
        if (self.nmass == 16):
            rhob = self.xmstar*self.xml/(np.pi*4*self.rjaf**3)
            self.xsc = ScreenTot(self.Y2,self.Y1*1e-5,self.r200,self.rs,self.rhostar,
                                 self.rstar, self.rhogas, self.rsgas, rhob, self.gamma, self.rjaf,
                                 H0 = self.H0, z = self.z,Om= self.Om, Ol = self.Ol)
         
            
    def model_info(self):
    #recap of all the model parameters specified in the class,
    #including cosmology, mass, number density and anisotropy     
        
#                  [1,2,3,4,5,6,7,8,9,10,11,12,13]
        stringmass=['NFW','Burkert','Hernquist','gNFW','DHOST NFW','Chameleon b-NFW', 
                    'Chameleon Burkert','Einasto', 'Chameleon Einasto', 'DHOST b-NFW',
                    'DHOST Burkert', 'DHOST Einasto', 'Chameleon gNFW','Chameleon Isothermal', 
                    'Refracted Gravity', "multicomponent (GR or Chameleon)"]
        if self.nmass>len(stringmass) or self.nmass<=0:
            print('ERROR: invalid mass number')
            return
#                  [1,2,3,4,5,6,7]        
        stringani=['C','ML','OM','T','O','gOM','gT']
        if self.nbeta>len(stringani) or self.nbeta<=0:
            print('ERROR: invalid anisotropy number')
            return        
        
        stringnum=['pNFW', 'pHernquist']
        
        print("*****************************************")
        print(" ")
        print("Cosmology:\n\n H0 = {:.2f} km/s/Mpc\t z = {:.2f}\n Omega_m = {:.2f}\t \t Omega_l = {:.2f}".format(
            self.H0,self.z,self.Om,self.Ol))
        print(" ")
        print("-----------------------------------------")
        print(" ")
        print("mass model: {:s} ".format(stringmass[self.nmass-1]))
        print(" ")
        if(self.nmass<4):
            print("with r200:\t {:.4f} Mpc\n and rs: \t {:.4f} Mpc".format(self.r200,self.rs))
        elif(self.nmass==4):
            print("with r200:\t {:.4f} Mpc\n \t rs:  \t {:.4f} Mpc\n and gamma:\t {:.4f}".format(self.r200,self.rs,self.gamma))
        elif(self.nmass==5):
            print("with r200:\t {:.4f} Mpc\n \t rs:  \t {:.4f} Mpc\n and Y1:\t {:.4f}".format(self.r200,self.rs,self.Y1))
        elif(self.nmass==6):
            print("with r200:\t {:.4f} Mpc\n \t rs:  \t {:.4f} Mpc\n \t Phi:\t {:.4f}e-5\n \t Q:   \t {:.4f}\n and b:   \t {:.4f}".format(self.r200,
                            self.rs,self.Y1,self.Y2,self.b))
        elif(self.nmass==7):
            #modified Burkert
            print("with r200:\t {:.4f} Mpc\n \t rs:  \t {:.4f} Mpc\n \t Phi:\t {:.4f}e-5\n and \t Q:   \t {:.4f}".format(self.r200,
                            self.rs,self.Y1,self.Y2))
            
        elif(self.nmass==8):
            print("with r200:\t {:.4f} Mpc\n \t rs:  \t {:.4f} Mpc\n and \t n:   \t {:.4f}".format(self.r200,
                 self.rs,self.b))            
        elif(self.nmass==9):
            print("with r200:\t {:.4f} Mpc\n \t rs:  \t {:.4f} Mpc\n \t b:\t {:.4f}\n \t Phi:\t {:.4f}e-5\n and \t Q:   \t {:.4f}".format(self.r200,
                 self.rs,self.b,self.Y1,self.Y2))      
        elif(self.nmass==10):
            print("with r200:\t {:.4f} Mpc\n \t rs:  \t {:.4f} Mpc\n \t Y1:\t {:.4f}\n \t Y2:   \t {:.4f}\n and b:   \t {:.4f}".format(self.r200,
                self.rs,self.Y1,self.Y2,self.b))
        elif(self.nmass==11):
            print("with r200:\t {:.4f} Mpc\n \t rs:  \t {:.4f} Mpc\n \t Y1:\t {:.4f}\n and \t Y2:   \t {:.4f}\n".format(self.r200,
                self.rs,self.Y1,self.Y2))         
        elif(self.nmass==12):
            print("with r200:\t {:.4f} Mpc\n \t rs:  \t {:.4f} Mpc\n \t Y1:\t {:.4f}\n \t Y2:   \t {:.4f}\n and n:   \t {:.4f}".format(self.r200,
                self.rs,self.Y1,self.Y2,self.b)) 
        elif(self.nmass==13):
            print(r"with r200:\t {:.4f} Mpc\n \t rs:  \t {:.4f} Mpc\n \t Phi:\t {:.4f}e-5\n \t Q:   \t {:.4f}\n and $\gamma$:   \t {:.4f}".format(self.r200,
                            self.rs,self.Y1,self.Y2,self.gamma))
        elif(self.nmass==14):
            print("with r200:\t {:.4f} Mpc\n \t rs:  \t {:.4f} Mpc\n \t Phi:\t {:.4f}e-5\n \t Q:   \t {:.4f}\n".format(self.r200,
                            self.rs,self.Y1,self.Y2))
            
        elif(self.nmass==15):
            print("with LBCG:\t {:.3e} L_sun\n \t xmstar:  \t {:.3f} \n \t rjaf:     \t {:.3f} Mpc\n \t epsilon0:\t {:.3f}\n \t log(Rhoc):\t{:.4f}\n \t Q:      \t {:.4f}".format(self.xml,
                            self.xmstar,self.rjaf,self.Y1,self.Y2,self.Y3))
        
        elif(self.nmass==16):
            print("with LBCG:\t \t {:.3e} L_sun\n \t xmstar:  \t {:.3f} \n \t rjaf:" \
                  "     \t {:.3f} Mpc\n \t rhogas:\t \t{:.3e} M_sun/Mpc**3 \n" \
                  "   \t rsgas:\t \t {:.4f} Mpc\n \t rhostar:\t {:.4e} Msun/Mpc**3 \n \t rstar: \t  \t {:.4f} Mpc\n" \
                  "\t r200:\t \t {:.3f} Mpc\n \t rs: \t \t {:.4f} Mpc\n \t gamma: \t \t {:.4f}\n"\
                      " Phi:\t \t {:.4f}e-5\n \t \t Q: \t {:.4f}".format(self.xml,
                            self.xmstar,self.rjaf,self.rhogas,self.rsgas, self.rhostar,
                            self.rstar,self.r200,self.rs,self.gamma,self.Y1,self.Y2))
            
        print(" ")
        print("-----------------------------------------")
        print(" ")
        print("anisotropy model: {:s} ".format(stringani[self.nbeta-1]))
        print(" ")
        if(self.nbeta==2 or self.nbeta==3):
            print("with r_beta:\t {:.4f} Mpc".format(self.beta))
        elif(self.nbeta<6):
            print("with beta_inf:\t {:.4f}".format(self.beta))
        else:
            print("with beta_inf:\t {:.4f}\n and beta_0:\t     {:.4f}".format(self.beta,self.beta0))    
        print(" ")
        print("-----------------------------------------")
        print(" ")
        print("Number density model: {:s} ".format(stringnum[self.nnu-1]))
        print(" ")
        print("with rnu:\t {:.4f} Mpc\n".format(self.rnu))
        print("*****************************************")
        
        
    def screen(self):
        if (self.nmass == 7):
            self.xsc = ScreenBur(self.r200,self.rs,self.Y2,self.Y1*1e-5,self.H0,self.z,self.Om,self.Ol)
            print("new screening S/rs = ", self.xsc)
        elif (self.nmass == 16):
            rhob = self.xmstar*self.xml/(np.pi*4*self.rjaf**3)
            self.xsc = ScreenTot(self.Y2,self.Y1*1e-5,self.r200,self.rs,self.rhostar,
                                 self.rstar, self.rhogas, self.rsgas, rhob, self.gamma, self.rjaf,
                                 H0 = self.H0, z = self.z,Om= self.Om, Ol = self.Ol)
            print("new screening S = ", self.xsc)
        

    def __numdens(self,r):
        
        r200 = self.r200
        rc = self.rnu
        t = r/rc
        
        if (self.nnu == 1):
         
             if (t > 1):
                cm1x=np.arccos(1.0/t)
                ft=(1.-1./np.sqrt(abs(t*t-1))*cm1x)/(t*t-1.)
             elif (t < 1):
                cm1x=np.arccosh(1./t)
                ft=(1.-1./np.sqrt(abs(t*t-1))*cm1x)/(t*t-1.)
             else:
                ft=1./3.
    
             c200=r200/rc
             gc=1./(np.log(c200+1)-c200/(c200+1))
             
             return ft*c200*c200/(2.*np.pi)*gc/(r200*r200)
        elif (self.nnu == 2):
            s=r/rc

            if (s > 1):
                xs=np.arccos(1./s)/np.sqrt(s*s-1.)
                fr=((2.+s*s)*xs-3.)/(2.*np.pi*rc*rc*(1.-s*s)**2)
            elif (s < 1):
                xs = np.log((1.+np.sqrt(1.-s*s))/s)/np.sqrt(1.-s*s)
                fr=((2.+s*s)*xs-3.)/(2.*np.pi*rc*rc*(1.-s*s)**2)
            else:
                xs=1.
                fr=2./(15.*np.pi*rc*rc)


            return fr
            
    
        
    #now define the mass profile according to the model
    def get_mass(self,r,expc=True):
        r200=self.r200
        rs=self.rs
        Y1=self.Y1
        Y2=self.Y2
        gamma=self.gamma
        H0=self.H0
        z=self.z
        Om=self.Om
        Ol=self.Ol
        b=self.b
        
        Y3 = self.Y3
        xml =self.xml 
        xmstar = self.xmstar 
        rjaf = self.rjaf
        rsgas = self.rsgas 
        rhogas = self.rhogas 
        rhostar = self.rhostar
        rstar = self.rstar
        
        if self.nmass==1:
            return M_nfw(r, r200, rs, z, H0, Om, Ol)
        
        if self.nmass==2:
            return M_Bur(r, r200, rs, z, H0, Om, Ol)
            
        if self.nmass==3:
            return M_Her(r, r200, rs, z, H0, Om, Ol)

        if self.nmass==4:
            return M_gNFW(r,r200,rs,gamma,z,H0,Om,Ol)

        if self.nmass==5:
            return M_nfw_BH(r, r200, rs, Y1, z, H0, Om, Ol)
        
        if self.nmass==6:
            return M_bnfw_CS(r, r200, rs, Y1, Y2, b, z, H0, Om, Ol,
                             expcutoff=expc)
        if self.nmass==7:
            return M_bur_CS(r, r200, rs, Y1, Y2, z, H0, Om, Ol,
                             xsc = self.xsc, expcutoff=expc)
        if self.nmass==8:
            return M_Ein(r, r200, rs, b, z, H0, Om, Ol)
                            
        if self.nmass==9:
            return M_Ein_CS(r, r200, rs, b, Y1, Y2, z, H0, Om, Ol,
                             expcutoff=expc)
        if self.nmass==10:
            return M_bnfw_BH(r,r200,rs,Y1,b,z,H0,Om,Ol)
        if self.nmass==11:
            return M_Bur_BH(r,r200,rs,Y1,z,H0,Om,Ol)
        if self.nmass==12:
            return M_Ein_BH(r,r200,rs,Y1,b,z,H0,Om,Ol)
        
        if self.nmass==13: #Chameleon gNFW
            return M_gnfw_CS(r, r200, rs, gamma, Y1, Y2, z, H0, Om, Ol,
                             expcutoff=expc)
        if self.nmass==14: #Chameleon gNFW
            return MIso_CS(r, r200, rs, Y1, Y2, z, H0, Om, Ol,
                             expcutoff=expc)
        if self.nmass==15: #Refracted gravity
        #M_RG(x, eps, ldens, Q, xml, xmstar, rjaf, rsgas = 0.16, rhogas = 60.80e13 )
            return M_RG(r, Y1, Y2, Y3, xml, xmstar, rjaf, rhogas, rhostar, 
                        rsgas, rstar)
        
        
        if self.nmass == 16: #TOtal multicomponent GR: gas+BCG+DM+gal
            #define the central density of the Jaffe profile
            rhob = xmstar*xml/(np.pi*4*rjaf**3)
            
            return Mtot(r,Y2,Y1*1e-5,r200,rs,rhostar,rstar,rhogas,
                 rsgas,rhob, gamma,rjaf,xsc =self.xsc,
                    z= z, H0= H0,  Om = Om, Ol  = Ol )
             
            
        print("No other model implemented")
        return 0
        
    
    def get_beta(self,r):
        beta=self.beta
        beta0=self.beta0
        if(self.nbeta==2 or self.nbeta==3):
            rbeta=beta
        else:
            rbeta=self.rs            # NFW or other mass profiles
            if (self.nmass==3 or ((self.nmass==6 or self.nmass==10) 
                                  and self.b==3)):
                rbeta=0.5*self.rs    # Hernquist
            if (self.nmass==2 or self.nmass==7 or self.nmass==11):
                rbeta=1.521*self.rs  # Burker
            
            if(self.nmass==4): #gNFW
                rbeta=(2-self.gamma)*self.rs
        return betaf(r,rbeta,beta,self.nbeta,beta0)
            
    
#function sigma^r to be integrated in the Maximum Likelihood procedure ++++++++
# the integral is made in x=log(r)
    def sr2int(self,r):
        
        rc=self.rnu
        r200=self.r200
           
        t=np.exp(r)
        G=4.302e-9
        #get the mass:
        xm=self.get_mass(t,expc=self.expc)    
        #get the number density profile
        
        if (self.nnu==1):
            c=r200/rc
            gc=1./(np.log(c+1)-c/(c+1))
            xnu=1.0/(t/rc*(1.+t/rc)**2)*c*c*c*gc/4./np.pi/(r200*r200*r200)
        elif (self.nnu==2): 
            xnu=rc/(2.*np.pi*t)/(t+rc)**3        
        
        #get the anisotropy kernel
        anis=self.nbeta
        bec=self.beta
        bec0=self.beta0
        if (anis==1):           
             app=1.e-2                                         
    #
    #     Radial anisotropy
    #
             if(abs(bec-1.0)<app): 
                 sr2int=xnu*xm
    #
    #     Isotropic case
    #
             elif (abs(bec)<app):
                sr2int=xnu*xm/(t*t)
    #
    #     General constant anisotropy
    #
             else:
                sr2int=xnu*xm*t**(2.*bec-2.)
             
    #
    #     Mamon Lokas anisotropy profile
    #
        elif (anis==2):
            sr2int=xnu*xm*(t+bec)/(t*t)
            
    #  OM profile        
        elif(anis==3):
            sr2int=xnu*xm*(t*t+bec*bec)/(t*t)
    #generalised OM profile
        elif(anis==6):
            rm2=self.rs            # NFW or other mass profiles
            if (self.nmass==3):
                rm2=0.5*self.rs    # Hernquist
            if (self.nmass==2 or self.nmass==7 or self.nmass==11):
                rm2=1.521*self.rs  # Burkert
            sr2int=xnu*xm*t**(2.*bec0-2.)*(t*t+rm2*rm2)**(bec-bec0)
            
    #     Simplified Tiret (modified ML profile)
    #
        elif (anis==4):
            rm2=self.rs            # NFW or other mass profiles
            if (self.nmass==3):
                rm2=0.5*self.rs    # Hernquist
            if (self.nmass==2 or self.nmass==7 or self.nmass==11):
                rm2=1.521*self.rs  # Burker
            sr2int=xnu*xm*(t+rm2)**(2.*bec)/(t*t)
            
    #
    #     Simplified Tiret (modified ML profile)
    #
        elif (anis==7):
            rm2=self.rs            # NFW or other mass profiles
            if (self.nmass==3):
                rm2=0.5*self.rs    # Hernquist
            if (self.nmass==2 or self.nmass==7 or self.nmass==11):
                rm2=1.521*self.rs  # Burker
            sr2int=xnu*xm*t**(2*bec0-2.)*(t+rm2)**(2.*(bec-bec0))
    #
    #     modified Tiret (non-zero central anisotropy)
    #
        elif (anis==5):
            rm2=self.rs            # NFW or other mass profiles
            if (self.nmass==3):
                rm2=0.5*self.rs    # Hernquist
            if (self.nmass==2 or self.nmass==7 or self.nmass==11):
                rm2=1.521*self.rs  # Burker
            sr2int=xnu*xm*(t+rm2)**(4.*bec)*t**(-2.-2.*bec)
            
        return G*t*sr2int #integral in dlog(t)    
    
#factor outside the integral of the velocity dispersion profile ***************

    def sr2out(self,t):
        #number density parameters
        rc=self.rnu
        r200=self.r200
        knfit=self.nnu 
        
        #get the anisotropy kernel
        anis=self.nbeta
        bec=self.beta
        bec0=self.beta0     

        if (knfit==1):
            c200=r200/rc
            gc=1./(np.log(c200+1)-c200/(c200+1))
            xnu=1.0/(t/rc*(1.+t/rc)**2)*c200**3*gc/4./np.pi/(r200**3)
        elif (knfit==2): 
            xnu=rc/(2.*np.pi*t)/(t+rc)**3
        #else:
#      BETA MODEL FOR n(r) TO BE IMPLEMENTED **********************************            
#            xnu=-1.d0/np.sqrt(np.pi*1.0)*gammarec(-al+0.5d0)/gammarec(-al+1.0)*
#            al/rc*(1.d0+(t/rc)*(t/rc))**(-0.5d0+al)

        if (anis==1):           
             app=1.0e-2                                         
    #
    #     Radial anisotropy
    #
             if(abs(bec-1.0)<app): 
                 sr2out=1./(t*t)
    #
    #     Isotropic case
    #
             elif (abs(bec)<app):
                sr2out=1.
    #
    #     General constant anisotropy
    #
             else:
                sr2out=t**(-2.*bec)
             
    #
    #     Mamon Lokas anisotropy profile
    #
        elif (anis==2):
            sr2out=1./(t+bec)
            
    #  OM profile        
        elif(anis==3):
            sr2out=1./(t*t+bec*bec)
    #generalised OM profile
        elif(anis==6):
            rm2=self.rs            # NFW or other mass profiles
            if (self.nmass==3):
                rm2=0.5*self.rs    # Hernquist
            if (self.nmass==2 or self.nmass==7 or self.nmass==11):
                rm2=1.521*self.rs  # Burkert
            sr2out=t**(-2.*bec0)*(t*t+rm2*rm2)**(-1.*(bec-bec0))
            
    #
    #     Simplified Tiret (modified ML profile)
    #
        elif (anis==4):
            rm2=self.rs            # NFW or other mass profiles
            if (self.nmass==3):
                rm2=0.5*self.rs    # Hernquist
            if (self.nmass==2 or self.nmass==7 or self.nmass==11):
                rm2=1.521*self.rs  # Burker
            sr2out=(t+rm2)**(-2.*bec)
            
    #
    #     Simplified Tiret (modified ML profile)
    #
        elif (anis==7):
            rm2=self.rs            # NFW or other mass profiles
            if (self.nmass==3):
                rm2=0.5*self.rs    # Hernquist
            if (self.nmass==2 or self.nmass==7 or self.nmass==11):
                rm2=1.521*self.rs  # Burker
            sr2out=t**(-2.*bec0)*(t+rm2)**(-2.*(bec-bec0))
    #
    #     modified Tiret (non-zero central anisotropy)
    #
        elif (anis==5):
            rm2=self.rs            # NFW or other mass profiles
            if (self.nmass==3):
                rm2=0.5*self.rs    # Hernquist
            if (self.nmass==2 or self.nmass==7 or self.nmass==11):
                rm2=1.521*self.rs  # Burker
            sr2out=(t+rm2)**(-4.*bec)*t**(2.*bec)
#
#     Divide by the density profile
#
        return sr2out/xnu
 
    
#Private module of the cluster class 
    def __get_sigma(self,r,rinfinity=70):
        
        #integrate to have sigma^2_r(r)
        
        sintegral=ints.quad(self.sr2int,np.log(r),np.log(2.0*rinfinity),
                        epsabs=1.49e-8, epsrel=1.49e-8)
        return sintegral[0]*self.sr2out(r)
 
#Public module of the cluster class ******************************************   
#    def get_sigma2(self,r,rinfinity=70):
#        #integrate to have sigma^2_r(r)
        
#        return self.__get_sigma(r,rinfinity)   
    

#************************************************+*****************************    
    
    def generate_data(self,lenout,title=None,lmin=-3, savefiles = False):
        r=np.logspace(lmin,1.3,lenout)
        
        sigr = self.get_sigmar(r)
        siglos = self.get_sigmaLOS(r)
        mass = self.get_mass(r)
        yb=np.column_stack((r,np.sqrt(sigr)))
        #xb=np.array(list(vars(self).values())[:-1]).astype(float)
        
        outsigmar = np.column_stack((r,siglos))
        outMass = np.column_stack((r,mass))
        #if title is not None:
            #d=np.array([xb,yb],dtype='object')
           # np.save(title+'.npy',d)
        if (savefiles == True):
            x = datetime.datetime.now()
            distdate='data'+str(x.day)+'_'+str(x.month)+'_'+str(x.year)
            titlesr = 'Sigmar'+str(distdate)
            titleSlos = 'SigmaLOS'+str(distdate)
            titleMass = 'Mass'+str(distdate)
            np.savetxt(titlesr, yb, fmt='%10.3e',header='r [Mpc] \t sigma_r [km/s]')
            np.savetxt(titleSlos, outsigmar, fmt='%10.3e',header='R [Mpc] \t sigma_LOS [km/s]')
            np.savetxt(titleMass, outMass, fmt='%10.3e',header='r [Mpc] \t M(r) [M_sun]')

            
        return yb, outsigmar,outMass
    
    def generate_WL_shear(self,sigmael=0.25,sigmalls=0.005,nrmax=2.9,nbin=10,ngal=30,zm=10):
        r200=self.r200
        rs=self.rs

        gamma=self.gamma
        H0=self.H0
        zc=self.z
        Om=self.Om
        Ol=self.Ol
        bb=self.b
        
        rmax=nrmax*r200
        
        errcheck = False
        if self.nmass==1 or self.nmass==6:
            
            if self.nmass==1:
                bb=2.0
            
            xmin=1.5*Einstein_radius(r200,rs,zc,H0=H0,Om=Om,Ol=Ol,zm=zm,deg=False,
                                model='bNFW',b=bb)
            
            if xmin<=1e-1*r200: #avoid too small values that can give a negative tangential shear
                xmin=1.01e-1*r200
            xl=np.logspace(np.log10(xmin),np.log10(rmax),nbin+1)
            
            angarcmin=xl/np.pi/Da(zc,0,H0,Om,Ol)*(180*60)
            errstat=(sigmael**2)/(np.pi*(angarcmin[1:]**2-angarcmin[:-1]**2)*ngal)
            err_gtl=np.sqrt(errstat+sigmalls**2)
            gtl=gt_ave_bNFW(xl[:-1],r200,rs,bb,zc,H0,Om,Ol,ngal,zm)
            
        elif self.nmass==4:
            xmin=1.5*Einstein_radius(r200,rs,zc,H0,Om,Ol,zm,deg=False,
                                model='gNFW',b=gamma)
            if xmin<=1e-1*r200: #avoid too small values that can give a negative tangential shear
                xmin=1.01e-1*r200
            
            xl=np.logspace(np.log10(xmin),np.log10(rmax),nbin+1)
            
            angarcmin=xl/np.pi/Da(zc,0,H0=H0,Om=Om,Ol=Ol)*(180*60)
            errstat=(sigmael**2)/(np.pi*(angarcmin[1:]**2-angarcmin[:-1]**2)*ngal)
            err_gtl=np.sqrt(errstat+sigmalls**2)
            gtl=gt_ave_bNFW(xl[:-1],r200,rs,bb,zc,H0,Om,Ol,ngal,zm)
            
        else:
            print("Weak Lensing not implemented with these models!")
            xl=np.zeros(nbin+1)
            gtl=xl[:-1]
            err_gtl=xl[:-1]
            errcheck = True
            
        return xl[:-1], gtl, err_gtl, errcheck
  





    def fsigmaLOS(self,tlog):
        """
        Integrand function for the LOS velocity dispersion
        (eq. A1, A16 from Mamon&Lokas 2014)

        Parameters
        ----------
        tlog : FLOAT
            Log of the distance in Mpc.

        Returns
        -------
        FLOAT: the integrand f(x) times x (since the integral is in log)

        """
        
        r200=self.r200
        rs=self.rs
        rc = self.rnu

        xmin = self.xmin
        G = 4.302e-9
        
#        rsvalues = np.array([0.01,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,
#         0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,
#        0.30,0.35,0.40,0.45,0.50,1.00])
        
#        r100values = np.array([1.288,1.307,1.311,1.314,1.318,1.321,1.324,1.326,
#                                1.329,1.332,1.334,1.337,1.339,1.342,1.344,1.347,1.349,1.351,
#                                1.353,1.356,1.358,1.36,1.37,1.38,1.389,1.398,1.406,1.476],dtype = float) 
        t=np.exp(tlog)

# check on value of t and xmin (t must always be > xmin)

        if (xmin >= t):
            return 0.0
# c
# c     physical units
# c
        #gm200=r200*v200*v200

#       Mass profile
        xm = self.get_mass(t)*G

# c     nu(r) from NFW, Hernquist - beta-model to be added
# c
        if (self.nnu==1):
            c=r200/rc
            gc=1./(np.log(c+1)-c/(c+1))
            xnu=1.0/(t/rc*(1.+t/rc)**2)*c*c*c*gc/4./np.pi/(r200*r200*r200)
        elif (self.nnu==2): 
            xnu=rc/(2.*np.pi*t)/(t+rc)**3   
        #else:
        #    xnu=-1.0/np.sqrt(np.pi*1.0)* \
        #    gammarec(-al+0.5d0)/gammarec(-al+1.d0)* \
        #   al/rc*(1.d0+(t/rc)*(t/rc))**(-0.5d0+al)
        
        #integration variable
        u=t/xmin
        
        #get the anisotropy kernel
        anis=self.nbeta
        bec=self.beta
        bec0=self.beta0
        
# c
# c     Constant anisotropy
# c
        if (anis == 1):
            app=1.0e-2
# c
# c     Radial
# c
            if(abs(bec-1.0) < app):
                xk = np.pi/4.0*u-0.5*np.sqrt(1.0-1.0/(u*u))-0.50*u* \
                np.arcsin(1.0/u)
# c
# c     Isotropic
# c
            elif (abs(bec)< app):
                xk=np.sqrt(1.0-1.0/(u*u))
# c
# c     Negative integer
# c
            elif (abs((bec-round(bec))/bec) < app and bec < 0):
                xk1=np.sqrt(u*u-1.0)/u**(1.0-2.0*bec)
                xk2=0.0
                xk3=0.0
                k = np.arange(0,round(-bec))
                xk2 = np.sum(binom(round(-bec),k)* \
                (u*u-1.0)**k/(2.0*k+1.0))

                k2 = np.arange(0,round(-bec-1))
                xk3= np.sum(binom(round(-bec-1),k2)* \
                    (u*u-1.0)**k2/(2.0*k2+1.0))

                xk=xk1*(xk2-bec*xk3)
# c
# c     +1/2 or -1/2
            elif ((abs(bec-0.5) < app) or (abs(bec+0.50) < app)):
                xk=u**(2.0*bec-1.0)*np.arccosh(u) - \
                bec*np.sqrt(1.0-1.0/(u*u))
# c
# c     Generic constant anisotropy, using
# c     expression in Mamon+Lokas (2006, MNRAS, 370, 1582)
# c
            else:
                xxx=1.0/(u*u)
                xk1=np.sqrt(1.0-xxx)/(1.0-2.0*bec)
                xk2=np.sqrt(np.pi)/2.0*gamma(bec-0.50)/gamma(bec)
                xk3=(1.50-bec)*u**(2.0*bec-1.0)
                xk4=1.0-betainc(bec+0.5,0.5, xxx)
                xk=xk1+xk2*xk3*xk4

            fa=2.*xk*xm*xnu/t
# c
# c     cbe is the radius in eq.(61) of Mamon & Lokas 2005 
# c     (Osipkov-Merritt anisotropy)
# c
        elif (anis == 3):
          
          ua=bec/xmin
          xk=(ua*ua+0.50)/(ua*ua+1.0)**(1.50)*(u*u+ua*ua)/u* \
              np.arctan(np.sqrt((u*u-1.0)/(ua*ua+1.0)))- \
              0.50/(ua*ua+1.0)*np.sqrt(1.0-1.0/(u*u))
 
          fa=2.*xk*xm*xnu/t
# c
# c     cbe is the radius ra in eq.(60) of Mamon & Lokas 2005
# c     in this case (mamon e Lokas profile)
        elif (anis == 2):
            ua=bec/xmin
            
            if (abs(ua-1.0) < 1.0e-5):
                xk=(1.0+1.0/u)*np.cosh(u)- \
                (8.0/u+7.0)/6.0*np.sqrt((u-1.0)/(u+1.0))
            else:
                if (ua > 1.0):
                    ccc=np.arccosh((ua*u+1.0)/(u+ua))
                    sss=1.0
                else:
                    ccc=np.arccos((ua*u+1.0)/(u+ua))
                    sss=-1.0

                fff=ua*(ua**2-0.50)/(abs(ua**2.0-1.0))**1.50* \
                (1.0+ua/u)*ccc

                xk=0.50/(ua*ua-1.0)*np.sqrt(1.0-1.0/(u*u))+ \
                (1.0+ua/u)*np.arccosh(u)-sss*fff

            fa=2.*xk*xm*xnu/t       

# c
# c     if the anisotropy profile is a simplified Wojtak,
# c     simplified Tiret (modified ML), generalized Tiret,
# c     generalized OM or modifed Tiret (Opposite), there are no 
# c     analytical solutions to provide K(r,ra) in eq.(A8) 
# c     of ML05b, so we use the expression that includes
# c     sigma_r for which we have a spline interpolation
# c
        else:

# c     choose between simplified Wojtak, simpl. Tiret,
# c     a model similar to Tiret (modified Tiret)
# c     with non-zero central anisotropy, and Hansen+Moore

            if (anis == 6): #generalized OM

                rm2=rs                     # NFW or other mass profiles
                if (self.nmass==3):
                    rm2=0.5*self.rs    # Hernquist
                if (self.nmass==2 or self.nmass == 11):
                    rm2=1.521*self.rs  # Burkert
                if (self.nmass == 4): 
                    rm2=(2-self.gamma)*rs   #gNFW
#             !avoid unphyiscal values 
                if (rm2 < 0):
                    rm2=1e-4
         
                banis=bec0+(bec-bec0)*t*t/(t*t+rm2*rm2)        
                if (bec >= 1.0):
                    bec=0.9990  # avoid unphysical values            
         
            elif (anis == 4):
                rm2=rs                     # NFW or other mass profiles
                if (self.nmass==3):
                    rm2=0.5*self.rs    # Hernquist
                if (self.nmass==2 or self.nmass == 11):
                    rm2=1.521*self.rs  # Burkert
                if (self.nmass == 4): 
                    rm2=(2-self.gamma)*rs   #gNFW
#             !avoid unphyiscal values 
                if (rm2 < 0):
                    rm2=1e-4

            
                banis=(bec)*t/(t+rm2)
            
            elif (anis == 7):
                rm2=rs                     # NFW or other mass profiles
                if (self.nmass==3):
                    rm2=0.5*self.rs    # Hernquist
                if (self.nmass==2 or self.nmass == 11):
                    rm2=1.521*self.rs  # Burkert
                if (self.nmass == 4): 
                    rm2=(2-self.gamma)*rs   #gNFW
#             !avoid unphyiscal values 
                if (rm2 < 0):
                    rm2=1e-4
         
                banis=bec0+(bec-bec0)*t/(t+rm2)        
#             if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values 
        
            elif (anis == 5):
                rm2=rs                     # NFW or other mass profiles
                if (self.nmass==3):
                    rm2=0.5*self.rs    # Hernquist
                if (self.nmass==2 or self.nmass == 11):
                    rm2=1.521*self.rs  # Burkert
                if (self.nmass == 4): 
                    rm2=(2-self.gamma)*rs   #gNFW
#             !avoid unphyiscal values 
                if (rm2 < 0):
                    rm2=1e-4
            
                banis=(bec)*(t-rm2)/(t+rm2)
#                else   ! Hansen+Moore
#             rhm=rc ! use nu(r)
# cc            rhm=rs ! use rho(r)
#             b=ahm-bhm*(rhm+3.*t)/(rhm+t)            

            sr = self.get_sigmar(t)

            fa=2.*xnu*sr*(1.-banis*(xmin/t)**2)*t/np.sqrt(t*t-xmin*xmin)


        return t*fa 
    
    #fsigmaLOS = np.vectorize(fsigmaLOS)
        
    def __get_sigmaLOS(self,r,rinfinity=70):
        x = np.log(r)
        xinf = np.log(2.0*rinfinity)
        self.xmin = r
        
        sintegral=ints.quad(self.fsigmaLOS,x,xinf,
            epsabs=1.49e-4, epsrel=1.49e-4)
        return np.sqrt(sintegral[0]/self.get_Np(r))
        



#===========================================================================
        

##### per integrare!
     #  errrel=0.005d0
     #  errabs=0.

     #  do i=1,50
     #     xx1 = dlog(rlow)+dfloat(i-1)*
     # &    dlog(1.001*rinfinity/rlow)/dfloat(50-1)
     #     xx2=dlog(2.*rinfinity)
     #     xmin=dexp(xx1)
         
     #     call dgaus8 (fa,xx1,xx2, errrel, ris2n, IERR)

     #     if ((knfit.eq.1).or.(nrc.eq.-1.and.kmp.eq.1)) then 
     #        pnr=sigmar1(xmin)
     #     elseif ((knfit.eq.2).or.(nrc.eq.-1.and.kmp.eq.2)) then
     #        pnr=sigmar2(xmin)
     #     else
     #        pnr=sigmar3(xmin)
     #     endif
     #     ris=dsqrt(ris2n/pnr)
     #     if (ris.eq.ris) write(ius,*) xmin,ris
     #  enddo
    


# !*********************** BCG *******************************************
#       fbarr=0.0
#       if (kbary.eq.1) fbarr=fbars(t)      
      
#       if (kbcg.eq.1) then

#        !GR setup
#        if(kmp.eq.10.or.kmp.le.6) then 
        
#         fjaffe=grav*(xmstar*xlumbcg*t/rjaf/(1.+t/rjaf)+fbarr)/gm200
#         xm=xm+fjaffe !inclusion of gas,galaxies and bcg
#         !Chameleon setup
#        elseif (kmp.eq.11.or.kmp.eq.12.or.kmp.eq.9) then 
#         dpj=dphiJ(t/rjaf) !Include the chameleon field of the BCG
#        fjaffe=grav*(xmstar*xlumbcg*(t/rjaf/(1.+t/rjaf)+dpj)+fbarr)/gm200
#         xm=xm+fjaffe !inclusion of gas,galaxies and bcg
        
#         !Refracted gravity setup
#         !Hernquist mass for the BCG
#        elseif (kmp.eq.13) then
#         !the mass is all given by the BCG and gas+baryons, no DM
        
#         !refracted gravity with the function of Valeria
#         xm= grav*rgmass(t)/gm200

#         !Covariant formulation. It is not good for the moment
        
#         !xm=grav*(xmstar*xlumbcg*t**2/(t+rjaf)**2+fbarr)/fmodRG(t)/gm200
        
#        elseif (kmp.eq.14) then
#         !Boson star + BCG
#         fjaffe=grav*(xmstar*xlumbcg*t/rjaf/(1.+t/rjaf)+fbarr)/gm200
        
#         xt = t/rs
#         if (xt.le.2) then 
#          xmer = dtanh((t/rs)**tmass) + 1.0e-3*dA3*(t/rs)**screen
#         else
#          xmer =  dtanh((2.0d0)**tmass) + 1.0e-3*dA3*(2.0d0)**screen
#         endif 
        
#          xm = xmer+fjaffe
         
#        endif  
        
#       endif     
# !***********************************************************************


# c

# c
# c     Constant anisotropy
# c
#       if (kani.eq.0) then
# c
# c     Convert from cbe=beta' to bec=beta 
# c     in the case of constant anisotropy
# c
#          bec=1.-1./(cbe*cbe)
#          if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values
#          app=1.d-2
# c
# c     Radial
# c
#          if(dabs(bec-1.d0).lt.app) then
#             xk=pig/4.d0*u-0.5d0*dsqrt(1.d0-1.d0/(u*u))-0.5d0*u*
#      *           dasin(1.d0/u)
# c
# c     Isotropic
# c
#          elseif (dabs(bec).lt.app) then
#             xk=dsqrt(1.d0-1.d0/(u*u))
# c
# c     Negative integer
# c
#          elseif (dabs((bec-nint(bec))/bec).lt.app.and.bec.lt.0) then
#             xk1=dsqrt(u*u-1.d0)/u**(1.d0-2.d0*bec)
#             xk2=0.d0
#             xk3=0.d0
#             do k=0,nint(-bec)
#                xk2 = xk2 + dbinom(nint(-bec),k)*
#      &              (u*u-1.d0)**k/(2.d0*k+1.d0)
#             enddo
#             do k=0,nint(-bec-1)
#                xk3=xk3+dbinom(nint(-bec-1),k)*
#      &              (u*u-1.d0)**k/(2.d0*k+1.d0)
#             enddo
#             xk=xk1*(xk2-bec*xk3)
# c
# c     +1/2 or -1/2
# c
#          elseif ((dabs(bec-0.5d0).lt.app).or.
#      &           (dabs(bec+0.5d0).lt.app)) then
#             xk=u**(2.d0*bec-1.d0)*dacosh(u)-
#      &           bec*dsqrt(1.d0-1.d0/(u*u))
# c
# c     Generic constant anisotropy, using
# c     expression in Mamon+Lokas (2006, MNRAS, 370, 1582)
# c
#          else
#             xxx=1.d0/(u*u)
#             xk1=dsqrt(1.d0-xxx)/(1.d0-2.d0*bec)
#             xk2=sqrt(pig)/2.d0*gammarec(bec-0.5d0)/gammarec(bec)
#             xk3=(1.5d0-bec)*u**(2.d0*bec-1.d0)
#             xk4=1.d0-betairec(xxx)
#             xk=xk1+xk2*xk3*xk4
#          endif
#          fa=2.*xk*xm*xnu/t
# c
# c     cbe is the radius in eq.(61) of Mamon & Lokas 2005 
# c     (Osipkov-Merritt anisotropy)
# c
#       elseif (kani.eq.2) then
#          ua=cbe/xmin
#          xk=(ua*ua+0.5d0)/(ua*ua+1.d0)**(1.5d0)*(u*u+ua*ua)/u*
#      &        datan(dsqrt((u*u-1.d0)/(ua*ua+1.d0)))-
#      &        0.5d0/(ua*ua+1.d0)*dsqrt(1.d0-1.d0/(u*u))
#          fa=2.*xk*xm*xnu/t
# c
# c     cbe is the radius ra in eq.(60) of Mamon & Lokas 2005
# c     in this case

#       elseif (kani.eq.1) then
#          ua=cbe/xmin
#          if (dabs(ua-1.d0).lt.1.d-5) then
#             xk=(1.d0+1.d0/u)*dcosh(u)-
#      &           (8.d0/u+7.d0)/6.d0*dsqrt((u-1.d0)/(u+1.d0))
#          else
#             if (ua.gt.1.d0) then
#                ccc=dacosh((ua*u+1.d0)/(u+ua))
#                sss=1.d0
#             else
#                ccc=dacos((ua*u+1.d0)/(u+ua))
#                sss=-1.d0
#             endif
#             fff=ua*(ua**2d0-0.5d0)/(dabs(ua**2.d0-1.d0))**1.5d0*
#      &           (1.d0+ua/u)*ccc
#             xk=0.5d0/(ua*ua-1.d0)*dsqrt(1.d0-1.d0/(u*u))+
#      &           (1.d0+ua/u)*dacosh(u)-sss*fff
#          endif
#          fa=2.*xk*xm*xnu/t
# c
# c     if the anisotropy profile is a simplified Wojtak,
# c     simplified Tiret (modified ML), generalized Tiret,
# c     generalized OM or modifed Tiret (Opposite), there are no 
# c     analytical solutions to provide K(r,ra) in eq.(A8) 
# c     of ML05b, so we use the expression that includes
# c     sigma_r for which we have a spline interpolation
# c
#       else

# c     choose between simplified Wojtak, simpl. Tiret,
# c     a model similar to Tiret (modified Tiret)
# c     with non-zero central anisotropy, and Hansen+Moore

#          if (kani.eq.3) then
#             b=t**1.5/(t**1.5+cbe**1.5)
            
#          elseif (kani.eq.21) then
#             rm2=rs                     ! NFW or other mass profiles
#             if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
#             if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
#             if (kmp.eq.10) rm2=(2-tmass)*rs !gNFW
#             !avoid unphyiscal values 
#             if (rm2.lt.0) rm2=1e-4
            
#             becin=1.-1./(cbe*cbe)
#             bec0=1.-1./(cbe0*cbe0)
#             if (becin.ge.1.) becin=0.999d0  ! avoid unphysical values
#             if (bec0.ge.1.) bec0=0.999d0  ! avoid unphysical values
         
#             b=bec0+(becin-bec0)*t*t/(t*t+rm2*rm2)        
#             if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values            
         
#          elseif (kani.eq.4) then
#             rm2=rs              ! NFW or other mass profiles
#             if (kmp.eq.2) rm2=0.5*rs ! Hernquist
#             if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
#             if (kmp.eq.10) rm2=(2-tmass)*rs !gNFW
#             !avoid unphyiscal values 
#             if (rm2.lt.0) rm2=1e-4
            
#             b=(1.-1./(cbe*cbe))*t/(t+rm2)
            
#          elseif (kani.eq.41) then
#             rm2=rs                     ! NFW or other mass profiles
#             if (kmp.eq.2) rm2=0.5*rs   ! Hernquist
#             if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
#             if (kmp.eq.10) rm2=(2-tmass)*rs !gNFW
#             !avoid unphyiscal values 
#             if (rm2.lt.0) rm2=1e-4
            
#             becin=1.-1./(cbe*cbe)
#             bec0=1.-1./(cbe0*cbe0)
#             if (becin.ge.1.) becin=0.999d0  ! avoid unphysical values
#             if (bec0.ge.1.) bec0=0.999d0  ! avoid unphysical values
         
#             b=bec0+(becin-bec0)*t/(t+rm2)        
#             if (bec.ge.1.) bec=0.999d0  ! avoid unphysical values 
        
#          elseif (kani.eq.5) then
#             rm2=rs              ! NFW or other mass profiles
#             if (kmp.eq.2) rm2=0.5*rs ! Hernquist
#             if (kmp.eq.4.or.kmp.eq.11) rm2=1.521*rs ! Burkert
#             if (kmp.eq.10) rm2=(2-tmass)*rs !gNFW
#             !avoid unphyiscal values 
#             if (rm2.lt.0) rm2=1e-4
            
#             b=(1.-1./(cbe*cbe))*(t-rm2)/(t+rm2)
#          else   ! Hansen+Moore
#             rhm=rc ! use nu(r)
# cc            rhm=rs ! use rho(r)
#             b=ahm-bhm*(rhm+3.*t)/(rhm+t)            
#          endif

# c     interpolate sigma_r

#          rjl(1)=tlog
# c         call icsevu(xris,yris,ninterp,csr,icsr,rjl,yrjl,1,ier)
#          call SPLINE_CUBIC_VAL(ninterp,xris,yris,ypp2,rjl,yrjl,
#      &  ydev1,ydev2)
#          sr=dexp(yrjl(1))
#          fa=2.*xnu*sr**2.*(1.-b*(xmin/t)**2)*t/dsqrt(t*t-xmin*xmin)
#       endif

# ccc checking constant beta solution using double integral rather than kernel

# ccc      fa1=fa
# ccc      if (kani.eq.0) then
# ccc         rjl(1)=tlog
# ccc         call icsevu(xris,yris,ninterp,csr,icsr,rjl,yrjl,1,ier)
# ccc         sr=dexp(yrjl(1))
# ccc         fa=2.*xnu*sr**2.*(1.-bec*(xmin/t)**2)*t/dsqrt(t*t-xmin*xmin)
# ccc         if (abs(1.-fa/fa1).gt.0.1)  write(*,915) xmin,t,fa1/fa,bec
# ccc 915     format(5(2x,f6.3))
# ccc      endif

# ccc remove the above part after check done

#       fa=t*fa   ! integral in dlog

#       return
#       end


#Auxiliary functions
#****************************************************************************** 
# def fgen(x,r200,rs,rhog,rg,M,rj,Q,phinf,H0=70,z=0.3,Om=0.3,Ol=0.7):
#     cl=3e5
#     #limit=0.001 #limiting value for the screening radius
#     G=4.302e-9
#     c200=r200/rs

#     Hz2=H0**2*(Om*(1+z)**3+Ol)
    
    
#     M200=100*Hz2/G*r200**3


#     fac200=np.log(1+c200)-c200/(1+c200) 

         
#     rhos=M200/(4*np.pi*rs**3*fac200)   
    
#     rhob=M/(np.pi*4*rj**3)

        
#    # res = 4*cl**2*phinf/(8*np.pi*G*Q) + (4*rg**2*rhog+(rhob*rj**3)/(np.pi*(rj+x))
#    #                                - (4*rho0DM*rs**3)/((rs+x))) \
#    #     + (rhob*rj**2*np.log(x/(rj+x))+2*np.pi*rg**2*rhog*np.log(rg**2+x**2))/np.pi
        
#     res = phinf+((Q * ((rhob * rj**3)/(rj+x)-(rhos* rs**3)/(rs+x)-(rg**3*rhog)/np.sqrt(rg**2+x**2))) \
#                  +(Q*rhob* rj**2*np.log(x/(rj+x)) ) )*8*np.pi*G/cl**2
    
#     return res

# # This function compute the screening radius for aLL THE PROFILES
# def ScreenGen(r200,rs,rhog,rg,M,rj,Q,tmass,H0=70,z=0.3,Om=0.3,Ol=0.7):
#     phinf = tmass*1e-5
#     def fgen2(y):
#         return fgen(y,r200,rs,rhog,rg,M,rj,Q,phinf,H0,z,Om,Ol)
#     limit=1e-6
#     xsc=limit
#     x0=limit   
    
#     if(fgen2(x0)<0):
#         for i in range(0,4000):
#             xx=x0+0.01*i
            
            
#             if (fgen2(x0)*fgen2(xx) < 0) :
#                 print(xx)
#                 sol=opt.root(fgen2,xx)
            
#                 if sol.success:
#                     xsc=sol.x[0]
#                 else:
#                     xsc=xx                
#                     break
#                 x0=xx
#         if (xsc==limit):
#               xsc=xx
#     return xsc 

# ScreenGen=np.vectorize(ScreenGen) 



# def dphiGen(x,r200,rs,rhog,rg,M,rj,Q,phinf,xsc = 0.001, expcutoff = False,
#             H0=70,z=0.3,Om=0.3,Ol=0.7):
#     clight=3e5
#     limit=1e-6 #limiting value for the screening radius
#     G=4.302e-9

#     c200=r200/rs

#     Hz2=H0**2*(Om*(1+z)**3+Ol)
    
    
#     M200=100*Hz2/G*r200**3


#     fac200=np.log(1+c200)-c200/(1+c200) 

         
#     rhos=M200/(4*np.pi*rs**3*fac200)  

    
#     Mpl2 = (clight**2)/(8*np.pi*G)
    
#     rhob=M/(np.pi*4*rj**3)    
#     #eponential cutoff:
#     dminf = np.exp(-x*rs*np.sqrt(1e-7*Q/phinf))
#     if expcutoff == False:
#         dminf = 1.0
    
#     if (xsc>limit):
    
#         Ctot = Q*((rhob*rj**4)/(rj + xsc) - (rhos*rs**4)/(rs + xsc) + 
#                (rg**3*rhog*xsc)/np.sqrt(rg**2 + xsc**2) - rg**3*rhog*np.arcsinh(xsc/rg) - 
#                rhos* rs**3*np.log(rs + xsc))/Mpl2


   

#         if (x<xsc):
#              out=0
#         else:
#              out= Ctot + (Q*(-((rhob* rj**4)/(rj + x)) + (rhos* rs**4)/(rs + x) - 
#                              (rg**3 * rhog * x)/np.sqrt(rg**2 + x**2)))/Mpl2 +  \
#             Q*(rg**3 * rhog *np.arcsinh(x/rg) + rhos* rs**3* np.log(rs + x))/Mpl2
             
             

#     else:
#         xsc =limit
#         Ctot = Q*((rhob*rj**4)/(rj + xsc) - (rhos*rs**4)/(rs + xsc) + 
#                (rg**3*rhog*xsc)/np.sqrt(rg**2 + xsc**2) - rg**3*rhog*np.arcsinh(xsc/rg) - 
#                rhos* rs**3*np.log(rs + xsc))/Mpl2

#         if (x<xsc):
#              out=0
#         else:
#              out= Ctot + (Q*(-((rhob* rj**4)/(rj + x)) + (rhos* rs**4)/(rs + x) - 
#                              (rg**3 * rhog * x)/np.sqrt(rg**2 + x**2)))/Mpl2 +  \
#             Q*(rg**3 * rhog *np.arcsinh(x/rg) + rhos* rs**3* np.log(rs + x))/Mpl2        
        
#          #test!
#     #     #Cb= 1/4*(-4*phinf + Facda*np.pi) NO! Viene negativa la derivata 
#     #     out=Facda/(4)*(-2*np.arctan(x) + 
#     #         2*np.log(1 + x)+np.log(1 + x**2)) #+ Cb #/x**2
        
#     return out*Q*clight**2/G*dminf
# dphiGen=np.vectorize(dphiGen,otypes=[float],excluded=[1,2,3,4]) 



# def PhiGen(x,r200,rs,rhog,rg,M,rj,Q,phinf,xsc = 0.001, expcutoff = False,
#             H0=70,z=0.3,Om=0.3,Ol=0.7):
#     clight=3e5
#     limit=1e-6 #limiting value for the screening radius
#     G=4.302e-9

#     c200=r200/rs

#     Hz2=H0**2*(Om*(1+z)**3+Ol)
    
    
#     M200=100*Hz2/G*r200**3


#     fac200=np.log(1+c200)-c200/(1+c200) 

         
#     rhos=M200/(4*np.pi*rs**3*fac200)  

    
#     Mpl2 = (clight**2)/(8*np.pi*G)
    
#     rhob=M/(np.pi*4*rj**3)    
#     #eponential cutoff:
#     dminf = np.exp(-x*rs*np.sqrt(1e-7*Q/phinf))
#     if expcutoff == False:
#         dminf = 1.0
    
#     if (xsc>limit):
    
#         Ctot = Q*((rhob*rj**4)/(rj + xsc) - (rhos*rs**4)/(rs + xsc) + 
#                (rg**3*rhog*xsc)/np.sqrt(rg**2 + xsc**2) - rg**3*rhog*np.arcsinh(xsc/rg) - 
#                rhos* rs**3*np.log(rs + xsc))/Mpl2


#         if (x<xsc):
#              out=0
#         else:
#              out= (-Ctot*Mpl2+ Q* rhob * rj**3-Q * rhos * rs**3+Mpl2* phinf* x-
#               Q* rg**3 *rhog*np.arcsinh(x/rg)+Q*rhob*rj**2*x*np.log(x/(rj+x))-
#               Q*rhos*rs**3*np.log(rs+x))/(Mpl2*x)
             
             
             

#     else:
#         xsc =limit
#         Ctot = Q*((rhob*rj**4)/(rj + xsc) - (rhos*rs**4)/(rs + xsc) + 
#                (rg**3*rhog*xsc)/np.sqrt(rg**2 + xsc**2) - rg**3*rhog*np.arcsinh(xsc/rg) - 
#                rhos* rs**3*np.log(rs + xsc))/Mpl2


   

#         if (x<xsc):
#              out=0
#         else:
#              out= (-Ctot*Mpl2+ Q* rhob * rj**3-Q * rhos * rs**3+Mpl2* phinf* x-
#               Q* rg**3 *rhog*np.arcsinh(x/rg)+Q*rhob*rj**2*x*np.log(x/(rj+x))-
#               Q*rhos*rs**3*np.log(rs+x))/(Mpl2*x)
        
#          #test!
#     #     #Cb= 1/4*(-4*phinf + Facda*np.pi) NO! Viene negativa la derivata 
#     #     out=Facda/(4)*(-2*np.arctan(x) + 
#     #         2*np.log(1 + x)+np.log(1 + x**2)) #+ Cb #/x**2
        
#     return out

# PhiGen=np.vectorize(PhiGen,otypes=[float],excluded=[1,2,3,4]) 

############ MODIFIED GRAVITY TOTAL MASS PROFILE ##############################


# def Masstot_GC(x,r200,rs,r2b,rb,M,rj,tmass,Q,H0=70,z=0.3,Om=0.3,Ol=0.7,expcutoff=False):
#     cl=3e5
#     limit=0.001 #limiting value for the screening radius
#     dminf=1

    
#     G=4.302e-9
#     c200=r200/rs
#     c200b=r2b/rb
#     Hz2=H0**2*(Om*(1+z)**3+Ol)
#     phinf=tmass*1e-5
    
#     if expcutoff==True:
#         dminf=np.exp(-x*np.sqrt(1e-7*Q/phinf))
    
#     M200=100*Hz2/G*r200**3
#     M200b=100*Hz2/G*r2b**3

#     fac200=np.log(1+c200)-c200/(1+c200) 
#     fac200b=np.log(1+c200b)-c200b/(1+c200b)
         
#     rho0DM=M200/(4*np.pi*rs**3*fac200)   
#     rho0bar=M200b/(4*np.pi*rb**3*fac200b) 
#     rho0BCG=M/(np.pi*4*rj**3)
    
#     #Find the screening radius
#     xsc=ScreenGen(r200,rs,r2b,rb,M,rj,Q,phinf,H0,z,Om,Ol)
#     Ctot=-8*G*np.pi*Q*(((rb**3*rho0bar)/(rb + xsc)) - (rho0BCG*rj**3)/(rj + xsc) + 
#         (rho0DM*rs**3)/(rs + xsc) +rb**3*rho0bar*np.log(rb+xsc) + rs**3*rho0DM*np.log(rs+xsc))/cl**2
    
#     #screened case
#     if (xsc>limit):
#         if (x<xsc):
#             out=0
#         else:
            
#             out=Ctot+8*G*np.pi*Q*(((rb**3*rho0bar)/(rb + x)) - (rho0BCG*rj**3)/(rj + x) + 
#                 (rho0DM*rs**3)/(rs + x) +rb**3*rho0bar*np.log(rb+x) + rs**3*rho0DM*np.log(rs+x))/cl**2 #/x**2

#     #un-screened case
#     else:
        
#         Cce=-8*G*np.pi*Q*(((rb**3*rho0bar)/(rb)) - (rho0BCG*rj**3)/(rj) + 
#             (rho0DM*rs**3)/(rs) +rb**3*rho0bar*np.log(rb) + rs**3*rho0DM*np.log(rs))/cl**2
        
#         out=Cce+8*G*np.pi*Q*(((rb**3*rho0bar)/(rb + x)) - (rho0BCG*rj**3)/(rj + x) + 
#             (rho0DM*rs**3)/(rs + x) +rb**3*rho0bar*np.log(rb+x) + rs**3*rho0DM*np.log(rs+x))/cl**2 #/x**2
    
    
#     xD=x/rs
#     xJ=x/rj
#     xB=x/rb
#     MassDM=M200/fac200*(np.log(1+xD)-xD/(1+xD))
#     MassB=M200b/fac200b*(np.log(1+xB)-xD/(1+xB))
#     MassJ=M*xJ/(1+xJ)
    
    
#     return out*Q*cl**2/G*dminf+ MassDM+MassB+MassJ
    

# Masstot_GC=np.vectorize(Masstot_GC)
