# from __future__ import division

#-------------------------------------
#
#        Started at 10/18/2017 (YuE)
# 
#-------------------------------------

import os, sys
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
from matplotlib import ticker
from matplotlib import markers
import matplotlib as mpl

from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D

import scipy.integrate as integrate
from scipy.integrate import quad, nquad, dblquad

from scipy.constants import pi

#
# All physical constants have its dimension in units in the system CI.
# This code uses units in the system CGS!
#
from scipy.constants import speed_of_light as clight
from scipy.constants import epsilon_0 as eps0
from scipy.constants import mu_0 as mu0
from scipy.constants import elementary_charge as qe
from scipy.constants import electron_mass as me
from scipy.constants import proton_mass as mp
from scipy.constants import Boltzmann as kB

#========================================================
#
# This code compairs two approaches: "classical" (from [1]) and
# "magnus" (from [2]).
#
# For "classical" approach the magnetized interaction between ion
# and electron is considered for ion velocities V_i > rmsTrnsvVe.
#
# References:
#
# [1] I. Meshkov, A. Sidorin, A. Smirnov, G. Trubnikov.
#    "Physics guide of BETACOOL code. Version 1.1". C-A/AP/#262, November
#    2006, Brookhaven National Laboratory, Upton, NY 11973.
# [2] 
#
#========================================================

pi=3.14159265358 

#
# Physical constants:
#
m_e=9.10938356e-28          # electron mass, g
m_p=1.672621898e-24         # electron mass, g
q_e=4.803204673e-10         # electron charge, sqrt(g*cm^3/sec^2)
cLight=2.99792458e10        # speed of light, cm/sec
eVtoErg=1.6021766208e-12    # 1 eV = 1.6...e-12 erg
CtoPart=2.99792458e9        # 1 C = 1 A*sec = 2.9...e9 particles

# 
# Electron beam parameters:
#
Ekin=3.0e4                             # kinetic energy, eV
curBeam=0.5                            # current density, A/cm^2
dBeam=3.0                              # beam diameter, cm
angSpread=3.0                          # angular spread, mrad
alpha=0.                               # angle between V_i and magnetic field
trnsvT=0.5                             # transversal temperature, eV
longT=0.0                              # longitudinal temperature, eV (2.0e-4)
nField=1                               # number ov values  of the magnetic field
fieldB=np.zeros(nField)                # 
fieldB[0]=3.e3                         # magnetic field, Gs
omega_p=1.0e9                          # plasma frequency, 1/sec
n_e=omega_p**2*m_e/(4.*pi*q_e**2)        # plasma density, 3.1421e+08 cm-3

n_e1=8.e7                              # plasma density, cm-3
omega_p1=np.sqrt(4.*pi*n_e1*q_e**2/m_e) # plasma frequency, 5.0459e+08 1/s  
#
# Cooling system parameter:
#
coolLength=150.0        # typical length of the coolong section, cm

# 
# Calculated parameters of the electron beam:
#
V0=np.sqrt(2.*Ekin*eVtoErg/m_e)                # average velocity, cm/s
rmsTrnsvVe=np.sqrt(2.*trnsvT*eVtoErg/m_e)      # RMS transversal velocity, cm/s
rmsLongVe=np.sqrt(2.*longT*eVtoErg/m_e)        # RMS longitudinal velocity, cm/s
# dens=curBeam*CtoPart/V0                      # density, 1/cm^3
# omega=np.sqrt(4.*pi*dens*q_e**2/m_e)          # plasma frequency, 1/s
cyclFreq=q_e*fieldB/(m_e*cLight)               # cyclotron frequency, 1/s
rmsRoLarm=rmsTrnsvVe*cyclFreq**(-1)            # RMS Larmor radius, cm
dens=omega_p**2*m_e/(4.*pi*q_e**2)               # density, 1/cm^3
likeDebyeR=(3./dens)**(1./3.)                  # "Debye" sphere with 3 electrons, cm

coolPassTime=coolLength/V0                     # time pass through cooling section, cm

powV0=round(np.log10(V0)) 
mantV0=V0/(10**powV0) 
# powV0=str(powV0) 
#
# Formfactor ffForm for friction force:
#
# ffForm = 2*pi*dens*q_e**4/(m_e*V0**2)=
#        = 0.5*omega_p**2*q_e**2/V0**2 
#
# Dimension of ffForm is force:  g*cm/sec**2=erg/cm 
#
# 1 MeV/m = 1.e6*eVtoErg/100. g*cm/sec**2 = 1.e4*eVtoErg erg/cm
MeV_mToErg_cm=1.e4*eVtoErg
ffForm=-.5*omega_p**2*q_e**2/V0**2/MeV_mToErg_cm      # MeV/m
eV_mToErg_m=100.*eVtoErg
ffForm=-.5*omega_p**2*q_e**2/V0**2/eV_mToErg_m        # =-6.8226e-12 eV/m
# eV_mInErg_cm=100.*eVtoErg
ffForm=-.5*omega_p**2*q_e**2/V0**2/eVtoErg            # =-6.8226e-10 eV/cm
ffForm=100.*ffForm                                 # =-6.8226e-08 eV/m

#
# Relative velocities of electrons:
#
relVeTrnsv=rmsTrnsvVe/V0 
relVeLong=rmsLongVe/V0

#
# Scanning over relative ion velocity 
#
#=====================================================================
#
# "Low" area (rmsLongVe < Vion < rmsTrnsvVe):
#
#=====================================================================
minRelVion_l=.1*relVeTrnsv 
maxRelVion_l=relVeTrnsv 
nVion_l=50 

minRelVion_l=.4e-2*relVeTrnsv 
maxRelVion_l=relVeTrnsv 
nVion_l=100 

stepRelVion_l=(maxRelVion_l-minRelVion_l)/(nVion_l-1)

relVion_l=np.zeros(nVion_l)        # relative ion velocity (related V0) 
debyeR_l=np.zeros(nVion_l)         # Debye radius (Betacool), cm
rhoMax_l=np.zeros(nVion_l)         # maximal impact parameter, cm
rhoMin_l=np.zeros(nVion_l)         # minimal impact parameter, cm
rhoPass_l=np.zeros(nVion_l) 
rhoFast_l=np.zeros((nVion_l,nField)) 

for n in range(nVion_l):
    relVion_l[n]=minRelVion_l+stepRelVion_l*n 
    absVion=V0*relVion_l[n] 
    rhoPass_l[n]=np.sqrt(absVion**2+rmsLongVe**2)*coolPassTime 
    debyeR_l[n]=np.sqrt(absVion**2+rmsLongVe**2+rmsTrnsvVe**2)/omega_p      
    for k in range(nField):
        rhoFast_l[n,k]=np.sqrt(absVion**2+rmsLongVe**2)/cyclFreq[k] 

#=====================================================================
#
# "High" area (rmsLongVe < rmsTrnsvVe < Vion):
#
#=====================================================================
minRelVion_h=relVeTrnsv 
maxRelVion_h=10.*relVeTrnsv 
nVion_h=50 

stepRelVion_h=(maxRelVion_h-minRelVion_h)/(nVion_h-1)

relVion_h=np.zeros(nVion_h)    #    # relative ion velocity (related V0) 
debyeR_h=np.zeros(nVion_h)         # Debye radius (Betacool), cm
rhoMax_h=np.zeros(nVion_h)         # maximal impact parameter, cm
rhoMin_h=np.zeros(nVion_h)         # minimal impact parameter, cm
rhoPass_h=np.zeros(nVion_h) 
rhoFast_h=np.zeros((nVion_h,nField)) 

for n in range(nVion_h):
    relVion_h[n]=minRelVion_h+stepRelVion_h*n 
    absVion=V0*relVion_h[n] 
    rhoPass_h[n]=np.sqrt(absVion**2+rmsLongVe**2)*coolPassTime 
    debyeR_h[n]=np.sqrt(absVion**2+rmsLongVe**2+rmsTrnsvVe**2)/omega_p      
    for k in range(nField):
        rhoFast_h[n,k]=np.sqrt(absVion**2+rmsLongVe**2)/cyclFreq[k] 

rhoMax_l=np.zeros(nVion_l)   
for n in range(nVion_l):
    help=max(debyeR_l[n],likeDebyeR)
    rhoMax_l[n]=min(rhoPass_l[n],help)

   
rhoMax_h=np.zeros(nVion_h)   
for n in range(nVion_h):
    help=max(debyeR_h[n],likeDebyeR)
    rhoMax_h[n]=min(rhoPass_h[n],help)
   
rhoCrit_l=np.zeros((nVion_l,nField))
for n in range(nVion_l):
    for k in range(nField):
        rhoCrit_l[n,k]=(q_e**2/m_e/cyclFreq[k]**2)**(1./3)

rhoCrit_h=np.zeros((nVion_h,nField))
for n in range(nVion_h):
    for k in range(nField):
        rhoCrit_h[n,k]=(q_e**2/m_e/cyclFreq[k]**2)**(1./3)

rhoMin_l=np.zeros(nVion_l)   
for n in range(nVion_l):
    rhoMin_l[n]=q_e**2/m_e/((relVion_l[n]*V0)**2+rmsTrnsvVe**2+rmsLongVe**2)    
   
rhoMin_h=np.zeros(nVion_h)   
for n in range(nVion_h):
    rhoMin_h[n]=q_e**2/m_e/((relVion_h[n]*V0)**2+rmsTrnsvVe**2+rmsLongVe**2)    

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Main calculations
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#=====================================================================
#
# "Low" area (rmsLongVe < Vion < rmsTrnsvVe):
#
#=====================================================================

rhoLarm_l=np.zeros((nVion_l,nField))

CL_l=np.zeros((nVion_l,nField))        # Coulomb logarithm

trnsvFF_l=np.zeros((nVion_l,nField))   # transverse friction force  
longFF_l=np.zeros((nVion_l,nField))    # longitufinal friction force          

for n in range(nVion_l):
    for k in range(nField):
        rhoLarm_l[n,k]=rmsTrnsvVe/cyclFreq[k]
        fctr=.5*rhoMax_l[n]/rhoLarm_l[n,k]
        if (fctr > 1.):
            CL_l[n,k]=math.log(fctr)

    relVion_l[n]=minRelVion_l+stepRelVion_l*n
    absVion=V0*relVion_l[n]
    trnsvVion=absVion*math.sin(alpha)  # transverse ion velocity
    longVion=absVion*math.cos(alpha)   # longitudinal ion velocity
    
    for k in range(nField):
# Transverse force:        
        trnsvFF_l[n,k]=ffForm*CL_l[n,k]/relVion_l[n]**2;
# For future: 
#         trnsvFF_l[n,k]=ffForm*CL_l[n,k]*trnsvVion/absVion**3* ...
#                        (relVeTrnsv**2-2.*relVeLong**2)
# Longitudinal force:        
        longFF_l[n,k]=ffForm*2./relVion_l[n]**2;
# For future: 
#         longFF_l[n,k]=ffForm*longVion/absVion**3* ...
#                      (3*CL_l[n,k]*relVeTrnsv**2+2.);

#=====================================================================
#
# "High" area (rmsLongVe < rmsTrnsvVe < Vion):
#
#=====================================================================

rhoLarm_h=np.zeros((nVion_h,nField))

CL_h=np.zeros((nVion_h,nField))        # Coulomb logarithm

trnsvFF_h=np.zeros((nVion_h,nField))   # transverse friction force  
longFF_h=np.zeros((nVion_h,nField))    # longitufinal friction force          

for n in range(nVion_h):
    for k in range(nField):
        rhoLarm_h[n,k]=rmsTrnsvVe/cyclFreq[k]
        fctr=.5*rhoMax_h[n]/rhoLarm_h[n,k]
        if (fctr > 1.):
            CL_h[n,k]=math.log(fctr)

    relVion_h[n]=minRelVion_h+stepRelVion_h*n
    absVion=V0*relVion_h[n]
    trnsvVion=absVion*math.sin(alpha)  # transverse ion velocity
    longVion=absVion*math.cos(alpha)   # longitudinal ion velocity
    
    for k in range(nField):
# Transverse force:        
        trnsvFF_h[n,k]=ffForm*CL_h[n,k]/relVion_h[n]**2
# For future: 
#         trnsvFF_h[n,k]=ffForm*CL_h[n,k]*trnsvVion/absVion**3* \
#                      (relVeTrnsv**2-2.*relVeLong**2)
# Longitudinal force:        
        longFF_h[n,k]=ffForm*2./relVion_h[n]**2
# For future: 
#         longFF_h[n,k]=ffForm*longVion/absVion**3* \
#                      (3*CL_h[n,k]*relVeTrnsv^2**2.);

# sys.exit()

xLimit=[.9*minRelVion_l,1.1*maxRelVion_h]

plt.figure(205)
plt.loglog(relVion_l,debyeR_l,'-b',linewidth=2)
plt.grid(True)
hold=True
plt.loglog(relVion_h,debyeR_h,'-m',linewidth=2)
plt.plot([relVion_l[0],relVion_h[nVion_h-1]],[likeDebyeR,likeDebyeR],'k',linewidth=2)
plt.xlabel('Relative Ion Velocity, $V_i/V_{e0}$',color='m',fontsize=16)
plt.ylabel('$R_{Debye}$ and $R_e$, cm',color='m',fontsize=16)
titleHeader='$R_{Debye}$ and $R_e$: $V_{e0}=%5.3f\cdot10^{%2d}$cm/s'
plt.title(titleHeader % (mantV0,powV0),color='m',fontsize=16)
plt.xlim(xLimit)
yLimit=[1.e-3,.5]
plt.ylim(yLimit)
plt.plot([relVeTrnsv,relVeTrnsv],yLimit,'--m',linewidth=1)
plt.text(2.e-3,7.5e-4,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
plt.text(2.e-4,2.4e-3,'$R_e=(3/n_e)^{1/3}$',color='k',fontsize=16)
plt.text(2.e-4,2.5e-2,'$R_{Debye}=$',color='k',fontsize=16)
plt.plot([5.8e-4,2.95e-3],[2.7e-2,2.7e-2],color='k')
plt.text(5.8e-4,3.2e-2,'$<|V_i-\Delta_{e||}|>$',color='k',fontsize=16)
plt.text(1.e-3,2.2e-2,'$\omega_{ep}$',color='k',fontsize=16)

plt.figure(207)
plt.loglog(relVion_l,rhoPass_l,'-b',linewidth=2)
plt.grid(True)
hold=True
plt.loglog(relVion_h,rhoPass_h,'-m',linewidth=2)
plt.xlabel('Relative Ion Velocity, $V_i/V_{e0}$',color='m',fontsize=16)
plt.ylabel('$R_{Pass}$, cm',color='m',fontsize=16)
titleHeader='$V_{e0}=%5.3f\cdot10^{%2d}$cm/s'
plt.title(titleHeader % (mantV0,powV0),color='m',fontsize=16)
plt.xlim(xLimit)
yLimit=[2.e-3,7.0]
plt.ylim(yLimit)
plt.plot([relVeTrnsv,relVeTrnsv],yLimit,'--m',linewidth=1)
plt.text(2.e-3,1.35e-3,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
plt.text(4.e-5,1.e-0,'$R_{Pass}=<|V_i-\Delta_{e||}|>\cdot$',color='k',fontsize=16)
plt.plot([6.1e-4,1.15e-3],[1.09,1.09],color='k',linewidth=1)
plt.text(6.1e-4,1.25,'$L_{Cool}$',color='k',fontsize=16)
plt.text(7.e-4,.8,'$V_{e0}$',color='k',fontsize=16)

plt.figure(209)
plt.loglog(relVion_l,debyeR_l,'-b',relVion_h,debyeR_h,'-m', \
           relVion_l,rhoPass_l,'-b',relVion_h,rhoPass_h,'-m',linewidth=2)
plt.grid(True)
hold=True
plt.plot ([relVion_l[0],relVion_h[nVion_h-1]],[likeDebyeR,likeDebyeR],color='k',linewidth=2)
plt.xlabel('Relative Ion Velocity, $V_i/V_{e0}$',color='m',fontsize=16)
plt.ylabel('$R_{Debye}$, $R_{Pass}$, $R_e$, cm',color='m',fontsize=16)
titleHeader='$R_{Debye}$, $R_{Pass}$, $R_e$: $V_{e0}=%5.3f\cdot10^{%2d}$cm/s'
plt.title(titleHeader % (mantV0,powV0),color='m',fontsize=16)
plt.xlim(xLimit)
yLimit=[1.e-3,10.]
plt.ylim(yLimit)
plt.plot([relVeTrnsv,relVeTrnsv],yLimit,'--m',linewidth=1)
plt.text(2.e-3,6.5e-4,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
plt.text(3.e-5,2.45e-3,'$R_e$',color='k',fontsize=16)
plt.text(3.e-5,5.e-2,'$R_{Debye}$',color='k',fontsize=16)
plt.text(3.e-5,1.e-2,'$R_{Pass}$',color='k',fontsize=16)

plt.figure(215)
plt.loglog(relVion_l,rhoMax_l,'-b',relVion_h,rhoMax_h,'-m',linewidth=2)
plt.grid(True)
hold=True
plt.xlabel('Relative Ion Velocity, $V_i/V_{e0}$',color='m',fontsize=16)
plt.ylabel('$R_{max}$, cm',color='m',fontsize=16)
titleHeader='$V_{e0}=%5.3f\cdot10^{%2d}$cm/s'
plt.title(titleHeader % (mantV0,powV0),color='m',fontsize=16)
plt.xlim(xLimit)
yLimit=[2.e-3,.5]
plt.ylim(yLimit)
plt.plot([relVeTrnsv,relVeTrnsv],yLimit,'--m',linewidth=1)
plt.text(2.e-3,1.5e-3,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
plt.text(6.e-5,0.2,'$R_{max}=min\{max\{R_{Debye},R_e\},R_{Pass}\}$',color='k',fontsize=16)

plt.figure(305)
plt.loglog(relVion_l,rhoFast_l[:,0],'-b',relVion_h,rhoFast_h[:,0],'-m', \
       relVion_l,rhoCrit_l[:,0],'-b',relVion_h,rhoCrit_h[:,0],'-m',linewidth=2)
plt.grid(True)
hold=True
plt.xlabel('Relative Ion Velocity, $V_i/V_{e0}$',color='m',fontsize=16)
plt.ylabel('$R_{Fast}$ and $R_{Crit}$, cm',color='m',fontsize=16)
titleHeader='$V_{e0}=%5.3f\cdot10^{%2d}$cm/s, $B=%6.1f$ Gs'
plt.title(titleHeader % (mantV0,powV0,fieldB[0]),color='m',fontsize=16)
plt.xlim(xLimit)
yLimit=[2.e-6,.01]
plt.ylim(yLimit)
plt.plot([relVeTrnsv,relVeTrnsv],yLimit,'--m',linewidth=1)
plt.text(2.e-3,1.325e-6,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
plt.text(1.e-3,1.5e-4,'$R_{Fast}=<|V_i-\Delta V_{e||}|>/\omega_{Larm}$',color='k',fontsize=16)
plt.text(1.e-3,2.8e-5,'$R_{Crit}=(q_e^2/m_e/\omega_{Larm}^2)^{1/3}$',color='k',fontsize=16)

plt.figure(307)
plt.loglog(relVion_l,rhoMin_l,'-b',relVion_h,rhoMin_h,'-m',linewidth=2)
plt.grid(True)
hold=True
plt.xlabel('Relative Ion Velocity, $V_i/V_{e0}$',color='m',fontsize=16)
plt.ylabel('$R_{min}$, cm',color='m',fontsize=16)
titleHeader='$V_{e0}=%5.3f\cdot10^{%2d}$cm/s'
plt.title(titleHeader % (mantV0,powV0),color='m',fontsize=16)
plt.xlim(xLimit)
yLimit=[1.e-9,2.e-7]
plt.ylim(yLimit)
plt.plot([relVeTrnsv,relVeTrnsv],yLimit,'--m',linewidth=1)
plt.text(2.e-3,7.e-10,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
plt.text(1.e-4,3.e-8,'$R_{min}=$',color='k',fontsize=16)
plt.plot([2.4e-4,1.75e-3],[3.2e-8,3.2e-8],color='k')
plt.text(4.e-4,3.5e-8,'$q_e^2/m_e$',color='k',fontsize=16)
plt.text(2.4e-4,2.3e-8,'$(<|V_i-V_{e\perp}|>)^2$',color='k',fontsize=16)

plt.figure(315)
plt.loglog(relVion_l,rhoMax_l,'-b',relVion_h,rhoMax_h,'-m', \
       relVion_l,2.*rhoLarm_l[:,0],'-b',relVion_h,2.*rhoLarm_h[:,0],'-m',linewidth=2)
plt.grid(True)
hold=True
plt.xlabel('Relative Ion Velocity, $V_i/V_{e0}$',color='m',fontsize=16)
# plt.ylabel('$R_{max}$ and $2\cdot<rho_\perp>$, cm',color='m',fontsize=16)
plt.ylabel('Impact Parameter, cm',color='m',fontsize=16)
titleHeader='$V_{e0}=%5.3f\cdot10^{%2d}$cm/s, $B=%6.1f$ Gs'
plt.title(titleHeader % (mantV0,powV0,fieldB[0]),color='m',fontsize=16)
titleHeader='$V_{e0}=%5.3f\cdot10^{%2d}$cm/s, $B=%6.1f$ Gs'
plt.title(titleHeader % (mantV0,powV0,fieldB[0]),color='m',fontsize=16)
plt.xlim(xLimit)
yLimit=[1.e-3,.6]
plt.ylim(yLimit)
plt.plot([relVeTrnsv,relVeTrnsv],yLimit,'--m',linewidth=1)
plt.text(2.e-3,7.25e-4,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
plt.text(7.e-4,1.75e-3,'$2\cdot<rho_\perp>$',color='k',fontsize=16)
plt.text(7.e-4,5.e-2,'$R_{max}$',color='k',fontsize=16)
plt.text(4.5e-3,.3*1.e-2,'$<rho_\perp>=$',color='k',fontsize=16)
plt.plot([2.e-2,4.1e-2],[.3*1.1e-2,.3*1.1e-2],color='k')   
plt.text(2.e-2,.3*1.25e-2,'$ \Delta V_{e\perp}$',color='k',fontsize=16)
plt.text(2.e-2,.3*9.e-3,'$\omega_{Larm}$',color='k',fontsize=16)
plt.text(1.5e-4,7.e-3,'Magnetized Collisions',color='r',fontsize=25)
plt.text(1.5e-4,1.05e-3,'Adiabatic Collisions',color='r',fontsize=25)
plt.text(2.8e-5,.2,'Collisions are Sreened',color='r',fontsize=25)

#----------------------------------------------
#
# Figures of Coulomb logarithms:
#
#----------------------------------------------
plt.figure(320)
plt.loglog(relVion_l,CL_l[:,0],'-b',relVion_h,CL_h[:,0],'-m',linewidth=2)
plt.grid(True)
hold=True
plt.xlabel('Relative Ion Velocity, $V_i/V_{e0}$',color='m',fontsize=16)
plt.ylabel('$CL_{magnetized}$',color='m',fontsize=16)
titleHeader='Coulomb Logarithm: $V_{e0}=%5.3f\cdot10^{%2d}$cm/s, $B=%6.1f$ Gs'
plt.title(titleHeader % (mantV0,powV0,fieldB[0]),color='m',fontsize=16)
plt.xlim(xLimit)
yLimit=[.4,6.]
plt.ylim(yLimit)
plt.plot([relVeTrnsv,relVeTrnsv],yLimit,'--m',linewidth=1)
plt.text(2.e-3,.35,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
plt.text(1.e-4,1.25,'$CL_{magnetized}=$',color='k',fontsize=16)
plt.plot([5.7e-4,2.4e-3],[1.29,1.29],color='k')
plt.text(8.e-4,1.35,'$R_{max}$',color='k',fontsize=16)
plt.text(5.7e-4,1.15,'$2\cdot<rho_\perp>$',color='k',fontsize=16)

plt.figure(3201)
plt.semilogx(relVion_l,CL_l[:,0],'-b',relVion_h,CL_h[:,0],'-m',linewidth=2)
plt.grid(True)
hold=True
plt.xlabel('Relative Ion Velocity, $V_i/V_{e0}$',color='m',fontsize=16)
plt.ylabel('$CL_{magnetized}$',color='m',fontsize=16)
titleHeader='Coulomb Logarithm: $V_{e0}=%5.3f\cdot10^{%2d}$cm/s, $B=%6.1f$ Gs'
plt.title(titleHeader % (mantV0,powV0,fieldB[0]),color='m',fontsize=16)
plt.xlim(xLimit)
yLimit=[0,6.]
plt.ylim(yLimit)
plt.plot([relVeTrnsv,relVeTrnsv],yLimit,'--m',linewidth=1)
plt.text(2.e-3,-.35,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
plt.text(1.e-4,4.25,'$CL_{magnetized}=$',color='k',fontsize=16)
plt.plot([5.7e-4,2.4e-3],[4.31,4.31],color='k')
plt.text(8.e-4,4.4,'$R_{max}$',color='k',fontsize=16)
plt.text(5.7e-4,4.05,'$2\cdot<rho_\perp>$',color='k',fontsize=16)

#
# Friction forces:
#       
plt.figure(395)
plt.loglog(relVion_l,abs(trnsvFF_l[:,0]),'-b',relVion_h,abs(trnsvFF_h[:,0]),'-m',linewidth=2)
plt.grid(True)
hold=True
plt.xlabel('Relative Ion Velocity, $V_i/V_{e0}$',color='m',fontsize=16)
plt.ylabel('$F_\perp$, eV/m',color='m',fontsize=16)
titleHeader='Transverse Friction Force: $V_{e0}=%5.3f\cdot10^{%2d}$cm/s, $B=%6.1f$ Gs'
plt.title(titleHeader % (mantV0,powV0,fieldB[0]),color='m',fontsize=16)
plt.xlim(xLimit)
yLimit=[2.e-4,2.e+2]
plt.ylim(yLimit)
plt.plot([relVeTrnsv,relVeTrnsv],yLimit,'--m',linewidth=1)
plt.text(2.e-3,1.e-4,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)

plt.figure(400)
plt.loglog(relVion_l,abs(longFF_l[:,0]),'-b',relVion_h,abs(longFF_h[:,0]),'-m',linewidth=2)
plt.grid(True)
hold=True
plt.xlabel('Relative Ion Velocity, $V_i/V_{e0}$',color='m',fontsize=16)
plt.ylabel('$F_\perp$, eV/m',color='m',fontsize=16)
titleHeader='Transverse Friction Force: $V_{e0}=%5.3f\cdot10^{%2d}$cm/s, $B=%6.1f$ Gs'
plt.title(titleHeader % (mantV0,powV0,fieldB[0]),color='m',fontsize=16)
# plt.xlim(xLimit)
# yLimit=[2.e-4,2.e+2]
# plt.ylim(yLimit)
# plt.plot([relVeTrnsv,relVeTrnsv],yLimit,'--m',linewidth=1)
# plt.text(2.e-3,1.e-4,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)

plt.show()

sys.exit()
