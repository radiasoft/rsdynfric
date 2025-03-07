# from __future__ import division

#-------------------------------------
#
#        Started at 06/08/2018 (YuE)
#
# This script based on the previous script
# threeApproachesComparison_v6.py
#
## Upgraded version of python (python3.4): script was rewritten to take into
# account some differences in the descriptions and using of some functions
# (version cma_v3 and more earlier scripts are written under python2).    
# 
# 07/24/2018: IT IS NOT FINISHED:
# 
# Which are still unsatisfactory: 
# 1) the absolute values of frictional forces for all methods of calculation,
# 2) their dependence on the ion velocity.
#
# But nevertheless, the dependences of the transmitted energy on the impact
# parameter are close to the inverse quadratic (as it should be!) at all velocities.
# 
# 07/27/2018: IT IS NOT FINISHED:
# 
# Which are still unsatisfactory: 
# 1) the absolute values of frictional forces for all methods of calculation,
# 2) their dependence on the ion velocity.
# The investigation of that is in progress.
#
# Some features were improved, some figures were corrected.
# 
#-------------------------------------

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
# [2] David L. Bruhwiler, Stephen D Webb. "New Algorithm for Dynamical Friction
#    of Ions in a Magnetized Electron Beam". AIP Conf. Proc. 1812, 05006 (2017). 
#
#========================================================

#########################################################
#
# Main issues of the calculations:
#
# 1) Friction force (FF) is calculated in the (P)article (R)est (F)rame,
#    i.e. in the frame moving together with both (cooled and cooling)
#    beams at a velocity V0;
# 2) Friction force is calculated  for each value of ion velocity
#    in the interval from .1*rmsTrnsvVe till 10*rmsTrnsvVe;
# 3) Initially assumped that all electrons have a logitudinal 
#    velocity rmsLongVe and transversal velocity rmsTrnsvVe;
# 4) For each ion velocity the minimal and maximal values of the 
#    impact parameter are defined. Radius of the shielding of the 
#    electric field of the ion equals to the value of the maximal
#    impact parameter;  
# 5) For each impact parameter in the interval from minimal till 
#    maximal values the transfered momenta deltap_x,y,z are
#    calculated;
# 6) Founded transfered momenta allow to calculate the transfered
#    energy delta_E =deltap^2/(2*m_e) and to integrate it over
#    impact parameter; then (expressions (3.4), (3.5) from [1]): 
#        FF =-2*pi*n_e*integral_rhoMin^rhoMax delta_E*rho*drho;
# 7) For taking into account the velocity distribution of the
#    electrons it is necessary to repeat these calculations for
#    each value of the electron's velocity and then integrate result
#    over distribution of the velocities.
#
# 10/26/2018:  
#
# 8) Item 6 is wrong and correct expression for transfered
#    energy delta_E will be used;
# 9) Method (my own) Least Squares Method - LSM is used to fit the
#    dependence of transferred momenta on impact parameter;
#
#
# 11/08/2018:  
#
# 10) Two functions ('fitting' and 'errFitAB' are defined to realize 
#     my LSM to find the parameters of the fitting end error of this
#     fitting;
#
# 11) Analys of different dependeces between values; graphical
#     presentation of these dependences;
#
#########################################################
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
from scipy import optimize
from statistics import mean
from array import array

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

pi=3.14159265358 

#
# Physical constants:
#
m_e=9.10938356e-28          # electron mass, g
m_elec=m_e                  # to keep variable from previous script
m_p=1.672621898e-24         # electron mass, g
M_ion = m_p                 # to keep variable from previous script  
q_e=4.803204673e-10         # electron charge, CGSE unit: sqrt(g*cm^3/sec^2)
q_elec=q_e                  # to keep variable from previous script
Z_ion = q_e                 # to keep variable from previous script
cLight=2.99792458e10        # speed of light, cm/sec
eVtoErg=1.6021766208e-12    # 1 eV = 1.6...e-12 erg
CtoPart=2.99792458e9        # 1 C = 1 A*sec = 2.9...e9 particles

# 
# Electron beam parameters:
#
Ekin=3.0e4                              # kinetic energy, eV
curBeam=0.5                             # current density, A/cm^2
dBeam=3.0                               # beam diameter, cm
angSpread=3.0                           # angular spread, mrad
trnsvT=0.5                              # transversal temperature, eV
longT=2.0e-4                            # longitudinal temperature, eV (was 2.0e-4)
eTempTran=trnsvT                        # to keep variable from previous script
eTempLong=longT                         # to keep variable from previous script
nField=1                                # number ov values  of the magnetic field
fieldB=np.zeros(nField)                 # magnetic field
fieldB[0]=3.e3                          # Gs
omega_p=1.0e9                           # plasma frequency, 1/sec
n_e=omega_p**2*m_e/(4.*pi*q_e**2)       # plasma density, 3.1421e+08 cm-3

n_e1=8.e7                               # plasma density, cm-3
omega_p1=np.sqrt(4.*pi*n_e1*q_e**2/m_e) # plasma frequency, 5.0459e+08 1/s  

#
# Cooling system parameter:
#
coolLength=150.0                        # typical length of the coolong section, cm

# 
# Calculated parameters of the electron beam:
#
V0=np.sqrt(2.*Ekin*eVtoErg/m_e)           # longitudinal velocity, cm/s
tetaV0=0.                                 # angle between V0 and magnetic field, rad
B_mag=fieldB[0]*np.cos(tetaV0)            # magnetic field acting on an electron, Gs
rmsTrnsvVe=np.sqrt(2.*trnsvT*eVtoErg/m_e) # RMS transversal velocity, cm/s
rmsLongVe=np.sqrt(2.*longT*eVtoErg/m_e)   # RMS longitudinal velocity, cm/s
# dens=curBeam*CtoPart/V0                 # density, 1/cm^3
# omega=np.sqrt(4.*pi*dens*q_e**2/m_e)    # plasma frequency, 1/s
cyclFreq=q_e*B_mag/(m_e*cLight)           # cyclotron frequency, 1/s
rmsRoLarm=rmsTrnsvVe*cyclFreq**(-1)       # RMS Larmor radius, cm
dens=omega_p**2*m_e/(4.*pi*q_e**2)        # density, 1/cm^3
likeDebyeR=(3./dens)**(1./3.)             # "Debye" sphere with 3 electrons, cm

coolPassTime=coolLength/V0                # time pass through cooling section, cm

thetaVi=0.                                # polar angle ion and cooled electron beams, rad
phiVi=0.                                  # azimuth angle ion and cooled electron beams, rad

powV0=round(np.log10(V0)) 
mantV0=V0/(10**powV0) 

pow_n_e=round(np.log10(n_e)) 
mant_n_e=n_e/(10**pow_n_e) 

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
# ffForm=-.5*omega_p**2*q_e**2/V0**2/MeV_mToErg_cm     # MeV/m
eV_mToErg_m=100.*eVtoErg
# ffForm=-.5*omega_p**2*q_e**2/V0**2/eV_mToErg_m       # =-6.8226e-12 eV/m
eV_mInErg_cm=100.*eVtoErg
ffForm=-.5*omega_p**2*q_e**2/V0**2/eVtoErg            # =-6.8226e-10 eV/cm
ffForm=100.*ffForm                                    # =-6.8226e-08 eV/m

ergToEV = 1./1.60218e-12

#
# Relative velocities of electrons:
#
relVeTrnsv=rmsTrnsvVe/V0 
relVeLong=rmsLongVe/V0

print ('V0=%e cm/s, rmsTrnsvVe=%e cm/s (rel = %e), rmsLongVe=%e cm/s (rel = %e)' % \
       (V0,rmsTrnsvVe,relVeTrnsv,rmsLongVe,relVeLong))

# Indices:
(Ix, Ipx, Iy, Ipy, Iz, Ipz) = range(6)

stepsNumberOnGyro = 25                                             # number of the steps on each Larmour period

'''
#
# Opening the input file: 
#
inputFile='areaOfImpactParameter_tAC-v6_fig110.data'
print ('Open input file "%s"...' % inputFile)
inpfileFlag=0
try:
   inpfile = open(inputFile,'r')
   inpfileFlag=1
except:
   print ('Problem to open input file "%s"' % inputFile)
if inpfileFlag == 1:
   print ('No problem to open input file "%s"' % inputFile)

lines=0                                                            # Number of current line from input file   
dataNumber=0                                                       # Number of current value of any types of Data
xAboundary=np.zeros(100)
xBboundary=np.zeros(100)
while True:
   lineData=inpfile.readline()
#    print ('line=%d: %s' % (lines,lineData))
   if not lineData:
      break
   lines += 1
   if lines > 4:
      words=lineData.split()
      nWords=len(words)
#      print ('Data from %d: words=%s, number of entries = %d' % (lines,words,nWords))
      xAboundary[dataNumber]=float(words[0])
      xBboundary[dataNumber]=float(words[1])
      dataNumber += 1

inpfile.close()
print ('Close input file "%s"' % inputFile)
'''

#====================================================================
#
#------------------ Begin of defined functions -----------------------

#
# Larmor frequency electron:
#
def omega_Larmor(mass,B_mag):
    return (q_elec)*B_mag/(mass*clight*1.e+2)                      # rad/sec

#
# Convertion from 6-vector of relectron's "coordinates" to 6-vector
# of guiding-center coordinates:
# z_e=(x_e,px_e,y_e,py_e,z_e,pz_e) --> zgc_e=(phi,p_phi,y_gc,p_gc,z_e,pz_e);
#
def toGuidingCenter(z_e):
    mOmega=m_elec*omega_L                                          # g/sec
    zgc_e=z_e.copy()                                               # 6-vector
    zgc_e[Ix] = np.arctan2(z_e[Ipx]+mOmega*z_e[Iy],z_e[Ipy])       # radians
    zgc_e[Ipx]= (((z_e[Ipx]+mOmega*z_e[Iy])**2+z_e[Ipy]**2)/(2.*mOmega)) # g*cm**2/sec
    zgc_e[Iy] =-z_e[Ipx]/mOmega                                    # cm
    zgc_e[Ipy]= z_e[Ipy]+mOmega*z_e[Ix]                            # g/sec
    return zgc_e

#
# Convertion from 6-vector of guiding-center coordinates to 6-vector
# of electron's "coordinates":
# zgc_e=(phi,p_phi,y_gc,p_gc,z_e,pz_e) --> z_e=(x_e,px_e,y_e,py_e,z_e,pz_e);
#
def fromGuidingCenter(zgc_e):
    mOmega=m_elec*omega_L                                          # g/sec
    rho_larm=np.sqrt(2.*zgc_e[Ipx]/mOmega)                         # cm
    z_e = zgc_e.copy()                                             # 6-vector
    z_e[Ix] = zgc_e[Ipy]/mOmega-rho_larm*np.cos(zgc_e[Ix])         # cm
    z_e[Ipx]=-mOmega*zgc_e[Iy]                                     # g*cm/sec
    z_e[Iy] = zgc_e[Iy]+rho_larm*np.sin(zgc_e[Ix])                 # cm
    z_e[Ipy]= mOmega*rho_larm*np.cos(zgc_e[Ix])                    # g*cm/sec
    return z_e

#
# Matrix to dragg electron through the solenoid with field 'B_mag' 
# during time interval 'deltaT':
#
def solenoid_eMatrix(B_mag,deltaT):
    slndMtrx=np.identity(6)
    omega_L=omega_Larmor(m_elec,B_mag)                             # rad/sec 
    mOmega= m_elec*omega_L                                         # g/sec
    phi=omega_L*deltaT                                             # phase, rad
    cosPhi=math.cos(phi)                                           # dimensionless
    sinPhi=math.sin(phi)                                           # dimensionless
    cosPhi_1=2.*math.sin(phi/2.)**2                                # dimensionless
    slndMtrx[Iy, Iy ]= cosPhi                                      # dimensionless
    slndMtrx[Ipy,Ipy]= cosPhi                                      # dimensionless
    slndMtrx[Iy, Ipy]= sinPhi/mOmega                               # sec/g
    slndMtrx[Ipy,Iy ]=-mOmega*sinPhi                               # g/sec
    slndMtrx[Iz, Ipz]= deltaT/m_elec                               # sec/g
    slndMtrx[Ix, Ipx]= sinPhi/mOmega                               # sec/g
    slndMtrx[Ix, Iy ]= sinPhi                                      # dimensionless
    slndMtrx[Ix, Ipy]= cosPhi_1/mOmega                             # sec/g 
    slndMtrx[Iy, Ipx]=-cosPhi_1/mOmega                             # sec/g
    slndMtrx[Ipy,Ipx]=-sinPhi                                      # dimensionless
    return slndMtrx

#
# Matrix to dragg particle through the drift during time interval 'deltaT':
#
def drift_Matrix(M_prtcl,deltaT):
    driftMtrx = np.identity(6)
    for i in (Ix,Iy,Iz):
        driftMtrx[i,i+1]=deltaT/M_prtcl                            # sec/g
    return driftMtrx

#
# Matrix to dragg electron in the "guiding center" system during time interval 'deltaT':
#
def guidingCenter_Matrix(deltaT):
    gcMtrx = np.identity(6)
    gcMtrx[Iz,Ipz]=deltaT/m_elec                                   # sec/g
    return gcMtrx

#
# Description of the collision during time interval 'deltaT'
# in the system coordinates of "guiding center" of electron
# input - 6-vectors for electron and ion before collision and time step deltaT; 
# output - transfered momenta to ion and electron: 
#
def guidingCenterCollision(vectrElec_gc,vectrIon,deltaT):

   dpIon=np.zeros(3)
   dpElec=np.zeros(3)
   mOmegaLarm=m_elec*omega_L                                 # g/sec
   dpFactor_gc=q_elec**2                                     # g*cm^3/sec^2
   rhoLarm_gc=np.sqrt(2.*vectrElec_gc[1]/mOmegaLarm)         # cm
   sinOmega_gc=math.sin(vectrElec_gc[0])
   cosOmega_gc=math.cos(vectrElec_gc[0])
   x_gc=vectrElec_gc[3]/mOmegaLarm                           # cm
   numer=(vectrIon[0]-x_gc)*cosOmega_gc- \
         (vectrIon[2]-vectrElec_gc[2])*sinOmega_gc           # cm
   denom=((vectrIon[0]-x_gc)**2+(vectrIon[2]-vectrElec_gc[2])**2+ \
          (vectrIon[4]-vectrElec_gc[4])**2+rhoLarm_gc**2)**(3/2)         # cm^3
   action=vectrElec_gc[1]+dpFactor_gc*numer*rhoLarm_gc/(omega_L*denom)   # g*cm^2/sec
   b_gc=np.sqrt((vectrIon[0]-x_gc)**2+ \
                (vectrIon[2]-vectrElec_gc[2])**2+ \
                (vectrIon[4]-vectrElec_gc[4])**2+2.*action/mOmegaLarm)   # cm
# Dimensions of dpIon, dpElec are g*cm/sec:   
   dpIon[0]=-dpFactor_gc*deltaT*(vectrIon[0]-x_gc)/b_gc**3               
   dpIon[1]=-dpFactor_gc*deltaT*(vectrIon[2]-vectrElec_gc[2])/b_gc**3  
   dpIon[2]=-dpFactor_gc*deltaT*(vectrIon[4]-vectrElec_gc[4])/b_gc**3
   dpElec[0]=-dpIon[0] 
   dpElec[1]=-dpIon[1] 
   dpElec[2]=-dpIon[2]
#    print ('dpIon[0]=%e, dpIon[1]=%e, dpIon[2]=%e' % \
#           (dpIon[0],dpIon[1],dpIon[2]))
   return dpIon,dpElec,action,b_gc                                      

#
# "Magnus expansion" description of the collision during time interval 'deltaT'
# in the system coordinates of "guiding center" of electron
# input - 6-vectors for electron and ion before collision and time step deltaT; 
# output - transfered momenta to ion and electron and electron y_gc coordinate
# as well calculated parameters C1,C2,C3,b,D1,D2,q for testing: 
#
def MagnusExpansionCollision(vectrElec_gc,vectrIon,deltaT):

#    print ('Ion: x=%e, y=%e, z=%e' % (vectrIon[0],vectrIon[2],vectrIon[4]))
#    print ('Electron: x=%e, y=%e, z=%e' % 
#          (vectrElec_gc[0],vectrElec_gc[4],vectrElec_gc[4]))
   dpIon=np.zeros(3)
   dpElec=np.zeros(3)
   mOmegaLarm=m_elec*omega_L                                 # g/sec
   dpFactor_gc=q_elec**2                                     # g*cm^3/sec^2
   rhoLarm_gc=np.sqrt(2.*vectrElec_gc[1]/mOmegaLarm)         # cm
   sinOmega_gc=math.sin(vectrElec_gc[0])
   cosOmega_gc=math.cos(vectrElec_gc[0])
   x_gc=vectrElec_gc[3]/mOmegaLarm                           # cm
   numer=(vectrIon[0]-x_gc)*cosOmega_gc- \
         (vectrIon[2]-vectrElec_gc[2])*sinOmega_gc           # cm
   denom=((vectrIon[0]-x_gc)**2+(vectrIon[2]-vectrElec_gc[2])**2+ \
          (vectrIon[4]-vectrElec_gc[4])**2+rhoLarm_gc**2)**(3./2.)     # cm^3
   action=vectrElec_gc[1]+dpFactor_gc*numer*rhoLarm_gc/(omega_L*denom) # g*cm^2/sec
#    C1=np.sqrt((vectrIon[0]-x_gc)**2+ \
#               (vectrIon[2]-vectrElec_gc[2])**2+ \
#               (vectrIon[4]-vectrElec_gc[4])**2+2.*action/mOmegaLarm)   # cm^2
   C1=(vectrIon[0]-x_gc)**2+(vectrIon[2]-vectrElec_gc[2])**2+ \
      (vectrIon[4]-vectrElec_gc[4])**2+2.*action/mOmegaLarm        # cm^2
   C2=2.*((vectrIon[0]-x_gc)*vectrIon[1]/M_ion+ \
          (vectrIon[2]-vectrElec_gc[2])*vectrIon[3]/M_ion+ \
          (vectrIon[4]-vectrElec_gc[4])* \
	  (vectrIon[5]/M_ion-vectrElec_gc[5]/m_elec))              # cm^2/sec
   C3=(vectrIon[1]/M_ion)**2+(vectrIon[3]/M_ion)**2+ \
      (vectrIon[5]/M_ion-vectrElec_gc[5]/m_elec)**2                # cm^2/sec^2
   b=np.sqrt(C1+C2*deltaT+C3*deltaT**2)                            # cm
   D1=(2.*C3*deltaT+C2)/b-C2/np.sqrt(C1)                           # cm/sec
   D2=(C2*deltaT+2.*C1)/b-2.*np.sqrt(C1)                           # cm
   q=4.*C1*C3-C2**2                                                # cm^4/sec^2
# Dimensions of dpIon, deElec are g*cm/sec:   
   dpIon[0]=-2.*dpFactor_gc/q*((vectrIon[0]-x_gc)*D1-vectrIon[1]/M_ion*D2)
   dpIon[1]=-2.*dpFactor_gc/q*((vectrIon[2]-vectrElec_gc[2])*D1- \
                               vectrIon[3]/M_ion*D2)
   dpIon[2]=-2.*dpFactor_gc/q*((vectrIon[4]-vectrElec_gc[4])*D1- \
                               (vectrIon[5]/M_ion-vectrElec_gc[5]/m_elec)*D2) 
   dpElec[0]=-dpIon[0] 
   dpElec[1]=-dpIon[1] 
   dpElec[2]=-dpIon[2] 
   dy_gc=dpIon[0]/mOmegaLarm                                        # cm
#    print ('dpIon[0]=%e, dpIon[1]=%e, dpIon[2]=%e' % \
#           (dpIon[0],dpIon[1],dpIon[2]))
   return dpIon,dpElec,action,dy_gc,C1,C2,C3,b,D1,D2,q                                      

#
# Minimized functional (my own Least Squares Method - LSM;
# Python has own routine for LSM - see site  
#     http://scipy-cookbook.readthedocs.io/items/FittingData.html):
#
# Funcional = {log10(funcY) - [fitB*log10(argX) + fitA]}^2
#

def fitting(nPar1,nPar2,argX,funcY):

   log10argX = np.zeros((nPar1,nPar2))
   log10funcY = np.zeros((nPar1,nPar2))
   
   for i in range(nVion):
      for n in range(nPar1):
         log10argX[n,i] = np.log10(argX[n,i])
         log10funcY[n,i] = np.log10(funcY[n,i])

   sumArgX = np.zeros(nPar2)
   sumArgX2 = np.zeros(nPar2)
   sumFuncY = np.zeros(nPar2) 
   sumArgXfuncY= np.zeros(nPar2)
   fitA = np.zeros(nPar2) 
   fitB = np.zeros(nPar2)

   for i in range(nPar2):
      for n in range(nPar1):
         sumArgX[i] += log10argX[n,i]
         sumArgX2[i] += log10argX[n,i]**2
         sumFuncY[i] += log10funcY[n,i]
         sumArgXfuncY[i] += log10argX[n,i]*log10funcY[n,i]

      delta = sumArgX[i]**2-nPar1*sumArgX2[i]
      fitA[i] = (sumArgX[i]*sumArgXfuncY[i]-sumArgX2[i]*sumFuncY[i])/delta
      fitB[i] = (sumArgX[i]*sumFuncY[i]-nPar1*sumArgXfuncY[i])/delta
#      print ('fitA(%d) = %e, fitB(%d) = %e' % (i,fitA[i],i,fitB[i]))

   argXfit = np.zeros((nPar1,nPar2))
   funcYfit = np.zeros((nPar1,nPar2))
   funcHi2 = np.zeros(nPar2)

   for i in range(nPar2):
      factorA = math.pow(10.,fitA[i])
      for n in range(nPar1):
         argXfit[n,i] = math.pow(10.,log10argX[n,i])
         funcYfit[n,i] = factorA*math.pow(argXfit[n,i],fitB[i])
         funcHi2[i] += (np.log10(abs(funcY[n,i])) - np.log10(abs(funcYfit[n,i])))**2  

   return fitA,fitB,funcHi2,argXfit,funcYfit

#
# +-Errors for fitied parameters fitA and fitB:
#
def errFitAB(nPar1,nPar2,argX,funcY,fitA,fitB,funcHi2,errVar,errType):
 
   log10argX = np.zeros((nPar1,nPar2))
   log10funcY = np.zeros((nPar1,nPar2))
   sumArgX = np.zeros(nPar2)
   sumArgX2 = np.zeros(nPar2)

   stepA = 5.e-4*mean(funcHi2)
   stepB = 1.e-4*mean(funcHi2)
#   print ('errFitAB: mean(funcHi2) = %e, stepA = %e, stepB = %e'  % (mean(funcHi2),stepA,stepB))
   for i in range(nPar2):
      for n in range(nPar1):
         log10argX[n,i] = np.log10(argX[n,i])
         log10funcY[n,i] = np.log10(funcY[n,i])
         sumArgX[i] += log10argX[n,i]
         sumArgX2[i] += log10argX[n,i]**2

   posErrFit = np.zeros(nPar2)
   for i in range(nPar2):
      k = 0
      deltaFuncHi2 = 0.
      while (deltaFuncHi2 < 1.):
         k += 1 
         if k > 2000:
            print ('Break in errFitAB (Fit funcY: case %d); positive error) for %d' % (errVar,i))
            break
#         print ('i=%d: fitParamtr = %e, funcHi2 = %e' % (i,fitParamtr[i], funcHi2[i]))
         curFitA = fitA[i]
         if (int(errVar) == 1):
            curFitA = fitA[i] + k*stepA
         curFuncHi2 = 0.
         factorA = math.pow(10.,curFitA)
         curFitB = fitB[i]
         if (int(errVar) == 2):
            curFitB = fitB[i] + k*stepB
         curFuncHi2 = 0.
         for n in range(nPar1):
            curArgX = math.pow(10.,log10argX[n,i])
            curFuncYfit = factorA*math.pow(curArgX,curFitB)
            curFuncHi2 += (np.log10(abs(curFuncYfit)) - log10funcY[n,i])**2 
         deltaFuncHi2 = curFuncHi2 - funcHi2[i]
      if (int(errVar) == 1):	  
         posErrFit[i] = abs(curFitA - fitA[i])
      else:
         posErrFit[i] = abs(curFitB - fitB[i])
      func1sigma2 = funcHi2[i]/(nPar2-3)
      if (int(errVar) == 1):	  
         fitSigma = np.sqrt(sumArgX2[i]/(nPar2*sumArgX2[i]-sumArgX[i]**2)*func1sigma2)
      else:
         fitSigma = np.sqrt(nPar2/(nPar2*sumArgX2[i]-sumArgX[i]**2)*func1sigma2) 
      if (int(errType) == 2):
         posErrFit[i] = fitSigma
#      if (int(errVar) == 1):	  
#         print ('i=%d: fitA  = %e + %e (%e), funcHi2 = %e (for %d steps curFuncHi2 = %e)' % \
#                (i,fitA[i],posErrFit[i],fitSigma,funcHi2[i],k,curFuncHi2))
#      else:
#         print ('i=%d: fitB  = %e + %e (%e), funcHi2 = %e (for %d steps curFuncHi2 = %e)' % \
#                (i,fitB[i],posErrFit[i],fitSigma,funcHi2[i],k,curFuncHi2))

   negErrFit = np.zeros(nPar2)
   for i in range(nPar2):
      k = 0
      deltaFuncHi2 = 0.
      while (deltaFuncHi2 < 1.):
         k += 1 
         if k > 2000:
            print ('Break in errFitAB (Fit funcY: case %d); negative error) for %d' % (errVar,i))
            break
         curFitA = fitA[i]
         if (int(errVar) == 1):
            curFitA = fitA[i] - k*stepA
         factorA = math.pow(10.,curFitA)
         curFitB = fitB[i]
         if (int(errVar) == 2):
            curFitB = fitB[i] - k*stepB
         curFuncHi2 = 0.
         for n in range(nPar1):
            curArgX = math.pow(10.,log10argX[n,i])
            curFuncYfit = factorA*math.pow(curArgX,curFitB)
            curFuncHi2 += (np.log10(abs(curFuncYfit)) - log10funcY[n,i])**2 
         deltaFuncHi2 = curFuncHi2 - funcHi2[i] 
      if (int(errVar) == 1):	  
         negErrFit[i] = abs(curFitA - fitA[i])
      else:
         negErrFit[i] = abs(curFitB - fitB[i])
      if (int(errType) == 2):
         negErrFit[i] = posErrFit[i]
#      if (errVar == 1):	  
#         print ('i=%d: fitA  = %e - %e, funcHi2 = %e (for %d steps curFuncHi2 = %e)' % \
#                (i,fitA[i],posErrFit[i],funcHi2[i],k,curFuncHi2))
#      else:
#         print ('i=%d: fitB  = %e - %e, funcHi2 = %e (for %d steps curFuncHi2 = %e)' % \
#                (i,fitB[i],negErrFit[i],funcHi2[i],k,curFuncHi2))
   return posErrFit,negErrFit


#----------------------------------------------------
#
#  "Gauss-Kronrod" method of integration (GK)
#
# Points (psi_i) and weigths (w_i) to integrate for interval from -1 to 1;
# These data are from William H. Beyer. "Handbook of Mathematical Science".
# 5th Edition, CRC Press, Inc, 1978 (pp. 661,662).
#
# To integrate for interval from 0 to 1 it is necessary to change points
# psi_i with points ksi_i=(1+psi_i)/2;
#
# For method with order N for function F(x):
#      int_(-1)^1 = sum_1^N [w_i* F(psi_i)];
#
# In case of integration over interval from a to b:  
#      int_(a)^b = (b-a)/2 * sum_1^N [w_i* F(x_i)], where
#            x_i = (b-a)*psi_i/2+(a+b)/2.
#
# See also script myTest_Gauss-Kronrod.py.
#
#----------------------------------------------------
# Data for GK:

nPoints_GK = 16
psi_16=np.array([-0.9894009, -0.9445750, -0.8656312, -0.7554044, -0.6178762, \
                 -0.4580168, -0.2816036, -0.0950125,  0.0950125,  0.2816036, \
		  0.4580168,  0.6178762,  0.7554044,  0.8656312,  0.9445750, \
		  0.9894009])
w_16  =np.array([ 0.0271525,  0.0622535,  0.0951585,  0.1246290,  0.1495960, \
                  0.1691565,  0.1826034,  0.1894506,  0.1894506,  0.1826034, \
		  0.1691565,  0.1495960,  0.1246290,  0.0951585,  0.0622535, \
		  0.0271525])

def fittedGKintegration(xMin,xMax,fitA,fitB):

   y = np.zeros(nPoints_GK)
   yIntegrated = 0.
   for n in range(nPoints_GK):
      xCrrnt = psi_16[n]*(xMax-xMin)/2 + (xMax+xMin)/2.
      factorA = math.pow(10.,fitA)
      y[n] = factorA*math.pow(xCrrnt,fitB)
      yIntegrated += (xMax-xMin)*w_16[n]*y[n]*xCrrnt
   return y,yIntegrated      

#------------------ End of defined functions -----------------------
#

#====================================================================
#
#--------- Input data, parameters and derived quantities-------------
#

omega_L = omega_Larmor(m_elec,B_mag)                               # rad/sec 
T_larm = 2*pi/omega_L                                              # sec
timeStep = T_larm/stepsNumberOnGyro                                # time step, sec
print ('omega_Larmor= %e rad/sec, T_larm = %e sec, timeStep = %e sec' % \
       (omega_L,T_larm,timeStep))

nLarmorAvrgng=10               # number of averaged Larmor rotations 
#
# Data to integrate transferred momemta over the track:
#
timeStep_c=nLarmorAvrgng*stepsNumberOnGyro*timeStep                # sec      
print ('timeStep_c = %e s' % timeStep_c)

eVrmsTran = np.sqrt(2.*eTempTran*eVtoErg/m_elec)                   # cm/sec
eVrmsLong = np.sqrt(2.*eTempLong*eVtoErg/m_elec)                   # cm/sec
kinEnergy = m_elec*(eVrmsTran**2+eVrmsLong**2)/2.     # kinetic energy; erg
print ('eVrmsTran = %e cm/sec, eVrmsLong = %e cm/sec, kinEnergy = %e eV' % \
       (eVrmsTran,eVrmsLong,ergToEV*kinEnergy))

ro_larmRMS = eVrmsTran/omega_L                                     # cm
print ('ro_larmRMS =%e mkm' % (1.e4*ro_larmRMS))

#
# Electrons are magnetized for impact parameter >> rhoCrit:
#
rhoCrit=math.pow(q_elec**2/(m_elec*omega_L**2),1./3)               # cm
print ('rhoCrit (mkm) = ' , 1.e+4*rhoCrit)

sphereNe=3.
R_e=math.pow(sphereNe/n_e,1./3)                                     # cm
print ('R_e (cm)=%e' % R_e)

ro_Larm = eVrmsTran/omega_L                                         # cm
print ('ro_Larm (cm)=%e' % ro_Larm)

impctPrmtrMin=2.*ro_Larm

#
# Electrons are magnetized for impact parameter >> rhoCrit:
#
rhoCrit=math.pow(q_elec**2/(m_elec*omega_L**2),1./3)               # cm
print ('rhoCrit (mkm) = ' , 1.e+4*rhoCrit)


#============ Important flags ===========================
#
# Taking into account the transfer of momenta for both particles
# (for "classical" only):
dpTransferFlag = 1        # no taking into account if = 0!
#
saveFilesFlag = 0         # no saving if = 0!
#
plotFigureFlag = 1        # plot if = 1!
#
#========================================================
nVion=50
Vion=np.zeros(nVion)
VionLong=np.zeros(nVion)
VionTrnsv=np.zeros(nVion)
VionRel=np.zeros(nVion)

vIonMin=4.e-3*eVrmsTran
vIonMax=10.*eVrmsTran

vIonMinRel=vIonMin/V0
vIonMaxRel=vIonMax/V0
print ('VionMin=%e (vIonMinRel=%e), vIonMax=%e (vIonMaxRel=%e)' % \
       (vIonMin,vIonMinRel,vIonMax,vIonMaxRel))   

vIonLogStep=math.log10(vIonMax/vIonMin)/(nVion-1)

R_debye=np.zeros(nVion)
R_pass=np.zeros(nVion)
R_pass_1=np.zeros(nVion)             # for longT=0. --> eVrmsLong=0.
impctPrmtrMax=np.zeros(nVion)
impctPrmtrMax_1=np.zeros(nVion)      # for longT=0. --> eVrmsLong=0.

for i in range(nVion):
   crrntLogVionRel=math.log10(vIonMinRel)+i*vIonLogStep
   VionRel[i]=math.pow(10.,crrntLogVionRel)
   Vion[i]=VionRel[i]*V0
   VionLong[i]=Vion[i]*np.cos(thetaVi)
   VionTrnsv[i]=Vion[i]*np.sin(thetaVi)
   R_debye[i]=np.sqrt(Vion[i]**2+eVrmsTran**2+eVrmsLong**2)/omega_p
   R_pass[i]=np.sqrt(Vion[i]**2+eVrmsLong**2)*coolPassTime
   R_pass_1[i]=np.sqrt(Vion[i]**2+0.*eVrmsLong**2)*coolPassTime
   help=max(R_debye[i],R_e)
   impctPrmtrMax[i]=min(help,R_pass[i])
   impctPrmtrMax_1[i]=min(help,R_pass_1[i])

#-----------------------------------------------------------------
# Checking of corection of the maximal impact parameter on depence
# of preset number of minimal Larmor turns
# 
larmorTurnsMin=[10,20,30,40]

impctPrmtrMaxCrrctd=np.zeros((nVion,4))
impctPrmtrMaxCrrctdRel=np.zeros((nVion,4))

for n in range (4):
   for i in range(nVion):
      impctPrmtrMaxCrrctd[i,n]=impctPrmtrMax[i]* \
      np.sqrt(1.- (pi*larmorTurnsMin[n]*eVrmsLong/omega_L/impctPrmtrMax[i])**2)
      impctPrmtrMaxCrrctdRel[i,n]=impctPrmtrMaxCrrctd[i,n]/impctPrmtrMax[i]

#      
# First plotting:      
#   
if (plotFigureFlag == 0):   
   fig10 = plt.figure(10)
   plt.semilogx(impctPrmtrMax,impctPrmtrMaxCrrctdRel[:,0],'-r', \
                  impctPrmtrMax,impctPrmtrMaxCrrctdRel[:,1],'-b', \
                  impctPrmtrMax,impctPrmtrMaxCrrctdRel[:,2],'-g', \
                  impctPrmtrMax,impctPrmtrMaxCrrctdRel[:,3],'-m',linewidth=2)
   plt.grid(True)
   hold=True
   plt.xlabel('Maximal Impact parameter $R_{max}$, cm',color='m',fontsize=16)
   plt.ylabel('$R_{max}^{Crrctd}/R_{Max}$',color='m',fontsize=16)
#   plt.xlim([.9*min(impctPrmtrMax),1.1*max(impctPrmtrMax)])
   plt.xlim([1.e-2,1.1*max(impctPrmtrMax)])
   plt.ylim([.986,1.001])
   titleHeader='$R_{max}^{Crrctd}=R_{Max} \cdot [1-(\pi\cdot N_{Larm} \cdot' 
   titleHeader += '\Delta_{e||}/(\omega_{Larm} \cdot R_{max})]^{1/2}$'
   plt.title(titleHeader,color='m',fontsize=16)
   plt.legend([('$N_{Larm}=$%2d' % larmorTurnsMin[0]), \
               ('$N_{Larm}=$%2d' % larmorTurnsMin[1]), \
               ('$N_{Larm}=$%2d' % larmorTurnsMin[2]), \
               ('$N_{Larm}=$%2d' % larmorTurnsMin[3])],loc='lower center',fontsize=14)
   if (saveFilesFlag == 1):
      fig10.savefig('picturesCMA/correctedRmax_fig10cma.png')  
      print ('File "picturesCMA/correctedRmax_fig10cma.png" is written')   

xLimit=[.9*VionRel[0],1.1*VionRel[nVion-1]]

#
# Types of collisions:
#
if (plotFigureFlag == 1):   
   fig3151=plt.figure (3151)
   plt.loglog(VionRel,impctPrmtrMax,'-r', VionRel,impctPrmtrMax_1,'--r', \
      [VionRel[0],VionRel[nVion-1]],[impctPrmtrMin,impctPrmtrMin],'-b',linewidth=2)
   plt.grid(True)
   hold=True
   plt.xlabel('Relative Ion Velocity, $V_i/V_{e0}$',color='m',fontsize=14)
   plt.ylabel('Impact Parameter, cm',color='m',fontsize=14)
   titleHeader= \
             'Types of Collisions: $V_{e0}=%4.2f\cdot10^{%2d}$ cm/s, $B=%6.1f$ Gs'
   plt.title(titleHeader % (mantV0,powV0,fieldB[0]),color='m',fontsize=16)
   plt.xlim(xLimit)
   yLimit=[8.e-4,.6]
   plt.ylim(yLimit)
   plt.plot([relVeTrnsv,relVeTrnsv],yLimit,'--m',linewidth=1)
   plt.text(1.6e-3,5.e-4,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([relVeLong,relVeLong],yLimit,'--m',linewidth=1)
   plt.text(4.4e-5,.0018,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   plt.text(3.e-4,1.75e-3,'$R_{min}=2\cdot<rho_\perp>$',color='k',fontsize=16)
   plt.text(7.e-4,5.e-2,'$R_{max}$',color='k',fontsize=16)
   plt.text(2.85e-5,3.3e-3,'$R_{max}$ $for$ $T_{e||}=0$',color='k',fontsize=16)
   plt.plot([VionRel[0],VionRel[nVion-1]],[20.*rhoCrit,20.*rhoCrit],color='k')   
   plt.text(1.e-4,7.e-3,'Magnetized Collisions',color='r',fontsize=20)
   plt.text(1.e-4,10.e-4,'Adiabatic or Fast Collisions',color='r',fontsize=20)
   plt.text(2.25e-5,.275,'Collisions are Screened',color='r',fontsize=20)
   plt.text(1.6e-5,1.e-3,'$ \cong 20\cdot R_{Crit}$',color='k',fontsize=16)
   if (saveFilesFlag == 1):
      fig3151.savefig('picturesCMA_v7/impctPrmtr_fig3151cma.png')    
      print ('File "picturesCMA_v7/impctPrmtr_fig3151cma.png" is written')   

#------- 11/19/2018 (velocities of electrons are distributed flattenedly) ------

eVlongMin = .01*eVrmsLong
eVlongMax = 5.*eVrmsLong
print ('eVlongMin = %e, eVlongMax = %e' %(eVlongMin,eVlongMax)) 

eVlongNumb = 16
eVlongStep=math.log10(eVlongMax/eVlongMin)/(eVlongNumb-1)

eVlong = np.zeros(eVlongNumb)

for n in range(eVlongNumb):
#   eVlongLog = math.log10(eVlongMin)+n*eVlongStep
#   eVlong[n] = math.pow(10.,eVlongLog)
# xCrrnt = psi_16_ccrnt*(xMax-xMin)/2. + (xMax+xMin)/2.
   eVlong[n] = psi_16[n]*(eVlongMax-eVlongMin)/2. + (eVlongMax+eVlongMin)/2.
#   print ('eVlong(%d) = %e ==> %e' % (n,eVlong[n],eVlong[n]/eVrmsLong))

eVtrnsvMin = .01*eVrmsTran
eVtrnsvMax = 5.*eVrmsTran
print ('eVtrnsvMin = %e, eVtrnsvMax = %e' %(eVtrnsvMin,eVtrnsvMax)) 

eVtrnsvNumb = 16
eVtrnsvStep=math.log10(eVtrnsvMax/eVtrnsvMin)/(eVtrnsvNumb-1)

eVtrnsv = np.zeros(eVtrnsvNumb)

for n in range(eVtrnsvNumb):
#     eVtrnsvLog = math.log10(eVtrnsvMin)+n*eVtrnsvStep
#     eVtrnsv[n] = math.pow(10.,eVtrnsvLog)
# xCrrnt = psi_16_ccrnt*(xMax-xMin)/2. + (xMax+xMin)/2.
     eVtrnsv[n] = psi_16[n]*(eVtrnsvMax-eVtrnsvMin)/2. + (eVtrnsvMax+eVtrnsvMin)/2.
#     print ('eVtrnsv(%d) = %e ==> %e' % (n,eVtrnsv[n],eVtrnsv[n]/eVrmsTran))

R_min = np.zeros((eVlongNumb,eVtrnsvNumb,nVion))
R_pass_map = np.zeros((eVlongNumb,nVion))
R_debye_map = np.zeros((eVlongNumb,eVtrnsvNumb,nVion))
R_min_map = np.zeros((eVtrnsvNumb,nVion))
R_max_map = np.zeros((eVlongNumb,eVtrnsvNumb,nVion))

for i in range(nVion):
  for n in range(eVlongNumb):
     R_pass_map[n,i] = np.sqrt(Vion[i]**2+eVlong[n]**2)*coolPassTime
     for j in range(eVtrnsvNumb):
        R_debye_map[n,j,i]=np.sqrt(Vion[i]**2+eVtrnsv[j]**2+eVlong[n]**2)/omega_p
        help=max(R_debye_map[n,j,i],R_e)
        R_max_map[n,j,i]=min(help,R_pass_map[n,i])

#
# Map for R_pass:
#
X = np.zeros((eVlongNumb,nVion))
Y = np.zeros((eVlongNumb,nVion))
Z = np.zeros((eVlongNumb,nVion))

for i in range(nVion):
   for n in range(eVlongNumb):
      X[n,i] = np.log10(VionRel[i])
      Y[n,i] = np.log10(eVlong[n]/eVrmsLong)
      Z[n,i] = np.log10(R_pass_map[n,i])                    

yLimit = [np.log10(min(eVlong)/eVrmsLong),np.log10(max(eVlong)/eVrmsLong)]
# print ('yLimMin = %e, yLimMax = %e' %(yLimit[0],yLimit[1]))

locs = np.zeros(10)
abels = np.zeros(10)
myYticks = np.zeros(10)
if (plotFigureFlag == 1): 
      figCrrnt = plt.figure(3152)
      ax = figCrrnt.add_subplot(111)                                       # for contours plotting
      mapCrrnt = ax.contourf(X,Y,Z,cmap='jet') 
#
# Curves don't make much sense on this graph because ordinate-axis for them is
# impact parameter and not electron longitudinal velocity! Nevertheless...
#      
#      plt.plot(np.log10(VionRel),np.log10(R_pass),'-r',linewidth=2)
#      plt.plot(np.log10(VionRel),np.log10(R_debye),'-w',linewidth=2)
#      plt.plot(np.log10(VionRel),1.02*np.log10(impctPrmtrMax),'-k',linewidth=2)
#      plt.text(-4.6,-1.3,'$R_{Debye}$',color='w',fontsize=16)
#      plt.text(-3.5,-.5,'$R_{Pass}$ for\n$V_{e||}=\Delta V_{e||}$',color='r',fontsize=16)
#      plt.text(-4.,-1.9,'$R_{max}$ = $min(R_{Debye},R_{Pass})$',color='k',fontsize=16)
#
      titleHeader = \
                 'Map for $R_{Pass}$ = $\sqrt{V_{ion}^2+V_{e||}^2}\cdot T_{cool}$, cm (Log Scale)'
      plt.title(titleHeader,color='m',fontsize=14)
      plt.xlabel('Ion Velocity,  $log_{10}(V_{ion}/V_0$)',color='m',fontsize=14)
      plt.ylabel('Longitudinal Velocity (Log Scale), $V_{e||}/\Delta V_{e||}$', \
                 color='m',fontsize=14)
      plt.ylim(yLimit)
      figCrrnt.colorbar(mapCrrnt)
      plt.grid(True)
      locs,labels = plt.yticks()
      for i in range(10):
         myYticks[i] = 10.**locs[i]
         myYticks[i] = "{:5.3f}".format(myYticks[i])
      ax.set_yticklabels(myYticks)
      if (saveFilesFlag == 1):
         fileName = 'picturesCMA_v7/mapRpass_fig3152cma.png'
         figCrrnt.savefig(fileName) 
         print ('File "',fileName,'" is written')

for j in range(eVtrnsvNumb):
   powVeTran=round(np.log10(eVtrnsv[j])) 
   mantVeTran=eVtrnsv[j]/(10**powVeTran) 
   for i in range(nVion):
      for n in range(eVlongNumb):
#         Z[n,i] = np.log10(R_debye_map[n,j,i])                    
         Z[n,i] = R_debye_map[n,j,i]                    
   if (plotFigureFlag == 1): 
      figCrrnt = plt.figure(4100+j)
      ax = figCrrnt.add_subplot(111)                   # for contours plotting
      mapCrrnt = ax.contourf(X,Y,Z,cmap='jet') 
#      mapCrrnt1 = ax.contour(X,Y,Z,7,colors='white') 
#      plt.clabel(mapCrrnt1,fmt='%4.2f',inline=True)
      titleHeader = '$R_{Debye}$ = $\sqrt{V_{ion}^2+V_{e||}^2+V_{e\perp}^2}/\omega_p$, cm'
      plt.title(titleHeader,color='m',fontsize=14)
      plt.xlabel('Ion Velocity,  $log_{10}(V_{ion}/V_0$)',color='m',fontsize=14)
      plt.ylabel(' Longitudinal Velocity,cm: $V_{e||}/\Delta V_{e||}$ (Log Scale)', \
	         color='m',fontsize=14)
      plt.text(-4,.5,('$V_{ion \perp}=%3.1f\cdot10^{%2d}$cm/s' % (mantVeTran,powVeTran)), \
               color='w',fontsize=14)
      plt.ylim(yLimit)
      figCrrnt.colorbar(mapCrrnt)
      plt.grid(True)
# Array myYticks is defined firstly for figure 3152 for map of R_pass:      
      ax.set_yticklabels(myYticks)
      if (saveFilesFlag == 1):
         fileName = 'picturesCMA_v7/mapRpass_fig'+str(4100+j)+'cma.png'
         figCrrnt.savefig(fileName) 
         print ('File "',fileName,'" is written')

for j in range(eVtrnsvNumb):
   powVeTran=round(np.log10(eVtrnsv[j])) 
   mantVeTran=eVtrnsv[j]/(10**powVeTran) 
   for i in range(nVion):
      for n in range(eVlongNumb):
#         Z[n,i] = np.log10(R_max_map[n,j,i])                    
         Z[n,i] = R_max_map[n,j,i]                    
   if (plotFigureFlag == 1): 
      figCrrnt = plt.figure(4300+j)
      ax = figCrrnt.add_subplot(111)                  # for contours plotting
      mapCrrnt = ax.contourf(X,Y,Z,cmap='jet') 
#      mapCrrnt1 = ax.contour(X,Y,Z,7,colors='white') 
#      plt.clabel(mapCrrnt1,fmt='%4.2f',inline=True)
      titleHeader = 'Maximal Impact Parameter $R_{max}$, cm' 
      plt.title(titleHeader,color='m',fontsize=14)
      plt.xlabel('Ion Velocity,  $log_{10}(V_{ion}/V_0$)',color='m',fontsize=14)
      plt.ylabel('Longitudinal Velocity, $V_{e||}/\Delta V_{e||}$ (Log Scale)', \
                 color='m',fontsize=14)
      plt.text(-4,.5,('$V_{e\perp}=%3.1f\cdot10^{%2d}$cm/s' % (mantVeTran,powVeTran)), \
               color='w',fontsize=14)
      plt.ylim(yLimit)
      figCrrnt.colorbar(mapCrrnt)
      plt.grid(True)
# Array myYticks is defined firstly for figure 3152 for map of R_pass:      
      ax.set_yticklabels(myYticks)
      if (saveFilesFlag == 1):
         fileName = 'picturesCMA_v7/mapRmax_fig'+str(4300+j)+'cma.png'
         figCrrnt.savefig(fileName) 
         print ('File "',fileName,'" is written')
	 
X = np.zeros((eVtrnsvNumb,eVlongNumb))		     
Y = np.zeros((eVtrnsvNumb,eVlongNumb))		     
Z = np.zeros((eVtrnsvNumb,eVlongNumb))		     

for i in range(eVtrnsvNumb):
   for n in range(eVlongNumb):
      X[i,n] = np.log10(eVtrnsv[i]/eVrmsTran)		     
      Y[i,n] = np.log10(eVlong[n]/eVrmsLong)		     

for j in range(0,nVion,3):
   powVion=round(np.log10(Vion[j])) 
   mantVion=Vion[j]/(10**powVion) 
   for i in range(eVtrnsvNumb):
      for n in range(eVlongNumb):
#         Z[n,i] = np.log10(R_debye_map[n,i,j])                    
         Z[i,n] = R_debye_map[n,i,j]                    
   if (plotFigureFlag == 0): 
      figCrrnt = plt.figure(4200+j)
      ax = figCrrnt.add_subplot(111)                        # for contours plotting
      mapCrrnt = ax.contourf(X,Y,Z,cmap='jet') 
#      mapCrrnt1 = ax.contour(X,Y,Z,7,colors='white') 
#      plt.clabel(mapCrrnt1,fmt='%4.2f',inline=True)
      titleHeader = '$R_{Debye}$ = $\sqrt{V_{ion}^2+V_{e||}^2+V_{e\perp}^2}/\omega_p$, cm'
      plt.title(titleHeader,color='m',fontsize=14)
      plt.title(titleHeader,color='m',fontsize=14)
#      plt.xlabel('Ion Velocity,  $log_{10}(V_{ion}/V_0$)',color='m',fontsize=14)
      plt.xlabel('Electron Transversal Velocity, $log_{10}(V_{e\perp}/\Delta V_{e\perp}$)',\
      color='m',fontsize=14)
      plt.ylabel('Longitudinal Velocity, $V_{e||}/\Delta V_{e||}$ (Log Scale)', \
                 color='m',fontsize=14)
      plt.text(-1,-1.25,('$V_{ion \perp}=%3.1f\cdot10^{%2d}$cm/s' % (mantVion,powVion)), \
               color='w',fontsize=14)
      plt.ylim(yLimit)
      figCrrnt.colorbar(mapCrrnt)
      plt.grid(True)
# Array myYticks is defined firstly for figure 3152 for map of R_pass:      
      ax.set_yticklabels(myYticks)
      if (saveFilesFlag == 1):
         fileName = 'picturesCMA_v7/mapRmax_fig'+str(4200+j)+'cma.png'
         figCrrnt.savefig(fileName) 
         print ('File "',fileName,'" is written')

for j in range(0,nVion,3):
   powVion=round(np.log10(Vion[j])) 
   mantVion=Vion[j]/(10**powVion) 
   for i in range(eVtrnsvNumb):
      for n in range(eVlongNumb):
#         Z[n,i] = np.log10(R_max_map[n,i,j])                    
         Z[i,n] = R_max_map[n,i,j]                    
   if (plotFigureFlag == 0): 
      figCrrnt = plt.figure(4400+j)
      ax = figCrrnt.add_subplot(111)                        # for contours plotting
      mapCrrnt = ax.contourf(X,Y,Z,cmap='jet') 
      mapCrrnt1 = ax.contour(X,Y,Z,7,colors='white') 
      plt.clabel(mapCrrnt1,fmt='%5.3f',inline=True)
      titleHeader = 'Maximal Impact Parameter $R_{max}$, cm' 
      plt.title(titleHeader,color='m',fontsize=14)
#      plt.xlabel('Ion Velocity,  $log_{10}(V_{ion}/V_0$)',color='m',fontsize=14)
      plt.xlabel('Electron Transversal Velocity, $log_{10}(V_{e\perp}/\Delta V_{e\perp}$)',\
      color='m',fontsize=14)
      plt.ylabel('Longitudinal Velocity, $V_{e||}/\Delta V_{e||}$ (Log Scale)', \
                 color='m',fontsize=14)
      plt.text(-1,-1.25,('$V_{ion \perp}=%3.1f\cdot10^{%2d}$cm/s' % (mantVion,powVion)), \
               color='w',fontsize=14)
      plt.ylim(yLimit)
      figCrrnt.colorbar(mapCrrnt)
      plt.grid(True)
# Array myYticks is defined firstly for figure 3152 for map of R_pass:      
      ax.set_yticklabels(myYticks)
      if (saveFilesFlag == 1):
         fileName = 'picturesCMA_v7/mapRmax_fig'+str(4400+j)+'cma.png'
         figCrrnt.savefig(fileName) 
         print ('File "',fileName,'" is written')
	 
gaussDnst = np.zeros(50)
gaussVelct = np.zeros(50)

for i in range(50):
   longVlog = math.log10(eVlongMin)+i*math.log10(eVlongMax/eVlongMin)/49
   longVcrrnt = math.pow(10.,longVlog)
   gaussVelct[i] = longVcrrnt/eVrmsLong
   gaussDnst[i] = math.exp(-(gaussVelct[i]-0.*V0/eVrmsLong)**2)	 

gaussVelctSlctd = np.zeros(eVlongNumb)
gaussDnstSlctd = np.zeros(eVlongNumb)

for i in range(eVlongNumb):
   gaussVelctSlctd[i] = eVlong[i]/eVrmsLong
   gaussDnstSlctd[i] = math.exp(-(gaussVelctSlctd[i])**2)	 

if (plotFigureFlag == 0): 
   plt.figure(2)
   plt.plot(gaussVelctSlctd,gaussDnstSlctd,'or',gaussVelct,gaussDnst,'-xb')
   plt.xlabel('Relative Velocity, $V_e/\Delta V_e$',color='m',fontsize=14)
   plt.ylabel('Relative Density',color='m',fontsize=14)
   plt.title('Gaussian Density Distribution',color='m',fontsize=14)
   plt.legend(['16 Points for GK Integration','Standard Set of 50 Points'],loc='upper right',fontsize=14)
   plt.grid(True)
   if (saveFilesFlag == 1):
      fileName = 'picturesCMA_v7/densityDistribution_fig2cma.png'
      figCrrnt.savefig(fileName) 
      print ('File "',fileName,'" is written')

eDenst = np.zeros((eVtrnsvNumb,eVlongNumb))

for i in range(eVtrnsvNumb):
   for n in range(eVlongNumb):
      eDenst[i,n] = math.exp(-(eVtrnsv[i]/eVrmsTran)**2)* \
                    math.exp(-(eVlong[n]/eVrmsLong)**2)  
		    
X = np.zeros((eVtrnsvNumb,eVlongNumb))		     
Y = np.zeros((eVtrnsvNumb,eVlongNumb))		     

for i in range(eVtrnsvNumb):
   for n in range(eVlongNumb):
      X[i,n] = eVtrnsv[i]/eVrmsTran		     
      Y[i,n] = eVlong[n]/eVrmsLong		     

if (plotFigureFlag == 0): 
   figCrrnt = plt.figure(3)
   ax = figCrrnt.add_subplot(111)                                       # for contours plotting
   mapCrrnt = ax.contourf(X,Y,np.log10(eDenst),cmap='jet') 
   mapLevels = ax.contour(X,Y,np.log10(eDenst),15,colors='black') 
   plt.clabel(mapLevels,fmt='%4.2f',inline=True)
   titleHeader = '"Flattened" Electron Density ($log_{10}$): $T_{\perp}$=%3.1f eV, $T_{||}$=%3.1f meV' 
   plt.title(titleHeader % (eTempTran,1.e3*eTempLong),color='m',fontsize=12)
   plt.xlabel('Transverse Velocity,  $V_{e\perp}/\Delta V_{e\perp}$',color='m',fontsize=14)
   plt.ylabel('Longitudinal Velocity,  $V_{e||}/\Delta V_{e||}$',color='m',fontsize=14)
   figCrrnt.colorbar(mapCrrnt)
   plt.grid(True)

plt.show()

sys.exit()

#
# Magnetized  collisions:
#
if (plotFigureFlag == 0):   
   fig209=plt.figure (209)
   plt.loglog(VionRel,R_debye,'-r',VionRel,R_pass,'-b', \
                                   VionRel,R_pass_1,'--b',linewidth=2)
   plt.grid(True)
   hold=True
   plt.plot([VionRel[0],VionRel[nVion-1]],[R_e,R_e],color='m',linewidth=2)   
   plt.xlabel('Relative Ion Velocity, $V_i/V_{e0}$',color='m',fontsize=16)
   plt.ylabel('$R_{Debye}$, $R_{Pass}$, $R_e$, cm',color='m',fontsize=16)
#   titleHeader='Magnetized Collision: $R_{Debye}$, $R_{Pass}$, $R_e$: $V_{e0}=%5.3f\cdot10^{%2d}$cm/s'
#   plt.title(titleHeader % (mantV0,powV0),color='m',fontsize=16)
   plt.title('Magnetized Collisions: $R_{Debye}$, $R_{Pass}$, $R_e$',color='m',fontsize=16)
   plt.xlim(xLimit)
   yLimit=[1.e-3,10.]
   plt.ylim(yLimit)
   plt.plot([relVeTrnsv,relVeTrnsv],yLimit,'--m',linewidth=1)
   plt.text(1.6e-3,5.5e-4,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([relVeLong,relVeLong],yLimit,'--m',linewidth=1)
   plt.text(4.4e-5,0.001175,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   plt.text(3.e-5,2.45e-3,'$R_e$',color='k',fontsize=16)
   plt.text(3.e-5,5.e-2,'$R_{Debye}$',color='k',fontsize=16)
   plt.text(3.e-5,1.8e-2,'$R_{Pass}$',color='k',fontsize=16)
   plt.text(4.5e-5,4.8e-3,'$R_{Pass}$ $for$ $T_{e||}=0$',color='k',fontsize=16)
   plt.text(8.3e-5,4.0,('$V_{e0}=%5.3f\cdot10^{%2d}$cm/s' % (mantV0,powV0)), \
            color='m',fontsize=16)
   if (saveFilesFlag == 1):
      fig209.savefig('picturesCMA/rDebye_rLikeDebye_rPass_fig209cma.png')    
      print ('File "picturesCMA/rDebye_rLikeDebye_rPass_fig209cma.png" is written')   

#
# Coulomb logarithm evaluation:
#

clmbLog = np.zeros(nVion)

for i in range(nVion):
   clmbLog[i] = math.log(impctPrmtrMax[i]/impctPrmtrMin)
#   clmbLog[i] = math.log(impctPrmtrMax_1[i]/impctPrmtrMin)

if (plotFigureFlag == 0):   
   fig3155=plt.figure (3155)
   plt.semilogx(VionRel,clmbLog,'-xr',linewidth=2)
   plt.xlabel('Relative Ion Velocity, $V_i/V_{e0}$',color='m',fontsize=14)
   plt.ylabel('Coulomb Logarithm $L_c$',color='m',fontsize=14)
   plt.title('Coulomb Logarithm: $L_c$ = $ln(R_{max}/R_{min})$',color='m',fontsize=16)
   yLimit=[min(clmbLog)-.1,max(clmbLog)+.1]
   plt.ylim(yLimit)
   plt.plot([relVeTrnsv,relVeTrnsv],yLimit,'--m',linewidth=1)
   plt.text(1.6e-3,5.,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([relVeLong,relVeLong],yLimit,'--m',linewidth=1)
   plt.text(3.4e-5,5.,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   if (saveFilesFlag == 1):
      fig3155.savefig('picturesCMA_v7/coulombLogrthm_fig3155cma.png')    
      print ('File "picturesCMA_v7/coulombLogrthm_fig3155cma.png" is written')   

# plt.show()

# sys.exit()

#
# matrix for electron with .5*timeStep_c:
#
matr_elec_c=guidingCenter_Matrix(.5*timeStep_c)        
#
# matrix for ion with mass M_ion  and .5*timeStep_c:
#
matr_ion_c=drift_Matrix(M_ion,.5*timeStep_c) 

larmorTurns = 10
nImpctPrmtr = 50

rhoMin = impctPrmtrMin
rhoMax = np.zeros(nVion)
log10rhoMin = math.log10(rhoMin)
crrntImpctPrmtr = np.zeros(nImpctPrmtr)
halfLintr = np.zeros((nImpctPrmtr,nVion))
pointAlongTrack = np.zeros((nImpctPrmtr,nVion))

totalPoints = 0
for i in range(nVion):
   rhoMax[i] = impctPrmtrMax[i]* \
   np.sqrt(1.- (pi*larmorTurns*eVrmsLong/omega_L/impctPrmtrMax[i])**2)
   rhoMax[i] = impctPrmtrMax[i]
#   rhoMax[i] = impctPrmtrMax_1[i]              # for checking!
#   print ('rhoMax(%d) = %e' % (i,rhoMax[i]))
   log10rhoMax = math.log10(rhoMax[i])
   log10rhoStep = (log10rhoMax-log10rhoMin)/(nImpctPrmtr)
#   print ('Vion(%d) = %e, rhoMax = %e' % (i,Vion[i],rhoMax[i]))
   for n in range(nImpctPrmtr):
      log10rhoCrrnt = log10rhoMin+(n+0.5)*log10rhoStep 
      rhoCrrnt = math.pow(10.,log10rhoCrrnt)
#      print ('    rhoCrrnt(%d) = %e' % (n,rhoCrrnt))
      halfLintr[n,i] = np.sqrt(rhoMax[i]**2-rhoCrrnt**2)   # half length of interaction; cm
      timeHalfPath = halfLintr[n,i]/eVrmsLong     # 0.5 time of interaction; sec
      numbLarmor = int(2.*timeHalfPath/T_larm)             
      pointAlongTrack[n,i] = int(2.*timeHalfPath/timeStep_c)
      totalPoints += pointAlongTrack[n,i]
#      print ('     %d: rhoCrrnt = %e, numbLarmor = %d, pointAlongTrack = %d' % \
#            (n,rhoCrrnt,numbLarmor,pointAlongTrack[n,i]))
# print ('totalPoints = %d' % totalPoints)

totalPoints = int(totalPoints)
nnTotalPoints=np.arange(0,2*totalPoints-1,1)

arrayA=np.zeros(2*totalPoints)      
arrayB=np.zeros(2*totalPoints)     
bCrrnt_c = np.zeros(2*totalPoints)

#
# Variables for different testing:
#
b_gc = np.zeros(totalPoints)
action_gc = np.zeros(totalPoints)
C1test = np.zeros(totalPoints)
C2test = np.zeros(totalPoints)
C3test = np.zeros(totalPoints)
b_ME = np.zeros(totalPoints)
D1test = np.zeros(totalPoints)
D2test = np.zeros(totalPoints)
qTest = np.zeros(totalPoints)                
action_ME = np.zeros(totalPoints)
actn_gc_ME_rel = np.zeros(totalPoints)
indxTest = 0

rhoInit = np.zeros((nImpctPrmtr,nVion))
#
# "Classical" approach:
#
deltaPx_c = np.zeros((nImpctPrmtr,nVion))  
deltaPy_c = np.zeros((nImpctPrmtr,nVion))  
deltaPz_c = np.zeros((nImpctPrmtr,nVion)) 
ionVx_c = np.zeros((nImpctPrmtr,nVion)) 
ionVy_c = np.zeros((nImpctPrmtr,nVion)) 
ionVz_c = np.zeros((nImpctPrmtr,nVion)) 
deltaEnrgIon_c = np.zeros((nImpctPrmtr,nVion))
#   
# "Magnus Expand" approach:
#
deltaPx_m = np.zeros((nImpctPrmtr,nVion))  
deltaPy_m = np.zeros((nImpctPrmtr,nVion))  
deltaPz_m = np.zeros((nImpctPrmtr,nVion))  
ionVx_m = np.zeros((nImpctPrmtr,nVion)) 
ionVy_m = np.zeros((nImpctPrmtr,nVion)) 
ionVz_m = np.zeros((nImpctPrmtr,nVion)) 
deltaEnrgIon_m = np.zeros((nImpctPrmtr,nVion))
#   
# Comparison of approaches (ratio deltaEnrgIon_c/deltaEnrgIon_m):
#
deltaPx_c_m = np.zeros((nImpctPrmtr,nVion)) 
deltaPy_c_m = np.zeros((nImpctPrmtr,nVion)) 
deltaPz_c_m = np.zeros((nImpctPrmtr,nVion)) 
dEion_c_m = np.zeros((nImpctPrmtr,nVion))

#
# Factor to calculate transferred energy to ion
# (the friction force is defined by this transfered energy): 
#
deFactor = 0.5/M_ion                               # 1/g

frctnForce_cSM = np.zeros(nVion)      # integration, using Simpson method
frctnForce_mSM = np.zeros(nVion)      # integration, using Simpson method

numberWrongSign_c=0
numberWrongSign_m=0

posSignDeltaEnrgIon_c=0
negSignDeltaEnrgIon_c=0
posSignDeltaEnrgIon_m=0
negSignDeltaEnrgIon_m=0

timeRun = np.zeros(nVion)     
totalTimeRun = 0.

indx = 0

###################   Main simulation   ###############
#
# Initial electron longitudinal and transversal velocities do not 
# distributted, but equel to eVrmsLong and eVrmsTran correspondingly
#
for i in range(nVion):
# Taking into account the corection of the maximal impact parameter
# on depence of preset number of minimal Larmor turns:
   rhoMax[i] = impctPrmtrMax[i]* \
   np.sqrt(1.- (pi*larmorTurns*eVrmsLong/omega_L/impctPrmtrMax[i])**2)
# Without taking into account the corection of the maximal impact parameter
# on depence of preset number of minimal Larmor turns:
   rhoMax[i] = impctPrmtrMax[i]
#   rhoMax[i] = impctPrmtrMax_1[i]              # for checking!
   log10rhoMax = math.log10(rhoMax[i])
   log10rhoStep = (log10rhoMax-log10rhoMin)/(nImpctPrmtr)
#   print ('Vion(%d) = %e, rhoMax = %e' % (i,Vion[i],rhoMax[i]))
   timeStart=os.times()
   for n in range(nImpctPrmtr):
      log10rhoCrrnt = log10rhoMin+(n+0.5)*log10rhoStep 
      rhoCrrnt = math.pow(10.,log10rhoCrrnt)
#      rhoInit[i*nImpctPrmtr+n] = rhoCrrnt
      rhoInit[n,i] = rhoCrrnt
      halfLintr[n,i] = np.sqrt(rhoMax[i]**2-rhoCrrnt**2)   # half length of interaction; cm
      z_ionCrrnt_c = np.zeros(6)      # Zeroing out of vector for ion ("GC"-approach) 
      z_elecCrrnt_c = np.zeros(6)     # Zeroing out of vector for electron ("GC"-approach)
      z_ionCrrnt_m = np.zeros(6)      # Zeroing out of vector for ion ("ME"-approach) 
      z_elecCrrnt_m = np.zeros(6)     # Zeroing out of vector for electron ("ME"-approach)
# Zeroing out of "guiding center" vector for electron (both approaches):
      z_elecCrrnt_gc_c = np.zeros(6)  
      z_elecCrrnt_gc_m = np.zeros(6)  
# Current values of transfered momemta 
# (second index numerates "Guiding Center", (if 0) and 
#                         "Magnus Expantion" (if 1) approaches: 
      dpCrrnt = np.zeros((3,2))
# Intermediate arrays:
      dpIon_c = np.zeros(3) 
      dpIon_m = np.zeros(3) 
      dpElec_c = np.zeros(3) 
      dpElec_m = np.zeros(3) 
# Current initial vector for electron:
      z_elecCrrnt_c[Ix] = rhoCrrnt                     # x, cm
      z_elecCrrnt_c[Iz] = -halfLintr[n,i]              # z, cm
      z_elecCrrnt_c[Ipy] = m_elec*eVrmsTran            # py, g*cm/sec
      z_elecCrrnt_c[Ipz] = m_elec*eVrmsLong            # pz, g*cm/sec
      z_elecCrrnt_m[Ix] = rhoCrrnt                     # x, cm
      z_elecCrrnt_m[Iz] = -halfLintr[n,i]              # z, cm
      z_elecCrrnt_m[Ipy] = m_elec*eVrmsTran            # py, g*cm/sec
      z_elecCrrnt_m[Ipz] = m_elec*eVrmsLong            # pz, g*cm/sec
# Current initial vector for ion velocity for both approaches:
      ionVx_c[n,i] = VionTrnsv[i]*np.cos(phiVi)
      ionVy_c[n,i] = VionTrnsv[i]*np.sin(phiVi) 
      ionVz_c[n,i] = VionLong[i]
      ionVx_m[n,i] = VionTrnsv[i]*np.cos(phiVi)
      ionVy_m[n,i] = VionTrnsv[i]*np.sin(phiVi) 
      ionVz_m[n,i] = VionLong[i]
# transfer to system of guiding center:
      z_elecCrrnt_gc_c=toGuidingCenter(z_elecCrrnt_c)  
      z_elecCrrnt_gc_m=toGuidingCenter(z_elecCrrnt_m)  
#
# Main loop along the each track:
#
      for k in range(int(pointAlongTrack[n,i])):
#
# Dragging both particles through first half of the step of the track:
#
         z_elecCrrnt_gc_c = np.dot(matr_elec_c,z_elecCrrnt_gc_c) # electron
         z_elecCrrnt_gc_m = np.dot(matr_elec_c,z_elecCrrnt_gc_m) # electron
         z_ionCrrnt_c = np.dot(matr_ion_c,z_ionCrrnt_c)          # ion
         z_ionCrrnt_m = np.dot(matr_ion_c,z_ionCrrnt_m)          # ion
# transfer from system of guiding center: 
         z_elecCrrnt_c=fromGuidingCenter(z_elecCrrnt_gc_c)     
         z_elecCrrnt_m=fromGuidingCenter(z_elecCrrnt_gc_m)     
# Current distance between ion and electron; cm:
         bCrrnt_c[indx]=np.sqrt((z_ionCrrnt_c[0]-z_elecCrrnt_c[0])**2+ \
                                (z_ionCrrnt_c[2]-z_elecCrrnt_c[2])**2+ \
                                (z_ionCrrnt_c[4]-z_elecCrrnt_c[4])**2)
# Current values of parameters A,B:  
         arrayA[indx] = math.log10(ro_Larm/bCrrnt_c[indx])     
         arrayB[indx] = math.log10((q_elec**2/bCrrnt_c[indx])/kinEnergy)
         indx += 1
#
# Dragging both particles through interaction during this step of track
# (for both approaches):
#
#    "Guiding Center":
         dpIon_c,dpElec_c,action,b_gc_c = \
                guidingCenterCollision(z_elecCrrnt_gc_c,z_ionCrrnt_c,timeStep_c) 
#    "Magnus Expantion":
         dpIon_m,dpElec_m,actionME,dy_gc_m,C1,C2,C3,b,D1,D2,q = \
                MagnusExpansionCollision(z_elecCrrnt_gc_m,z_ionCrrnt_m,timeStep_c) 
# Save data for testing:
         b_gc[indxTest] = b_gc_c           # "Guiding Center" approach
         action_gc[indxTest] = action      # -"- -"- -"- -"- -"- -"- 
         C1test[indxTest] = C1             # "Magnus expansion" approach
         C2test[indxTest] = abs(C2)        # -"- -"- -"- -"- -"- -"-
         C3test[indxTest] = C3             # -"- -"- -"- -"- -"- -"-
         b_ME[indxTest] = b                # -"- -"- -"- -"- -"- -"-
         D1test[indxTest] = D1             # -"- -"- -"- -"- -"- -"-
         D2test[indxTest] = D2             # -"- -"- -"- -"- -"- -"-
         qTest[indxTest] = q               #-"- -"- -"- -"- -"- -"- 
         action_ME[indxTest] = actionME    #-"- -"- -"- -"- -"- -"- 
         indxTest += 1  
         indxTestMax = indxTest           
#
# Taking into account transfer of momentum for both particles:
#
         if (dpTransferFlag == 1):
            for ic in range(3):
               z_ionCrrnt_c[2*ic+1] += dpIon_c[ic]   
               z_elecCrrnt_c[2*ic+1] += dpElec_c[ic]
               z_ionCrrnt_m[2*ic+1] += dpIon_m[ic]   
               z_elecCrrnt_m[2*ic+1] += dpElec_m[ic]
# transfer to system of guiding center:
         z_elecCrrnt_gc_c=toGuidingCenter(z_elecCrrnt_c)  
         z_elecCrrnt_gc_m=toGuidingCenter(z_elecCrrnt_m)  
# Accumulation of the transfered momenta to ion along the track for both approaches:  
         for ic in range(3):
#	    if i == 0:
#	       print ('dpIon_c[%2d] = %20.14e, dpIon_m[%2d] = %20.14e' % \
#	             (ic,dpIon_c[ic],ic,dpIon_m[ic]))
            dpCrrnt[ic,0] += dpIon_c[ic]       # "Guiding Center", g*cm/sec  
            dpCrrnt[ic,1] += dpIon_m[ic]       # "Magnus Expansion", g*cm/sec  
#
# Ion's elocity change along the track - both approaches: 
#
         ionVx_c[n,i] += dpCrrnt[0,0]/M_ion                    # cm/sec
         ionVy_c[n,i] += dpCrrnt[1,0]/M_ion                    # cm/sec
         ionVz_c[n,i] += dpCrrnt[2,0]/M_ion                    # cm/sec
         ionVx_m[n,i] += dpCrrnt[0,1]/M_ion                    # cm/sec
         ionVy_m[n,i] += dpCrrnt[1,1]/M_ion                    # cm/sec
         ionVz_m[n,i] += dpCrrnt[2,1]/M_ion                    # cm/sec
#
# Dragging both particles through second half of the step of the track:
#
         z_elecCrrnt_gc_c = np.dot(matr_elec_c,z_elecCrrnt_gc_c)     # electron
         z_ionCrrnt_c = np.dot(matr_ion_c,z_ionCrrnt_c)              # ion
         z_elecCrrnt_gc_m = np.dot(matr_elec_c,z_elecCrrnt_gc_m)     # electron
         z_ionCrrnt_m = np.dot(matr_ion_c,z_ionCrrnt_m)              # ion
# transfer from system of guiding center: 
         z_elecCrrnt_c=fromGuidingCenter(z_elecCrrnt_gc_c)     
         z_elecCrrnt_m=fromGuidingCenter(z_elecCrrnt_gc_m)     
# Current distance between ion and electron; cm:
         bCrrnt_c[indx]=np.sqrt((z_ionCrrnt_c[0]-z_elecCrrnt_c[0])**2+ \
                                (z_ionCrrnt_c[2]-z_elecCrrnt_c[2])**2+ \
                                (z_ionCrrnt_c[4]-z_elecCrrnt_c[4])**2)
# Current values of parameters A,B:  
         arrayA[indx] = math.log10(ro_Larm/bCrrnt_c[indx])     
         arrayB[indx] = math.log10((q_elec**2/bCrrnt_c[indx])/kinEnergy)
         indx += 1
#
# Transferred momenta along the track - "Guiding Center" approach:
#
      deltaPx_c[n,i] = dpCrrnt[0,0]                         # dpx, g*cm/sec
#      if deltaPx_c[n,i] <= 0.: 
#         print ('deltaPx_c[%2d,%2d] = %e, dpCrrnt[%2d,%2d] = %e' % \
#	       (n,i,deltaPx_c[n,i],n,i,dpCrrnt[0,0]))
      deltaPy_c[n,i] = dpCrrnt[1,0]                        # dpy, g*cm/sec 
#      if deltaPy_c[n,i] <= 0.: 
#         print ('deltaPy_c[%2d,%2d] = %e' % (n,i,deltaPy_c[n,i]))
      deltaPz_c[n,i] = dpCrrnt[2,0]                         # dpz, g*cm/sec 
#      if deltaPz_c[n,i] <= 0.: 
#         print ('deltaPz_c[%2d,%2d] = %e' % (n,i,deltaPz_c[n,i]))

# Incorrect value:
#      deltaEnrgIon_c[n,i] = (dpCrrnt[0,0]**2+dpCrrnt[1,0]**2+dpCrrnt[2,0]**2)* \
#                             deFactor/eVtoErg                     # eV
# Correct value:
      crrntDeltaEnrg = (dpCrrnt[0,0]*ionVx_c[n,i]+ \
                         dpCrrnt[1,0]*ionVy_c[n,i]+ \
			 dpCrrnt[2,0]*ionVz_c[n,i])*deFactor/eVtoErg   # eV
      absDeltaEnrgIon_c = abs(crrntDeltaEnrg)
      if (crrntDeltaEnrg != 0.):
         signDeltaEnrgIon_c = crrntDeltaEnrg/abs(crrntDeltaEnrg)
      deltaEnrgIon_c[n,i] = crrntDeltaEnrg
      if (deltaEnrgIon_c[n,i] > 0.):
         posSignDeltaEnrgIon_c += 1
      else:
         negSignDeltaEnrgIon_c += 1
#
# Transferred momenta along the track - "Magnus expansion" approach:
#
      deltaPx_m[n,i] = dpCrrnt[0,1]                # dpx, g*cm/sec
#      if deltaPx_m[n,i] <= 0.: 
#         print ('deltaPx_m[%2d,%2d] = %e' % (n,i,deltaPx_m[n,i]))
      deltaPy_m[n,i] = dpCrrnt[1,1]
#      if deltaPy_m[n,i] <= 0.: 
#         print ('deltaPy_m[%2d,%2d] = %e' % (n,i,deltaPy_m[n,i]))
      deltaPz_m[n,i] = dpCrrnt[2,1] 
#      if deltaPz_m[n,i] <= 0.: 
#         print ('deltaPz_m[%2d,%2d] = %e' % (n,i,deltaPz_m[n,i]))

# Incorrect value:
#      deltaEnrgIon_m[n,i] = (dpCrrnt[0,1]**2+dpCrrnt[1,1]**2+dpCrrnt[2,1]**2)* \
#                            deFactor/eVtoErg                     # eV
# Correct value absolute value):
      crrntDeltaEnrg = (dpCrrnt[0,1]*ionVx_m[n,i]+ \
                         dpCrrnt[1,1]*ionVy_m[n,i]+ \
			 dpCrrnt[2,1]*ionVz_m[n,i])*deFactor/eVtoErg   # eV
      absDeltaEnrgIon_m = abs(crrntDeltaEnrg)
      if (crrntDeltaEnrg != 0.):
         signDeltaEnrgIon_m = crrntDeltaEnrg/abs(crrntDeltaEnrg)
      deltaEnrgIon_m[n,i] = crrntDeltaEnrg
      if (deltaEnrgIon_m[n,i] > 0.):
         posSignDeltaEnrgIon_m += 1
      else:
         negSignDeltaEnrgIon_m += 1
#
# Comparison of the approaches (%):
#
      if (deltaPx_m[n,i] != 0.):
         deltaPx_c_m[n,i] = 100.*(deltaPx_c[n,i]/deltaPx_m[n,i]-1.)
      else:
         print ('Bad value (=0.) of deltaPx_m[%d,%d] = ' % (n,i))

      if (deltaPy_m[n,i] != 0.):
         deltaPy_c_m[n,i] = 100.*(deltaPy_c[n,i]/deltaPy_m[n,i]-1.)
      else:
         print ('Bad value (=0.) of deltaPy_m[%d,%d] = ' % (n,i))

      if (deltaPz_m[n,i] != 0.):
         deltaPz_c_m[n,i] = 100.*(deltaPz_c[n,i]/deltaPz_m[n,i]-1.)
      else:
         print ('Bad value (=0.) of deltaPz_m[%d,%d] = ' % (n,i))

      if (deltaEnrgIon_m[n,i] != 0.):
         dEion_c_m[n,i] = 100.*(deltaEnrgIon_c[n,i]/deltaEnrgIon_m[n,i]-1.)
      else:
         print ('Bad value (=0.) of deltaEnrgIon_m[%d,%d] = ' % (n,i))
#
# Integration using Simpson method:
#
      if (n > 0):
         frctnForce_cSM[i] +=  pi*n_e*100.*(deltaEnrgIon_c[n,i]+deltaEnrgIon_c[n-1,i])* \
                               .5*(rhoInit[n,i]+rhoInit[n-1,i])* \
                               (rhoInit[n,i]-rhoInit[n-1,i])                 # eV/m 
         frctnForce_mSM[i] +=  pi*n_e*100.*(deltaEnrgIon_m[n,i]+deltaEnrgIon_m[n-1,i])* \
                               .5*(rhoInit[n,i]+rhoInit[n-1,i])* \
                               (rhoInit[n,i]-rhoInit[n-1,i])                 # eV/m 
   timeEnd = os.times()
   timeRun[i] = float(timeEnd[0])-float(timeStart[0])  # CPU time , sec
   totalTimeRun += timeRun[i]
   print ('timeRun(%2d) = %6.3f seconds' % (i,timeRun[i]))

print ('Total time (icluding Simpson integration) = %6.3f seconds' % totalTimeRun)

print ('deltaEnrgIon_c: nPos=%d, nNeg=%d;  deltaEnrgIon_m: nPos=%d, nNeg=%d' % \
       (posSignDeltaEnrgIon_c,negSignDeltaEnrgIon_c, \
        posSignDeltaEnrgIon_m,negSignDeltaEnrgIon_m))

#
# Output for checking:
#
# print \
# ('n     Px_c          Px_m         Py_c          Py_m         Pz_c        Pz_m         Pz_c_m')  
# for i in range(10,11,1):     
#    for n in range(nImpctPrmtr):
#       print ('%d: %e %e %e %e %e %e %e' % \
#              (n,deltaPx_c[n,i],deltaPx_m[n,i],deltaPy_c[n,i], \
# 	        deltaPy_m[n,i],deltaPz_c[n,i],deltaPz_m[n,i],deltaPz_c_m[n,i]))
# print ('n     dEion_c      dEion_m')  
# for i in range(10,11,1):     
#    for n in range(nImpctPrmtr):
#       print ('%d: %e %e ' % (n,deltaEnrgIon_c[n,i],deltaEnrgIon_m[n,i]))

# print ('indxTestMax = %d' % indxTestMax)

#
# Plotting of the tests:
#
nn=np.arange(0,indxTestMax-1,1)

#
# C1:
#
if (plotFigureFlag == 0):   
   fig2020=plt.figure (2020)
   plt.plot(nn,C1test[0:indxTestMax-1],'.r')
   plt.xlabel('Points of Tracks',color='m',fontsize=16)
   plt.ylabel('$C1$, $cm^2$',color='m',fontsize=16)
   plt.title('$C1=[x_{gc}^2+y_{gc}^2+z_e^2+2J/(m_e \cdot \Omega_e)]^{0.5}$', \
             color='m',fontsize=16)
   plt.xlim([-5000,indxTestMax+5000])
   plt.grid(True)
if (saveFilesFlag == 1):
   fig2020.savefig('picturesCMA_v7/magnusExpansion_C1_fig2020cma.png')    
   print ('File "picturesCMA_v7/magnusExpansion_C1_fig2020cma.png" is written')   

#
# C2:
#
if (plotFigureFlag == 0):   
   fig2030=plt.figure (2030)
   plt.plot(nn,1.e-5*C2test[0:indxTestMax-1],'.r')
   plt.xlabel('Points of Tracks',color='m',fontsize=16)
   plt.ylabel('$C2$, $\cdot 10^5$ $cm^2/s$',color='m',fontsize=16)
   plt.title('$C2=2\cdot[V_{ix}\cdot(x_i-x_{gc})+V_{iy}\cdot(y_i-y_{gc})+(V_{iz}-V_{ez})\cdot(z_i-z_e)]$', \
             color='m',fontsize=14)
   plt.xlim([-5000,indxTestMax+5000])
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig2030.savefig('picturesCMA_v7/magnusExpansion_C2_fig2030cma.png')    
      print ('File "picturesCMA_v7/magnusExpansion_C2_fig2030cma.png" is written')   

#
# C3:
#
if (plotFigureFlag == 0):   
   fig2040=plt.figure (2040)
   plt.plot(nn,1e-11*C3test[0:indxTestMax-1],'.r')
   plt.xlabel('Points of Tracks',color='m',fontsize=16)
   plt.ylabel('$C3$, $\cdot 10^{11}$ $cm^2/s^2$',color='m',fontsize=16)
   plt.title('$C3=V_{ix}^2+V_{iy}^2+(V_{iz}-V_{ez})^2$',color='m',fontsize=16)
   plt.xlim([-5000,indxTestMax+5000])
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig2040.savefig('picturesCMA_v7/magnusExpansion_C3_fig2040cma.png')    
      print ('File "picturesCMA_v7/magnusExpansion_C3_fig2040cma.png" is written')   

#
# D1:
#
if (plotFigureFlag == 0):   
   fig2025=plt.figure (2025)
   plt.plot(nn,1.e-5*D1test[0:indxTestMax-1],'.r')
   plt.xlabel('Points of Tracks',color='m',fontsize=16)
   plt.ylabel('$10^{-5}\cdot D1$, $cm/s$',color='m',fontsize=16)
   plt.title('$D1=(2C_3\cdot \Delta t+C_2)/b_{ME}$ $-$ $C_2/C_1^{0.5}$',color='m',fontsize=16)
   plt.xlim([-5000,indxTestMax+5000])
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig2025.savefig('picturesCMA_v7/magnusExpansion_D1_fig2025cma.png')    
      print ('File "picturesCMA_v7/magnusExpansion_D1_fig2025cma.png" is written')   

#
# D2:
#
if (plotFigureFlag == 0):   
   fig2035=plt.figure (2035)
   plt.plot(nn,1.e4*D2test[0:indxTestMax-1],'.r')
   plt.xlabel('Points of Tracks',color='m',fontsize=16)
   plt.ylabel('$10^4\cdot D2$, $cm$',color='m',fontsize=16)
   plt.title('$D2=(2C_1+C_2\cdot \Delta t)/b_{ME}$ $-$ $2C_1^{0.5}$',color='m',fontsize=16)
   plt.xlim([-5000,indxTestMax+5000])
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig2035.savefig('picturesCMA_v7/magnusExpansion_D2_fig2035cma.png')    
      print ('File "picturesCMA_v7/magnusExpansion_D2_fig2035cma.png" is written')   

#
# Distance b_ME between particles for "ME" approach:
#
if (plotFigureFlag == 0):   
   fig2050=plt.figure (2050)
   plt.plot(nn,b_ME[0:indxTestMax-1],'.r')
   plt.xlabel('Points of Tracks',color='m',fontsize=16)
   plt.ylabel('$b_{ME}$, $cm$',color='m',fontsize=16)
   plt.title('Distance $b_{ME}$ between Particles for "ME" Approach', color='m',fontsize=16)
   plt.text(3500,.4,'$b_{ME}=[C1+C2\cdot \Delta t +C3 \cdot \Delta t^2]^{0.5}$', \
            color='m',fontsize=16)
   plt.text(33000,.36,('$(\Delta t=%8.2e$ $s)$' % timeStep_c),color='m',fontsize=16)
   plt.xlim([-5000,indxTestMax+5000])
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig2050.savefig('picturesCMA_v7/particleDistance_me_fig2050cma.png')    
      print ('File "picturesCMA_v7/particleDistance_me_fig2050cma.png" is written')   

#
# Distance b_gc between particles for "GC" approach:
#
if (plotFigureFlag == 0):   
   fig2055=plt.figure (2055)
   plt.plot(nn,b_gc[0:indxTestMax-1],'.r')
   plt.xlabel('Points of Tracks',color='m',fontsize=16)
   plt.ylabel('$b_{GC}$, $cm$',color='m',fontsize=16)
   plt.title('Distance $b_{GC}$ between Particles for "GC" Approach', color='m',fontsize=16)
   plt.text(0,.4,'$b_{GC}=[(x_i-x_{gc})^2+(y_i-y_{gc})^2+$',color='m',fontsize=16)
   plt.text(55500,.36,'$+(z_i-z_e)^2+2J/(m_e \cdot \Omega_e)]^{0.5}$', \
            color='m',fontsize=16)
   plt.xlim([-5000,indxTestMax+5000])
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig2055.savefig('picturesCMA/particleDistance_gc_fig2055cma.png')    
      print ('File "picturesCMA/particleDistance_gc_fig2055cma.png" is written')   

#
# Comparison of bCrrnt_c from "Guiding Center" with bTest from
# "Magnus expansion" approaches:
#
bCrrnt_cTest = np.zeros(indxTestMax)
bCrrnt_cTestRel = np.zeros(indxTestMax)
b_gc_ME_rel = np.zeros(indxTestMax)

for k in range(indxTestMax):
   bCrrnt_cTest[k] = .5*(bCrrnt_c[2*k]+bCrrnt_c[2*k+1])
#   bCrrnt_cTestRel[k] = bCrrnt_cTest[k]/b_ME[k]
   b_gc_ME_rel[k] = b_gc[k]/b_ME[k]
   actn_gc_ME_rel[k] = 1.e7*(action_gc[k]/action_ME[k]-1.)
   
if (plotFigureFlag == 0):   
   fig2060=plt.figure (2060)
#   plt.semilogy(nn,bCrrnt_cTest[0:indxTestMax-1],'.r')
   plt.plot(nn,bCrrnt_cTest[0:indxTestMax-1],'.r')
   plt.xlabel('Points of Tracks',color='m',fontsize=16)
   plt.ylabel('Test $b_{crrntTest}$, $cm$',color='m',fontsize=16)
   plt.title('Test $b_{crrntTest} = .5 \cdot [b_{crrnt}(k)+b_{crrnt}(k+1)]$',color='m', \
             fontsize=16)
   plt.xlim([-5000,indxTestMax+5000])
   # plt.ylim([.9*min(bCrrnt_cTest),1.1*max(bCrrnt_cTest)])
   plt.grid(True)

#
# Ratio b_gc/b_ME (absolute value):
#
if (plotFigureFlag == 0):   
   fig2070=plt.figure (2070)
#   plt.semilogy(nn,b_gc_ME_rel[0:indxTestMax-1],'.r')
   plt.plot(nn,b_gc_ME_rel[0:indxTestMax-1],'.r')
   plt.xlabel('Points of Tracks',color='m',fontsize=16)
   plt.ylabel('$b_{GC}/b_{ME}$',color='m',fontsize=16)
   plt.title('Comparison of Distances $b_{GC}$ and $b_{ME}$ between Particles',color='m',fontsize=16)
   plt.xlim([-5000,indxTestMax+5000])
   # plt.ylim([.9*min(b_gc_ME_rel),1.1*max(b_gc_ME_rel)])
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig2070.savefig('picturesCMA_v7/particleDistanceComprsn_gc_me_fig2070cma.png')    
      print ('File "picturesCMA_v7/particleDistanceComprsn_gc_me_fig2070cma.png" is written')   

#
# Ratio b_gc/b_ME (relative value):
#
if (plotFigureFlag == 0):   
   fig2080=plt.figure (2080)
   # plt.semilogy(nn,actn_gc_ME_rel[0:indxTestMax-1],'.r')
   plt.plot(nn,actn_gc_ME_rel[0:indxTestMax-1],'.r')
   plt.xlabel('Points of Tracks',color='m',fontsize=16)
   plt.ylabel('$10^7\cdot (J_{GC}/J_{ME}$ $-$ $1)$',color='m',fontsize=16)
   plt.title('Comparison of Actions $J_{GC}$ and $J_{ME}$',color='m',fontsize=16)
   plt.xlim([-5000,indxTestMax+5000])
   plt.ylim([.99*min(actn_gc_ME_rel),1.01*max(actn_gc_ME_rel)])
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig2080.savefig('picturesCMA_v7/actionComprsn_gc_me_fig2080cma.png')    
      print ('File "picturesCMA_v7/actionComprsn_gc_me_fig2080cma.png" is written')   

#
# Total length of interaction (1/2 of value):
#
nn=np.arange(0,nVion*nImpctPrmtr,1)
halfLintrTest = np.zeros(nVion*nImpctPrmtr)
for i in range(nVion):
   for n in range(nImpctPrmtr):
     halfLintrTest[nVion*i+n] = halfLintr[i,n] 

if (plotFigureFlag == 0):   
   fig2090=plt.figure (2090)
   plt.semilogy(nn,halfLintrTest,'.r')
   plt.xlabel('Points of Tracks',color='m',fontsize=16)
   plt.ylabel('$0.5 \cdot L_{Intrctn}$, $cm$',color='m',fontsize=16)
   plt.title('Total Length of Interaction: $L_{Intrctn}=2 \cdot [R_{max}^2-rho_{Init}^2)]^{0.5}$', \
             color='m',fontsize=16)
   plt.xlim([-100,nVion*nImpctPrmtr+100])
   plt.ylim([.9*min(halfLintrTest),1.1*max(halfLintrTest)])
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig2090.savefig('picturesCMA/totalLengthIntrsctn_fig2090cma.png')    
      print ('File "picturesCMA/totalLengthIntrsctn_fig2090cma.png" is written')   

#===================================================
#
# There is fitting for correct values of deltaEnrgIon_m
#
#===================================================
#
# Fitting for figures with deltaEnrgIon_m (my own Least Squares Method - LSM;
# Python has own routine for LSM - see site  
#     http://scipy-cookbook.readthedocs.io/items/FittingData.html):
#
#
# Fittied function: 
#   
#  |deltaEnrgIon| = 10^fitA * rho^fitB,
#  so that
#
# log10(|deltaEnrgIon|) = fitB*log10(rho) + fitA
#
# So, the dimension of expression (10^fitA * rho^fitB) is the same
# as deltaEnrgIon,  i.e. eV
#

timeStart = os.times()

fitA_dEion = np.zeros(nVion)            # dimensionless 
fitB_dEion = np.zeros(nVion)            # dimensionless
rhoInitFit_dEion = np.zeros((nImpctPrmtr,nVion))
deltaEnrgIon_m_fit = np.zeros((nImpctPrmtr,nVion))

funcHi2_dEion = np.zeros(nVion)
fitA_dEion,fitB_dEion,funcHi2_dEion,rhoInitFit_dEion, deltaEnrgIon_m_fit = \
fitting(nImpctPrmtr,nVion,rhoInit,deltaEnrgIon_m)

dPosA_dEion = np.zeros(nVion)
dNegA_dEion = np.zeros(nVion)

dPosA_dEion,dNegA_dEion = \
errFitAB(nImpctPrmtr,nVion,rhoInit,deltaEnrgIon_m_fit,fitA_dEion,fitB_dEion,funcHi2_dEion,1,2)

dPosB_dEion = np.zeros(nVion)
dNegB_dEion = np.zeros(nVion)

dPosB_dEion,dNegB_dEion = \
errFitAB(nImpctPrmtr,nVion,rhoInit,deltaEnrgIon_m_fit,fitA_dEion,fitB_dEion,funcHi2_dEion,2,2)
# print ('Fitting for deltaEion:')
# for i in range(nVion):
#    print ('i=%2d: fitA_dEion = %e (+%e,-%e), fitB_dEion = %e (+%e,-%e), hi2_1 = %e'  % \
#           (i,fitA_dEion[i],dPosA_dEion[i],dNegA_dEion[i], \
# 	     fitB_dEion[i],dPosB_dEion[i],dNegB_dEion[i],funcHi2_dEion[i]))

#
# Analytical Integration of the fitted dependence 10**A*rho**B.
#
# For this dependece on rho:
#
# Friction force = 10**A*n_e*integral_rhoMin^rhoMax (rho**B*rho)*dRho = 
# = 10**A*n_e/(B+2)*[rhoMax**(B+2)-rhoMax**(B+2)] (dimension=eV/cm):
# 
frctnForce_AI = np.zeros(nVion)

for i in range(nVion):
   factorA1 = math.pow(10.,fitA_dEion[i])
   factorB1 = 2.+fitB_dEion[i]
   frctnForce_AI[i] = 2.*pi*n_e*100.*factorA1/factorB1* \
                      (math.pow(impctPrmtrMax[i],factorB1)- \
                       math.pow(impctPrmtrMin,factorB1))             # eV/m

timeEnd = os.times()
timeFitting = float(timeEnd[0])-float(timeStart[0])  # CPU time , sec
print ('Time of integration = %6.3f seconds' % timeFitting)

#
# Dependences of transferred energy to ion on ion velocity for 
# different initial impact parameters:
#
rhoSlctd = [.004,.02,.06,.1]
nRhoSlctd = len(rhoSlctd)

deltaEnrgIon_dpnd_Vi = np.zeros((nRhoSlctd,nVion))
npStart = np.zeros((nRhoSlctd,), dtype=int)

for k in range(nRhoSlctd):
   slctdFlag = 0
   for i in range(nVion):
      if (slctdFlag == 0):
         for n in range(nImpctPrmtr):
            if (rhoInit[n,i] >= rhoSlctd[k]):
               npStart[k] = i
               slctdFlag = 1
               break

for k in range(nRhoSlctd):
   for i in range(npStart[k],nVion,1):
      factorA = math.pow(10.,fitA_dEion[i])
      deltaEnrgIon_dpnd_Vi[k,i] = factorA*math.pow(rhoSlctd[k],fitB_dEion[i])
#      print ('deltaEnrgIon_dpnd_Vi[%d,%d] = %e' %(k,i,deltaEnrgIon_dpnd_Vi[k,i]))

#===================================================
#
# There is fitting of deltaPz_m (these values > 0 always) !!!
#
#===================================================
#
# Fitting for figures with deltaPz_m (my own Least Squares Method - LSM;
# Python has own routine for LSM - see site  
#     http://scipy-cookbook.readthedocs.io/items/FittingData.html):
#
#
# Fittied function: 
#   
#  deltaPz_m = 10^fitA_pz * rho^fitB_pz,
#  so that
#
# log10(deltaPz_m) = fitB_pz*log10(rho) + fitA_pz
#
# So, the dimension of expression (10^fitA_pz * rho^fitB_pz) is the same
# as deltaPz_m,  i.e. eV
#

fitA_pz = np.zeros(nVion)            # dimensionless 
fitB_pz = np.zeros(nVion)            # dimensionless
rhoInitFit_pz = np.zeros((nImpctPrmtr,nVion))
deltaPz_m_fit = np.zeros((nImpctPrmtr,nVion))
fitA_pz,fitB_pz,funcHi2_pz,rhoInitFit_pz, deltaPz_m_fit = \
fitting(nImpctPrmtr,nVion,rhoInit,deltaPz_m)

dPosA_pz = np.zeros(nVion)
dNegA_pz = np.zeros(nVion)
dPosA_pz,dNegA_pz = \
errFitAB(nImpctPrmtr,nVion,rhoInit,deltaPz_m_fit,fitA_pz,fitB_pz,funcHi2_pz,1,2)
dPosB_pz = np.zeros(nVion)
dNegB_pz = np.zeros(nVion)
dPosB_pz,dNegB_pz = \
errFitAB(nImpctPrmtr,nVion,rhoInit,deltaPz_m_fit,fitA_pz,fitB_pz,funcHi2_pz,2,2)
# print ('Fitting fordeltaPz_m:')
# for i in range(nVion):
#    print ('i=%2d: fitA_pz = %e (+%e,-%e), fitB_pz = %e (+%e,-%e), hi2_1 = %e'  % \
#           (i,fitA_pz[i],dPosA_pz[i],dNegA_pz[i], \
# 	     fitB_pz[i],dPosB_pz[i],dNegB_pz[i],funcHi2_pz[i]))
# print ('<fitA_pz> = %e +- %e' % (mean(fitA_pz),mean(dNegA_pz)))
# print ('<fitB_pz> = %e +- %e' % (mean(fitB_pz),mean(dNegB_pz)))

#===================================================
#
# There is fitting of deltaPx_m (these values > 0 always) !!!
#
#===================================================
#

rhoInitFit_px = np.zeros((nImpctPrmtr,nVion))
deltaPx_m_fit = np.zeros((nImpctPrmtr,nVion))
funcHi2__px = np.zeros(nVion)
fitA_px = np.zeros(nVion)            # dimensionless 
fitB_px = np.zeros(nVion)            # dimensionless

fitA_px,fitB_px,funcHi2_px,rhoInitFit_px, deltaPx_m_fit = \
fitting(nImpctPrmtr,nVion,rhoInit,deltaPx_m)

dPosA_px = np.zeros(nVion)
dNegA_px = np.zeros(nVion)
dPosA_px,dNegA_px = \
errFitAB(nImpctPrmtr,nVion,rhoInit,deltaPx_m_fit,fitA_px,fitB_px,funcHi2_px,1,2)
dPosB_px = np.zeros(nVion)
dNegB_px = np.zeros(nVion)
dPosB_px,dNegB_px = \
errFitAB(nImpctPrmtr,nVion,rhoInit,deltaPx_m_fit,fitA_px,fitB_px,funcHi2_px,2,2)
# print ('Fitting for deltaPx_m:')
# for i in range(nVion):
#    print ('i=%2d: fitA_px = %e (+%e,-%e), fitB_px = %e (+%e,-%e), hi2_1 = %e'  % \
#           (i,fitA_px[i],dPosA_px[i],dNegA_px[i], \
# 	     fitB_px[i],dPosB_px[i],dNegB_px[i],funcHi2_px[i]))

xLimit = [1.015*np.log10(VionRel[0]),.95*np.log10(VionRel[nVion-1])]

yLimMin = 0.
yLimMax = 10.*min(fitA_pz)
if (min(fitA_pz) > 0):
   yLimMin = 10.*max(fitA_pz)
   yLimMax = 0.
for i in range(nVion):
   if (fitA_pz[i] - dNegA_pz[i]) < yLimMin:
      yLimMin = fitA_pz[i] - dNegA_pz[i]
   if (fitA_pz[i] + dPosA_pz[i]) > yLimMax:
      yLimMax = fitA_pz[i] + dPosA_pz[i]
# print ('Exponent A (pz): yLimMin = %e, yLimMax = %e' % (yLimMin,yLimMax))

yLimit = [yLimMin-.25,yLimMax+.25]

if (plotFigureFlag == 0):   
   fig3000=plt.figure (3000)
   plt.errorbar(np.log10(VionRel),fitA_pz,yerr=[dNegA_pz,dPosA_pz],fmt='-ro', \
                ecolor='b',capsize=5,capthick=1)
   plt.xlabel('Relative Ion Velocity, $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
   plt.ylabel('Exponent $A$', color='m',fontsize=14)
   titleHeader = 'Dependence of Transferred Momenta to Single Ion: '
   titleHeader += '$\Delta P_z$ = $10^A\cdot rho^B$'
   plt.title(titleHeader,color='m',fontsize=12)
   plt.text(-3.75,-26.0,('$V_{e0}=%5.3f\cdot10^{%2d}$cm/s' % (mantV0,powV0)), \
            color='m',fontsize=16)
   plt.text(-4.0,-28.,('<A>=%7.3f $\pm$ %5.3f' % (mean(fitA_pz),mean(dNegA_pz))), \
            color='r',fontsize=16)
#   plt.text(-3.25,-29.65,('$-$%5.3f' % (mean(dNegA_pz))),color='r',fontsize=12)
#   plt.text(-3.25,-29.15,('$+$%5.3f' % (mean(dPosA_pz))),color='r',fontsize=12)
   plt.xlim(xLimit)
   plt.ylim(yLimit)
   plt.grid(True)
   plt.plot([np.log10(relVeTrnsv),np.log10(relVeTrnsv)],yLimit,'--m',linewidth=1)
   plt.text(-2.55,-28.25,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([np.log10(relVeLong),np.log10(relVeLong)],yLimit,'--m',linewidth=1)
   plt.text(-4.24,-28.25,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   if (saveFilesFlag == 1):
      fig3000.savefig('picturesCMA_v7/fitA_dPz_fig3000cma.png')    
      print ('File "picturesCMA_v7/fitA_dPz_fig3000cma.png" is written')   

yLimMin = 0.
yLimMax = 10.*min(fitB_pz)
if (min(fitB_pz) > 0):
   yLimMin = 10.*max(fitB_pz)
   yLimMax = 0.
for i in range(nVion):
   if (fitB_pz[i] - dNegB_pz[i]) < yLimMin:
      yLimMin = fitB_pz[i] - dNegB_pz[i]
   if (fitB_pz[i] + dPosB_pz[i]) > yLimMax:
      yLimMax = fitB_pz[i] + dPosB_pz[i]
# print ('Exponent B (pz): yLimMin = %e, yLimMax = %e' % (yLimMin,yLimMax))

yLimit = [yLimMin-.1,yLimMax+.1]
if (plotFigureFlag == 0):   
   fig3010=plt.figure (3010)
   plt.errorbar(np.log10(VionRel),fitB_pz,yerr=[dNegB_pz,dPosB_pz],fmt='-ro', \
                ecolor='b',capsize=5,capthick=1)
   plt.xlabel('Relative Ion Velocity, $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
   plt.ylabel('Exponent $B$', color='m',fontsize=14)
   titleHeader = 'Dependence of Transferred Momenta to Single Ion: '
   titleHeader += '$\Delta P_z$ = $10^A\cdot rho^B$'
   plt.title(titleHeader,color='m',fontsize=12)
   plt.text(-3.75,-.87,('$V_{e0}=%5.3f\cdot10^{%2d}$cm/s' % (mantV0,powV0)), \
            color='m',fontsize=16)
   plt.text(-3.9,-1.55,('<B>=%6.3f $\pm$ %5.3f' % (mean(fitB_pz),mean(dNegB_pz))), \
            color='r',fontsize=16)
#   plt.text(-2.85,-2.25,('$-$%5.3f' % (mean(dNegB_pz))),color='r',fontsize=12)
#   plt.text(-2.85,-1.75,('$+$%5.3f' % (mean(dPosB_pz))),color='r',fontsize=12)
   plt.xlim(xLimit)
   plt.ylim(yLimit)
   plt.grid(True)
   plt.plot([np.log10(relVeTrnsv),np.log10(relVeTrnsv)],yLimit,'--m',linewidth=1)
   plt.text(-2.55,-1.74,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([np.log10(relVeLong),np.log10(relVeLong)],yLimit,'--m',linewidth=1)
   plt.text(-4.24,-1.74,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   if (saveFilesFlag == 1):
      fig3010.savefig('picturesCMA_v7/fitB_dPz_fig3010cma.png')    
      print ('File "picturesCMA_v7/fitB_dPz_fig3010cma.png" is written')   

yLimMin = 0.
yLimMax = 10.*min(fitA_px)
if (min(fitA_px) > 0):
   yLimMin = 10.*max(fitA_px)
   yLimMax = 0.
for i in range(nVion):
   if (fitA_px[i] - dNegA_px[i]) < yLimMin:
      yLimMin = fitA_px[i] - dNegA_px[i]
   if (fitA_px[i] + dPosA_px[i]) > yLimMax:
      yLimMax = fitA_px[i] + dPosA_px[i]
# print ('Exponent A (px): yLimMin = %e, yLimMax = %e' % (yLimMin,yLimMax))

yLimit = [yLimMin-.15,yLimMax+.15]

if (plotFigureFlag == 0):   
   fig3020=plt.figure (3020)
   plt.errorbar(np.log10(VionRel),fitA_px,yerr=[dNegA_px,dPosA_px],fmt='-ro', \
                ecolor='b',capsize=5,capthick=1)
   plt.xlabel('Relative Ion Velocity, $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
   plt.ylabel('Exponent $A$', color='m',fontsize=14)
   titleHeader = 'Dependence of Transferred Momenta to Single Ion: '
   titleHeader += '$\Delta P_x$ = $10^A\cdot rho^B$'
   plt.title(titleHeader,color='m',fontsize=12)
   plt.text(-3.75,-24.2,('$V_{e0}=%5.3f\cdot10^{%2d}$cm/s' % (mantV0,powV0)), \
            color='m',fontsize=16)
   plt.text(-3.9,-24.8,('<A>=%6.3f $\pm$ %5.3f' % (mean(fitA_px),mean(dNegA_px))), \
            color='r',fontsize=16)
   plt.xlim(xLimit)
   plt.ylim(yLimit)
   plt.grid(True)
   plt.plot([np.log10(relVeTrnsv),np.log10(relVeTrnsv)],yLimit,'--m',linewidth=1)
   plt.text(-2.55,-25.05,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([np.log10(relVeLong),np.log10(relVeLong)],yLimit,'--m',linewidth=1)
   plt.text(-4.24,-25.05,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   if (saveFilesFlag == 1):
      fig3020.savefig('picturesCMA_v7/fitA_dPx_fig3020cma.png')    
      print ('File "picturesCMA_v7/fitA_dPx_fig3020cma.png" is written')   

yLimMin = 0.
yLimMax = 10.*min(fitB_px)
if (min(fitB_px) > 0):
   yLimMin = 10.*max(fitB_px)
   yLimMax = 0.
for i in range(nVion):
   if (fitB_px[i] - dNegB_px[i]) < yLimMin:
      yLimMin = fitB_px[i] - dNegB_px[i]
   if (fitB_px[i] + dPosB_px[i]) > yLimMax:
      yLimMax = fitB_px[i] + dPosB_px[i]
# print ('Exponent B (px): yLimMin = %e, yLimMax = %e' % (yLimMin,yLimMax))

yLimit = [yLimMin-.05,yLimMax+.05]

if (plotFigureFlag == 0):   
   fig3030=plt.figure (3030)
   plt.errorbar(np.log10(VionRel),fitB_px,yerr=[dNegB_px,dPosB_px],fmt='-ro', \
                ecolor='b',capsize=5,capthick=1)
   plt.xlabel('Relative Ion Velocity, $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
   plt.ylabel('Exponent $B$', color='m',fontsize=14)
   titleHeader = 'Dependence of Transferred Momenta to Single Ion: '
   titleHeader += '$\Delta P_x$ = $10^A\cdot rho^B$'
   plt.title(titleHeader,color='m',fontsize=12)
   plt.text(-3.75,-.95,('$V_{e0}=%5.3f\cdot10^{%2d}$cm/s' % (mantV0,powV0)), \
            color='m',fontsize=16)
   plt.text(-3.9,-1.15,('<B>=%6.3f $\pm$ %5.3f' % (mean(fitB_px),mean(dNegB_px))), \
            color='r',fontsize=16)
   plt.xlim(xLimit)
   plt.ylim(yLimit)
   plt.grid(True)
   plt.plot([np.log10(relVeTrnsv),np.log10(relVeTrnsv)],yLimit,'--m',linewidth=1)
   plt.text(-2.55,-1.22,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([np.log10(relVeLong),np.log10(relVeLong)],yLimit,'--m',linewidth=1)
   plt.text(-4.24,-1.22,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   if (saveFilesFlag == 1):
      fig3030.savefig('picturesCMA_v7/fitB_dPx_fig3030cma.png')    
      print ('File "picturesCMA/_v7/fitB_dPx_fig3030cma.png" is written')   

# plt.show()

# sys.exit()

#
#=======================================================
#      
# Main plotting:      
#      
if (plotFigureFlag == 0):   
   fig110=plt.figure (110)
   plt.plot(arrayA,arrayB,'.r')
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   plt.title('Map of Parameters A,B', color='m',fontsize=16)
   # plt.xlim([minA,maxA])
   # plt.ylim([minB,maxB])
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig110.savefig('picturesCMA/mapA-B_fig110cma.png')    
      print ('File "picturesCMA/mapA-B_fig110cma.png" is written')   

if (plotFigureFlag == 0):   
   fig20=plt.figure (20)
   plt.plot(nnTotalPoints,bCrrnt_c[0:2*totalPoints-1],'.r')
   # plt.semilogy(nn,bCrrnt_c[0:2*totalPoints-1],'.r')
   plt.xlabel('Points of Tracks',color='m',fontsize=16)
   plt.ylabel('$b_{Lab.Sys}$, $cm$',color='m',fontsize=16)
   plt.title('Distance $b_{Lab.Sys}$ between Particles in Lab.System', color='m',fontsize=16)
   plt.xlim([-5000,2*totalPoints+5000])
   # plt.xlim([0,2000])
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig20.savefig('picturesCMA/particleDistance_ls_fig20cma.png')    
      print ('File "picturesCMA/particleDistance_ls_fig20cma.png" is written')   

if (plotFigureFlag == 0):   
   fig30=plt.figure (30)
   plt.plot(nnTotalPoints,arrayA[0:2*totalPoints-1],'.r', \
            nnTotalPoints,arrayB[0:2*totalPoints-1],'.b')
   plt.xlabel('Points of Tracks',color='m',fontsize=16)
   plt.ylabel('$A$, $B$',color='m',fontsize=16)
   plt.title('$A=log_{10}(q_e^2/b/E_{kin})$, $B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   plt.xlim([-5000,2*totalPoints+5000])
   # plt.ylim([minB,maxB])
   plt.grid(True)
   plt.legend(['A','B'],loc='lower left',fontsize=14)
   if (saveFilesFlag == 1):
      fig30.savefig('picturesCMA/parametersA-B_fig30cma.png')    
      print ('File "picturesCMA/parametersA-B_fig30cma.png" is written')   

xVionRel = np.zeros((nImpctPrmtr,nVion))
for i in range(nVion):
   for n in range(nImpctPrmtr):
       xVionRel[n,i] = VionRel[i]

if (plotFigureFlag == 0):   
   fig40=plt.figure (40)
   for i in range(nVion):
       plt.semilogx(xVionRel[0:nImpctPrmtr,i],rhoInit[0:nImpctPrmtr,i],'.r')
   plt.xlabel('Relative Ion Velocity, $V_i/V_{e0}$',color='m',fontsize=16)
   plt.ylabel('$rho_{Init}$, cm',color='m',fontsize=16)
   plt.title('Subdivisions for $rho_{Init}$ for Integration: Simpson Method', \
             color='m',fontsize=16)
   plt.grid(True)
   yLimit=[0.,.405]
   plt.ylim(yLimit)
   plt.plot([relVeTrnsv,relVeTrnsv],yLimit,'--m',linewidth=1)
   plt.text(1.6e-3,-.026,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([relVeLong,relVeLong],yLimit,'--m',linewidth=1)
   plt.text(3.9e-5,.05,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   if (saveFilesFlag == 1):
      fig40.savefig('picturesCMA/initialImpactParameter_SM_fig40cma.png')    
      print ('File "picturesCMA/initialImpactParameter_SM_fig40cma.png" is written')   

if (plotFigureFlag == 0):   
   fig45=plt.figure (45)
   for i in range(nVion):
       plt.loglog(xVionRel[0:nImpctPrmtr,i],rhoInit[0:nImpctPrmtr,i],'.r')
   plt.xlabel('Relative Ion Velocity, $V_i/V_{e0}$',color='m',fontsize=16)
   plt.ylabel('$rho_{Init}$, cm',color='m',fontsize=16)
   plt.title('Subdivisions for $rho_{Init}$ for Integration: Simpson Method', \
             color='m',fontsize=16)
   plt.grid(True)
   yLimit=[1.3e-3,.45]
   plt.ylim(yLimit)
   plt.plot([relVeTrnsv,relVeTrnsv],yLimit,'--m',linewidth=1)
   plt.text(1.6e-3,.15,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([relVeLong,relVeLong],yLimit,'--m',linewidth=1)
   plt.text(3.9e-5,.15,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   if (saveFilesFlag == 1):
      fig45.savefig('picturesCMA/initialImpactParameter_SM_fig45cma.png')    
      print ('File "picturesCMA/initialImpactParameter_SM_fig45cma.png" is written')   

'''
#
# Figure compares calculated values of of deltaEnrgIon (their dependences 
# on impact parameter for different ion velocities) for two approaches
# (figure numbrFigures[0]+1 is the same and taking into account positive and 
#  negative values of the deltaEnrgIon_c for guiding center approach): 
#

VionCrrnt = V0*VionRel[0]
powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
plt.figure (50)
plt.loglog(rhoInit[0:nImpctPrmtr-1,0],deltaEnrgIon_c[0:nImpctPrmtr-1,0],'-xr', \
           rhoInit[0:nImpctPrmtr-1,0],deltaEnrgIon_m[0:nImpctPrmtr-1,0],'--or', \
           linewidth=1)
plt.xlabel('Track Initial Impact Parameter $rho_{Init}$', color='m',fontsize=16)
plt.ylabel('$\Delta E_{ion}$, $erg$', color='m',fontsize=16)
titleHeader = 'Transferred Energy $\Delta E_{ion}$ to Ion for '
titleHeader += ' $V_{ion}=%5.3f\cdot10^{%2d}$ $cm/s$'
plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=16)
plt.xlim([.95*rhoInit[0,0],1.05*rhoInit[nImpctPrmtr-1,0]])
plt.grid(True)

xRhoInitPx_c = np.zeros(nImpctPrmtr*nVion)
xRhoInitPx_m = np.zeros(nImpctPrmtr*nVion)
yDeltaPx_c = np.zeros(nImpctPrmtr*nVion)
yDeltaPx_m = np.zeros(nImpctPrmtr*nVion)
indx_c = 0
indx_m = 0
for n in range(nImpctPrmtr):
   if deltaPx_c[n,0] > 0.:
      xRhoInitPx_c[indx_c] = rhoInit[n,0]
      yDeltaPx_c[indx_c] = deltaPx_c[n,0]
#      print ('n_c=%2d: xRhoInitPx_c = %e, yDeltaPx_c = %e' % \
#            (indx_c,xRhoInitPx_c[indx_c],yDeltaPx_c[indx_c])) 
      indx_c += 1
   if deltaPx_m[n,0] > 0.:
      xRhoInitPx_m[indx_c] = rhoInit[n,0]
      yDeltaPx_m[indx_c] = deltaPx_m[n,0]
#      print ('n_m=%2d: xRhoInitPx_m = %e, yDeltaPx_m = %e' % \
#            (indx_m,xRhoInitPx_m[indx_m],yDeltaPx_m[indx_m])) 
      indx_m += 1
maxIndx_c = indx_c-1
maxIndx_m = indx_m-1
# print ('maxIndx_c = %d, maxIndx_m = %d' % (maxIndx_c,maxIndx_m))

#
# Figure compares calculated values of of deltaPx (their dependences 
# on impact parameter for different ion velocities) for two approaches
# (figure numbrFigures[0]+2 is the same and taking into account positive and 
#  negative values of the deltaPx_c for guiding center approach): 
#
plt.figure (51)
plt.loglog(xRhoInitPx_c[0:maxIndx_c],yDeltaPx_c[0:maxIndx_c],'-xr', \
           xRhoInitPx_m[0:maxIndx_m],yDeltaPx_m[0:maxIndx_m],'--or', \
           linewidth=1)
plt.xlabel('Track Initial Impact Parameter $rho_{Init}$, $cm$', \
           color='m',fontsize=16)
plt.ylabel('$\Delta P_{ix}$, $g\cdot cm/s$', color='m',fontsize=16)
titleHeader = 'Transferred Momenta $\Delta P_{ix}$ to Ion for '
titleHeader += ' $V_{ion}=%5.3f\cdot10^{%2d}$ $cm/s$'
plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=16)
plt.xlim([.95*min(xRhoInitPx_c[0],xRhoInitPx_m[0]), \
         1.05*max(xRhoInitPx_c[maxIndx_c],xRhoInitPx_m[maxIndx_m])])
plt.grid(True)

xRhoInitPz_c = np.zeros(nImpctPrmtr*nVion)
xRhoInitPz_m = np.zeros(nImpctPrmtr*nVion)
yDeltaPz_c = np.zeros(nImpctPrmtr*nVion)
yDeltaPz_m = np.zeros(nImpctPrmtr*nVion)
indx_c = 0
indx_m = 0
for n in range(nImpctPrmtr):
   if deltaPz_c[n,0] > 0.:
      xRhoInitPz_c[indx_c] = rhoInit[n,0]
      yDeltaPz_c[indx_c] = deltaPz_c[n,0]
#      print ('n_c=%2d: xRhoInitPz_c = %e, yDeltaPz_c = %e' % \
#            (indx_c,xRhoInitPz_c[indx_c],yDeltaPz_c[indx_c])) 
      indx_c += 1
   if deltaPz_m[n,0] > 0.:
      xRhoInitPz_m[indx_c] = rhoInit[n,0]
      yDeltaPz_m[indx_c] = deltaPz_m[n,0]
#      print ('n_m=%2d: xRhoInitPz_m = %e, yDeltaPz_m = %e' % \
#            (indx_m,xRhoInitPz_m[indx_m],yDeltaPz_m[indx_m])) 
      indx_m += 1
maxIndx_c = indx_c-1
maxIndx_m = indx_m-1
# print ('maxIndx_c = %d, maxIndx_m = %d' % (maxIndx_c,maxIndx_m))

#
# Figure compares calculated values of of deltaPz (their dependences 
# on impact parameter for different ion velocities) for two approaches
# (figure numbrFigures[0]+5): 
#
plt.figure (53)
plt.loglog(xRhoInitPz_c[0:maxIndx_c],yDeltaPz_c[0:maxIndx_c],'-xr', \
           xRhoInitPz_m[0:maxIndx_m],yDeltaPz_m[0:maxIndx_m],'--or', \
           linewidth=1)
plt.xlabel('Track Initial Impact Parameter $rho_{Init}$, $cm$', \
           color='m',fontsize=16)
plt.ylabel('$\Delta P_{iz}$, $g\cdot cm/s$', color='m',fontsize=16)
titleHeader = 'Transferred Momenta $\Delta P_{iz}$ to Ion for '
titleHeader += ' $V_{ion}=%5.3f\cdot10^{%2d}$ $cm/s$'
plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=16)
plt.xlim([.95*min(xRhoInitPz_c[0],xRhoInitPz_m[0]), \
         1.05*max(xRhoInitPz_c[maxIndx_c],xRhoInitPz_m[maxIndx_m])])
plt.grid(True)
'''

#
# Figures 60,70,80, and 90 compare calculated values of of deltaEnrgIon 
# (their dependences on impact parameter for first values of ion velocities)
# for two approaches (figure numbrFigures[*]+1 is the same and taking into 
# account positive and negative values of the deltaEnrgIon_c for guiding center approach): 
#
'''
VionCrrnt = V0*VionRel[1]
powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
plt.figure (60)
plt.loglog(rhoInit[0:nImpctPrmtr-1,1],deltaEnrgIon_c[0:nImpctPrmtr-1,1],'-xb', \
           rhoInit[0:nImpctPrmtr-1,1],deltaEnrgIon_m[0:nImpctPrmtr-1,1],'--ob', \
           linewidth=1)
plt.xlabel('Track Initial Impact Parameter $rho_{Init}$', \
           color='m',fontsize=16)
plt.ylabel('$\Delta E_{ion}$, $erg$', color='m',fontsize=16)
titleHeader = 'Transferred Energy $\Delta E_{ion}$ to Ion for '
titleHeader += ' $V_{ion}=%5.3f\cdot10^{%2d}$ $cm/s$'
plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=16)
plt.xlim([.95*rhoInit[0,1],1.05*rhoInit[nImpctPrmtr-1,1]])
plt.grid(True)

VionCrrnt = V0*VionRel[2]
powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
plt.figure (70)
plt.loglog(rhoInit[0:nImpctPrmtr-1,2],deltaEnrgIon_c[0:nImpctPrmtr-1,2],'-xg', \
           rhoInit[0:nImpctPrmtr-1,2],deltaEnrgIon_m[0:nImpctPrmtr-1,2],'--og', \
           linewidth=1)
plt.xlabel('Track Initial Impact Parameter $rho_{Init}$', \
           color='m',fontsize=16)
plt.ylabel('$\Delta E_{ion}$, $erg$', color='m',fontsize=16)
titleHeader = 'Transferred Energy $\Delta E_{ion}$ to Ion for '
titleHeader += ' $V_{ion}=%5.3f\cdot10^{%2d}$ $cm/s$'
plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=16)
plt.xlim([.95*rhoInit[0,2],1.05*rhoInit[nImpctPrmtr-1,2]])
plt.grid(True)

VionCrrnt = V0*VionRel[3]
powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
plt.figure (80)
plt.loglog(rhoInit[0:nImpctPrmtr-1,3],deltaEnrgIon_c[0:nImpctPrmtr-1,3],'-xk', \
           rhoInit[0:nImpctPrmtr-1,3],deltaEnrgIon_m[0:nImpctPrmtr-1,3],'--ok', \
           linewidth=1)
plt.xlabel('Track Initial Impact Parameter $rho_{Init}$', \
           color='m',fontsize=16)
plt.ylabel('$\Delta E_{ion}$, $erg$', color='m',fontsize=16)
titleHeader = 'Transferred Energy $\Delta E_{ion}$ to Ion for '
titleHeader += ' $V_{ion}=%5.3f\cdot10^{%2d}$ $cm/s$'
plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=16)
plt.xlim([.95*rhoInit[0,3],1.05*rhoInit[nImpctPrmtr-2,3]])
plt.grid(True)

VionCrrnt = V0*VionRel[4]
powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
plt.figure (90)
plt.loglog(rhoInit[0:nImpctPrmtr-1,4],deltaEnrgIon_c[0:nImpctPrmtr-1,4],'-xm', \
           rhoInit[0:nImpctPrmtr-1,4],deltaEnrgIon_m[0:nImpctPrmtr-1,4],'--om', \
           linewidth=1)
plt.xlabel('Track Initial Impact Parameter $rho_{Init}$', \
           color='m',fontsize=16)
plt.ylabel('$\Delta E_{ion}$, $erg$', color='m',fontsize=16)
titleHeader = 'Transferred Energy $\Delta E_{ion}$ to Ion for '
titleHeader += ' $V_{ion}=%5.3f\cdot10^{%2d}$ $cm/s$'
plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=16)
plt.xlim([.95*rhoInit[0,4],1.05*rhoInit[nImpctPrmtr-1,4]])
plt.grid(True)
'''

#
# Dependences of transferred energy to ion and different momenta on initial 
# impact parameter for different ion velocity (calculated and fitted values):
#
indxFigures = [0,9,12,18,19,23,27,29,31,34,39,49]
numbrFigures = [500,600,630,660,700,730,760,800,830,860,900,1000]
xPos = [.00218,.0022,.0024,.0027,.0026,.00265,.00265,.00265,.00265,.0028,.0029,.0035]
yPos = [6.4e-9,6.7e-9,6.4e-9,5.9e-9,6.2e-9,5.6e-9,5.8e-9,6.3e-9,5.8e-9,5.9e-9,5.8e-9,4.7e-9]

if (plotFigureFlag == 0):   
   for i in range(12):
      VionCrrnt = V0*VionRel[indxFigures[i]]
      powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
      mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 

#
# Pz:
#
posPz_c = np.zeros((12,nImpctPrmtr))
rhoPosPz_c = np.zeros((12,nImpctPrmtr)) 
negPz_c = np.zeros((12,nImpctPrmtr)) 
rhoNegPz_c = np.zeros((12,nImpctPrmtr)) 
posPz_m = np.zeros((12,nImpctPrmtr)) 
rhoPosPz_m = np.zeros((12,nImpctPrmtr)) 
negPz_m = np.zeros((12,nImpctPrmtr)) 
rhoNegPz_m = np.zeros((12,nImpctPrmtr)) 
nPosPz_c = array('i',[0]*12)
nNegPz_c = array('i',[0]*12)
nPosPz_m = array('i',[0]*12)
nNegPz_m = array('i',[0]*12)
for i in range(12):
   nPosPz_c[i] = -1
   nNegPz_c[i] = -1
   nPosPz_m[i] = -1
   nNegPz_m[i] = -1
   for k in range(nImpctPrmtr):
      if (deltaPz_c[k,indxFigures[i]] > 0):
         nPosPz_c[i] += 1
         rhoPosPz_c[i,nPosPz_c[i]] = rhoInit[k,indxFigures[i]]
         posPz_c[i,nPosPz_c[i]] = deltaPz_c[k,indxFigures[i]] 
      if (deltaPz_c[k,indxFigures[i]] <= 0):
         nNegPz_c[i] += 1 
         rhoNegPz_c[i,nNegPz_c[i]] = rhoInit[k,indxFigures[i]]
         negPz_c[i,nNegPz_c[i]] = abs(deltaPz_c[k,indxFigures[i]]) 
      if (deltaPz_m[k,indxFigures[i]] > 0):
         nPosPz_m[i] += 1 
         rhoPosPz_m[i,nPosPz_m[i]] = rhoInit[k,indxFigures[i]]
         posPz_m[i,nPosPz_m[i]] = deltaPz_m[k,indxFigures[i]] 
      if (deltaPz_m[k,indxFigures[i]] <= 0):
         nNegPz_m[i] += 1 
         rhoNegPz_m[i,nNegPz_m[i]] = rhoInit[k,indxFigures[i]]
         negPz_m[i,nNegPz_m[i]] = abs(deltaPz_m[k,indxFigures[i]]) 
#    print ('i=%d: nPosPz_c=%d, nNegPz_c=%d, nPosPz_m=%d, nNegPz_m=%d' % \
#              (i,nPosPz_c[i],nNegPz_c[i],nPosPz_m[i],nNegPz_m[i]))
	      
#
# Figures to compare calculated values of of deltaPz (their dependences 
# on impact parameter for different ion velocities) for two approaches
#
if (plotFigureFlag == 0):   
   for i in range(12):
      VionCrrnt = V0*VionRel[indxFigures[i]]
      powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
      mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
      helpValue_c =  int(nPosPz_c[i])      
      helpValue_m =  int(nPosPz_m[i])     
      figCrrnt = plt.figure(numbrFigures[i]+1)
      plt.loglog(rhoPosPz_c[i,0:nPosPz_c[i]],posPz_c[i,0:nPosPz_c[i]] ,'xb', \
                 rhoNegPz_c[i,0:nPosPz_c[i]],negPz_c[i,0:nPosPz_c[i]] ,'ob', \
                 rhoPosPz_m[i,0:nPosPz_m[i]],posPz_m[i,0:nPosPz_m[i]] ,'xr', \
                 rhoNegPz_m[i,0:nPosPz_m[i]],negPz_m[i,0:nPosPz_m[i]] ,'or',linewidth=2)
      plt.ylabel('$|\Delta P_z|$, $eV$', color='m',fontsize=14)
      plt.legend(['$\Delta P_z > 0$ (CG)','$\Delta P_z < 0$ (CG)', \
                  '$\Delta P_z > 0$ (ME)','$\Delta P_z < 0$ (ME)'], \
		  loc='lower left',fontsize=10)
      plt.xlabel('Initial Impact Parameter $rho_{Init}$, $cm$',color='m',fontsize=14)
      titleHeader = 'Transferred Momenta $\Delta P_z$ to Single Ion:'
      titleHeader += ' $V_{ion}=%3.1f\cdot10^{%2d}$ $cm/s$'
      plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=14)
      plt.xlim([.95*rhoInit[0,indxFigures[i]],1.05*rhoInit[nImpctPrmtr-1,indxFigures[i]]])
#      plt.text(xPos[i],yPos[i],'Fitted $\Delta E_{ion}$ are proportional to $rho_{Init}^{-B}$', \
#	       color='m',fontsize=16)      	
      plt.grid(True)
      if (saveFilesFlag == 1):
         fileName = 'picturesCMA/deltaEtransf_indxPlot-'+str(indxFigures[i])+'_fig'
         fileName += str(numbrFigures[i])+'cma.png' 
         figCrrnt.savefig(fileName) 
         print ('File "',fileName,'" is written')

#
# Px:
#
posPx_c = np.zeros((12,nImpctPrmtr)) 
rhoPosPx_c = np.zeros((12,nImpctPrmtr)) 
negPx_c = np.zeros((12,nImpctPrmtr)) 
rhoNegPx_c = np.zeros((12,nImpctPrmtr)) 
posPx_m = np.zeros((12,nImpctPrmtr)) 
rhoPosPx_m = np.zeros((12,nImpctPrmtr))
negPx_m = np.zeros((12,nImpctPrmtr)) 
rhoNegPx_m = np.zeros((12,nImpctPrmtr)) 
nPosPx_c = array('i',[0]*12)
nNegPx_c = array('i',[0]*12)
nPosPx_m = array('i',[0]*12)
nNegPx_m = array('i',[0]*12)
for i in range(12):
   nPosPx_c[i] = -1
   nNegPx_c[i] = -1
   nPosPx_m[i] = -1
   nNegPx_m[i] = -1
   for k in range(nImpctPrmtr):
      if (deltaPx_c[k,indxFigures[i]] > 0):
         nPosPx_c[i] += 1
         rhoPosPx_c[i,nPosPx_c[i]] = rhoInit[k,indxFigures[i]]
         posPx_c[i,nPosPx_c[i]] = deltaPx_c[k,indxFigures[i]] 
      if (deltaPx_c[k,indxFigures[i]] <= 0):
         nNegPx_c[i] += 1 
         rhoNegPx_c[i,nRegPx_c[i]] = rhoInit[k,indxFigures[i]]
         negPx_c[i,nRegPx_c[i]] = abs(deltaPx_c[k,indxFigures[i]]) 
      if (deltaPx_m[k,indxFigures[i]] > 0):
         nPosPx_m[i] += 1 
         rhoPosPx_m[i,nPosPx_m[i]] = rhoInit[k,indxFigures[i]]
         posPx_m[i,nPosPx_m[i]] = deltaPx_m[k,indxFigures[i]] 
      if (deltaPx_m[k,indxFigures[i]] <= 0):
         nNegPx_m[i] += 1 
         rhoNegPx_m[i,nNegPx_m[i]] = rhoInit[k,indxFigures[i]]
         negPx_m[i,nNegPx_m[i]] = abs(deltaPx_m[k,indxFigures[i]]) 
#       print ('nPosPx_c=%d, nNegPx_c=%d, nPosPx_m=%d, nNegPx_m=%d' % \
#              (nPosPx_c,nNegPx_c,nPosPx_m,nNegPx_m))
	      
#
# Comparison of calculated values of deltaPx (their dependences 
# on impact parameter for different ion velocities) for two approaches:
#
if (plotFigureFlag == 0):   
   for i in range(12):
      VionCrrnt = V0*VionRel[indxFigures[i]]
      powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
      mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
      figCrrnt = plt.figure(numbrFigures[i]+2)
      plt.loglog(rhoPosPx_c[i,0:nPosPx_c[i]],.99*posPx_c[i,0:nPosPx_c[i]] ,'xb', \
                 rhoPosPx_m[i,0:nPosPx_m[i]],1.01*posPx_m[i,0:nPosPx_m[i]] ,'xr',linewidth=2)
#      plt.loglog(rhoPosPx_c[i,0:nPosPx_c[i]],.99*posPx_c[i,0:nPosPx_c[i]] ,'xb', \
#                 rhoNegPx_c[i,0:nNegPx_c[i]],.99*negPx_c[i,0:nNegPx_c[i]] ,'ob', \
#                 rhoPosPx_m[i,0:nPosPx_m[i]],1.01*posPx_m[i,0:nPosPx_m[i]] ,'xr', \
#                 rhoNegPx_m[i,0:nNegPx_m[i]],1.01*negPx_m[i,0:nNegPx_m[i]] ,'or',linewidth=2)
      plt.ylabel('$\Delta P_x$, $eV$', color='m',fontsize=14)
#      plt.ylabel('$|\Delta P_x|$, $eV$', color='m',fontsize=14)
      plt.legend(['$0.99\cdot\Delta P_x$: CG - Center Guide', \
                  '$1.01\cdot\Delta P_x$: ME - Magnus Expansion'],loc='lower left',fontsize=10)
#      plt.legend(['$\Delta P_x > 0$ (CG)','$\Delta P_x < 0$ (CG)', \
#                  '$\Delta P_x > 0$ (ME)','$\Delta P_x < 0$ (ME)'], \
#		 loc='lower left',fontsize=11)
      plt.xlabel('Initial Impact Parameter $rho_{Init}$, $cm$',color='m',fontsize=14)
      titleHeader = 'Transferred Momenta $\Delta P_x$ to Single Ion:'
      titleHeader += ' $V_{ion}=%3.1f\cdot10^{%2d}$ $cm/s$'
      plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=14)
      plt.xlim([.95*rhoInit[0,indxFigures[i]],1.05*rhoInit[nImpctPrmtr-1,indxFigures[i]]])
#      plt.text(xPos[i],yPos[i],'Fitted $\Delta E_{ion}$ are proportional to $rho_{Init}^{-B}$', \
#	       color='m',fontsize=16)      	
      plt.grid(True)
      if (saveFilesFlag == 1):
         fileName = 'picturesCMA/deltaEtransf_indxPlot-'+str(indxFigures[i])+'_fig'
         fileName += str(numbrFigures[i])+'cma.png' 
         figCrrnt.savefig(fileName) 
         print ('File "',fileName,'" is written')

timeEnd = os.times()
timeIntgrtn = float(timeEnd[0])-float(timeStart[0])  # CPU time , sec
print ('Time of plotting = %6.3f seconds' % timeIntgrtn)

yPosText = [-2.12,-2.12,-2.12,-2.20,-2.12,-2.12,-2.12,-2.20,-2.12,-2.12,-2.12,-2.12]

#
# Dependence of calculated and fitted values of deltaPz on impact parameter
# for different ion velocities:
#
if (plotFigureFlag == 0):   
   for i in range(12):
      VionCrrnt = V0*VionRel[indxFigures[i]]
      powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
      mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
      figCrrnt = plt.figure(numbrFigures[i]+5)
      plt.loglog(rhoInit[0:nImpctPrmtr,indxFigures[i]], \
                 1.e24*deltaPz_m[0:nImpctPrmtr,indxFigures[i]],'xr', \
                 rhoInitFit_pz[0:nImpctPrmtr,indxFigures[i]], \
	         1.e24*deltaPz_m_fit[0:nImpctPrmtr,indxFigures[i]],'ob',linewidth=2)
      plt.xlabel('Initial Impact Parameter $rho_{Init}$, $cm$',color='m',fontsize=14)
      plt.ylabel('$10^{24} \cdot \Delta P_z$, $eV$', color='m',fontsize=14)
      titleHeader = 'Transferred Momenta $\Delta P_z$ to Single Ion:'
      titleHeader += ' $V_{ion}=%3.1f\cdot10^{%2d}$ $cm/s$'
      plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=14)
      plt.xlim([.95*rhoInit[0,indxFigures[i]],1.05*rhoInit[nImpctPrmtr-1,indxFigures[i]]])
      plt.ylim([ .9e24*deltaPz_m[nImpctPrmtr-1,indxFigures[i]], \
                1.1e24*deltaPz_m_fit[0,indxFigures[i]]])
      plt.legend(['Calculated Data', \
                  ('Fitting: $\Delta P_z=10^A\cdot$rho$_{init}^B$; B = %5.3f $\pm$ %5.3f' % \
                   (fitB_pz[indxFigures[i]],dNegB_pz[indxFigures[i]]))],loc='lower left',fontsize=11)
#      plt.text(xPos[i],yPos[i],'Fitted $\Delta P_z$ are proportional to $rho_{Init}^{-B}$', \
#	       color='m',fontsize=16)
      plt.grid(True)
      if (saveFilesFlag == 1):
         fileName = 'picturesCMA_v7/dPz_withFit_indxPlot'+str(indxFigures[i])+'_fig'
         fileName += str(numbrFigures[i]+5)+'cma.png' 
         figCrrnt.savefig(fileName) 
         print ('File "',fileName,'" is written')

#
# Dependence of calculated and fitted values of deltaPx on impact parameter
# for different ion velocities:
#
if (plotFigureFlag == 0):   
   for i in range(12):
      VionCrrnt = V0*VionRel[indxFigures[i]]
      powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
      mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
      figCrrnt = plt.figure(numbrFigures[i]+7)
      plt.loglog(rhoInit[0:nImpctPrmtr,indxFigures[i]], \
                 deltaPx_m[0:nImpctPrmtr,indxFigures[i]],'xr', \
                 rhoInitFit_px[0:nImpctPrmtr,indxFigures[i]], \
	         deltaPx_m_fit[0:nImpctPrmtr,indxFigures[i]],'ob',linewidth=2)
      plt.xlabel('Initial Impact Parameter $rho_{Init}$, $cm$',color='m',fontsize=14)
      plt.ylabel('$\Delta P_x$, $eV$', color='m',fontsize=14)
      titleHeader = 'Transferred Momenta $\Delta P_x$ to Single Ion:'
      titleHeader += ' $V_{ion}=%3.1f\cdot10^{%2d}$ $cm/s$'
      plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=14)
      plt.xlim([.95*rhoInit[0,indxFigures[i]],1.05*rhoInit[nImpctPrmtr-1,indxFigures[i]]])
      plt.ylim([ .9*deltaPx_m[nImpctPrmtr-1,indxFigures[i]], \
                1.1*deltaPx_m_fit[0,indxFigures[i]]])
      plt.legend(['Calculated Data', \
                  ('Fitting: $\Delta P_x=10^A\cdot rho_{init}^B$; B = %5.3f $\pm$ %5.3f' % \
                  (fitB_px[indxFigures[i]],dNegB_px[indxFigures[i]]))],loc='lower left',fontsize=11)
#      plt.text(xPos[i],yPos[i],'Fitted $\Delta P_x$ are proportional to $rho_{Init}^{-B}$', \
#	       color='m',fontsize=16)
      plt.grid(True)
      if (saveFilesFlag == 1):
         fileName = 'picturesCMA_v7/dPx_withFit_indxPlot-'+str(indxFigures[i])+'_fig'
         fileName += str(numbrFigures[i]+7)+'cma.png' 
         figCrrnt.savefig(fileName) 
         print ('File "',fileName,'" is written')

#
# Dependence of calculated and fitted values of deltaEnrgIon
# on impact parameter for different ion velocities:
#
if (plotFigureFlag == 0):   
   for i in range(12):
      VionCrrnt = V0*VionRel[indxFigures[i]]
      powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
      mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
      figCrrnt = plt.figure(numbrFigures[i]+6)
      plt.loglog(rhoInit[0:nImpctPrmtr,indxFigures[i]], \
                 1.e-18*deltaEnrgIon_m[0:nImpctPrmtr,indxFigures[i]],'xr', \
                 rhoInitFit_dEion[0:nImpctPrmtr,indxFigures[i]], \
	         1.e-18*deltaEnrgIon_m_fit[0:nImpctPrmtr,indxFigures[i]],'ob',linewidth=2)
      plt.xlabel('Initial Impact Parameter $rho_{Init}$, $cm$',color='m',fontsize=14)
      plt.ylabel('$10^{-18} \cdot \Delta E_{ion}$, $eV$', color='m',fontsize=14)
      titleHeader = 'Transferred Energy $\Delta E_{ion}$ to Single Ion:'
      titleHeader += ' $V_{ion}=%3.1f\cdot10^{%2d}$ $cm/s$'
      plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=14)
      plt.xlim([.95*rhoInit[0,indxFigures[i]],1.05*rhoInit[nImpctPrmtr-1,indxFigures[i]]])
      plt.ylim([ .9e-18*min(deltaEnrgIon_m[nImpctPrmtr-1,indxFigures[i]], \
                            deltaEnrgIon_m_fit[nImpctPrmtr-1,indxFigures[i]]), \
                1.1e-18*max(deltaEnrgIon_m[0,indxFigures[i]], \
		            deltaEnrgIon_m_fit[0,indxFigures[i]])])
      plt.legend(['Calculated Data', \
                  ('Fitting: $\Delta E_{ion}=10^A\cdot rho_{init}^B$; B = %5.3f $\pm$ %5.3f' % \
                   (fitB_dEion[indxFigures[i]],dNegB_dEion[indxFigures[i]]))],loc='lower left',fontsize=11)
      plt.grid(True)
      if (saveFilesFlag == 1):
         fileName = 'picturesCMA_v7/dEion_withFit_indxPlot'+str(indxFigures[i])+'_fig'
         fileName += str(numbrFigures[i]+6)+'cma.png' 
         figCrrnt.savefig(fileName) 
         print ('File "',fileName,'" is written')

#
# Dependence of fitted values of deltaEnrgIon on ion velocity
# for different impact parameters:
#
yPos = [1.45,6.1,3.3,2.01]
viewFctr = [1.,1.,1.,1.]
viewFctr = [1.e-19,1.e-18,1.e-17,1.e-17]
if (plotFigureFlag == 0):   
   for k in range(nRhoSlctd):
      powViewFctr = math.floor(np.log10(viewFctr[k])) 
      mantViewFctr = viewFctr[k]/(10**powViewFctr) 
      figCrrnt=plt.figure (7000+100*k)
      plt.semilogx(VionRel[npStart[k]:nVion-1], \
                   viewFctr[k]*deltaEnrgIon_dpnd_Vi[k,npStart[k]:nVion-1],'.r')
      plt.xlabel('Relative Ion Velocity, $V_i/V_{e0}$',color='m',fontsize=16)
      plt.ylabel('$10^{%2d} \cdot \Delta E_{ion}$, eV' % powViewFctr,color='m',fontsize=16)
      titleHeader = 'Transferred Energy $\Delta E_{ion}$ to Single Ion: $rho=%5.3f$ cm'
      titleHeader += '\n$\Delta E_{ion}=10^A \cdot rho^{-B}$ with $A=%7.3f$, $B=%5.3f$' 
      plt.title(titleHeader % (rhoSlctd[k],fitA_dEion[k],fitB_dEion[k]),color='m',fontsize=14)
      xLimit = [.95*VionRel[npStart[k]],1.05*VionRel[nVion-1]]
      plt.xlim(xLimit)
      yLimit = [0.99*viewFctr[k]*deltaEnrgIon_dpnd_Vi[k,npStart[k]], \
                1.01*viewFctr[k]*deltaEnrgIon_dpnd_Vi[k,nVion-1]]
      plt.ylim(yLimit)
      if ((relVeLong >= xLimit[0]) and (relVeLong <= xLimit[1])):
         plt.plot([relVeLong,relVeLong],yLimit,'--m',linewidth=1)
         plt.text(3.8e-5,yPos[k],'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
      if ((relVeTrnsv >= xLimit[0]) and (relVeTrnsv <= xLimit[1])):
         plt.plot([relVeTrnsv,relVeTrnsv],yLimit,'--m',linewidth=1)
         plt.text(1.6e-3,yPos[k],'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
      plt.grid(True)
      if (saveFilesFlag == 1):
         fileName = 'picturesCMA_v7/deltaEtransfOnVion_rhoIndx-'+str(k)+'_fig'
         fileName += str(7000+100*k)+'cma.png' 
         figCrrnt.savefig(fileName) 
         print ('File "',fileName,'" is written')

#----------------------------------------------------
#
# Integration of the transferred energy to ion from electron beam
#
# "Gauss-Kronrod" (GK) method is used
# Technic of integration see above (function fittedGKintegration)
#----------------------------------------------------
#
# Data for GK:
#
#----------------------------------------------------

timeStart = os.times()
nPointsGK = 16

frctnForce_GK = np.zeros(nVion)   # integration using "Gauss-Kronrod" method 

psi16=np.array([-0.9894009, -0.9445750, -0.8656312, -0.7554044, -0.6178762, \
                -0.4580168, -0.2816036, -0.0950125,  0.0950125,  0.2816036, \
		 0.4580168,  0.6178762,  0.7554044,  0.8656312,  0.9445750, \
		 0.9894009])
w16  =np.array([ 0.0271525,  0.0622535,  0.0951585,  0.1246290,  0.1495960, \
                 0.1691565,  0.1826034,  0.1894506,  0.1894506,  0.1826034, \
		 0.1691565,  0.1495960,  0.1246290,  0.0951585,  0.0622535, \
		 0.0271525])

rhoCrrntGK = np.zeros((nPointsGK,nVion))
deltaEnrgIon_cGK = np.zeros((nPointsGK,nVion))
deltaPx_cGK = np.zeros((nPointsGK,nVion))  
deltaPy_cGK = np.zeros((nPointsGK,nVion))  
deltaPz_cGK = np.zeros((nPointsGK,nVion))  
deltaEnrgIon_GK = np.zeros((nPointsGK,nVion))
deltaEnrgIon_GKfit = np.zeros((nPointsGK,nVion))
deltaPx_mGK = np.zeros((nPointsGK,nVion))  
deltaPy_mGK = np.zeros((nPointsGK,nVion))  
deltaPz_mGK = np.zeros((nPointsGK,nVion))  
ionVx_cGK = np.zeros((nPointsGK,nVion))
ionVy_cGK = np.zeros((nPointsGK,nVion))
ionVz_cGK = np.zeros((nPointsGK,nVion))
ionVx_mGK = np.zeros((nPointsGK,nVion))
ionVy_mGK = np.zeros((nPointsGK,nVion))
ionVz_mGK = np.zeros((nPointsGK,nVion))

indx = 0
totalPointsIntgrtn = 0
for i in range(nVion):
#
# Some possible correction of maximal values of the impact parameter (on depence
# of preset number of minimal Larmor turns) for each value of ion velocity
# (see for example figure 10):
#
   rhoMaxCrrnt = impctPrmtrMax[i]* \
   np.sqrt(1.- (pi*larmorTurns*eVrmsLong/omega_L/impctPrmtrMax[i])**2)
#
# The above correction is not taken into account:
   rhoMaxCrrnt = impctPrmtrMax[i]
#
   for n in range(nPointsGK):
      rhoCrrntGK[n,i] = psi16[n]*(rhoMaxCrrnt-rhoMin)/2 + \
                       (rhoMaxCrrnt+rhoMin)/2
#      print ('    rhoCrrntGK[%2d,%2d] = %e' % (n,i,rhoCrrntGK[n,i]))   
# Half length of interaction (cm):
      halfLintrCrrnt = np.sqrt(rhoMaxCrrnt**2-rhoCrrntGK[n,i]**2)
# 0.5 time of interaction (sec):
      timeHalfPath = halfLintrCrrnt/eVrmsLong     
      numbLarmor = int(2.*timeHalfPath/T_larm)             
      pointAlongTrackCrrnt = int(2.*timeHalfPath/timeStep_c)
      totalPointsIntgrtn += pointAlongTrackCrrnt
      z_ionCrrnt_c = np.zeros(6)      # Zeroing out of vector for ion
      z_elecCrrnt_c = np.zeros(6)     # Zeroing out of vector for electron
      z_ionCrrnt_m = np.zeros(6)      # Zeroing out of vector for ion
      z_elecCrrnt_m = np.zeros(6)     # Zeroing out of vector for electron
# Zeroing out of "guiding center" vector for electron:
      z_elecCrrnt_gc_c = np.zeros(6)  
# Zeroing out of "Magnus expansion" vector for electron:
      z_elecCrrnt_gc_m = np.zeros(6)  
# Current values to transfered momemta 
# (second index numerates "Guiding Center" approach, if = 0, and 
#                         "Magnus Expansion" approach  if = 1: 
      dpCrrnt = np.zeros((3,2))
# Intermediate arrays:
      dpIon_c = np.zeros(3) 
      dpIon_m = np.zeros(3) 
      dpElec_c = np.zeros(3) 
      dpElec_m = np.zeros(3) 
# Current initial vector for electron:
      z_elecCrrnt_c[Ix] = rhoCrrntGK[n,i]              # x, cm
      z_elecCrrnt_c[Iz] = -halfLintrCrrnt              # z, cm
      z_elecCrrnt_c[Ipy] = m_elec*eVrmsTran            # py, g*cm/sec
      z_elecCrrnt_c[Ipz] = m_elec*eVrmsLong            # pz, g*cm/sec
      z_elecCrrnt_m[Ix] = rhoCrrntGK[n,i]              # x, cm
      z_elecCrrnt_m[Iz] = -halfLintrCrrnt              # z, cm
      z_elecCrrnt_m[Ipy] = m_elec*eVrmsTran            # py, g*cm/sec
      z_elecCrrnt_m[Ipz] = m_elec*eVrmsLong            # pz, g*cm/sec
# Current vector for ion velocities for both approaches:
      ionVx_cGK[n,i] = VionTrnsv[i]*np.cos(phiVi)
      ionVy_cGK[n,i] = VionTrnsv[i]*np.sin(phiVi) 
      ionVz_cGK[n,i] = VionLong[i]
      ionVx_mGK[n,i] = VionTrnsv[i]*np.cos(phiVi)
      ionVy_mGK[n,i] = VionTrnsv[i]*np.sin(phiVi) 
      ionVz_mGK[n,i] = VionLong[i]
# transfer to system of guiding center:
      z_elecCrrnt_gc_c=toGuidingCenter(z_elecCrrnt_c)  
      z_elecCrrnt_gc_m=toGuidingCenter(z_elecCrrnt_m)  
#
# Main loop along the each track:
#
      for k in range(int(pointAlongTrackCrrnt)):
#
# Dragging both electrons through first half of the step of the track:
#
         z_elecCrrnt_gc_c = np.dot(matr_elec_c,z_elecCrrnt_gc_c) # electron
         z_ionCrrnt_c = np.dot(matr_ion_c,z_ionCrrnt_c)          # ion
         z_elecCrrnt_gc_m = np.dot(matr_elec_c,z_elecCrrnt_gc_m) # electron
         z_ionCrrnt_m = np.dot(matr_ion_c,z_ionCrrnt_m)          # ion
# transfer from system of guiding center: 
         z_elecCrrnt_c=fromGuidingCenter(z_elecCrrnt_gc_c)     
         z_elecCrrnt_m=fromGuidingCenter(z_elecCrrnt_gc_m)     
# Current distance between ion and electron; cm:
         bCrrnt_c[indx]=np.sqrt((z_ionCrrnt_c[0]-z_elecCrrnt_c[0])**2+ \
                          (z_ionCrrnt_c[2]-z_elecCrrnt_c[2])**2+ \
                          (z_ionCrrnt_c[4]-z_elecCrrnt_c[4])**2)
         indx += 1
#
# Dragging electrons through interaction during this step of track
# (for both approaches):
#
#    "Guiding Center":
         dpIon_c,dpElec_c,action,b_gc_c = \
                 guidingCenterCollision(z_elecCrrnt_gc_c,z_ionCrrnt_c,timeStep_c) 
#    "Magnus Expansion":
         dpIon_m,dpElec_m,action,dy_gc_m,C1,C2,C3,b,D1,D2,q = \
                 MagnusExpansionCollision(z_elecCrrnt_gc_m,z_ionCrrnt_m,timeStep_c) 
#
# Taking into account transfer of momentum for both particles
# and both approaches:
#
         if (dpTransferFlag == 1):
            for ic in range(3):
               z_ionCrrnt_c[2*ic+1] += dpIon_c[ic]   
               z_elecCrrnt_c[2*ic+1] += dpElec_c[ic]
               z_ionCrrnt_m[2*ic+1] += dpIon_m[ic]   
               z_elecCrrnt_m[2*ic+1] += dpElec_m[ic]
# Dragging ion velocities for both approaches
         ionVx_cGK[n,i] += dpIon_c[0]/M_ion                    # cm/sec
         ionVy_cGK[n,i] += dpIon_c[1]/M_ion                    # cm/sec
         ionVz_cGK[n,i] += dpIon_c[2]/M_ion                    # cm/sec
         ionVx_mGK[n,i] += dpIon_m[0]/M_ion                    # cm/sec
         ionVy_mGK[n,i] += dpIon_m[1]/M_ion                    # cm/sec
         ionVz_mGK[n,i] += dpIon_m[2]/M_ion                    # cm/sec
# transfer to system of guiding center:
         z_elecCrrnt_gc_c=toGuidingCenter(z_elecCrrnt_c)  
         z_elecCrrnt_gc_m=toGuidingCenter(z_elecCrrnt_m)  
# Accumulation of the transfered momenta along the track:  
         for ic in range(3):
#	    if i == 0:
#	       print ('dpIon_c[%2d] = %20.14e, dpIon_m[%2d] = %20.14e' % \
#	             (ic,dpIon_c[ic],ic,dpIon_m[ic]))
            dpCrrnt[ic,0] += dpIon_c[ic]                      # g*cm/csec  
            dpCrrnt[ic,1] += dpIon_m[ic]                      # g*cm/csec  
# Dragging ion velocities for both approaches
         ionVx_cGK[n,i] += dpIon_c[0]/M_ion                    # cm/sec
         ionVy_cGK[n,i] += dpIon_c[1]/M_ion                    # cm/sec
         ionVz_cGK[n,i] += dpIon_c[2]/M_ion                    # cm/sec
         ionVx_mGK[n,i] += dpIon_m[0]/M_ion                    # cm/sec
         ionVy_mGK[n,i] += dpIon_m[1]/M_ion                    # cm/sec
         ionVz_mGK[n,i] += dpIon_m[2]/M_ion                    # cm/sec
#
# Dragging both particles through second half of the step of the track:
#
         z_elecCrrnt_gc_c = np.dot(matr_elec_c,z_elecCrrnt_gc_c)     # electron
         z_ionCrrnt_c = np.dot(matr_ion_c,z_ionCrrnt_c)              # ion
         z_elecCrrnt_gc_m = np.dot(matr_elec_c,z_elecCrrnt_gc_m)     # electron
         z_ionCrrnt_m = np.dot(matr_ion_c,z_ionCrrnt_m)              # ion
# transfer from system of guiding center: 
         z_elecCrrnt_c=fromGuidingCenter(z_elecCrrnt_gc_c)     
         z_elecCrrnt_m=fromGuidingCenter(z_elecCrrnt_gc_m)     
# Current distance between ion and electron; cm:
         bCrrnt_c[indx]=np.sqrt((z_ionCrrnt_c[0]-z_elecCrrnt_c[0])**2+ \
                          (z_ionCrrnt_c[2]-z_elecCrrnt_c[2])**2+ \
                          (z_ionCrrnt_c[4]-z_elecCrrnt_c[4])**2)
         indx += 1
#	 
# Transferred momenta and energy along the entire length of each track
# for both approaches:
#
      deltaPx_cGK[n,i] = dpCrrnt[0,0] 
      deltaPx_mGK[n,i] = dpCrrnt[0,1] 
#      if deltaPx_cGK[n,i] <= 0.: 
#         print ('deltaPx_cGK[%2d,%2d] = %e, dpCrrnt[%2d,%2d] = %e' % \
#	       (n,i,deltaPx_cGK[n,i],n,i,dpCrrnt[0,0]))
#      if deltaPx_mGK[n,i] <= 0.: 
#         print ('deltaPx_mGK[%2d,%2d] = %e, dpCrrnt[%2d,%2d] = %e' % \
#	       (n,i,deltaPx_mGK[n,i],n,i,dpCrrnt[0,1]))

      deltaPy_cGK[n,i] = dpCrrnt[1,0] 
      deltaPy_mGK[n,i] = dpCrrnt[1,1] 
#      if deltaPy_cGK[n,i] <= 0.: 
#         print ('deltaPy_cGK[%2d,%2d] = %e' % (n,i,deltaPy_cGK[n,i]))
#      if deltaPy_mGK[n,i] <= 0.: 
#         print ('deltaPy_mGK[%2d,%2d] = %e' % (n,i,deltaPy_mGK[n,i]))

      deltaPz_cGK[n,i] = dpCrrnt[2,0] 
      deltaPz_mGK[n,i] = dpCrrnt[2,1] 
#      if deltaPz_cGK[n,i] <= 0.: 
#         print ('deltaPz_cGK[%2d,%2d] = %e' % (n,i,deltaPz_cGK[n,i]))
#      if deltaPz_mGK[n,i] <= 0.: 
#         print ('deltaPz_mGK[%2d,%2d] = %e' % (n,i,deltaPz_MGK[n,i]))

# Incorrect value:
#      deltaEnrgIon_GK[n,i] = (dpCrrnt[0,0]**2+dpCrrnt[1,0]**2+dpCrrnt[2,0]**2)* \
#                            deFactor/eVtoErg                      # eV
# Correct value absolute value):
      deltaEnrgIon_GK[n,i] = (dpCrrnt[0,0]*ionVx_cGK[n,i]+ \
                               dpCrrnt[1,0]*ionVy_cGK[n,i]+ \
			       dpCrrnt[2,0]*ionVz_cGK[n,i])* deFactor/eVtoErg   # eV


      deltaPx_mGK[n,i] = dpCrrnt[0,1]
#      if deltaPx_mGK[n,i] <= 0.: 
#         print ('deltaPx_mGK[%2d,%2d] = %e' % (n,i,deltaPx_mGK[n,i]))

      deltaPy_mGK[n,i] = dpCrrnt[1,1]
#      if deltaPy_mGK[n,i] <= 0.: 
#         print ('deltaPy_mGK[%2d,%2d] = %e' % (n,i,deltaPy_mGK[n,i]))

      deltaPz_mGK[n,i] = dpCrrnt[2,1]
#      if deltaPz_mGK[n,i] <= 0.: 
#         print ('deltaPz_mGK[%2d,%2d] = %e' % (n,i,deltaPz_mGK[n,i]))

#
# Integration using "Gauss-Kronrod" method:
#
      frctnForce_GK[i] += pi*n_e*100.*(rhoMaxCrrnt-rhoMin)*w16[n]* \
                          deltaEnrgIon_GK[n,i]*rhoCrrntGK[n,i]                # eV/m

timeEnd = os.times()
timeIntgrtn = float(timeEnd[0])-float(timeStart[0])  # CPU time , sec
print ('Time of GK-Integration = %6.3f seconds' % timeIntgrtn)

xVionRelIntgr = np.zeros((nPointsGK,nVion))
for i in range(nVion):
   for n in range(nPointsGK):
       xVionRelIntgr[n,i] = VionRel[i]

if (plotFigureFlag == 0):   
   fig41=plt.figure (41)
   for i in range(nVion):
      plt.semilogx(xVionRelIntgr[0:nPointsGK,i],rhoCrrntGK[0:nPointsGK,i],'.r')
   plt.xlabel('Relative Ion Velocity, $V_i/V_{e0}$',color='m',fontsize=16)
   plt.ylabel('$rho_{Init}$, cm',color='m',fontsize=16)
   plt.title('Subdivisions for $rho_{Init}$: Gauss-Kronrod Method Integration', \
             color='m',fontsize=14)
   plt.grid(True)
   yLimit=[0.,max(rhoCrrntGK[0:nPointsGK,nVion-1])+.01]
   plt.ylim(yLimit)
   plt.plot([relVeTrnsv,relVeTrnsv],yLimit,'--m',linewidth=1)
   plt.text(1.6e-3,-.03,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([relVeLong,relVeLong],yLimit,'--m',linewidth=1)
   plt.text(4.4e-5,.05,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   if (saveFilesFlag == 1):
      fig41.savefig('picturesCMA/initialImpactParameter_GK_fig41cma.png') 
      print ('File "picturesCMA/initialImpactParameter_GK_fig41cma.png" is written')


frctnForce_GKfit = np.zeros(nVion)
# frctnForce_GKfit1 = np.zeros(nVion)
for i in range(nVion):
   rhoMaxCrrnt = impctPrmtrMax[i]
   arrayTemp,valueTemp = fittedGKintegration(rhoMin,rhoMaxCrrnt,fitA_dEion[i],fitB_dEion[i])
   for n in range(nPointsGK):
      deltaEnrgIon_GKfit[n,i] = arrayTemp[n]  
   frctnForce_GKfit[i] = pi*n_e*100.*valueTemp                      # eV/m
#    for n in range(nPointsGK):
#       rhoCrrntGK[n,i] = psi16[n]*(rhoMaxCrrnt-rhoMin)/2 + \
#                        (rhoMaxCrrnt+rhoMin)/2
#       factorA = math.pow(10.,fitA_dEion[i])
#       deltaEnrgIon_GKfit[n,i] = factorA*math.pow(rhoCrrntGK[n,i],fitB_dEion[i])
#
# Figures show dependence of calculated and fitted values of deltaEnrgIon
# on impact parameter for different ion velocities:
#
if (plotFigureFlag == 0):   
   for i in range(12):
      VionCrrnt = V0*VionRel[indxFigures[i]]
      powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
      mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
      figCrrnt = plt.figure(numbrFigures[i]+3)
      plt.loglog(rhoInit[0:nImpctPrmtr,indxFigures[i]], \
                 1.e-18*deltaEnrgIon_m[0:nImpctPrmtr,indxFigures[i]],'xr', \
                 rhoInitFit_dEion[0:nImpctPrmtr,indxFigures[i]], \
	         1.e-18*deltaEnrgIon_m_fit[0:nImpctPrmtr,indxFigures[i]],'ob', \
		 rhoCrrntGK[0:int(nPointsGK),indxFigures[i]], \
		 1.e-18*deltaEnrgIon_GKfit[0:int(nPointsGK),indxFigures[i]],'om',linewidth=2)
      plt.xlabel('Initial Impact Parameter $rho_{Init}$, $cm$',color='m',fontsize=14)
      plt.ylabel('$10^{-18} \cdot \Delta E_{ion}$, $eV$', color='m',fontsize=14)
      titleHeader = 'Transferred Energy $\Delta E_{ion}$ to Single Ion:'
      titleHeader += ' $V_{ion}=%3.1f\cdot10^{%2d}$ $cm/s$'
      plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=14)
      plt.xlim([.95*rhoInit[0,indxFigures[i]],1.05*rhoInit[nImpctPrmtr-1,indxFigures[i]]])
      plt.ylim([ .9e-18*min(deltaEnrgIon_m[nImpctPrmtr-1,indxFigures[i]], \
                            deltaEnrgIon_m_fit[nImpctPrmtr-1,indxFigures[i]], \
			    deltaEnrgIon_GK[int(nPointsGK)-1,indxFigures[i]]), \
                1.1e-18*max(deltaEnrgIon_m[0,indxFigures[i]], \
		            deltaEnrgIon_m_fit[0,indxFigures[i]], \
			    deltaEnrgIon_GKfit[0,indxFigures[i]])])
      plt.legend(['Calculated Data', \
                  ('Fitting: $\Delta E_{ion}=10^A\cdot rho_{init}^B$; B = %5.3f $\pm$ %5.3f' % \
                   (fitB_dEion[indxFigures[i]],dNegB_dEion[indxFigures[i]])), \
		   'Calculated Data for GK Integtation'],loc='lower left',fontsize=11)
      plt.grid(True)
      if (saveFilesFlag == 1):
         fileName = 'picturesCMA_v7/dPz_withFit_indxPlot'+str(indxFigures[i])+'_fig'
         fileName += str(numbrFigures[i]+3)+'cma.png' 
         figCrrnt.savefig(fileName) 
         print ('File "',fileName,'" is written')


indxFigures = [0,9,12,18,19,23,27,29,31,34,39,49]

if (plotFigureFlag == 0):   
   for i in range(2,10,1):
      VionCrrnt = V0*VionRel[indxFigures[i]]
      powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
      mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
      figCrrnt = plt.figure (3000+indxFigures[i])
      plt.semilogx(rhoInit[:,indxFigures[i]],dEion_c_m[:,indxFigures[i]],'-xr',linewidth=2)
      plt.xlabel('Initial Impact Parameter $rho_{Init}$, $cm$',color='m',fontsize=14)
      plt.ylabel('Difference of $\Delta E_{ion}$: "GC"/"ME" $-$ 1, %',color='m',fontsize=14)
      titleHeader = 'Comparison of Approaches for $\Delta E_{ion}$:'
      titleHeader += ' $V_{ion}=%3.1f\cdot10^{%2d}$ $cm/s$'
      plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=14)
      plt.xlim([.95*rhoInit[0,indxFigures[i]],1.1*rhoInit[nImpctPrmtr-1,indxFigures[i]]])
   #   plt.ylim([min(dEion_c_m[:,indxFigures[i]])-.0005,max(dEion_c_m[:,indxFigures[i]])+.0005])
      plt.grid(True)
      if (saveFilesFlag == 1):
         fileName = 'picturesCMA/deltaEcomprsn_indxPlot-'+str(indxFigures[i])+'_fig30'
         fileName += str(indxFigures[i])+'cma.png' 
         figCrrnt.savefig(fileName) 
         print ('File "',fileName,'" is written')

#
# Maps for relative differences of for energy dEion and momenta dPx, dPy, dPz:
#
X = np.zeros((nVion,nImpctPrmtr))
Y = np.zeros((nImpctPrmtr,nVion))
Z = np.zeros((nVion,nImpctPrmtr,5))

for i in range(nVion):
   for n in range(nImpctPrmtr):
      X[n,i] = np.log10(VionRel[i])
      Y[n,i] = np.log10(rhoInit[n,i])
      Z[i,n,0] = dEion_c_m[i,n]                    # dEion
      Z[i,n,1] = deltaPx_c_m[i,n]                  # dPx
      Z[i,n,2] = deltaPy_c_m[i,n]                  # dPy
      Z[i,n,3] = deltaPz_c_m[i,n]                  # dPz
      Z[i,n,4] = np.log10(abs(.01*deltaPz_c_m[i,n]+1.))                 # dPz
yLimit = [-2.92,-0.26]

if (plotFigureFlag == 0): 
   for k in range(1,5,1):  
      figCrrnt = plt.figure(1245+100*k)
      ax = figCrrnt.add_subplot(111)                                       # for contours plotting
      mapCrrnt = ax.contourf(X,Y,Z[:,:,k],cmap='jet') 
      plt.plot(np.log10(VionRel),np.log10(impctPrmtrMax),'-r',linewidth=2)
      plt.plot([np.log10(VionRel[0]),np.log10(VionRel[nVion-1])], \
               [np.log10(impctPrmtrMin),np.log10(impctPrmtrMin)],'-r',linewidth=2)
      plt.plot([np.log10(relVeTrnsv),np.log10(relVeTrnsv)],yLimit,'--m',linewidth=1)
      plt.plot([np.log10(relVeLong), np.log10(relVeLong)], yLimit,'--m',linewidth=1)
      plt.xlabel('Relative Ion Velocity,  $log_{10}(V_{ion}/V_0$)',color='m',fontsize=16)
      plt.ylabel('Initial Impact Parameter, $log_{10}(rho_{Init})$',color='m',fontsize=16)
      titleHeader = 'Difference'
      fileName = 'picturesCMA_v7/delta'
      if (k == 0):
         titleHeader += ' of $\Delta E_{ion}$: "GC"/"ME" $-$ 1, %' 
         fileName += 'EcomprsnMap_fig'
      if (k == 1):
         titleHeader += ' of $\Delta p_x$: "GC"/"ME" $-$ 1, %' 
         fileName += 'PXcomprsnMap_fig'
      if (k == 2):
         titleHeader += ' of $\Delta p_y$: "GC"/"ME" $-$ 1, %' 
         fileName += 'PYcomprsnMap_fig'
      if (k == 3):
         titleHeader += ' of $\Delta p_z$: "GC"/"ME" $-$ 1, %' 
         fileName += 'PZcomprsnMap_fig'
      if (k == 4):
         titleHeader += ' of $\Delta p_z$: $log_{10}($"GC"/"ME"$)$' 
         fileName += 'PZcomprsnMap_fig'
      plt.title(titleHeader,color='m',fontsize=16)
      plt.text(-4.14,-.6,'Screened Collisions',color='r',fontsize=16)
      plt.text(-4.14,-2.2,'Magnetized Collisions',color='k',fontsize=16)
      plt.text(-3.25,-1.3,'$R_{max}$',color='k',fontsize=16)
      plt.text(-3.25,-2.89,'$R_{min}$',color='k',fontsize=16)
      plt.text(-2.85,-0.95,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
      plt.text(-4.5,-0.95,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
      plt.ylim(yLimit)
      figCrrnt.colorbar(mapCrrnt)
      plt.grid(True)
      fileName += str(1245+100*k)
      fileName += 'cma.png'
      if (saveFilesFlag == 1):
         figCrrnt.savefig(fileName) 
         print ('File "',fileName,'" is written')

#--------------------------------------------------------
#
# cutoff approach is a bad decision!
#

# cutoffLevel = [.05,.05,.75,-60.]
# cutoffZ = np.zeros((nVion,nImpctPrmtr,4))
# 
# for i in range(nVion):
#    for n in range(nImpctPrmtr):
#       cutoffZ[i,n,0] = dEion_c_m[i,n]                    # dEion
#       cutoffZ[i,n,1] = deltaPx_c_m[i,n]                  # dPx
#       cutoffZ[i,n,2] = deltaPy_c_m[i,n]                  # dPy
#       cutoffZ[i,n,3] = deltaPz_c_m[i,n]                  # dPz
#       for k in range(3):
#          if (cutoffZ[i,n,k] > cutoffLevel[k]):
#             cutoffZ[i,n,k] = cutoffLevel[k]
#       if (cutoffZ[i,n,3] < cutoffLevel[3]):
#          cutoffZ[i,n,3] = cutoffLevel[3]
#         
# if (plotFigureFlag == 0): 
#    for k in range(4):  
#       figCrrnt = plt.figure(246+100*k)
#       ax = figCrrnt.add_subplot(111)                      # for contours plotting
#       mapCrrnt_co = ax.contourf(X,Y,cutoffZ[:,:,k],cmap='jet') 
#       mapCrrnt_cl = ax.contour(X,Y,cutoffZ[:,:,k],8,colors='black') 
#       plt.clabel(mapCrrnt_cl,fmt='%4.2f',inline=True)
#       plt.plot(np.log10(VionRel),np.log10(impctPrmtrMax),'-r',linewidth=2)
#       plt.plot([np.log10(VionRel[0]),np.log10(VionRel[nVion-1])], \
#                [np.log10(impctPrmtrMin),np.log10(impctPrmtrMin)],'-r',linewidth=2)
#       plt.plot([np.log10(relVeTrnsv),np.log10(relVeTrnsv)],yLimit,'--m',linewidth=1)
#       plt.plot([np.log10(relVeLong), np.log10(relVeLong)], yLimit,'--m',linewidth=1)
#       plt.xlabel('Relative Ion Velocity,  $log_{10}(V_{ion}/V_0$)',color='m',fontsize=16)
#       plt.ylabel('Initial Impact Parameter, $log_{10}(rho_{Init})$',color='m',fontsize=16)
#       titleHeader = 'Difference of $\Delta'
#       fileName = 'picturesCMA/delta'
#       if (k == 0):
#          titleHeader += ' E_{ion}$: $\widetilde{\Delta E}_{ion}$ = "GC"/"ME" $-$ 1, %'
#          titleHeader += '\n($\widetilde{\Delta E}_{ion}$ = .05, ' 
#          titleHeader += 'if $\widetilde{\Delta E}_{ion}$ > .05)'
#          fileName += 'EcomprsnMap_fig'
#       if (k == 1):
#          titleHeader += ' p_x$: $\widetilde{\Delta p}_x$ = "GC"/"ME" $-$ 1, %' 
#          titleHeader += '\n($\widetilde{\Delta p}_x$ = .05, if $\widetilde{\Delta p}_x$ > .5)'
#          fileName += 'PXcomprsnMap_fig'
#       if (k == 2):
#          titleHeader += ' p_y$: $\widetilde{\Delta p}_y$ = "GC"/"ME" $-$ 1, %' 
#          titleHeader += '\n($\widetilde{\Delta p}_y$ = .75, if $\widetilde{\Delta p}_y$ > .75)' 
#          fileName += 'PYcomprsnMap_fig'
#       if (k == 3):
#          titleHeader += ' p_z$: $\widetilde{\Delta p}_z$ = "GC"/"ME" $-$ 1, %' 
#          titleHeader += '\n($\widetilde{\Delta p}_z$ = -60, if $\widetilde\Delta {p}_z$ > -60)' 
#          fileName += 'PZcomprsnMap_fig'
#       plt.title(titleHeader,color='m',fontsize=16)
#       plt.text(-4.14,-.45,'Screened Collisions',color='r',fontsize=16)
#       plt.text(-3.25,-1.3,'$R_{max}$',color='k',fontsize=16)
#       plt.text(-3.25,-2.89,'$R_{min}$',color='k',fontsize=16)
#       plt.text(-2.85,-0.95,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
#       plt.text(-4.5,-0.95,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
#       plt.ylim(yLimit)
#       figCrrnt.colorbar(mapCrrnt_co)
#       plt.grid(True)
#       fileName += str(246+100*k)
#       fileName += 'cma.png'
#       if (saveFilesFlag == 1):
#          figCrrnt.savefig(fileName) 
#          print ('File "',fileName,'" is written')
#
#--------------------------------------------------------


viewFctr = 1.e-27   
powViewrFctr=round(np.log10(viewFctr)) 
mantViewrFctr=viewFctr/(10**powViewrFctr) 

if (plotFigureFlag == 0):   
   fig5000=plt.figure (5000)
   plt.semilogx(VionRel,viewFctr*frctnForce_GKfit,'-xr',linewidth=2)
   plt.xlabel('Relative Ion Velocity  $V_{ion}/V_{e0}$',color='m',fontsize=14)
   if (mantViewrFctr == 1.):
      plt.ylabel(('$10^{%2d} \cdot F_{ion}$, eV/m' % powViewrFctr), color='m',fontsize=14)
   else:
      plt.ylabel(('$%3.1f \cdot 10^{%2d} \cdot F_{ion}$, eV/m' % \
                  (mantViewrFctr,powViewrFctr)), color='m',fontsize=14)
   plt.title('Friction Force $F_{ion}$: "Gauss-Kronrod" Integration',color='m',fontsize=14)
#   plt.text(1.2e-5,7.5, \
#            ('$V_{e0}=%4.2f\cdot10^{%2d}$cm/s, $n_e=%4.2f\cdot10^{%2d}$cm$^3$, $B=%4d$kG' % \
#             (mantV0,powV0,mant_n_e,pow_n_e,int(fieldB[0]))),color='m',fontsize=14)
   plt.xlim([.9*VionRel[0],1.1*VionRel[nVion-1]])
#   yLimit=[min(.9*viewFctr*frctnForce_GKfit),max(1.1*viewFctr*frctnForce_GKfit)]
#    print ('ylim[0] = %e, ylim[1] = %e' % (yLimit[0],yLimit[1]))
#   plt.ylim(yLimit)
   plt.plot([relVeTrnsv,relVeTrnsv],yLimit,'--m',linewidth=1)
#   plt.text(1.6e-3,yLimit[0]-.3,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([relVeLong,relVeLong],yLimit,'--m',linewidth=1)
#   plt.text(4.4e-5,yLimit[0]+.1,'$ \Delta V_{e||}/ sV_{e0}$',color='m',fontsize=14)
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig5000.savefig('picturesCMA/frctnForce_GKfit_fig5000cma.png')    
      print ('File "picturesCMA/frctnForce_GKfit_fig5000cma.png" is written')

if (plotFigureFlag == 0):   
   fig5010=plt.figure (5010)
   plt.semilogx(VionRel,viewFctr*frctnForce_AI,'-xr',linewidth=2)
   plt.xlabel('Relative Ion Velocity  $V_{ion}/V_{e}0$',color='m',fontsize=14)
   if (mantViewrFctr == 1.):
      plt.ylabel(('$10^{%2d} \cdot F_{ion}$, eV/m' % powViewrFctr), color='m',fontsize=14)
   else:
      plt.ylabel(('$%3.1f \cdot 10^{%2d} \cdot F_{ion}$, eV/m' % \
                  (mantViewrFctr,powViewrFctr)), color='m',fontsize=14)
   plt.title('Friction Force $F_{ion}$: "Analytical" Integration for "ME" Approach', \
             color='m',fontsize=14)
#   plt.text(1.2e-5,7.15, \
#            ('$V_{e0}=%4.2f\cdot10^{%2d}$cm/s, $n_e=%4.2f\cdot10^{%2d}$cm$^3$, $B=%4d$kG' % \
#             (mantV0,powV0,mant_n_e,pow_n_e,int(fieldB[0]))),color='m',fontsize=14)
   plt.xlim([.9*VionRel[0],1.1*VionRel[nVion-1]])
#   yLimit=[min(.9*viewFctr*frctnForce_AI),max(1.1*viewFctr*frctnForce_AI)]
#   plt.ylim(yLimit)
   plt.plot([relVeTrnsv,relVeTrnsv],yLimit,'--m',linewidth=1)
#   plt.text(1.6e-3,yLimit[0]-.3,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([relVeLong,relVeLong],yLimit,'--m',linewidth=1)
#   plt.text(4.4e-5,yLimit[0]+.1,'$ \Delta V_{e||}/ sV_{e0}$',color='m',fontsize=14)
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig5000.savefig('picturesCMA/frctnForce_AI_fig5010cma.png')    
      print ('File "picturesCMA/frctnForce_AI_fig5010cma.png" is written')

if (plotFigureFlag == 0):   
   fig5020=plt.figure (5020)
   plt.semilogx(VionRel,1.05*viewFctr*frctnForce_GKfit,'-xr', \
                VionRel,0.95*viewFctr*frctnForce_AI,'-xb',linewidth=2)
   plt.xlabel('Relative Ion Velocity  $V_{ion}/V_{e0}$',color='m',fontsize=14)
   if (mantViewrFctr == 1.):
      plt.ylabel(('$10^{%2d} \cdot F_{ion}$, eV/m' % powViewrFctr), color='m',fontsize=14)
   else:
      plt.ylabel(('$%3.1f \cdot 10^{%2d} \cdot F_{ion}$, eV/m' % \
                  (mantViewrFctr,powViewrFctr)), color='m',fontsize=14)
   plt.title('Comparison Two Methods to Calculate Friction Force $F_{ion}$', \
             color='m',fontsize=14)
#            ('$V_{e0}=%4.2f\cdot10^{%2d}$cm/s, $n_e=%4.2f\cdot10^{%2d}$cm$^3$, $B=%4d$kG' % \
#             (mantV0,powV0,mant_n_e,pow_n_e,int(fieldB[0]))),color='m',fontsize=14)
   plt.xlim([.9*VionRel[0],1.1*VionRel[nVion-1]])
   yLimit=[0.,max(1.1*viewFctr*frctnForce_GKfit)]
   plt.ylim(yLimit)
   plt.plot([relVeTrnsv,relVeTrnsv],yLimit,'--m',linewidth=1)
   plt.text(1.6e-3,.5,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([relVeLong,relVeLong],yLimit,'--m',linewidth=1)
   plt.text(4.4e-5,.5,'$ \Delta V_{e||}/ sV_{e0}$',color='m',fontsize=14)
   plt.legend(['GK-Integration (x 1.05)','Analytical Integration (x 0.95)'],loc='lower left',fontsize=11)
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig5020.savefig('picturesCMA_v7/frctnForce_GKfit_fig5020cma.png')    
      print ('File "picturesCMA_v7/frctnForce_GKfit_fig5020cma.png" is written')

if (plotFigureFlag == 0):   
   fig5030=plt.figure (5030)
   plt.semilogx(VionRel,np.log10(viewFctr*frctnForce_GKfit),'-r', \
                VionRel,np.log10(viewFctr*frctnForce_AI),'xb',linewidth=2)
   plt.xlabel('Relative Ion Velocity  $V_{ion}/V_{e0}$',color='m',fontsize=14)
   if (mantViewrFctr == 1.):
      plt.ylabel(('$log_{10}(10^{%2d} \cdot F_{ion})$' % powViewrFctr), color='m',fontsize=14)
   else:
      plt.ylabel(('$log_{10}(%3.1f \cdot 10^{%2d} \cdot F_{ion})$' % \
                  (mantViewrFctr,powViewrFctr)), color='m',fontsize=14)
   plt.title('Comparison Two Methods to Calculate Friction Force $F_{ion}$', \
             color='m',fontsize=14)
#            ('$V_{e0}=%4.2f\cdot10^{%2d}$cm/s, $n_e=%4.2f\cdot10^{%2d}$cm$^3$, $B=%4d$kG' % \
#             (mantV0,powV0,mant_n_e,pow_n_e,int(fieldB[0]))),color='m',fontsize=14)
   plt.xlim([.9*VionRel[0],1.1*VionRel[nVion-1]])
   yLimit=[min(np.log10(.9*viewFctr*frctnForce_GKfit)), \
           max(np.log10(1.1*viewFctr*frctnForce_GKfit))]
   plt.ylim(yLimit)
   plt.plot([relVeTrnsv,relVeTrnsv],yLimit,'--m',linewidth=1)
   plt.text(1.6e-3,yLimit[1]-.3,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([relVeLong,relVeLong],yLimit,'--m',linewidth=1)
   plt.text(4.4e-5,yLimit[1]-.3,'$ \Delta V_{e||}/ sV_{e0}$',color='m',fontsize=14)
   plt.legend(['GK-Integration','Analytical Integration'],loc='lower right',fontsize=11)
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig5030.savefig('picturesCMA_v7/frctnForce_GKfit_fig5030cma.png')    
      print ('File "picturesCMA_v7/frctnForce_GKfit_fig5030cma.png" is written')

yLimit = [-2.8,-0.35]
log10relVeTrnsv = np.log10(relVeTrnsv)
log10relVeLong = np.log10(relVeLong)

#
# Figure (5100) incorrect due to incorrect using of arrays X,Y.
# Correct figures are (6302) and (6301)
#
# if (plotFigureFlag == 0):   
#    fig5100=plt.figure (5100)
#    ax = fig5100.add_subplot(111)                    # for contours plotting
#    mapDenrgF = ax.contourf(X,Y,1.e-19*deltaEnrgIon_m,cmap='jet') 
#    mapDenrg = ax.contour(X,Y,1.e9*deltaEnrgIon_m,levels=range(0,10,1),colors='black') 
#    plt.clabel(mapDenrg,fmt='%3.1f',inline=True)
#    plt.xlabel('Relative Ion Velocity  $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
#    plt.ylabel('Initial Impact Parameter, $log_{10}(rho_{Init})$',color='m',fontsize=14)
#    titleHeader = 'Transferred Energy $\Delta E_{ion}$ per Ion:'
#    titleHeader += '\n$10^{-19} \cdot \Delta E_{ion}$, eV'
#    plt.title(titleHeader,color='m',fontsize=16)
#    plt.plot([log10relVeTrnsv,log10relVeTrnsv],yLimit,'--m',linewidth=1)
#    plt.text(-2.85,-0.95,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
#    plt.plot([log10relVeLong,log10relVeLong],yLimit,'--m',linewidth=1)
#    plt.text(-4.5,-0.95,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
#    plt.ylim(yLimit)
#    fig5100.colorbar(mapDenrgF)
#    plt.grid(True)
#    if (saveFilesFlag == 1):
#       fig5100.savefig('picturesCMA/mapEion_m_fig5100cma.png')    
#       print ('File "picturesCMA/mapEion_m_fig5100cma.png" is written')

#
# Additional output for checking:
#
# print ('\n                       rhoInit')
# for i in range(nVion):
#    print ('    Vion(%d) = %e' % (i,Vion[i]))
#    nn = 0
#    for n in range (5):
#       print ('%e %e %e %e %e %e %e %e %e %e ' % \
#              (rhoInit[nn,i],  rhoInit[nn+1,i],rhoInit[nn+2,i], \
# 	      rhoInit[nn+3,i],rhoInit[nn+4,i],rhoInit[nn+5,i], \
# 	      rhoInit[nn+6,i],rhoInit[nn+7,i],rhoInit[nn+8,i], \
# 	      rhoInit[nn+9,i]))
#       nn += 10


# print ('\n                       deltaPx_m')
# for i in range(nVion):
#    print ('    Vion(%d) = %e' % (i,Vion[i]))
#    nn = 0
#    for n in range (5):
#       print ('%e %e %e %e %e %e %e %e %e %e ' % \
#              (deltaPx_m[nn,i],  deltaPx_m[nn+1,i],deltaPx_m[nn+2,i], \
# 	      deltaPx_m[nn+3,i],deltaPx_m[nn+4,i],deltaPx_m[nn+5,i], \
# 	      deltaPx_m[nn+6,i],deltaPx_m[nn+7,i],deltaPx_m[nn+8,i], \
# 	      deltaPx_m[nn+9,i]))
#       nn += 10

# print ('\n                       deltaPx_m')
# for n in range(nImpctPrmtr):
#    print ('    rhoInit(%d) = %e' % (n,rhoInit[n,0]))
#    nn = 0
#    for i in range (5):
#       print ('%e %e %e %e %e %e %e %e %e %e ' % \
#              (deltaPx_m[n,nn],  deltaPx_m[n,nn+1],deltaPx_m[n,nn+2], \
# 	      deltaPx_m[n,nn+3],deltaPx_m[n,nn+4],deltaPx_m[n,nn+5], \
# 	      deltaPx_m[n,nn+6],deltaPx_m[n,nn+7],deltaPx_m[n,nn+8], \
# 	      deltaPx_m[n,nn+9]))
#       nn += 10

# print ('\n                       deltaPz_m')
# for i in range(nVion):
#    print ('    Vion(%d) = %e' % (i,Vion[i]))
#    nn = 0
#    for n in range (5):
#       print ('%e %e %e %e %e %e %e %e %e %e ' % \
#              (deltaPz_m[nn,i],  deltaPz_m[nn+1,i],deltaPz_m[nn+2,i], \
# 	      deltaPz_m[nn+3,i],deltaPz_m[nn+4,i],deltaPz_m[nn+5,i], \
# 	      deltaPz_m[nn+6,i],deltaPz_m[nn+7,i],deltaPz_m[nn+8,i], \
# 	      deltaPz_m[nn+9,i]))
#       nn += 10

# print ('\n                       deltaPz_m')
# for n in range(nImpctPrmtr):
#    print ('    rhoInit(%d) = %e' % (n,rhoInit[n,0]))
#    nn = 0
#    for i in range (5):
#       print ('%e %e %e %e %e %e %e %e %e %e ' % \
#              (deltaPz_m[n,nn],  deltaPz_m[n,nn+1],deltaPz_m[n,nn+2], \
# 	      deltaPz_m[n,nn+3],deltaPz_m[n,nn+4],deltaPz_m[n,nn+5], \
# 	      deltaPz_m[n,nn+6],deltaPz_m[n,nn+7],deltaPz_m[n,nn+8], \
# 	      deltaPz_m[n,nn+9]))
#       nn += 10

nVion_c = 50
nImpctPrmtr_c = 50
X_c = np.zeros((nImpctPrmtr_c,nVion_c))
Y_c = np.zeros((nImpctPrmtr_c,nVion_c))
log10deltaPx_m = np.zeros((nImpctPrmtr_c,nVion_c))
log10deltaPz_m = np.zeros((nImpctPrmtr_c,nVion_c))

for i in range(nVion_c):
   for n in range(nImpctPrmtr_c):
      X_c[n,i] = np.log10(VionRel[i])
      Y_c[n,i] = np.log10(rhoInit[n,i])
      log10deltaPx_m[n,i] = np.log10(1.e22*deltaPx_m[n,i])
      log10deltaPz_m[n,i] = np.log10(1.e24*deltaPz_m[n,i])
      
if (plotFigureFlag == 0):   
   fig5201=plt.figure (5201)
   ax = fig5201.add_subplot(111)                    # for contours plotting
   mapDpxF1 = ax.contourf(X_c,Y_c,1.e22*deltaPx_m,cmap='jet') 
#   mapDpx1 = ax.contour(X_c,Y_c,1.e22*deltaPx_m,levels=range(0,2,1),colors='black') 
   mapDpx1 = ax.contour(X_c,Y_c,1.e22*deltaPx_m,7,colors='black') 
   plt.clabel(mapDpx1,fmt='%4.2f',inline=True)
   plt.xlabel('Relative Ion Velocity  $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
   plt.ylabel('Initial Impact Parameter, $log_{10}(rho_{Init})$',color='m',fontsize=14)
   titleHeader = 'Transferred Momenta $\Delta p_x$ per Ion'
   titleHeader += '\n$10^{22} \cdot \Delta p_x$, g$\cdot$cm/s'
   plt.title(titleHeader,color='m',fontsize=16)
   plt.plot([log10relVeTrnsv,log10relVeTrnsv],yLimit,'--m',linewidth=1)
   plt.text(-2.85,-0.95,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([log10relVeLong,log10relVeLong],yLimit,'--m',linewidth=1)
   plt.text(-4.5,-0.95,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-4.14,-.6,'Screened Collisions',color='r',fontsize=16)
   plt.text(-4.14,-2.2,'Magnetized Collisions',color='k',fontsize=16)
   plt.ylim(yLimit)
   fig5201.colorbar(mapDpxF1)
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig5201.savefig('picturesCMA_v7/mapDeltaPx_m_fig5201cma.png')    
      print ('File "picturesCMA_v7/mapDeltaPx_m_fig5201cma.png" is written')

if (plotFigureFlag == 0):   
   fig5202=plt.figure (5202)
   ax = fig5202.add_subplot(111)                    # for contours plotting
   mapDpxF2 = ax.contourf(X_c,Y_c,log10deltaPx_m,cmap='jet') 
#   mapDpx2 = ax.contour(X_c,Y_c,log10deltaPx_m,levels=range(0,2,1),colors='black') 
   mapDpx2 = ax.contour(X_c,Y_c,log10deltaPx_m,7,colors='black') 
   plt.clabel(mapDpx2,fmt='%4.2f',inline=True)
   plt.xlabel('Relative Ion Velocity  $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
   plt.ylabel('Initial Impact Parameter, $log_{10}(rho_{Init})$',color='m',fontsize=14)
   titleHeader = 'Transferred Momenta $\Delta p_x$ per Ion'
   titleHeader += '\n$log_{10}(10^{22} \cdot \Delta p_x)$'
   plt.title(titleHeader,color='m',fontsize=16)
   plt.plot([log10relVeTrnsv,log10relVeTrnsv],yLimit,'--m',linewidth=1)
   plt.text(-2.85,-0.95,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([log10relVeLong,log10relVeLong],yLimit,'--m',linewidth=1)
   plt.text(-4.5,-0.95,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-4.14,-.6,'Screened Collisions',color='r',fontsize=16)
   plt.text(-4.14,-2.2,'Magnetized Collisions',color='k',fontsize=16)
   plt.ylim(yLimit)
   fig5202.colorbar(mapDpxF2)
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig5202.savefig('picturesCMA_v7/mapLog10deltaPx_m_fig5202cma.png')    
      print ('File "picturesCMA_v7/mapLog10deltaPx_m_fig5202cma.png" is written')

if (plotFigureFlag == 0):   
   fig5401=plt.figure (5401)
   ax = fig5401.add_subplot(111)                    # for contours plotting
   mapDpzF1 = ax.contourf(X_c,Y_c,1.e24*deltaPz_m,cmap='jet') 
#   mapDpz1 = ax.contour(X_c,Y_c,1.e24*deltaPz_m,levels=range(0,5,1),colors='black') 
   mapDpz1 = ax.contour(X_c,Y_c,1.e24*deltaPz_m,7,colors='black') 
   plt.clabel(mapDpz1,fmt='%4.2f',inline=True)
   plt.xlabel('Relative Ion Velocity  $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
   plt.ylabel('Initial Impact Parameter, $log_{10}(rho_{Init})$',color='m',fontsize=14)
   titleHeader = 'Transferred Momenta $\Delta p_z$ per Ion'
   titleHeader += '\n$10^{24} \cdot \Delta p_z$, g$\cdot$cm/s'
   plt.title(titleHeader,color='m',fontsize=16)
   plt.plot([log10relVeTrnsv,log10relVeTrnsv],yLimit,'--m',linewidth=1)
   plt.text(-2.85,-0.95,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([log10relVeLong,log10relVeLong],yLimit,'--m',linewidth=1)
   plt.text(-4.5,-0.95,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-4.14,-.6,'Screened Collisions',color='r',fontsize=16)
   plt.text(-4.14,-2.2,'Magnetized Collisions',color='w',fontsize=16)
   plt.ylim(yLimit)
   fig5401.colorbar(mapDpzF1)
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig5401.savefig('picturesCMA_v7/mapDeltaPz_m_fig5401cma.png')    
      print ('File "picturesCMA_v7/mapDeltaPz_m_fig5401cma.png" is written')

if (plotFigureFlag == 0):
   fig5402=plt.figure (5402)
   ax = fig5402.add_subplot(111)                    # for contours plotting
   mapDpzF2 = ax.contourf(X_c,Y_c,log10deltaPz_m,cmap='jet') 
#   mapDpz2 = ax.contour(X_c,Y_c,log10deltaPz_m,levels=range(0,5,1),colors='black') 
   mapDpz2 = ax.contour(X_c,Y_c,log10deltaPz_m,7,colors='black') 
   plt.clabel(mapDpz2,fmt='%4.2f',inline=True)
   plt.xlabel('Relative Ion Velocity  $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
   plt.ylabel('Initial Impact Parameter, $log_{10}(rho_{Init})$',color='m',fontsize=14)
   titleHeader = 'Transferred Momenta $\Delta p_z$ per Ion'
   titleHeader += '\n$log_{10}(10^{24} \cdot \Delta p_z)$'
   plt.title(titleHeader,color='m',fontsize=16)
   plt.plot([log10relVeTrnsv,log10relVeTrnsv],yLimit,'--m',linewidth=1)
   plt.text(-2.85,-0.95,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([log10relVeLong,log10relVeLong],yLimit,'--m',linewidth=1)
   plt.text(-4.5,-0.95,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-4.14,-.6,'Screened Collisions',color='r',fontsize=16)
   plt.text(-4.14,-2.2,'Magnetized Collisions',color='k',fontsize=16)
   plt.ylim(yLimit)
   fig5402.colorbar(mapDpzF2)
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig5402.savefig('picturesCMA_v7/mapLog10deltaPz_m_fig5402cma.png')    
      print ('File "picturesCMA_v7/mapLog10deltaPz_m_fig5402cma.png" is written')

if (plotFigureFlag == 0):   
   fig5250=plt.figure (5250)
   ax = fig5250.add_subplot(111)                    # for contours plotting
   mapVixm = ax.contourf(X_c,Y_c,1.e-4*ionVx_m,cmap='jet') 
#   mapVix = ax.contour(X_c,Y_c,1.e-4*ionVx_m,levels=range(0,2,1),colors='black') 
   mapVix = ax.contour(X_c,Y_c,1.e-4*ionVx_m,7,colors='black') 
   plt.clabel(mapVix,fmt='%4.2f',inline=True)
   plt.xlabel('Relative Ion Velocity,  $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
   plt.ylabel('Initial Impact Parameter $log_{10}(rho_{Init})$',color='m',fontsize=14)
   titleHeader = 'x-Component of Ion Velocity'
   titleHeader += '\n$10^{-4} \cdot Vx_{ion}$, cm/s'
   plt.title(titleHeader,color='m',fontsize=16)
   plt.plot([log10relVeTrnsv,log10relVeTrnsv],yLimit,'--m',linewidth=1)
#   plt.text(-2.85,-0.95,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([log10relVeLong,log10relVeLong],yLimit,'--m',linewidth=1)
#   plt.text(-4.5,-0.95,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-2.85,-0.95,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-4.5,-0.95,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-4.14,-.6,'Screened Collisions',color='r',fontsize=16)
   plt.text(-4.14,-2.2,'Magnetized Collisions',color='w',fontsize=16)
   plt.ylim(yLimit)
   fig5250.colorbar(mapVixm)
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig5250.savefig('picturesCMA_v7/mapVix_m_fig5250cma.png')    
      print ('File "picturesCMA_v7/mapVix_m_fig5250cma.png" is written')

if (plotFigureFlag == 0):   
   fig5350=plt.figure (5350)
   ax = fig5350.add_subplot(111)                    # for contours plotting
   mapViym = ax.contourf(X_c,Y_c,ionVy_m,cmap='jet') 
#   mapVyx = ax.contour(X_c,Y_c,ionVy_m,levels=range(0,2,1),colors='black') 
   mapVyx = ax.contour(X_c,Y_c,ionVy_m,7,colors='black') 
   plt.clabel(mapVyx,fmt='%4.2f',inline=True)
   plt.xlabel('Relative Ion Velocity,  $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
   plt.ylabel('Initial Impact Parameter $log_{10}(rho_{Init})$',color='m',fontsize=14)
   titleHeader = 'y-Component of Ion Velocity'
   titleHeader += '\n$Vy_{ion}$, cm/s'
   plt.title(titleHeader,color='m',fontsize=16)
   plt.plot([log10relVeTrnsv,log10relVeTrnsv],yLimit,'--m',linewidth=1)
#   plt.text(-2.85,-0.95,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([log10relVeLong,log10relVeLong],yLimit,'--m',linewidth=1)
#   plt.text(-4.5,-0.95,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-2.85,-0.95,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-4.5,-0.95,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-4.14,-.6,'Screened Collisions',color='r',fontsize=16)
   plt.text(-4.14,-2.2,'Magnetized Collisions',color='w',fontsize=16)
   plt.ylim(yLimit)
   fig5350.colorbar(mapViym)
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig5350.savefig('picturesCMA_v7/mapViy_m_fig5350cma.png')    
      print ('File "picturesCMA_v7/mapViy_m_fig5350cma.png" is written')

if (plotFigureFlag == 0):   
   fig5450=plt.figure (5450)
   ax = fig5450.add_subplot(111)                    # for contours plotting
   mapVizm = ax.contourf(X_c,Y_c,1.e-8*ionVz_m,cmap='jet') 
#   mapViz = ax.contour(X_c,Y_c,1.e-8*ionVz_m,levels=range(0,2,1),colors='black') 
   mapViz = ax.contour(X_c,Y_c,1.e-8*ionVz_m,7,colors='black') 
   plt.clabel(mapViz,fmt='%4.2f',inline=True)
   plt.xlabel('Relative Ion Velocity,  $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
   plt.ylabel('Initial Impact Parameter $log_{10}(rho_{Init})$',color='m',fontsize=14)
   titleHeader = 'z-Component of Ion Velocity'
   titleHeader += '\n$10^{-8} \cdot Vz_{ion}$, cm/s'
   plt.title(titleHeader,color='m',fontsize=16)
   plt.plot([log10relVeTrnsv,log10relVeTrnsv],yLimit,'--m',linewidth=1)
#   plt.text(-2.85,-0.95,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([log10relVeLong,log10relVeLong],yLimit,'--m',linewidth=1)
#   plt.text(-4.5,-0.95,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-2.85,-0.95,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-4.5,-0.95,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-4.14,-.6,'Screened Collisions',color='r',fontsize=16)
   plt.text(-4.14,-2.2,'Magnetized Collisions',color='w',fontsize=16)
   plt.ylim(yLimit)
   fig5450.colorbar(mapVizm)
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig5450.savefig('picturesCMA_v7/mapViz_m_fig5450cma.png')    
      print ('File "picturesCMA_v7/mapViz_m_fig5450cma.png" is written')

log10ionVx_m = np.zeros((nImpctPrmtr_c,nVion_c))
log10ionVy_m = np.zeros((nImpctPrmtr_c,nVion_c))
log10ionVz_m = np.zeros((nImpctPrmtr_c,nVion_c))

for i in range(nVion_c):
   for n in range(nImpctPrmtr_c):
      log10ionVx_m[n,i] = np.log10(1.e-4*ionVx_m[n,i])
      log10ionVy_m[n,i] = np.log10(ionVy_m[n,i])
      log10ionVz_m[n,i] = np.log10(1.e-8*ionVz_m[n,i])

if (plotFigureFlag == 0):   
   fig5251=plt.figure (5251)
   ax = fig5251.add_subplot(111)                    # for contours plotting
   mapVixm1 = ax.contourf(X_c,Y_c,log10ionVx_m,cmap='jet') 
#   mapVix1 = ax.contour(X_c,Y_c,log10ionVx_m,levels=range(0,2,1),colors='black') 
   mapVix1 = ax.contour(X_c,Y_c,log10ionVx_m,7,colors='black') 
   plt.clabel(mapVix1,fmt='%4.2f',inline=True)
   plt.xlabel('Relative Ion Velocity,  $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
   plt.ylabel('Initial Impact Parameter $log_{10}(rho_{Init})$',color='m',fontsize=14)
   titleHeader = 'x-Component of Ion Velocity'
   titleHeader += '\n$log_{10}(10^{-4} \cdot Vx_{ion})$'
   plt.title(titleHeader,color='m',fontsize=16)
   plt.plot([log10relVeTrnsv,log10relVeTrnsv],yLimit,'--m',linewidth=1)
#   plt.text(-2.85,-0.95,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([log10relVeLong,log10relVeLong],yLimit,'--m',linewidth=1)
#   plt.text(-4.5,-0.95,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-2.85,-0.95,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-4.5,-0.95,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-4.14,-.6,'Screened Collisions',color='r',fontsize=16)
   plt.text(-4.14,-2.2,'Magnetized Collisions',color='k',fontsize=16)
   plt.ylim(yLimit)
   fig5251.colorbar(mapVixm1)
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig5251.savefig('picturesCMA_v7/mapLog10Vix_m_fig5251cma.png')    
      print ('File "picturesCMA_v7/mapLog10Vix_m_fig5251cma.png" is written')

if (plotFigureFlag == 0):   
   fig5351=plt.figure (5351)
   ax = fig5351.add_subplot(111)                    # for contours plotting
   mapViym1 = ax.contourf(X_c,Y_c,log10ionVy_m,cmap='jet') 
#   mapViy1 = ax.contour(X_c,Y_c,log10ionVy_m,levels=range(0,2,1),colors='black') 
   mapViy1 = ax.contour(X_c,Y_c,log10ionVy_m,7,colors='black') 
   plt.clabel(mapViy1,fmt='%4.2f',inline=True)
   plt.xlabel('Relative Ion Velocity,  $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
   plt.ylabel('Initial Impact Parameter $log_{10}(rho_{Init})$',color='m',fontsize=14)
   titleHeader = 'y-Component of Ion Velocity'
   titleHeader += '\n$log_{10}(Vy_{ion})$'
   plt.title(titleHeader,color='m',fontsize=16)
   plt.plot([log10relVeTrnsv,log10relVeTrnsv],yLimit,'--m',linewidth=1)
#   plt.text(-2.85,-0.95,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([log10relVeLong,log10relVeLong],yLimit,'--m',linewidth=1)
#   plt.text(-4.5,-0.95,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-2.85,-0.95,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-4.5,-0.95,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-4.14,-.6,'Screened Collisions',color='r',fontsize=16)
   plt.text(-4.14,-2.2,'Magnetized Collisions',color='k',fontsize=16)
   plt.ylim(yLimit)
   fig5351.colorbar(mapViym1)
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig5351.savefig('picturesCMA_v7/mapLog10Viy_m_fig5351cma.png')    
      print ('File "picturesCMA_v7/mapLog10Viy_m_fig5351cma.png" is written')

if (plotFigureFlag == 0):   
   fig5451=plt.figure (5451)
   ax = fig5451.add_subplot(111)                    # for contours plotting
   mapVizm1 = ax.contourf(X_c,Y_c,log10ionVz_m,cmap='jet') 
#   mapViz1 = ax.contour(X_c,Y_c,log10ionVz_m,levels=range(0,2,1),colors='black') 
   mapViz1 = ax.contour(X_c,Y_c,log10ionVz_m,7,colors='black') 
   plt.clabel(mapViz1,fmt='%4.2f',inline=True)
   plt.xlabel('Relative Ion Velocity,  $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
   plt.ylabel('Initial Impact Parameter $log_{10}(rho_{Init})$',color='m',fontsize=14)
   titleHeader = 'z-Component of Ion Velocity'
   titleHeader += '\n$log_{10}(10^{-8} \cdot Vz_{ion})$'
   plt.title(titleHeader,color='m',fontsize=16)
   plt.plot([log10relVeTrnsv,log10relVeTrnsv],yLimit,'--m',linewidth=1)
#   plt.text(-2.85,-0.95,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([log10relVeLong,log10relVeLong],yLimit,'--m',linewidth=1)
#   plt.text(-4.5,-0.95,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-2.85,-0.95,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-4.5,-0.95,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-4.14,-.6,'Screened Collisions',color='r',fontsize=16)
   plt.text(-4.14,-2.2,'Magnetized Collisions',color='k',fontsize=16)
   plt.ylim(yLimit)
   fig5451.colorbar(mapVizm1)
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig5451.savefig('picturesCMA_v7/mapLog10Viz_m_fig5451cma.png')    
      print ('File "picturesCMA_v7/mapLog10Viz_m_fig5451cma.png" is written')

#
# Some checkings:
#

# for i in range(nVion):
#    print ('\n                   i=%d' % i)
#    nn = 0
#    for n in range(5):
#       print ('%e %e %e %e %e %e %e %e %e %e ' % \
#              (ionVz_m[nn,i],  ionVz_m[nn+1,i],ionVz_m[nn+2,i], \
# 	      ionVz_m[nn+3,i],ionVz_m[nn+4,i],ionVz_m[nn+5,i], \
# 	      ionVz_m[nn+6,i],ionVz_m[nn+7,i],ionVz_m[nn+8,i], \
# 	      ionVz_m[nn+9,i]))
#       nn += 10

ionVx_dPx_m = np.zeros((nImpctPrmtr,nVion))
ionVy_dPy_m = np.zeros((nImpctPrmtr,nVion))
ionVz_dPz_m = np.zeros((nImpctPrmtr,nVion))
for i in range(nVion):
   for n in range(nImpctPrmtr):
      ionVx_dPx_m[n,i] = ionVx_m[n,i]*deltaPx_m[n,i]
      ionVy_dPy_m[n,i] = ionVy_m[n,i]*deltaPy_m[n,i]
      ionVz_dPz_m[n,i] = ionVz_m[n,i]*deltaPz_m[n,i]

yLimit = [-2.9,-0.3]

#--------------------------------------------
# 6000,6001 - wrong pictures and bad idea to use the xLabel and yLabel approach 
#             to better recognize the dependence of a function on parameters!!!
#      
# if (plotFigureFlag == 0):   
#    fig6000=plt.figure (6000)
#    ax = fig6000.add_subplot(111)                    # for contours plotting
#    mapVixPxm = ax.contourf(X,Y,1.e18*ionVx_dPx_m,cmap='jet') 
# #   mapVixPx = ax.contour(X,Y,1.e18*ionVx_dPx_m,levels=range(0,2,1),colors='black') 
#    mapVixPx = ax.contour(X,Y,1.e18*ionVx_dPx_m,7,colors='black') 
#    plt.clabel(mapVixPx,fmt='%4.2f',inline=True)
#    plt.xlabel('Relative Ion Velocity,  $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
#    plt.ylabel('Initial Impact Parameter $log_{10}(rho_{Init})$',color='m',fontsize=14)
#    titleHeader = '$10^{18}\cdot Vz_{ion}\cdot \Delta P_z$,  g$\cdot $cm$^2$/s$^2$'
#    plt.title(titleHeader,color='m',fontsize=12)
# #   titleHeader = 'z-Component of Ion Velocity $Vz_{ion}$, cm/s'
# #   titleHeader += '\n$\widetilde{Vz_{ion}}$ = $10^{-8} \cdot Vz_{ion}$ '
# #   titleHeader += '($\widetilde{Vz_{ion}}$ = 0.5, if $\widetilde{V_{ion}}$ > 0.5)'
# #   plt.title(titleHeader,color='m',fontsize=12)
#    plt.plot([log10relVeTrnsv,log10relVeTrnsv],yLimit,'--m',linewidth=1)
#    plt.text(-2.85,-0.8,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
#    plt.plot([log10relVeLong,log10relVeLong],yLimit,'--m',linewidth=1)
#    plt.text(-4.5,-0.8,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
#    plt.ylim(yLimit)
#    fig6000.colorbar(mapVixPxm)
#    plt.grid(True)
# if (saveFilesFlag == 1):
#    fig6000.savefig('picturesCMA/mapVix_m_fig6000cma.png')    
#    print ('File "picturesCMA/mapVix_m_fig6000cma.png" is written')
# 
# if (plotFigureFlag == 0):   
#    fig6001=plt.figure (6001)
#    ax = fig6001.add_subplot(111)                    # for contours plotting
#    mapVixPxm_1 = ax.contourf(X,Y,1.e18*ionVx_dPx_m,cmap='jet') 
# #   mapVixPx_1 = ax.contour(X,Y,1.e18*ionVx_dPx_m,levels=range(0,2,1),colors='black') 
#    mapVixPx_1 = ax.contour(X,Y,1.e18*ionVx_dPx_m,7,colors='black') 
#    plt.clabel(mapVixPx_1,fmt='%4.2f',inline=True)
#    plt.xlabel('Relative Ion Velocity,  $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
#    plt.ylabel('Initial Impact Parameter $log_{10}(rho_{Init})$',color='m',fontsize=14)
#    titleHeader = '$10^{18}\cdot Vz_{ion}\cdot \Delta P_z$,  g$\cdot $cm$^2$/s$^2$'
#    plt.title(titleHeader,color='m',fontsize=12)
# #   titleHeader = 'z-Component of Ion Velocity $Vz_{ion}$, cm/s'
# #   titleHeader += '\n$\widetilde{Vz_{ion}}$ = $10^{-8} \cdot Vz_{ion}$ '
# #   titleHeader += '($\widetilde{Vz_{ion}}$ = 0.5, if $\widetilde{V_{ion}}$ > 0.5)'
# #   plt.title(titleHeader,color='m',fontsize=12)
# #    plt.plot([log10relVeTrnsv,log10relVeTrnsv],yLimit,'--m',linewidth=1)
# #    plt.text(-2.85,-0.8,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
# #    plt.plot([log10relVeLong,log10relVeLong],yLimit,'--m',linewidth=1)
# #    plt.text(-4.5,-0.8,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
#    plt.xlim([-4.75,-4.])
#    plt.ylim([-2.5,-1.7,])
#    fig6001.colorbar(mapVixPxm_1)
#    plt.grid(True)
# if (saveFilesFlag == 1):
#    fig6001.savefig('picturesCMA/mapVix_m_fig6001cma.png')    
#    print ('File "picturesCMA/mapVix_m_fig6001cma.png" is written')
#--------------------------------------------

if (plotFigureFlag == 0):   
   fig6002=plt.figure (6002)
   ax = fig6002.add_subplot(111)                    # for contours plotting
   mapVixPxm2 = ax.contourf(X_c,Y_c,1.e18*ionVx_dPx_m,cmap='jet') 
#   mapVixPx2 = ax.contour(X_c,Y_c,1.e18*ionVx_dPx_m,levels=range(0,2,1),colors='black') 
   mapVixPx2 = ax.contour(X_c,Y_c,1.e18*ionVx_dPx_m,7,colors='black') 
   plt.clabel(mapVixPx2,fmt='%4.2f',inline=True)
   plt.xlabel('Relative Ion Velocity,  $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
   plt.ylabel('Initial Impact Parameter $log_{10}(rho_{Init})$',color='m',fontsize=14)
   titleHeader = '$10^{18}\cdot Vx_{ion}\cdot \Delta P_x$,  g$\cdot $cm$^2$/s$^2$'
   plt.title(titleHeader,color='m',fontsize=12)
#   titleHeader = 'z-Component of Ion Velocity $Vz_{ion}$, cm/s'
#   titleHeader += '\n$\widetilde{Vz_{ion}}$ = $10^{-8} \cdot Vz_{ion}$ '
#   titleHeader += '($\widetilde{Vz_{ion}}$ = 0.5, if $\widetilde{V_{ion}}$ > 0.5)'
#   plt.title(titleHeader,color='m',fontsize=12)
   plt.plot([log10relVeTrnsv,log10relVeTrnsv],yLimit,'--m',linewidth=1)
   plt.text(-2.85,-0.8,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([log10relVeLong,log10relVeLong],yLimit,'--m',linewidth=1)
   plt.text(-4.5,-0.8,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-4.14,-.6,'Screened Collisions',color='r',fontsize=16)
   plt.text(-4.14,-2.2,'Magnetized Collisions',color='w',fontsize=16)
   plt.ylim(yLimit)
   fig6002.colorbar(mapVixPxm2)
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig6002.savefig('picturesCMA_v7/mapVix_dPx_m_fig6002cma.png')    
      print ('File "picturesCMA_v7/mapVix_dPx_m_fig6002cma.png" is written')


yLimit = [-2.9,-0.3]
#--------------------------------------------
# 6100,6101 - wrong pictures and bad idea to use xLabel and yLabel approach 
#             to better recognize the dependence of a function on parameters!!!
#      
# if (plotFigureFlag == 0):   
#    fig6100=plt.figure (6100)
#    ax = fig6100.add_subplot(111)                    # for contours plotting
#    mapVizPzm = ax.contourf(X,Y,1.e16*ionVz_dPz_m,cmap='jet') 
# #   mapVizPz = ax.contour(X,Y,1.e16*ionVz_dPz_m,levels=range(0,2,1),colors='black') 
#    mapVizPz = ax.contour(X,Y,1.e16*ionVz_dPz_m,7,colors='black') 
#    plt.clabel(mapVizPz,fmt='%4.2f',inline=True)
#    plt.xlabel('Relative Ion Velocity,  $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
#    plt.ylabel('Initial Impact Parameter $log_{10}(rho_{Init})$',color='m',fontsize=14)
#    titleHeader = '$10^{16}\cdot Vx_{ion}\cdot \Delta P_x$,  g$\cdot $cm$^2$/s$^2$'
#    plt.title(titleHeader,color='m',fontsize=12)
#    plt.plot([log10relVeTrnsv,log10relVeTrnsv],yLimit,'--m',linewidth=1)
#    plt.text(-2.85,-0.8,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
#    plt.plot([log10relVeLong,log10relVeLong],yLimit,'--m',linewidth=1)
#    plt.text(-4.5,-0.8,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
#    plt.ylim(yLimit)
#    fig6100.colorbar(mapVizPzm)
#    plt.grid(True)
# if (saveFilesFlag == 1):
#    fig6100.savefig('picturesCMA/mapVix_m_fig6100cma.png')    
#    print ('File "picturesCMA/mapVix_m_fig6100cma.png" is written')
# 
# if (plotFigureFlag == 0):   
#    fig6101=plt.figure (6101)
#    ax = fig6101.add_subplot(111)                    # for contours plotting
#    mapVizPzm_1 = ax.contourf(X,Y,1.e16*ionVz_dPz_m,cmap='jet') 
# #   mapVizPz_1 = ax.contour(X,Y,1.e16*ionVz_dPz_m,levels=range(0,2,1),colors='black') 
#    mapVizPz_1 = ax.contour(X,Y,1.e16*ionVz_dPz_m,7,colors='black') 
#    plt.clabel(mapVizPz_1,fmt='%4.2f',inline=True)
#    plt.xlabel('Relative Ion Velocity,  $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
#    plt.ylabel('Initial Impact Parameter $log_{10}(rho_{Init})$',color='m',fontsize=14)
#    titleHeader = '$10^{16}\cdot Vx_{ion}\cdot \Delta P_x$,  g$\cdot $cm$^2$/s$^2$'
#    plt.title(titleHeader,color='m',fontsize=12)
#    plt.plot([log10relVeTrnsv,log10relVeTrnsv],yLimit,'--m',linewidth=1)
#    plt.text(-2.85,-0.8,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
#    plt.plot([log10relVeLong,log10relVeLong],yLimit,'--m',linewidth=1)
#    plt.text(-4.5,-0.8,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
#    plt.ylim(yLimit)
#    plt.xlim([-4.75,-4.25])
#    plt.ylim([-2.25,-1.75])
#    fig6101.colorbar(mapVizPzm_1)
#    plt.grid(True)
# if (saveFilesFlag == 1):
#    fig6101.savefig('picturesCMA/mapVix_m_fig6101cma.png')    
#    print ('File "picturesCMA/mapVix_m_fig6101cma.png" is written')
#--------------------------------------------

if (plotFigureFlag == 0):   
   fig6102=plt.figure (6102)
   ax = fig6102.add_subplot(111)                    # for contours plotting
   mapVizPzm2 = ax.contourf(X_c,Y_c,1.e16*ionVz_dPz_m,cmap='jet') 
#   mapVizPz2 = ax.contour(X_c,Y_c,1.e16*ionVz_dPz_m,levels=range(0,2,1),colors='black') 
   mapVizPz2 = ax.contour(X_c,Y_c,1.e16*ionVz_dPz_m,7,colors='black') 
   plt.clabel(mapVizPz2,fmt='%4.2f',inline=True)
   plt.xlabel('Relative Ion Velocity,  $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
   plt.ylabel('Initial Impact Parameter $log_{10}(rho_{Init})$',color='m',fontsize=14)
   titleHeader = 'Transferred Value: $10^{16}\cdot Vz_{ion}\cdot \Delta P_z$,  g$\cdot $cm$^2$/s$^2$'
   plt.title(titleHeader,color='m',fontsize=12)
   plt.plot([log10relVeTrnsv,log10relVeTrnsv],yLimit,'--m',linewidth=1)
   plt.text(-2.85,-0.8,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([log10relVeLong,log10relVeLong],yLimit,'--m',linewidth=1)
   plt.text(-4.5,-0.8,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-4.14,-.6,'Screened Collisions',color='r',fontsize=16)
   plt.text(-4.14,-2.2,'Magnetized Collisions',color='w',fontsize=16)
   plt.ylim(yLimit)
   fig6102.colorbar(mapVizPzm2)
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig6102.savefig('picturesCMA_v7/mapViz_dPz_m_fig6102cma.png')    
      print ('File "picturesCMA_v7/mapViz_dPz_m_fig6102cma.png" is written')

if (plotFigureFlag == 0):   
   fig6202=plt.figure (6202)
   ax = fig6202.add_subplot(111)                    # for contours plotting
   mapViyPym2 = ax.contourf(X_c,Y_c,1.e24*ionVy_dPy_m,cmap='jet') 
#   mapViyPy2 = ax.contour(X_c,Y_c,1.e24*ionVy_dPy_m,levels=range(0,2,1),colors='black') 
   mapViyPy2 = ax.contour(X_c,Y_c,1.e24*ionVy_dPy_m,7,colors='black') 
   plt.clabel(mapViyPy2,fmt='%4.2f',inline=True)
   plt.xlabel('Relative Ion Velocity,  $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
   plt.ylabel('Initial Impact Parameter $log_{10}(rho_{Init})$',color='m',fontsize=14)
   titleHeader = 'Transferred Value: $10^{24}\cdot Vy_{ion}\cdot \Delta P_y$,  g$\cdot $cm$^2$/s$^2$'
   plt.title(titleHeader,color='m',fontsize=12)
   plt.plot([log10relVeTrnsv,log10relVeTrnsv],yLimit,'--m',linewidth=1)
   plt.text(-2.85,-0.8,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([log10relVeLong,log10relVeLong],yLimit,'--m',linewidth=1)
   plt.text(-4.5,-0.8,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-4.14,-.6,'Screened Collisions',color='r',fontsize=16)
   plt.text(-4.14,-2.2,'Magnetized Collisions',color='w',fontsize=16)
   plt.ylim(yLimit)
   fig6202.colorbar(mapViyPym2)
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig6202.savefig('picturesCMA_v7/mapViy_dPy_m_fig6202cma.png')    
      print ('File "picturesCMA_v7/mapViy_dPy_m_fig6202cma.png" is written')

if (plotFigureFlag == 0):   
   fig6302=plt.figure (6302)
   ax = fig6302.add_subplot(111)                    # for contours plotting
   mapEnrgIon_m = ax.contourf(X_c,Y_c,1.e-19*deltaEnrgIon_m,cmap='jet') 
#   mapEnrgIon = ax.contour(X_c,Y_c,1.e-18*deltaEnrgIon_m,levels=range(0,2,1),colors='black') 
   mapEnrgIon = ax.contour(X_c,Y_c,1.e-18*deltaEnrgIon_m,7,colors='black') 
   plt.clabel(mapEnrgIon,fmt='%4.2f',inline=True)
   plt.xlabel('Relative Ion Velocity,  $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
   plt.ylabel('Initial Impact Parameter $log_{10}(rho_{Init})$',color='m',fontsize=14)
   titleHeader = 'Transferred Energy: $10^{-19}\cdot \Delta E_{ion}$, eV'
   plt.title(titleHeader,color='m',fontsize=12)
   plt.plot([log10relVeTrnsv,log10relVeTrnsv],yLimit,'--m',linewidth=1)
   plt.text(-2.85,-0.8,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([log10relVeLong,log10relVeLong],yLimit,'--m',linewidth=1)
   plt.text(-4.5,-0.8,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-4.14,-.6,'Screened Collisions',color='r',fontsize=16)
   plt.text(-4.14,-2.2,'Magnetized Collisions',color='w',fontsize=16)
   plt.ylim(yLimit)
   fig6302.colorbar(mapEnrgIon_m)
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig6302.savefig('picturesCMA_v7/mapEnrgIon_m_fig6302cma.png')    
      print ('File "picturesCMA_v7/mapEnrgIon_m_fig6302cma.png" is written')

#---------------------------------------
# Figure 6402 see after Figure 6400
#---------------------------------------

log10ionVx_dPx_m = np.zeros((nImpctPrmtr_c,nVion_c))
log10ionVy_dPy_m = np.zeros((nImpctPrmtr_c,nVion_c))
log10ionVz_dPz_m = np.zeros((nImpctPrmtr_c,nVion_c))
log10deltaEnrgIon_m = np.zeros((nImpctPrmtr_c,nVion_c))

for i in range(nVion_c):
   for n in range(nImpctPrmtr_c):
      log10ionVx_dPx_m[n,i] = np.log10(1.e18*ionVx_dPx_m[n,i])
      log10ionVy_dPy_m[n,i] = np.log10(1.e24*ionVy_dPy_m[n,i])
      log10ionVz_dPz_m[n,i] = np.log10(1.e16*ionVz_dPz_m[n,i])
      log10deltaEnrgIon_m[n,i] = np.log10(1.e-19*deltaEnrgIon_m[n,i])

if (plotFigureFlag == 0):   
   fig6000=plt.figure (6000)
   ax = fig6000.add_subplot(111)                    # for contours plotting
   mapVixPxm_c = ax.contourf(X_c,Y_c,log10ionVx_dPx_m[0:nImpctPrmtr_c,0:nVion_c],cmap='jet') 
#   mapVixPx_c = ax.contour(X_c,Y_c,log10ionVx_dPx_m[0:nImpctPrmtr_c,0:nVion_c],levels=range(0,2,1),colors='black') 
   mapVixPx_c = ax.contour(X_c,Y_c,log10ionVx_dPx_m[0:nImpctPrmtr_c,0:nVion_c],7,colors='black') 
   plt.clabel(mapVixPx_c,fmt='%4.2f',inline=True)
   plt.xlabel('Relative Ion Velocity,  $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
   plt.ylabel('Initial Impact Parameter $log_{10}(rho_{Init})$',color='m',fontsize=14)
   titleHeader = '$log_{10}(10^{18}\cdot Vx_{ion}\cdot \Delta P_x)$'
   plt.title(titleHeader,color='m',fontsize=12)
   plt.plot([log10relVeTrnsv,log10relVeTrnsv],yLimit,'--m',linewidth=1)
   plt.text(-2.85,-0.8,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([log10relVeLong,log10relVeLong],yLimit,'--m',linewidth=1)
   plt.text(-4.5,-0.8,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-4.14,-.6,'Screened Collisions',color='r',fontsize=16)
   plt.text(-4.14,-2.2,'Magnetized Collisions',color='k',fontsize=16)
# y   plt.xlim([-4.75,-3.25])
# y   plt.ylim([-2.75,-2.0])
   plt.ylim(yLimit)
   fig6000.colorbar(mapVixPxm_c)
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig6000.savefig('picturesCMA_v7/mapLog10Vix_dPx_m_fig6000cma.png')    
      print ('File "picturesCMA_v7/mapLog10Vix_dPx_m_fig6000cma.png" is written')

# yLimit = [-2.9,-0.3]
if (plotFigureFlag == 0):   
   fig6100=plt.figure (6100)
   ax = fig6100.add_subplot(111)                    # for contours plotting
   mapVizPzm_c = ax.contourf(X_c,Y_c,log10ionVz_dPz_m[0:nImpctPrmtr_c,0:nVion_c],cmap='jet') 
#   mapVizPz_c = ax.contour(X_c,Y_c,log10ionVz_dPz_m[0:nImpctPrmtr_c,0:nVion_c],levels=range(0,2,1),colors='black') 
   mapVizPz_c = ax.contour(X_c,Y_c,log10ionVz_dPz_m[0:nImpctPrmtr_c,0:nVion_c],7,colors='black') 
   plt.clabel(mapVizPz_c,fmt='%4.2f',inline=True)
   plt.xlabel('Relative Ion Velocity,  $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
   plt.ylabel('Initial Impact Parameter $log_{10}(rho_{Init})$',color='m',fontsize=14)
   titleHeader = 'Transferred Value: $log_{10}(10^{16}\cdot Vz_{ion}\cdot \Delta P_z)$'
   plt.title(titleHeader,color='m',fontsize=12)
   plt.plot([log10relVeTrnsv,log10relVeTrnsv],yLimit,'--m',linewidth=1)
   plt.text(-2.85,-0.8,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([log10relVeLong,log10relVeLong],yLimit,'--m',linewidth=1)
   plt.text(-4.5,-0.8,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-4.14,-.6,'Screened Collisions',color='r',fontsize=16)
   plt.text(-4.14,-2.2,'Magnetized Collisions',color='k',fontsize=16)
# y   plt.xlim([-4.75,-3.25])
# y   plt.ylim([-2.75,-2.0])
   plt.ylim(yLimit)
   fig6100.colorbar(mapVizPzm_c)
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig6100.savefig('picturesCMA_v7/mapLog10Viz_dPz_m_fig6100cma.png')    
      print ('File "picturesCMA_v7/mapLog10Viz_dPz_m_fig6100cma.png" is written')

if (plotFigureFlag == 0):   
   fig6200=plt.figure (6200)
   ax = fig6200.add_subplot(111)                    # for contours plotting
   mapViyPym_c = ax.contourf(X_c,Y_c,log10ionVy_dPy_m[0:nImpctPrmtr_c,0:nVion_c],cmap='jet') 
#   mapViyPy_c = ax.contour(X_c,Y_c,log10ionVy_dPy_m[0:nImpctPrmtr_c,0:nVion_c],levels=range(0,2,1),colors='black') 
   mapViyPy_c = ax.contour(X_c,Y_c,log10ionVy_dPy_m[0:nImpctPrmtr_c,0:nVion_c],7,colors='black') 
   plt.clabel(mapViyPy_c,fmt='%4.2f',inline=True)
   plt.xlabel('Relative Ion Velocity,  $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
   plt.ylabel('Initial Impact Parameter $log_{10}(rho_{Init})$',color='m',fontsize=14)
   titleHeader = 'Transferred Value: $log_{10}(10^{24}\cdot Vy_{ion}\cdot \Delta P_y)$'
   plt.title(titleHeader,color='m',fontsize=12)
   plt.plot([log10relVeTrnsv,log10relVeTrnsv],yLimit,'--m',linewidth=1)
   plt.text(-2.85,-0.8,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([log10relVeLong,log10relVeLong],yLimit,'--m',linewidth=1)
   plt.text(-4.5,-0.8,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-4.14,-.6,'Screened Collisions',color='r',fontsize=16)
   plt.text(-4.14,-2.2,'Magnetized Collisions',color='k',fontsize=16)
# y   plt.xlim([-4.75,-3.25])
# y   plt.ylim([-2.75,-2.0])
   plt.ylim(yLimit)
   fig6200.colorbar(mapViyPym_c)
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig6200.savefig('picturesCMA_v7/mapLog10Viy_dPy_m_fig6200cma.png')    
      print ('File "picturesCMA_v7/mapLog10Viy_dPy_m_fig6200cma.png" is written')

#
# Renumeration 6300 to 6301 because 6300 "is occupied"
#

if (plotFigureFlag == 0):   
   fig6301=plt.figure (6301)
   ax = fig6301.add_subplot(111)                    # for contours plotting
   mapIonEnrg_m1 = ax.contourf(X_c,Y_c,log10deltaEnrgIon_m[0:nImpctPrmtr_c,0:nVion_c],cmap='jet') 
#   mapIonEnrg1 = ax.contour(X_c,Y_c,log10deltaEnrgIon_m[0:nImpctPrmtr_c,0:nVion_c],levels=range(0,2,1),colors='black') 
   mapIonEnrg1 = ax.contour(X_c,Y_c,log10deltaEnrgIon_m[0:nImpctPrmtr_c,0:nVion_c],7,colors='black') 
   plt.clabel(mapIonEnrg1,fmt='%4.2f',inline=True)
   plt.xlabel('Relative Ion Velocity,  $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
   plt.ylabel('Initial Impact Parameter $log_{10}(rho_{Init})$',color='m',fontsize=14)
   titleHeader = 'Transferred Energy: $log_{10}(10^{-19}\cdot \Delta E_{ion})$'
   plt.title(titleHeader,color='m',fontsize=12)
   plt.plot([log10relVeTrnsv,log10relVeTrnsv],yLimit,'--m',linewidth=1)
   plt.text(-2.85,-0.8,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([log10relVeLong,log10relVeLong],yLimit,'--m',linewidth=1)
   plt.text(-4.5,-0.8,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-4.14,-.6,'Screened Collisions',color='r',fontsize=16)
   plt.text(-4.14,-2.2,'Magnetized Collisions',color='k',fontsize=16)
# y   plt.xlim([-4.75,-3.25])
# y   plt.ylim([-2.75,-2.0])
   plt.ylim(yLimit)
   fig6301.colorbar(mapIonEnrg_m1)
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig6301.savefig('picturesCMA_v7/mapLog10enrgIon_m_fig6301cma.png')    
      print ('File "picturesCMA_v7/mapLog10enrgIon_m_fig6301cma.png" is written')

X_f = np.zeros((nImpctPrmtr_c,nVion_c))
Y_f = np.zeros((nImpctPrmtr_c,nVion_c))
log10deltaEnrgIon_m_fit = np.zeros((nImpctPrmtr_c,nVion_c))

for i in range(nVion_c):
   for n in range(nImpctPrmtr_c):
      X_f[n,i] = np.log10(VionRel[i])
      Y_f[n,i] = np.log10(rhoInitFit_dEion[n,i])
      log10deltaEnrgIon_m_fit[n,i] = np.log10(1.e-19*deltaEnrgIon_m_fit[n,i])

if (plotFigureFlag == 0):   
   fig6400=plt.figure (6400)
   ax = fig6400.add_subplot(111)                    # for contours plotting
   mapIonEnrg_m_fit = ax.contourf(X_f,Y_f,log10deltaEnrgIon_m_fit[0:nImpctPrmtr_c,0:nVion_c],cmap='jet') 
#   mapIonEnrg = ax.contour(X_c,Y_c,log10deltaEnrgIon_m_fit[0:nImpctPrmtr_c,0:nVion_c],levels=range(0,2,1),colors='black') 
   mapIonEnrg_fit = ax.contour(X_c,Y_c,log10deltaEnrgIon_m_fit[0:nImpctPrmtr_c,0:nVion_c],7,colors='black') 
   plt.clabel(mapIonEnrg_fit,fmt='%4.2f',inline=True)
   plt.xlabel('Relative Ion Velocity,  $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
   plt.ylabel('Initial Impact Parameter $log_{10}(rho_{Init})$',color='m',fontsize=14)
   titleHeader = 'Transferred Energy (Fitted Data): $log_{10}(10^{-19}\cdot \Delta E_{ion})$'
   plt.title(titleHeader,color='m',fontsize=12)
   plt.plot([log10relVeTrnsv,log10relVeTrnsv],yLimit,'--m',linewidth=1)
   plt.text(-2.85,-0.8,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([log10relVeLong,log10relVeLong],yLimit,'--m',linewidth=1)
   plt.text(-4.5,-0.8,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-4.14,-.6,'Screened Collisions',color='r',fontsize=16)
   plt.text(-4.14,-2.2,'Magnetized Collisions',color='k',fontsize=16)
# y   plt.xlim([-4.75,-3.25])
# y   plt.ylim([-2.75,-2.0])
   plt.ylim(yLimit)
   fig6400.colorbar(mapIonEnrg_m_fit)
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig6400.savefig('picturesCMA_v7/mapLogDeltaEnrgFit_fig6400cma.png')    
      print ('File "picturesCMA_v7/mapLog10DeltaEnrgFit_fig6400cma.png" is written')

if (plotFigureFlag == 0):   
   fig6402=plt.figure (6402)
   ax = fig6402.add_subplot(111)                    # for contours plotting
   mapEnrgIon_m_fit = ax.contourf(X_f,Y_f,1.e-19*deltaEnrgIon_m_fit,cmap='jet') 
#   mapEnrgIon_fit = ax.contour(X_c,Y_c,1.e-18*deltaEnrgIon_m_fit,levels=range(0,2,1),colors='black') 
   mapEnrgIon_fit = ax.contour(X_c,Y_c,1.e-18*deltaEnrgIon_m_fit,7,colors='black') 
   plt.clabel(mapEnrgIon_fit,fmt='%4.2f',inline=True)
   plt.xlabel('Relative Ion Velocity,  $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
   plt.ylabel('Initial Impact Parameter $log_{10}(rho_{Init})$',color='m',fontsize=14)
   titleHeader = 'Transferred Energy (Fitted Data): $10^{-19}\cdot \Delta E_{ion}$, eV'
   plt.title(titleHeader,color='m',fontsize=12)
   plt.plot([log10relVeTrnsv,log10relVeTrnsv],yLimit,'--m',linewidth=1)
   plt.text(-2.85,-0.8,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([log10relVeLong,log10relVeLong],yLimit,'--m',linewidth=1)
   plt.text(-4.5,-0.8,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   plt.text(-4.14,-.6,'Screened Collisions',color='r',fontsize=16)
   plt.text(-4.14,-2.2,'Magnetized Collisions',color='w',fontsize=16)
   plt.ylim(yLimit)
   fig6402.colorbar(mapEnrgIon_m_fit)
   plt.grid(True)
   if (saveFilesFlag == 1):
      fig6302.savefig('picturesCMA_v7/mapDeltaEnrgFit_fig6402cma.png')    
      print ('File "picturesCMA_v7/mapDeltaEnrgFit_fig6402cma.png" is written')


#
# Dependecies of ionVz_dPz_m on impact parameter rhoInit for different
# ion velocities Vion:
#
if (plotFigureFlag == 0):   
   for i in range(12):
      VionCrrnt = V0*VionRel[indxFigures[i]]
      powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
      mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
      figCrrnt = plt.figure(numbrFigures[i]+8)
      plt.loglog(rhoInit[0:nImpctPrmtr,indxFigures[i]], \
                 1.e18*ionVz_dPz_m[0:nImpctPrmtr,indxFigures[i]],'xr',linewidth=2)
      plt.xlabel('Initial Impact Parameter $rho_{Init}$, $cm$',color='m',fontsize=14)
      plt.ylabel('$10^{18}\cdot V_z \cdot \Delta P_z$, g$\cdot$cm$^2$/c$^2$', color='m',fontsize=14)
      titleHeader = 'Transferred Value $V_z\cdot \Delta P_z$ to Single Ion:'
      titleHeader += ' $V_{ion}=%3.1f\cdot10^{%2d}$ $cm/s$'
      plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=14)
#      plt.xlim([.95*rhoInit[0,indxFigures[i]],1.05*rhoInit[nImpctPrmtr-1,indxFigures[i]]])
#       plt.ylim([ .9*[nImpctPrmtr-1,indxFigures[i]], \
#                 1.1*[0,indxFigures[i]]])
#       plt.legend(['Calculated Data',('Fitted Data (Func1): B = %5.3f' % \
#                   abs(fitB_pz[indxFigures[i]])), \
#                  ('Fitted Data (Func2): B = %5.3f'% abs(fitB2_pz[indxFigures[i]]))], \
#                 loc='lower left',fontsize=11)
#       plt.text(xPos[i],yPos[i],'Fitted $\Delta P_z$ are proportional to $rho_{Init}^{-B}$', \
# 	       color='m',fontsize=16)
      plt.grid(True)
      if (saveFilesFlag == 1):
         fileName = 'picturesCMA_v7/ionVz_dPz_indxPlot-'+str(indxFigures[i])+'_fig'
         fileName += str(numbrFigures[i]+8)+'cma.png' 
         figCrrnt.savefig(fileName) 
         print ('File "',fileName,'" is written')

#===================================================
#
# There is fitting of ionVz_dPz_m = Vz_ion*deltaPz_m (these values > 0 always) !!!
#
#===================================================
#
fitA_vz_pz = np.zeros(nVion)            # dimensionless 
fitB_vz_pz = np.zeros(nVion)            # dimensionless
dPosA_vz_pz = np.zeros(nVion)
dNegA_vz_pz = np.zeros(nVion)
dPosB_vz_pz = np.zeros(nVion)
dNegB_vz_pz = np.zeros(nVion)

funcHi2_vz_pz = np.zeros(nVion)
rhoInitFit_vz_pz = np.zeros((nImpctPrmtr,nVion))
ionVz_dPz_m_fit = np.zeros((nImpctPrmtr,nVion))

fitA_vz_pz,fitB_vz_pz,funcHi2_vz_pz,rhoInitFit_vz_pz, ionVz_dPz_m_fit = \
fitting(nImpctPrmtr,nVion,rhoInit,ionVz_dPz_m)
dPosA_vz_pz,dNegA_vz_pz = \
errFitAB(nImpctPrmtr,nVion,rhoInit,ionVz_dPz_m,fitA_vz_pz,fitB_vz_pz,funcHi2_vz_pz,1,2)
dPosB_vz_pz,dNegB_vz_pz = \
errFitAB(nImpctPrmtr,nVion,rhoInit,ionVz_dPz_m,fitA_vz_pz,fitB_vz_pz,funcHi2_vz_pz,2,2)
# print ('Fitting for ionVz_dPz_m:')
# for i in range(nVion):
#    print ('i=%2d: fitA_vz_pz = %e (+%e,-%e), fitB_vz_pz = %e (+%e,-%e), hi2_1 = %e'  % \
#           (i,fitA_vz_pz[i],dPosA_vz_pz[i],dNegA_vz_pz[i], \
# 	     fitB_vz_pz[i],dPosB_vz_pz[i],dNegB_vz_pz[i],funcHi2_vz_pz[i]))
# print ('<fitA_vz_pz> = %e +- %e' % (mean(fitA_vz_pz),mean(dNegA_vz_pz)))
# print ('<fitB_vz_pz> = %e +- %e' % (mean(fitB_vz_pz),mean(dNegB_vz_pz)))


xLimit = [1.015*np.log10(VionRel[0]),.95*np.log10(VionRel[nVion-1])]

yLimMin = 0.
yLimMax = 10.*min(fitA_vz_pz)
if (min(fitA_vz_pz) > 0):
   yLimMin = 10.*max(fitA_vz_pz)
   yLimMax = 0.
for i in range(nVion):
   if (fitA_vz_pz[i] - dNegA_vz_pz[i]) < yLimMin:
      yLimMin = fitA_vz_pz[i] - dNegA_vz_pz[i]
   if (fitA_vz_pz[i] + dPosA_vz_pz[i]) > yLimMax:
      yLimMax = fitA_vz_pz[i] + dPosA_vz_pz[i]
# print ('Exponent A (VzPz): yLimMin = %e, yLimMax = %e' % (yLimMin,yLimMax))

yLimit = [yLimMin-.5,yLimMax+.5]
if (plotFigureFlag == 0):   
   fig3040=plt.figure (3040)
   plt.errorbar(np.log10(VionRel),fitA_vz_pz,yerr=[dNegA_vz_pz,dPosA_vz_pz],fmt='-ro', \
                ecolor='b',capsize=5,capthick=1)
   plt.xlabel('Relative Ion Velocity, $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
   plt.ylabel('Exponent $A$', color='m',fontsize=14)
   titleHeader = 'Dependence of Transferred Value to Single Ion: '
   titleHeader += '$V_z \cdot \Delta P_z$ = $10^A\cdot rho^B$'
   plt.title(titleHeader,color='m',fontsize=12)
   plt.text(-3.75,-19.1,('$V_{e0}=%5.3f\cdot10^{%2d}$cm/s' % (mantV0,powV0)), \
            color='m',fontsize=16)
   plt.text(-4.0,-22.0,('<A>=%7.3f $\pm$ %5.3f' % (mean(fitA_vz_pz),mean(dNegA_vz_pz))), \
            color='r',fontsize=16)
#   plt.text(-2.77,-24.25,('$-$%5.3f' % (mean(dNegA_vz_pz))),color='r',fontsize=12)
#   plt.text(-2.77,-23.75,('$+$%5.3f' % (mean(dPosA_vz_pz))),color='r',fontsize=12)
   plt.xlim(xLimit)
   plt.ylim(yLimit)
   plt.grid(True)
   plt.plot([np.log10(relVeTrnsv),np.log10(relVeTrnsv)],yLimit,'--m',linewidth=1)
   plt.text(-2.55,-22.5,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([np.log10(relVeLong),np.log10(relVeLong)],yLimit,'--m',linewidth=1)
   plt.text(-4.24,-22.5,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   if (saveFilesFlag == 1):
      fig3040.savefig('picturesCMA_v7/fitA_ionVz_dPz_fig3040cma.png')    
      print ('File "picturesCMA_v7/fitA_ionVz_dPz_fig3040cma.png" is written')   

yLimMin = 0.
yLimMax = 10.*min(fitB_vz_pz)
if (min(fitB_vz_pz) > 0):
   yLimMin = 10.*max(fitB_vz_pz)
   yLimMax = 0.
for i in range(nVion):
   if (fitB_vz_pz[i] - dNegB_vz_pz[i]) < yLimMin:
      yLimMin = fitB_vz_pz[i] - dNegB_vz_pz[i]
   if (fitB_vz_pz[i] +  dPosB_vz_pz[i]) > yLimMax:
      yLimMax = fitB_vz_pz[i] +  dPosB_vz_pz[i]
# print ('Exponent B (VzPz): yLimMin = %e, yLimMax = %e' % (yLimMin,yLimMax))

yLimit = [yLimMin-.1,yLimMax+.1]
if (plotFigureFlag == 0):   
   fig3050=plt.figure (3050)
   plt.errorbar(np.log10(VionRel),fitB_vz_pz,yerr=[dNegB_vz_pz,dPosB_vz_pz],fmt='-ro', \
                ecolor='b',capsize=5,capthick=1)
   plt.xlabel('Relative Ion Velocity, $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
   plt.ylabel('Exponent $B$', color='m',fontsize=14)
   titleHeader = 'Dependence of Transferred Value to Single Ion: '
   titleHeader += '$V_z \cdot \Delta P_z$ = $10^A\cdot rho^B$'
   plt.title(titleHeader,color='m',fontsize=12)
   plt.text(-3.75,-.875,('$V_{e0}=%5.3f\cdot10^{%2d}$cm/s' % (mantV0,powV0)), \
            color='m',fontsize=16)
   plt.text(-3.9,-1.6,('<B>=%6.3f $\pm$ %5.3f' % (mean(fitB_vz_pz),mean(dNegB_vz_pz))), \
            color='r',fontsize=16)
#   plt.text(-2.87,-2.55,('$-$%5.3f' % (mean(dNegB_vz_pz))),color='r',fontsize=12)
#   plt.text(-2.87,-2.45,('$+$%5.3f' % (mean(dPosB_vz_pz))),color='r',fontsize=12)
   plt.xlim(xLimit)
   plt.ylim(yLimit)
   plt.grid(True)
   plt.plot([np.log10(relVeTrnsv),np.log10(relVeTrnsv)],yLimit,'--m',linewidth=1)
   plt.text(-2.55,-1.75,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([np.log10(relVeLong),np.log10(relVeLong)],yLimit,'--m',linewidth=1)
   plt.text(-4.24,-1.75,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   if (saveFilesFlag == 1):
      fig3050.savefig('picturesCMA_v7/fitB_ionVz_dPz_fig3050cma.png')    
      print ('File "picturesCMA_v7/fitB_ionVz_dPz_fig3050cma.png" is written')   

yLimMin = 0.
yLimMax = 10.*min(fitA_dEion)
if (min(fitA_dEion) > 0):
   yLimMin = 10.*max(fitA_dEion)
   yLimMax = 0.
for i in range(nVion):
   if (fitA_dEion[i] - dNegA_dEion[i]) < yLimMin:
      yLimMin = fitA_dEion[i] - dNegA_dEion[i]
   if (fitA_dEion[i] + dPosA_dEion[i]) > yLimMax:
      yLimMax = fitA_dEion[i] + dPosA_dEion[i]

yLimit = [yLimMin-.5,yLimMax+.5]
if (plotFigureFlag == 0):   
   fig3060=plt.figure (3060)
   plt.errorbar(np.log10(VionRel),fitA_dEion,yerr=[dNegA_dEion,dPosA_dEion],fmt='-ro', \
                ecolor='b',capsize=5,capthick=1)
   plt.xlabel('Relative Ion Velocity, $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
   plt.ylabel('Exponent $A$', color='m',fontsize=14)
   titleHeader = 'Dependence of Transferred Energy to Single Ion: '
   titleHeader += '$\Delta E_{ion}$ = $10^A\cdot rho^B$'
   plt.title(titleHeader,color='m',fontsize=12)
   plt.text(-3.75,16.25,('$V_{e0}=%5.3f\cdot10^{%2d}$cm/s' % (mantV0,powV0)), \
            color='m',fontsize=16)
   plt.text(-4.0,13.5,('<A>=%7.3f $\pm$ %5.3f' % (mean(fitA_dEion),mean(dNegA_dEion))), \
            color='r',fontsize=16)
#   plt.text(-2.77,-24.25,('$-$%5.3f' % (mean(dNegA_dEion))),color='r',fontsize=12)
#   plt.text(-2.77,-23.75,('$+$%5.3f' % (mean(dPosA_dEion))),color='r',fontsize=12)
   plt.xlim(xLimit)
   plt.ylim(yLimit)
   plt.grid(True)
   plt.plot([np.log10(relVeTrnsv),np.log10(relVeTrnsv)],yLimit,'--m',linewidth=1)
   plt.text(-2.55,12.5,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([np.log10(relVeLong),np.log10(relVeLong)],yLimit,'--m',linewidth=1)
   plt.text(-4.24,12.5,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   if (saveFilesFlag == 1):
      fig3060.savefig('picturesCMA_v7/fitA_dEion_fig3060cma.png')    
      print ('File "picturesCMA_v7/fitA_dEion_fig3060cma.png" is written')   

yLimMin = 0.
yLimMax = 10.*min(fitB_dEion)
if (min(fitB_dEion) > 0):
   yLimMin = 10.*max(fitB_dEion)
   yLimMax = 0.
for i in range(nVion):
   if (fitB_dEion[i] - dNegB_dEion[i]) < yLimMin:
      yLimMin = fitB_dEion[i] - dNegB_dEion[i]
   if (fitB_dEion[i] + dPosB_dEion[i]) > yLimMax:
      yLimMax = fitB_dEion[i] + dPosB_dEion[i]
# print ('Exponent B (dEion): yLimMin = %e, yLimMax = %e' % (yLimMin,yLimMax))

yLimit = [yLimMin-.1,yLimMax+.1]
if (plotFigureFlag == 0):   
   fig3070=plt.figure (3070)
   plt.errorbar(np.log10(VionRel),fitB_dEion,yerr=[dNegB_dEion,dPosB_dEion],fmt='-ro', \
                ecolor='b',capsize=5,capthick=1)
   plt.xlabel('Relative Ion Velocity, $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
   plt.ylabel('Exponent $B$', color='m',fontsize=14)
   titleHeader = 'Dependence of Transferred Energy to Single Ion: '
   titleHeader += '$\Delta E_{ion}$ = $10^A\cdot rho^B$'
   plt.title(titleHeader,color='m',fontsize=12)
   plt.text(-3.75,-.975,('$V_{e0}=%5.3f\cdot10^{%2d}$cm/s' % (mantV0,powV0)), \
            color='m',fontsize=16)
   plt.text(-4.0,-1.5,('<B>=%7.3f $\pm$ %5.3f' % (mean(fitB_dEion),mean(dNegB_dEion))), \
            color='r',fontsize=16)
#   plt.text(-2.77,-24.25,('$-$%5.3f' % (mean(dNegA_dEion))),color='r',fontsize=12)
#   plt.text(-2.77,-23.75,('$+$%5.3f' % (mean(dPosA_dEion))),color='r',fontsize=12)
   plt.xlim(xLimit)
   plt.ylim(yLimit)
   plt.grid(True)
   plt.plot([np.log10(relVeTrnsv),np.log10(relVeTrnsv)],yLimit,'--m',linewidth=1)
   plt.text(-2.55,-1.8,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
   plt.plot([np.log10(relVeLong),np.log10(relVeLong)],yLimit,'--m',linewidth=1)
   plt.text(-4.24,-1.8,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
   if (saveFilesFlag == 1):
      fig3070.savefig('picturesCMA_v7/fitB_dEion_fig3070cma.png')    
      print ('File "picturesCMA_v7/fitB_dEion_fig3070cma.png" is written')   

if (plotFigureFlag == 0):   
   for i in range(12):
      VionCrrnt = V0*VionRel[indxFigures[i]]
      powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
      mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
      figCrrnt = plt.figure(numbrFigures[i]+9)
      plt.loglog(rhoInit[0:nImpctPrmtr,indxFigures[i]], \
                 1.e16*ionVz_dPz_m[0:nImpctPrmtr,indxFigures[i]],'xr', \
                 rhoInitFit_vz_pz[0:nImpctPrmtr,indxFigures[i]], \
	         1.e16*ionVz_dPz_m_fit[0:nImpctPrmtr,indxFigures[i]],'ob',linewidth=2)
      plt.xlabel('Initial Impact Parameter $rho_{Init}$, $cm$',color='m',fontsize=14)
      plt.ylabel('$10^{16} \cdot Vi_z\cdot\Delta P_z$, $eV$', color='m',fontsize=14)
      titleHeader = 'Transferred Value $VI_z\cdot\Delta P_z$ to Single Ion:'
      titleHeader += ' $V_{ion}=%3.1f\cdot10^{%2d}$ $cm/s$'
      plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=14)
      plt.xlim([.95*rhoInit[0,indxFigures[i]],1.05*rhoInit[nImpctPrmtr-1,indxFigures[i]]])
#      plt.ylim([ .9e24*deltaPz_m[nImpctPrmtr-1,indxFigures[i]], \
#                1.1e24*deltaPz_m_fit[0,indxFigures[i]]])
      plt.legend(['Calculated Data', \
                  ('Fitting: $Vi_z\cdot\Delta P_z=10^A\cdot$rho$_{init}^B$; B = %5.3f $\pm$ %5.3f' % \
                   (fitB_vz_pz[indxFigures[i]],dNegB_vz_pz[indxFigures[i]]))],loc='lower left',fontsize=11)
      plt.grid(True)
      if (saveFilesFlag == 1):
         fileName = 'picturesCMA_v7/ionVz_dPz_withFit_indxPlot-'+str(indxFigures[i])+'_fig'
         fileName += str(numbrFigures[i]+9)+'cma.png' 
         figCrrnt.savefig(fileName) 
         print ('File "',fileName,'" is written')

#
# Dependecies of ionVx_dPx_m on ion velocity Vion for different 
# impact parameters rhoInit:
#
# if (plotFigureFlag == 0):   
#    for i in range(12):
#       VionCrrnt = V0*VionRel[indxFigures[i]]
#       powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
#       mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
#       figCrrnt = plt.figure(numbrFigures[i]+9)
#       plt.loglog(rhoInit[0:nImpctPrmtr,indxFigures[i]], \
#                  deltaPx_m[0:nImpctPrmtr,indxFigures[i]],'xr', \
#                  rhoInitFit_px[0:nImpctPrmtr,indxFigures[i]], \
# 	         deltaPx_m_fit[0:nImpctPrmtr,indxFigures[i]],'ob', \
#                  rhoInitFit2_px[0:nImpctPrmtr,indxFigures[i]], \
# 	         deltaPx_m_Fit2[0:nImpctPrmtr,indxFigures[i]], \
# 	         'or',linewidth=2)
#       plt.xlabel('Initial Impact Parameter $rho_{Init}$, $cm$',color='m',fontsize=14)
#       plt.ylabel('$\Delta P_x$, $eV$', color='m',fontsize=14)
#       titleHeader = 'Transferred Momenta $\Delta P_x$ to Single Ion:'
#       titleHeader += ' $V_{ion}=%3.1f\cdot10^{%2d}$ $cm/s$'
#       plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=14)
#       plt.xlim([.95*rhoInit[0,indxFigures[i]],1.05*rhoInit[nImpctPrmtr-1,indxFigures[i]]])
#       plt.ylim([ .9*deltaPx_m[nImpctPrmtr-1,indxFigures[i]], \
#                 1.1*deltaPx_m_Fit2[0,indxFigures[i]]])
#       plt.legend(['Calculated Data',('Fitted Data (Func1): B = %5.3f' % \
#                   abs(fitB_px[indxFigures[i]])), \
#                  ('Fitted Data (Func2): B = %5.3f'% abs(fitB2_px[indxFigures[i]]))], \
#                 loc='lower left',fontsize=11)
#       plt.text(xPos[i],yPos[i],'Fitted $\Delta P_x$ are proportional to $rho_{Init}^{-B}$', \
# 	       color='m',fontsize=16)
#       plt.grid(True)
#       if (saveFilesFlag == 1):
#          fileName = 'picturesCMA/deltaEtransf_indxPlot-'+str(indxFigures[i])+'_fig'
#          fileName += str(numbrFigures[i])+'cma.png' 
#          figCrrnt.savefig(fileName) 
#          print ('File "',fileName,'" is written')

plt.show()

sys.exit()

# if (plotFigureFlag == 0):   
#    fig6005=plt.figure (6005)
#    ax = fig6005.add_subplot(111)                    # for contours plotting
#    mapVixPxm1 = ax.contourf(X,Y,1.e18*ionVx_dPx_m-.6,cmap='jet') 
# #   mapVixPx1 = ax.contour(X,Y,1.e18*(ionVx_dPx_m-.6),levels=range(0,2,1),colors='black') 
#    mapVixPx1 = ax.contour(X,Y,1.e18*(ionVx_dPx_m-.6),7,colors='black') 
#    plt.clabel(mapVixPx1,fmt='%4.2f',inline=True)
#    plt.xlabel('Relative Ion Velocity,  $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
#    plt.ylabel('Initial Impact Parameter $log_{10}(rho_{Init})$',color='m',fontsize=14)
#    titleHeader = '$10^{18}\cdot Vz_{ion}\cdot \Delta P_z$,  g$\cdot $cm$^2$/s$^2$'
#    plt.title(titleHeader,color='m',fontsize=12)
# #   titleHeader = 'z-Component of Ion Velocity $Vz_{ion}$, cm/s'
# #   titleHeader += '\n$\widetilde{Vz_{ion}}$ = $10^{-8} \cdot Vz_{ion}$ '
# #   titleHeader += '($\widetilde{Vz_{ion}}$ = 0.5, if $\widetilde{V_{ion}}$ > 0.5)'
# #   plt.title(titleHeader,color='m',fontsize=12)
#    plt.plot([log10relVeTrnsv,log10relVeTrnsv],yLimit,'--m',linewidth=1)
#    plt.text(-2.85,-0.8,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
#    plt.plot([log10relVeLong,log10relVeLong],yLimit,'--m',linewidth=1)
#    plt.text(-4.5,-0.8,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
#    plt.ylim(yLimit)
#    fig6005.colorbar(mapVixPxm1)
#    plt.grid(True)
# if (saveFilesFlag == 1):
#    fig6005.savefig('picturesCMA/mapVix_m_fig6005cma.png')    
#    print ('File "picturesCMA/mapVix_m_fig6005cma.png" is written')

# if (plotFigureFlag == 0):   
#    fig6105=plt.figure (6105)
#    ax = fig6105.add_subplot(111)                    # for contours plotting
#    mapVizPzm1 = ax.contourf(X,Y,1.e16*ionVz_dPz_m-1.5,cmap='jet') 
#   # mapVizPz1 = ax.contour(X,Y,1.e16*(ionVz_dPz_m-1.5),levels=range(0,2,1),colors='black') 
#    mapVizPz1 = ax.contour(X,Y,1.e16*(ionVz_dPz_m-1.5),7,colors='black') 
#    plt.clabel(mapVizPz1,fmt='%4.2f',inline=True)
#    plt.xlabel('Relative Ion Velocity,  $log_{10}(V_{ion}/V_0)$',color='m',fontsize=14)
#    plt.ylabel('Initial Impact Parameter $log_{10}(rho_{Init})$',color='m',fontsize=14)
#    titleHeader = '$10^{16}\cdot Vx_{ion}\cdot \Delta P_x$,  g$\cdot $cm$^2$/s$^2$'
#    plt.title(titleHeader,color='m',fontsize=12)
#    plt.plot([log10relVeTrnsv,log10relVeTrnsv],yLimit,'--m',linewidth=1)
#    plt.text(-2.85,-0.8,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
#    plt.plot([log10relVeLong,log10relVeLong],yLimit,'--m',linewidth=1)
#    plt.text(-4.5,-0.8,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
#    plt.ylim(yLimit)
#    fig6105.colorbar(mapVizPzm1)
#    plt.grid(True)
# if (saveFilesFlag == 1):
#    fig6105.savefig('picturesCMA/mapVix_m_fig6105cma.png')    
#    print ('File "picturesCMA/mapVix_m_fig6105cma.png" is written')


'''
#
# Opening the output file: 
#
apprchClsscl_file='resultsClassicalApproach.dat'
print ('Open output file "%s"...' % apprchClsscl_file)
apprchClsscl_flag=0
try:
   outfile = open(apprchClsscl_file,'w')
   apprchClsscl_file_flag=1
except:
   print ('Problem to open output file "%s"' % apprchClsscl_file)

#
# Writing the results to output file: 
#
outfile.write ('\n                       Initial data\n\n')
outfile.write ('eVrmsTran=%5.3e cm/s, eVrmsLong=%5.3e cm/s, Ekin =%5.3f eV' \
               % (eVrmsTran,eVrmsLong,ergToEV*kinEnergy))
outfile.write ('\nBfield=%6.1f Gs, ro_Larm= %5.3f mkm, rho_min=%6.3f mkm' \
               % (fieldB[0],1.e4*ro_Larm,1.e4*rhoMin))
outfile.write ('\nV_0=%5.3e cm/s, VionMin=%5.3e cm/s, VionMax= %5.3e cm/s' \
               % (V0,vIonMin,vIonMax))
outfile.write ('\n    Vion/V_0:         min  =%5.3e,        max  = %5.3e' \
               % (vIonMin/V0,vIonMax/V0))
outfile.write ('\n    Vion/eVtran:      min  =%5.3e,        max  = %5.3e' \
               % (vIonMin/eVrmsTran,vIonMax/eVrmsTran))
outfile.write \
('\n\nVion, cm/s     Vion/V_0    Vion/eVtran  ro_max, mkm   Lintr_max, mkm      fitA      fitB\n')

for i in range(nVion):
   outfile.write \
    ('\n %5.3e    %5.3e     %5.3e    %8.2f       %8.2f        %8.4f   %5.3f' \
     % (Vion[i],Vion[i]/V0,Vion[i]/eVrmsTran,1.e4*rhoMax[i], \
        2.e4*halfLintr[0,i],fitA[i],fitB[i]))
'''
   
# sys.exit()

