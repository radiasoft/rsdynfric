# from __future__ import division

#-------------------------------------
#
#        Started at 06/08/2018 (YuE)
#
# This script based on the previous script
# threeApproachesComparison_v6.py
#
# 06/08/2018: IT IS NOT FINISHED!
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
Ekin=3.0e4                             # kinetic energy, eV
curBeam=0.5                            # current density, A/cm^2
dBeam=3.0                              # beam diameter, cm
angSpread=3.0                          # angular spread, mrad
alpha=0.                               # angle between V_i and magnetic field
trnsvT=0.5                             # transversal temperature, eV
longT=2.0e-4                           # longitudinal temperature, eV (was 2.0e-4)
eTempTran=trnsvT                       # to keep variable from previous script
eTempLong=longT                        # to keep variable from previous script
nField=1                               # number ov values  of the magnetic field
fieldB=np.zeros(nField)                 
fieldB[0]=3.e3                         # magnetic field, Gs
B_mag=fieldB[0]                        # to keep variable from previous script
omega_p=1.0e9                          # plasma frequency, 1/sec
n_e=omega_p**2*m_e/(4.*pi*q_e**2)      # plasma density, 3.1421e+08 cm-3

n_e1=8.e7                              # plasma density, cm-3
omega_p1=np.sqrt(4.*pi*n_e1*q_e**2/m_e) # plasma frequency, 5.0459e+08 1/s  
#
# Cooling system parameter:
#
coolLength=150.0        # typical length of the coolong section, cm

# 
# Calculated parameters of the electron beam:
#
V0=np.sqrt(2.*Ekin*eVtoErg/m_e)                # longitudinal velocity, cm/s
rmsTrnsvVe=np.sqrt(2.*trnsvT*eVtoErg/m_e)      # RMS transversal velocity, cm/s
rmsLongVe=np.sqrt(2.*longT*eVtoErg/m_e)        # RMS longitudinal velocity, cm/s
# dens=curBeam*CtoPart/V0                       # density, 1/cm^3
# omega=np.sqrt(4.*pi*dens*q_e**2/m_e)          # plasma frequency, 1/s
cyclFreq=q_e*fieldB/(m_e*cLight)               # cyclotron frequency, 1/s
rmsRoLarm=rmsTrnsvVe*cyclFreq**(-1)            # RMS Larmor radius, cm
dens=omega_p**2*m_e/(4.*pi*q_e**2)             # density, 1/cm^3
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
# ffForm=-.5*omega_p**2*q_e**2/V0**2/MeV_mToErg_cm     # MeV/m
eV_mToErg_m=100.*eVtoErg
# ffForm=-.5*omega_p**2*q_e**2/V0**2/eV_mToErg_m       # =-6.8226e-12 eV/m
eV_mInErg_cm=100.*eVtoErg
ffForm=-.5*omega_p**2*q_e**2/V0**2/eVtoErg            # =-6.8226e-10 eV/cm
ffForm=100.*ffForm                                    # =-6.8226e-08 eV/m

print 'V0=%e cm/s, rmsTrnsvVe=%e cm/s, rmsLongVe=%e cm/s' % (V0,rmsTrnsvVe,rmsLongVe)
rhoCrit=(m_e*(cLight/fieldB[0])**2)**(1./3)
print 'rhoCrit=%e cm' % rhoCrit

#
# Relative velocities of electrons:
#
relVeTrnsv=rmsTrnsvVe/V0 
relVeLong=rmsLongVe/V0

# Indices:
(Ix, Ipx, Iy, Ipy, Iz, Ipz) = range(6)

stepsNumberOnGyro = 25                                             # number of the steps on each Larmour period

#
# Opening the input file: 
#
inputFile='areaOfImpactParameter_tAC-v6_fig110.data'
print 'Open input file "%s"...' % inputFile
inpfileFlag=0
try:
   inpfile = open(inputFile,'r')
   inpfileFlag=1
except:
   print 'Problem to open input file "%s"' % inputFile
if inpfileFlag == 1:
   print 'No problem to open input file "%s"' % inputFile

lines=0                                                            # Number of current line from input file   
dataNumber=0                                                       # Number of current value of any types of Data
xAboundary=np.zeros(100)
xBboundary=np.zeros(100)
while True:
   lineData=inpfile.readline()
#    print 'line=%d: %s' % (lines,lineData)
   if not lineData:
      break
   lines += 1
   if lines > 4:
      words=lineData.split()
      nWords=len(words)
#      print 'Data from %d: words=%s, number of entries = %d' % (lines,words,nWords)
      xAboundary[dataNumber]=float(words[0])
      xBboundary[dataNumber]=float(words[1])
      dataNumber += 1

inpfile.close()
print 'Close input file "%s"' % inputFile

#
# Larmor frequency electron:
#
def omega_Larmor(mass,B_mag):
    return (q_elec)*B_mag/(mass*clight*1.e+2)                      # rad/sec

#
# Derived quantities:
#

omega_L = omega_Larmor(m_elec, B_mag)                              # rad/sec 
T_larm = 2*pi/omega_L                                              # sec
timeStep = T_larm/stepsNumberOnGyro                                # time step, sec
print 'omega_Larmor= %e rad/sec, T_larm = %e sec, timeStep = %e sec' % (omega_L,T_larm,timeStep)

eVrmsTran = np.sqrt(2.*eTempTran*eVtoErg/m_elec)                   # cm/sec
eVrmsLong = np.sqrt(2.*eTempLong*eVtoErg/m_elec)                   # cm/sec
print 'eVrmsTran = %e cm/sec, eVrmsLong = %e cm/sec' % (eVrmsTran,eVrmsLong)

ro_larmRMS = eVrmsTran/omega_L                                     # cm
print 'ro_larmRMS =%e mkm = ', 1.e4*ro_larmRMS
#
# Electrons are magnetized for impact parameter >> rhoCrit:
#
rhoCrit=math.pow(q_elec**2/(m_elec*omega_L**2),1./3)               # cm
print 'rhoCrit (mkm) = ' , 1.e+4*rhoCrit

#
# Convertion from 6-vector of relectron's "coordinates" to 6-vector of guiding-center coordinates:
# z_e=(x_e,px_e,y_e,py_e,z_e,pz_e) --> zgc_e=(phi,p_phi,y_gc,p_gc,z_e,pz_e);
#
def toGuidingCenter(z_e):
    mOmega=m_elec*omega_L                                          # g/sec
    zgc_e=z_e.copy()                                               # 6-vector
    zgc_e[Ix] = np.arctan2(z_e[Ipx]+mOmega*z_e[Iy],z_e[Ipy])       # radians
    zgc_e[Ipx]= (((z_e[Ipx]+mOmega*z_e[Iy])**2+z_e[Ipy]**2)/(2.*mOmega))     # g*cm**2/sec
    zgc_e[Iy] =-z_e[Ipx]/mOmega                                    # cm
    zgc_e[Ipy]= z_e[Ipy]+mOmega*z_e[Ix]                            # g/sec
    return zgc_e

#
# Convertion from 6-vector of guiding-center coordinates to 6-vector of electron's "coordinates":
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
# Matrix to dragg electron through the solenoid with field 'B_mag' during time interval 'deltaT':
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
   mOmegaLarm=m_elec*omega_L                                       # g/sec
   dpFactor_gc=q_elec**2                                           # g*cm^3/sec^2
   rhoLarm_gc=np.sqrt(2.*vectrElec_gc[1]/mOmegaLarm)               # cm
   sinOmega_gc=math.sin(vectrElec_gc[0])
   cosOmega_gc=math.cos(vectrElec_gc[0])
   x_gc=vectrElec_gc[3]/mOmegaLarm                                 # cm
   numer=(vectrIon[0]-x_gc)*cosOmega_gc- \
         (vectrIon[2]-vectrElec_gc[2])*sinOmega_gc                 # cm
   denom=((vectrIon[0]-x_gc)**2+(vectrIon[2]-vectrElec_gc[2])**2+ \
          (vectrIon[4]-vectrElec_gc[4])**2+rhoLarm_gc**2)**(3/2)                # cm^3
   action=vectrElec_gc[1]+dpFactor_gc*numer*rhoLarm_gc/(omega_L*denom)          # g*cm^2/sec
   b_gc=np.sqrt((vectrIon[0]-x_gc)**2+ \
                (vectrIon[2]-vectrElec_gc[2])**2+ \
                (vectrIon[4]-vectrElec_gc[4])**2+2.*action/mOmegaLarm)          # cm
   dpIon[0]=-dpFactor_gc*deltaT*(vectrIon[0]-x_gc)/b_gc**3                      # g*cm/sec
   dpIon[1]=-dpFactor_gc*deltaT*(vectrIon[2]-vectrElec_gc[2])/b_gc**3           # g*cm/sec
   dpIon[2]=-dpFactor_gc*deltaT*(vectrIon[4]-vectrElec_gc[4])/b_gc**3           # g*cm/sec
   dpElec[1]=-dpIon[1]                                                          # g*cm/sec
   dpElec[2]=-dpIon[2]                                                          # g*cm/sec
#    print 'dpIon[0]=%e, dpIon[1]=%e, dpIon[2]=%e' % (dpIon[0],dpIon[1],dpIon[2])
   return dpIon,dpElec,action,b_gc                                      

#
# "Magnus expansion" description of the collision during time interval 'deltaT'
# in the system coordinates of "guiding center" of electron
# input - 6-vectors for electron and ion before collision and time step deltaT; 
# output - transfered momenta to ion and electron and electron y_gc coordinate: 
#
def MagnusExpansionCollision(vectrElec_gc,vectrIon,deltaT):

#    print 'Ion: x=%e, y=%e, z=%e' % (vectrIon[0],vectrIon[2],vectrIon[4])
#    print 'Electron: x=%e, y=%e, z=%e' % (vectrElec_gc[0],vectrElec_gc[4],vectrElec_gc[4])
   dpIon=np.zeros(3)
   dpElec=np.zeros(3)
   mOmegaLarm=m_elec*omega_L                                       # g/sec
   dpFactor_gc=q_elec**2                                           # g*cm^3/sec^2
   rhoLarm_gc=np.sqrt(2.*vectrElec_gc[1]/mOmegaLarm)               # cm
   sinOmega_gc=math.sin(vectrElec_gc[0])
   cosOmega_gc=math.cos(vectrElec_gc[0])
   x_gc=vectrElec_gc[3]/mOmegaLarm                                 # cm
   numer=(vectrIon[0]-x_gc)*cosOmega_gc- \
         (vectrIon[2]-vectrElec_gc[2])*sinOmega_gc                 # cm
   denom=((vectrIon[0]-x_gc)**2+(vectrIon[2]-vectrElec_gc[2])**2+ \
          (vectrIon[4]-vectrElec_gc[4])**2+rhoLarm_gc**2)**(3./2.)              # cm^3
   action=vectrElec_gc[1]+dpFactor_gc*numer*rhoLarm_gc/(omega_L*denom)          # g*cm^2/sec
   C1=np.sqrt((vectrIon[0]-x_gc)**2+ \
              (vectrIon[2]-vectrElec_gc[2])**2+ \
              (vectrIon[4]-vectrElec_gc[4])**2+2.*action/mOmegaLarm)            # cm^2
   C2=2.*((vectrIon[0]-x_gc)*vectrIon[1]/M_ion+(vectrIon[2]-vectrElec_gc[2])*vectrIon[3]/M_ion+ \
          (vectrIon[4]-vectrElec_gc[4])*(vectrIon[5]/M_ion-vectrElec_gc[5]/m_elec))  # cm^2/sec
   C3=(vectrIon[1]/M_ion)**2+(vectrIon[3]/M_ion)**2+ \
      (vectrIon[5]/M_ion-vectrElec_gc[5]/m_elec)**2                # cm^2/sec^2
   b=np.sqrt(C1+C2*deltaT+C3*deltaT**2)                            # cm
   D1=(2.*C3*deltaT+C2)/b-C2/np.sqrt(C1)                           # cm/sec
   D2=(C2*deltaT+2.*C1)/b-2.*np.sqrt(C1)                              # cm
   q=4.*C1*C3-C2**2                                                # cm^4/sec^2
   dpIon[0]=-2.*dpFactor_gc/q*((vectrIon[0]-x_gc)*D1-vectrIon[1]/M_ion*D2)             # g*cm/sec
   dpIon[1]=-2.*dpFactor_gc/q*((vectrIon[2]-vectrElec_gc[2])*D1-vectrIon[3]/M_ion*D2)  # g*cm/sec
   dpIon[2]=-2.*dpFactor_gc/q*((vectrIon[4]-vectrElec_gc[4])*D1- \
                               (vectrIon[5]/M_ion-vectrElec_gc[5]/m_elec)*D2)          # g*cm/sec
   dpElec[1]=-dpIon[1]                                                                 # g*cm/sec
   dpElec[2]=-dpIon[2]                                                                 # g*cm/sec
   dy_gc=dpIon[0]/mOmegaLarm                                                           # cm
#    print 'dpIon[0]=%e, dpIon[1]=%e, dpIon[2]=%e' % (dpIon[0],dpIon[1],dpIon[2])
   return dpIon,dpElec,action,dy_gc                                      


# sys.exit()

z_elecCrrnt_c=np.zeros(6)                              # 6-vector for electron ("classic")
z_ionCrrnt_c=np.zeros(6)                               # 6-vector for ion ("classic")
z_elecCrrnt_gc=np.zeros(6)  # 6-vector for electron in #guiding center" system (both approaches)
z_elecCrrnt_m=np.zeros(6)                              # 6-vector for electron ("magnus")
z_ionCrrnt_m=np.zeros(6)                               # 6-vector for ion ("magnus")

#===============================================================================
#
# Case: the same longidudinal (eVrmsLong).
# The selection of main parameters is sescribed in the
# document "notesToChoiceCodeParameters.docx"
#

rhoMin=1.e-4                                                       # R_crit=1 mkm; cm
Rshield=[100.e-4,150.e-4,200.e-4]                                  # max value of the R_shield; cm
rhoMax=Rshield[0]                                                     # cm

bMin=rhoMin                                                        # cm
bMax=2*rhoMax                                                      # cm (In case with L_int=R_shield sqrt(2) instead 2 is  needed)

nTotal=500
bEdges=np.zeros(nTotal+1)
bStep=(bMax-bMin)/nTotal                                           # cm

print 'bMin(mkm)=%e, bMax(mkm)=%e, bStep(mkm)=%e' % (1.e+4*bMin,1.e+4*bMax,1.e+4*bStep)

for i in range(nTotal+1):
   bEdges[i]=bMin+bStep*i                                          # cm
   
print 'bEdges[0]=%e, bEdges[nTotal]=%e (all sizes in mkm)' % (1.e+4*bEdges[0],1.e+4*bEdges[nTotal])

#
# Transfered momenta for different approaches:
#      deltaP=q_e^2*timeStep*abs(dpApprch_NTot|).
# First index numerates the x,y and z-component of deltaP:
#
dpApprch_cTot=np.zeros((3,nTotal))                                 # 1/cm^2
dpApprch_mTot=np.zeros((3,nTotal))                                 # 1/cm^2

sumPoints_c=0
sumPoints_m=0

#
# Electrons are magnetized for impact parameter >> rhoCrit:
#
rhoCrit=math.pow(q_elec**2/(m_elec*omega_L**2),1./3)               # cm
alpha=2*q_elec**2*omega_L/(m_elec*eVrmsLong**3)                    # dimensionless
beta=2*q_elec**2/(m_elec*eVrmsLong**2)                             # cm
print 'Rshield(mkm)=%e, rhoCrit(mkm)=%f, alpha=%f' % (1.e+4*Rshield[0],1.e+4*rhoCrit,alpha)

#
# Array A=log10(Upot/Ekin)=log10([q_e^2/rho]/[m_e((V_transverse^2+V_longitudinal^2)/2]): 
#
nA=25
crrntA=np.zeros(nA)
minA=-5.
maxA=0.
stepA=(maxA-minA)/(nA-1)

#
# Array B=log10(R_larm/rho): 
#
nB=25
crrntB=np.zeros(nB)
minB=-3.
maxB=-.5
stepB=(maxB-minB)/(nB-1)

evTran=np.zeros((nA,nB))                                         # cm/sec
rho_larm=np.zeros((nA,nB))                                       # cm
kinEnergy=np.zeros((nA,nB))                                      # erg
rhoCrrnt=np.zeros((nA,nB))                                       # cm
#-----
# To check the solution of the  equation for V_transversal, using comparison 
# two formulas to calculate the impact parameter:
rhoCheck=np.zeros((nA,nB))                                       # cm
#-----
halfLintr=np.zeros((nA,nB))                                      # cm
timePath=np.zeros((nA,nB))                                       # sec
numbLarmor=np.zeros((nA,nB))                                     # dimensionless

# To plot area of the impact parameter b; 
# First index =0 for b < Rshield[0], 
#             =1 for Rshield[0] < b < Rshield[1], 
#             =2 for Rshield[1] < b < Rshield[2]: 

mapAimpctParmtr=np.zeros((3,nA*nB))
mapBimpctParmtr=np.zeros((3,nA*nB))
rhoCrrntMax=0.
rhoCrrntMax_A=0
rhoCrrntMax_B=0
trackNumb=-1                                                     # Tracks are enumerated from 0!  
track_1=-1                                                       # Tracks are enumerated from 0!
track_c=-1                                                       # Tracks are enumerated from 0!

for iA in range(nA):
   crrntA[iA]=minA+stepA*iA                                      # log10(Upot/Ekin)
   for iB in range(nB):
      crrntB[iB]=minB+stepB*iB                                   # log10(Rlarm/rho)
#
# Transverse electron velocity V_transverse=V_longidudinal * y, where
# y is the root of equation y^3+y+q=0 for each selected values of A and B:
#
      q=-alpha*math.pow(10.,crrntB[iB]-crrntA[iA])
      D=1./27+(q/2)**2
      root=math.pow(np.sqrt(D)-q/2,1./3)-math.pow(np.sqrt(D)+q/2,1./3)	 
# Checking, that root^3+root+q=0:
      leftPart=root**3+root+q  
# Corrections of the root:
      m=0
      while (abs(leftPart) > 1.e-8) and (m  < 5):
         deltaRoot=-leftPart/(root*(3.*root+1))
         root += deltaRoot
         leftPart=root**3+root+q  
	 m += 1
      evTran[iA,iB]=eVrmsLong*root                                 # transverse velocity; cm/sec
      rho_larm[iA,iB]=evTran[iA,iB]/omega_L                        # Larmor radius; cm
      kinEnergy[iA,iB]=m_elec*(evTran[iA,iB]**2+eVrmsLong**2)/2.   # kinetic energy; erg
      rhoCrrnt[iA,iB]=rho_larm[iA,iB]/math.pow(10.,crrntB[iB])     # current impact parameter; cm
#-----
# Comparison of the two formulas to calculate the impact parameter
# ("bad" solutions exist when abs(rhoCrrnt-rhoCheck) > 1.e-12):      
      rhoCheck[iA,iB]=beta/(1.+root**2)/math.pow(10.,crrntA[iA])   # impact parameter for checking; cm
      if abs(rhoCrrnt[iA,iB]-rhoCheck[iA,iB]) > 1.e-11:
         print 'A=%e, B=%e: abs(rhoCrrnt-rhoCheck)=%e' % \
	       (crrntA[iA],crrntB[iB],abs(rhoCrrnt[iA,iB]-rhoCheck[iA,iB]))
#-----
      if rhoCrrntMax < rhoCrrnt[iA,iB]:
         rhoCrrntMax=rhoCrrnt[iA,iB]
	 rhoCrrntMax_A=crrntA[iA]
	 rhoCrrntMax_B=crrntB[iB]
      if rhoCrrnt[iA,iB] < Rshield[0]:
         trackNumb += 1  
	 mapAimpctParmtr[0,trackNumb]=crrntA[iA]                   # Tracks are enumerated from 0!  
	 mapBimpctParmtr[0,trackNumb]=crrntB[iB]                   # Tracks are enumerated from 0!  
         halfLintr[iA,iB]=np.sqrt(Rshield[0]**2-rhoCrrnt[iA,iB]**2)   # half length of interaction; cm
         timePath[iA,iB]=halfLintr[iA,iB]/eVrmsLong                # time of interaction; sec
         numbLarmor[iA,iB]=int(timePath[iA,iB]/T_larm)             # dimensionless
###	 if trackNumb < 100:
###            print 'iA=%d, iB=%d: track=%d, numbLarmor[iA,iB]=%d: Lintr=%e, rhoLarm=%e, rhoCrrnt=%e' % \
###	          (iA,iB,trackNumb,numbLarmor[iA,iB],1.e4*halfLintr[iA,iB], \
###	           1.e4*rho_larm[iA,iB],1.e4*rhoCrrnt[iA,iB])  
###      else: 
###         print 'iA=%d, iB=%d, rhoLarm=%e, rhoCrrnt=%e' % \
###	       (iA,iB,1.e4*rho_larm[iA,iB],1.e4*rhoCrrnt[iA,iB]) 
#-----
# To plot area of the impact parameter b:
      if Rshield[0] <= rhoCrrnt[iA,iB] < Rshield[1]:
         track_1 += 1  
	 mapAimpctParmtr[1,track_1]=crrntA[iA]                     # Tracks are enumerated from 0!  
	 mapBimpctParmtr[1,track_1]=crrntB[iB]                     # Tracks are enumerated from 0!  
      if Rshield[1] <= rhoCrrnt[iA,iB] < Rshield[2]:
         track_c += 1  
	 mapAimpctParmtr[2,track_c]=crrntA[iA]                     # Tracks are enumerated from 0!  
	 mapBimpctParmtr[2,track_c]=crrntB[iB]                     # Tracks are enumerated from 0!  
#-----

lastTrackNumber=trackNumb
print 'lastTrackNumber=%d' % lastTrackNumber	  
print 'rhoCrrntMax=%e for iA=%e, iB=%e' % (1.e4*rhoCrrntMax,rhoCrrntMax_A,rhoCrrntMax_B)

#-----------------------------------------
#
# Plotting of working arrays:
#
#-----------------------------------------

plotFlagWorkingArrays=1         # Plot working arrays if = 1, orherwise don't plot

if plotFlagWorkingArrays ==1:

#
# Area of initial impact parameter:
#
   plt.figure (10)
   plt.plot(mapAimpctParmtr[0,0:trackNumb],mapBimpctParmtr[0,0:trackNumb],'.r',markersize=10)
   plt.plot(mapAimpctParmtr[1,0:track_1],mapBimpctParmtr[1,0:track_1],'.b',markersize=10)
   plt.plot(mapAimpctParmtr[2,0:track_c],mapBimpctParmtr[2,0:track_c],'.m',markersize=10)
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   Rshield_mkm=1.e4*Rshield[0]
   Rshield_150_mkm=1.e4*Rshield[1]
   Rshield_200_mkm=1.e4*Rshield[2]
   plt.title(('Area of Initial Impact Parameter $b$ (%d Tracks)' % lastTrackNumber), \
             color='m',fontsize=16)
   plt.xlim([minA,maxA])
   plt.ylim([minB,maxB])
   plt.legend([('$b$ < %6.1f $\mu$m' % Rshield_mkm), \
               ('%6.1f $\mu$m < $b$ < %6.1f $\mu$m' % (Rshield_mkm,Rshield_150_mkm)), \
      	       ('%6.1f $\mu$m < $b$ < %6.1f $\mu$m' % (Rshield_150_mkm,Rshield_200_mkm))], \
   	      loc='lower left',fontsize=16)
   plt.grid(True)

   X,Y=np.meshgrid(crrntA,crrntB)      

#
# Half length of interaction:
#
   fig20=plt.figure(20)
   ax20=fig20.gca(projection='3d')
   surf=ax20.plot_surface(X,Y,1.e+4*halfLintr,cmap=cm.coolwarm,linewidth=0,antialiased=False)
   plt.title('Half Length of Interaction $L_{interaction}$, $\mu$m', color='m',fontsize=16)
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   ax20.set_zlabel('$1/2 \cdot L_{interaction}$',color='m',fontsize=16)
   fig20.colorbar(surf, shrink=1., aspect=10)
   plt.grid(True)

#
# Number of Larmor turns:
#
   fig30=plt.figure(30)
   ax30=fig30.gca(projection='3d')
   surf=ax30.plot_surface(X,Y,numbLarmor,cmap=cm.coolwarm,linewidth=0,antialiased=False)
   plt.title('Number of Larmor Turns $N_{Larmor}$', color='m',fontsize=16)
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   ax30.set_zlabel('$N_{Larmor}$',color='m',fontsize=16)
   fig30.colorbar(surf, shrink=1., aspect=10)
   plt.grid(True)

# plt.show()

#-----------------------------------------
#
# Definitions:
#
#-----------------------------------------
larmorNumber=np.zeros(lastTrackNumber+1)
larmorRadius=np.zeros(lastTrackNumber+1)                           # cm

'''
# 
# For approach_1 (without averaging over larmor rotation):
#
matr_elec=solenoid_eMatrix(B_mag,timeStep)                         # matrix for electron for timeStep in magnetic field
matr_ion=drift_Matrix(M_ion,timeStep)                              # matrix for ion (with mass M_ion) for timeStep
'''
trackNumb_1=-1                                                     # Tracks will be enumerated from 0!
'''
dpFactor=q_elec**2*timeStep                                        # g*cm^3/sec
minLarmR_b_1=1.e8                                                  # dimensionless     
maxLarmR_b_1=-1.e8                                                 # dimensionless
minUpot_enrgKin_1=1.e8                                             # dimensionless      
maxUpot_enrgKin_1=-1.e8                                            # dimensionless      
timePoints_1=np.zeros(lastTrackNumber+1)                           # for approach_1; dimensionless
pointTrack_1=np.zeros(lastTrackNumber+1)                           # for approach_1; dimensionless  
'''

# 
# For approach_c (with averaging over nLarmorAvrgng larmor rotation) +
# "Magnus extansion" method to calculate the transferred momenmta):
#
nLarmorAvrgng=1                                              # number of averaged Larmor rotations 
timeStep_c=nLarmorAvrgng*stepsNumberOnGyro*timeStep          # time step for approach_c
matr_elec_c=guidingCenter_Matrix(.5*timeStep_c)       # matrix for electron for .5*timeStep_c 
matr_ion_c=drift_Matrix(M_ion,.5*timeStep_c)  # matrix for ion (with mass M_ion) for .5*timeStep_c
trackNumb_c=-1                                               # Tracks will be enumerated from 0!
minLarmR_b_c=1.e8                                            # dimensionless     
maxLarmR_b_c=-1.e8                                           # dimensionless
minUpot_enrgKin_c=1.e8                                       # dimensionless      
maxUpot_enrgKin_c=-1.e8                                      # dimensionless      
timePoints_c=np.zeros(lastTrackNumber+1)                     # for approach_c; dimensionless
print 'lastTrackNumber=%d' % lastTrackNumber
pointTrack_c=np.zeros(lastTrackNumber+1)                     # for approach_c; dimensionless  

# 
# For approach_m (with averaging over nLarmorAvrgng larmor rotation +
# "Guiding center" method to calculate the transferred momenmta):
#
nLarmorAvrgng=1                                              # number of averaged Larmor rotations 
timeStep_m=nLarmorAvrgng*stepsNumberOnGyro*timeStep          # time step for approach_m
matr_elec_m=guidingCenter_Matrix(.5*timeStep_m)       # matrix for electron for .5*timeStep_m 
matr_ion_m=drift_Matrix(M_ion,.5*timeStep_m)  # matrix for ion (with mass M_ion) for .5*timeStep_m
trackNumb_m=-1                                               # Tracks will be enumerated from 0!
minLarmR_b_m=1.e8                                            # dimensionless     
maxLarmR_b_m=-1.e8                                           # dimensionless
minUpot_enrgKin_m=1.e8                                       # dimensionless      
maxUpot_enrgKin_m=-1.e8                                      # dimensionless      
timePoints_m=np.zeros(lastTrackNumber+1)                     # for approach_m; dimensionless
pointTrack_m=np.zeros(lastTrackNumber+1)                     # for approach_m; dimensionless  

print 'timeStep=%e, timeStep_c=%e, timeStep_m=%e' % (timeStep,timeStep_c,timeStep_m)

# 
# To draw track with maximal transfer of px:
#
'''
maxAbsDpxApprch_1=0.
'''
maxAbsDpxApprch_c=0.
maxAbsDpxApprch_m=0.

'''
totAbsDpxApprch_1=0.
'''
totAbsDpxApprch_c=0.
totAbsDpxApprch_m=0.

'''
indxAmaxAbsDpxApprch_1=0
'''
indxAmaxAbsDpxApprch_c=0
indxAmaxAbsDpxApprch_m=0

'''
indxBmaxAbsDpxApprch_1=0
'''
indxBmaxAbsDpxApprch_c=0
indxBmaxAbsDpxApprch_m=0

'''
trackNumbMaxAbsDpxApprch_1=0
'''
trackNumbMaxAbsDpxApprch_c=0
trackNumbMaxAbsDpxApprch_m=0

#
# Working arrays:
#
larmorNumber=np.zeros(nA*nB)
larmorRadius=np.zeros(nA*nB)

'''
cpuTime_1=np.zeros(nA*nB)
'''

cpuTime_c=np.zeros(nA*nB)

pointTrack_m=np.zeros(nA*nB)  
timePoints_m=np.zeros(nA*nB)
cpuTime_m=np.zeros(nA*nB)

evTran=np.zeros((nA,nB))                                           # cm/sec
rho_larm=np.zeros((nA,nB))                                         # cm
kinEnergy=np.zeros((nA,nB))                                        # erg
rhoCrrnt=np.zeros((nA,nB))                                         # cm
halfLintr=np.zeros((nA,nB))                                        # cm
timePath=np.zeros((nA,nB))                                         # sec
numbLarmor=np.zeros((nA,nB))                                       # dimensionless

'''
dpxTotal_1=np.zeros((nA,nB))                                       # g*cm/sec
dpyTotal_1=np.zeros((nA,nB))                                       # g*cm/sec
dpzTotal_1=np.zeros((nA,nB))                                       # g*cm/sec

dpxTotalMax_1=0.                                                   # g*cm/sec
dpyTotalMax_1=0.                                                   # g*cm/sec
dpzTotalMax_1=0.                                                   # g*cm/sec
'''

dpxTotal_c=np.zeros((nA,nB))                                       # g*cm/sec
dpyTotal_c=np.zeros((nA,nB))                                       # g*cm/sec
dpzTotal_c=np.zeros((nA,nB))                                       # g*cm/sec

dpxTotalMax_c=0.                                                   # g*cm/sec
dpyTotalMax_c=0.                                                   # g*cm/sec
dpzTotalMax_c=0.                                                   # g*cm/sec

dpxTotal_m=np.zeros((nA,nB))                                       # g*cm/sec
dpyTotal_m=np.zeros((nA,nB))                                       # g*cm/sec
dpzTotal_m=np.zeros((nA,nB))                                       # g*cm/sec

dpxTotalMax_m=0.                                                   # g*cm/sec
dpyTotalMax_m=0.                                                   # g*cm/sec
dpzTotalMax_m=0.                                                   # g*cm/sec

'''
maxYcoorElec_1=0.
maxYcoorIon_1=0.
'''
maxYcoorElec_c=0.
maxYcoorIon_c=0.
maxYcoorElec_m=0.
maxYcoorIon_m=0.

cpuTimeTotal=0

'''
runFlagApproach_1=1             # run approach1 if = 1, otherwise don't run
'''
runFlagApproach_c=1             # run approach2 if = 1, otherwise don't run
runFlagApproach_m=1             # run approach3 if = 1, otherwise don't run

plotFlagTracks=1                # Plot everything about tracks if = 1, orherwise don't plot
plotFlagDpTransf=1              # Plot everything about dp transf. if = 1, orherwise don't plot

print 'nA=%d, nB=%d' % (nA,nB)

for iA in range(nA):
   crrntA[iA]=minA+stepA*iA                                      # log10(Upot/Ekin)
   for iB in range(nB):
      crrntB[iB]=minB+stepB*iB                                   # log10(Rlarm/rho)
#
# Transverse electron velocity V_transverse=V_longidudinal * y, where
# y is the root of equation y^3+y+q=0 for each selected values of A and B:
#
      q=-alpha*math.pow(10.,crrntB[iB]-crrntA[iA])
      D=1./27+(q/2)**2
      root=math.pow(np.sqrt(D)-q/2,1./3)-math.pow(np.sqrt(D)+q/2,1./3)	 
# Checking, that root^3+root+q=0:
      leftPart=root**3+root+q  
# Corrections of the root:
      m=0
      while (abs(leftPart) > 1.e-8) and (m  < 5):
         deltaRoot=-leftPart/(root*(3.*root+1))
         root += deltaRoot
         leftPart=root**3+root+q  
	 m += 1
      evTran[iA,iB]=eVrmsLong*root                                 # transverse velocity; cm/sec
      rho_larm[iA,iB]=evTran[iA,iB]/omega_L                        # Larmor radius; cm
      kinEnergy[iA,iB]=m_elec*(evTran[iA,iB]**2+eVrmsLong**2)/2.   # kinetic energy; erg
      rhoCrrnt[iA,iB]=rho_larm[iA,iB]/math.pow(10.,crrntB[iB])     # current impact parameter; cm
###      print 'A(%d)=%e, B(%d)=%e: rhoCrrnt=%e, Rshield[0]=%e' % \
###            (iA,crrntA[iA],iB,crrntB[iB],rhoCrrnt[iA,iB],Rshield[0])
      if rhoCrrnt[iA,iB] < Rshield[0]:
         print 'Track: %d' % trackNumb_1
         trackNumb_1 += 1  
         trackNumb_c += 1  
         trackNumb_m += 1  
         halfLintr[iA,iB]=np.sqrt(Rshield[0]**2-rhoCrrnt[iA,iB]**2)   # half length of interaction; cm
         timePath[iA,iB]=halfLintr[iA,iB]/eVrmsLong                # time of interaction; sec
         numbLarmor[iA,iB]=int(timePath[iA,iB]/T_larm)             # dimensionless
###         print 'iA=%d, iB=%d: track=%d, numbLarmor[iA,iB]=%d: Lintr=%e, rhoLarm=%e, rhoCrrnt=%e' % \
###	       (iA,iB,trackNumb_1,numbLarmor[iA,iB],1.e4*halfLintr[iA,iB], \
###	        1.e4*rho_larm[iA,iB],1.e4*rhoCrrnt[iA,iB])  
         larmorNumber[trackNumb_1]=numbLarmor[iA,iB]
         larmorRadius[trackNumb_1]=rho_larm[iA,iB]                 # cm
#--          timePoints_1[trackNumb_1]=int(numbLarmor[iA,iB]*stepsNumberOnGyro) # for approach_1; dimensionless
         timePoints_c[trackNumb_c]=int(numbLarmor[iA,iB]/nLarmorAvrgng)     # for approach_c; dimensionless
         timePoints_m[trackNumb_m]=int(numbLarmor[iA,iB]/nLarmorAvrgng)     # for approach_m; dimensionless
###	 print 'timePoints_1=%d, timePoints_c=%d, timePoints_m=%d, numbLarmor[iA,iB]=%d' % \
###	       (timePoints_1[trackNumb_1],timePoints_c[trackNumb_c],timePoints_m[trackNumb_m], \
###	        larmorNumber[trackNumb_1])   
	 print 'timePoints_c=%d, timePoints_m=%d, numbLarmor[iA,iB]=%d' % \
	       (timePoints_c[trackNumb_c],timePoints_m[trackNumb_m], \
	        larmorNumber[trackNumb_1])   

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Approach_1: dragging without averaging over larmor rotation
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# To draw only first trajectories (for checking only):
#
#------- Start of calculations for approach_1 --------------
#
#--          if runFlagApproach_1 == 1:
#--             timeStart=os.times()
#--             if trackNumb_1 == 0:
#-- 	       rhoFirstTurn=rhoCrrnt[iA,iB]
#-- 	       rhoLarmorFirstTurn=rho_larm[iA,iB]
#-- #
#-- # 6-vectors for ion and electron and distance 'b' between them for the first trajectory and 
#-- # (T)rajectory with (M)aximal (T)ransferred (dpx) ( trajectory TMTdpx);
#-- # (for checking only; indices 0-5 for electron, indices 6-11 for ion and index=12 for 'b'):
#-- #
#--                prtclCoorFirst_1=np.zeros((13,timePoints_1[trackNumb_1])) # first trajectory                 
#--             prtclCoorCrrnt_1=np.zeros((13,timePoints_1[trackNumb_1])) # current trajectory                
#--             prtclCoorMaxAbsDpx_1=np.zeros((13,timePoints_1[trackNumb_1]))      # trajectory TMTdpx
#-- # Current distance from origin of the coordinate system to electron along the trajectory; cm
#--             bCrrnt_1=np.zeros(timePoints_1[trackNumb_1])           # cm
#-- # Current log10 of two important ratios:
#-- # First - ratio R_larmor/b; dimensionless
#--             larmR_bCrrnt_1=np.zeros(timePoints_1[trackNumb_1])       
#-- # Second - ratio potential_energy/kinetic_energy; dimensionless
#--             uPot_enrgKinCrrnt_1=np.zeros(timePoints_1[trackNumb_1])  
#-- # deltaPapprch_1=dpFactor*dpApprch_1Crrnt; dpFactor=q_e^2*timeStep; 1/cm^2 
#--             dpApprch_1Crrnt=np.zeros((3,timePoints_1[trackNumb_1]))
#--             for m in range(6): 
#--                z_ionCrrnt_1[m]=0.                     # Current initial zero-vector for ion
#--                z_elecCrrnt_1[m]=0.                    # Zeroing out of vector for electron
#-- # Current initial vector for electron:
#--             z_elecCrrnt_1[Ix]=rhoCrrnt[iA,iB]                      # x, cm
#--             z_elecCrrnt_1[Iz]=-halfLintr[iA,iB]                    # z, cm
#--             z_elecCrrnt_1[Ipy]=m_elec*evTran[iA,iB]                # py, g*cm/sec
#--             z_elecCrrnt_1[Ipz]=m_elec*eVrmsLong                    # pz, g*cm/sec
#-- #            print 'Track %d: z=%e mkm' % (trackNumb_1,1.e4*z_elecCrrnt_1[Iz])
#-- #### To draw only first trajectory (for checking only):
#-- ###            for ic in range(6):
#-- ###               prtclCoorFirst_1[ic,0]=z_elecCrrnt_1[ic]  # 6-vector for electron
#-- ###               prtclCoorFirst_1[6+ic,0]=z_ionCrrnt_1[ic]       # 6-vector for ion 
#-- ###	    prtclCoorFirst_1[12,0]=np.sqrt(halfLintr[iA,iB]**2+rhoFirstTurn**2)   # cm
#-- #-----------------------------------------------
#-- # Main action - dragging of the current trajectories (for given i and j)
#-- #
#--             for k in range(int(timePoints_1[trackNumb_1])):
#--  	       z_elecCrrnt_1=np.dot(matr_elec,z_elecCrrnt_1)       # electron's dragging
#--  	       z_ionCrrnt_1=np.dot(matr_ion,z_ionCrrnt_1)          # ion's dragging
#-- # Current distance between ion and electron; cm:   
#--  	       bCrrnt_1[k]=np.sqrt((z_ionCrrnt_1[0]-z_elecCrrnt_1[0])**2+ \
#-- 	                           (z_ionCrrnt_1[2]-z_elecCrrnt_1[2])**2+ \
#-- 			           (z_ionCrrnt_1[4]-z_elecCrrnt_1[4])**2)
#-- # Current log10 of two important ratios:  
#-- 	       larmR_bCrrnt_1[k]=math.log10(rho_larm[iA,iB]/bCrrnt_1[k])      # dimensionless 
#-- 	       if maxLarmR_b_1 < larmR_bCrrnt_1[k]:
#-- 	          maxLarmR_b_1=larmR_bCrrnt_1[k]
#-- 	       if minLarmR_b_1 > larmR_bCrrnt_1[k]:
#-- 	          minLarmR_b_1=larmR_bCrrnt_1[k]
#-- 	       uPot_enrgKinCrrnt_1[k]=math.log10((q_elec**2/bCrrnt_1[k])/kinEnergy[iA,iB])
#-- 	       if maxUpot_enrgKin_1 < uPot_enrgKinCrrnt_1[k]:
#-- 	          maxUpot_enrgKin_1=uPot_enrgKinCrrnt_1[k]
#-- 	       if minUpot_enrgKin_1 > uPot_enrgKinCrrnt_1[k]:
#-- 	          minUpot_enrgKin_1=uPot_enrgKinCrrnt_1[k]
#-- # Current values to calculate deltaPapprch_1=dpFactor*dpApprch_1Crrnt; dpFactor=q_e^2*timeStep:  
#--                for ic in range(3):
#-- 	          dpApprch_1Crrnt[ic,k]=-(z_ionCrrnt_1[2*ic]-z_elecCrrnt_1[2*ic])/ \
#-- 		                         bCrrnt_1[k]**3        # 1/cm^2;  
#-- # To be prepared to draw future TMTdpx trajectory (for checking only):
#--                for ic in range(6):
#--                   prtclCoorCrrnt_1[ic,pointTrack_1[trackNumb_1]]=z_elecCrrnt_1[ic]  # 6-vector for electron
#--                   prtclCoorCrrnt_1[6+ic,pointTrack_1[trackNumb_1]]=z_ionCrrnt_1[ic] # 6-vector for ion 
#-- 	       prtclCoorCrrnt_1[12,pointTrack_1[trackNumb_1]]=bCrrnt_1[k]           # cm
#-- # To draw only first trajectory (for checking only):
#--                if trackNumb_1 == 0:
#--                   for ic in range(6):
#--                      prtclCoorFirst_1[ic,pointTrack_1[trackNumb_1]]=z_elecCrrnt_1[ic]  # 6-vector for electron
#--                      prtclCoorFirst_1[6+ic,pointTrack_1[trackNumb_1]]=z_ionCrrnt_1[ic]       # 6-vector for ion 
#-- 	          prtclCoorFirst_1[12,pointTrack_1[trackNumb_1]]=bCrrnt_1[k]                 # cm
#-- 	          if maxYcoorElec_1 < abs(z_elecCrrnt_1[2]):
#-- 	             maxYcoorElec_1=abs(z_elecCrrnt_1[2])
#-- 	          if maxYcoorIon_1 < abs(z_ionCrrnt_1[2]):
#-- 	             maxYcoorIon_1=abs(z_ionCrrnt_1[2])
#-- #	          crrntPoint=pointTrack_1[trackNumb_1]
#-- #	          if crrntPoint < 100 or crrntPoint > 45400:
#-- #                     print 'Point %d: x=%e, y=%e, z=%e' %  (crrntPoint,\
#-- #	                     1.e+7*prtclCoorFirst_1[6,crrntPoint], \
#-- #                            1.e+7*prtclCoorFirst_1[8,crrntPoint], \
#-- #                            1.e+7*prtclCoorFirst_1[10,crrntPoint])
#-- ####
#-- #### Taking into account transfer of momentum for both particles:
#-- ####
#-- ###               for ic in range(3):
#-- ###                  z_ionCrrnt_1[2*ic+1] += -dpFactor*dpApprch_1Crrnt[ic,k]
#-- ###                  z_elecCrrnt_1[2*ic+1] += dpFactor*dpApprch_1Crrnt[ic,k]
#--                pointTrack_1[trackNumb_1] += 1
#-- # To draw the distribution of transferred dp (g*cm/sec):
#--                dpxTotal_1[iA,iB] +=-dpFactor*dpApprch_1Crrnt[0,k] 
#--                dpyTotal_1[iA,iB] +=-dpFactor*dpApprch_1Crrnt[1,k] 
#--                dpzTotal_1[iA,iB] +=-dpFactor*dpApprch_1Crrnt[2,k] 
#-- #
#-- # Total transferred dpx  for current track:
#-- #
#-- 	       totAbsDpxApprch_1 += abs(dpFactor*dpApprch_1Crrnt[0,k])
#-- # 
#-- # End of dragging of the current trajectory	  
#-- #-----------------------------------------------
#-- #
#-- # Accumulate transferred momenta and other data: 
#-- #
#--             if trackNumb_1 == 0: 
#-- # First definition of the total distance from origin of the 
#-- # coordinate system to electron along the trajectory; cm:
#--  	       b_1=bCrrnt_1                                        # cm
#-- # First definition of the log10 of two important ratios; dimensionless:  
#-- 	       uPot_enrgKin_1=uPot_enrgKinCrrnt_1                  # ratio A 
#-- 	       larmR_b_1=larmR_bCrrnt_1                            # ratio B
#-- #### First definition of the values to calculate deltaPapprch_1 (g*cm/sec):
#-- ###               dpxApprch_1=dpFactor*dpApprch_1Crrnt[0,:] 
#-- ###               dpyApprch_1=dpFactor*dpApprch_1Crrnt[1,:]  
#-- ###               dpzApprch_1=dpFactor*dpApprch_1Crrnt[2,:] 
#--             else:  
#-- #### Total distance from origin of the coordinate system to electron along the trajectory:
#--  	       b_1=np.concatenate((b_1,bCrrnt_1),axis=0)           # cm
#-- # Total log10 of two important ratios; dimensionless:  
#-- 	       uPot_enrgKin_1=np.concatenate((uPot_enrgKin_1,uPot_enrgKinCrrnt_1),axis=0)        
#-- 	       larmR_b_1=np.concatenate((larmR_b_1,larmR_bCrrnt_1),axis=0)                  
#-- #### Total values to calculate deltaPapprch_1 (g*cm/sec):
#-- ###               dpxApprch_1=np.concatenate((dpxApprch_1,dpFactor*dpApprch_1Crrnt[0,:]),axis=0)
#-- ###               dpyApprch_1=np.concatenate((dpyApprch_1,dpFactor*dpApprch_1Crrnt[1,:]),axis=0)
#-- ###               dpzApprch_1=np.concatenate((dpzApprch_1,dpFactor*dpApprch_1Crrnt[2,:]),axis=0)
#-- # To draw TMTdpx trajectory (for checking only):
#-- #
#--             if maxAbsDpxApprch_1 < totAbsDpxApprch_1:
#-- 	       maxAbsDpxApprch_1=totAbsDpxApprch_1
#-- 	       indxAmaxAbsDpxApprch_1=iA
#-- 	       indxBmaxAbsDpxApprch_1=iB
#-- 	       trackNumbMaxAbsDpxApprch_1=trackNumb_1
#-- 	       prtclCoorMaxAbsDpx_1=prtclCoorCrrnt_1.copy()   
#-- 	       rhoMaxAbsDpxTurn_1=rhoCrrnt[iA,iB]
#-- 	       rhoLarmorMaxAbsDpxTurn_1=rho_larm[iA,iB]
#-- ###	       print 'iA=%d, iB=%d: TMT-track %d, points %d' % \
#-- ###	             (indxAmaxAbsDpxApprch_1,indxBmaxAbsDpxApprch_1, \
#-- ###		     trackNumbMaxAbsDpxApprch_1,pointTrack_1[trackNumbMaxAbsDpxApprch_1]) 
#-- ###	       print 'timePoints.shape: ', \
#-- ###	             (prtclCoorMaxAbsDpx_1.shape[0],prtclCoorMaxAbsDpx_1.shape[1])
#-- # End of all calculations for approach_1
#--             sumPoints_1 += pointTrack_1[trackNumb_1]
#-- # 
#-- # To draw the distribution of transferred dp:
#-- #
#-- 	    if dpxTotalMax_1 < abs(dpxTotal_1[iA,iB]):
#-- 	       dpxTotalMax_1=abs(dpxTotal_1[iA,iB])
#-- 	    if dpyTotalMax_1 < abs(dpyTotal_1[iA,iB]):
#-- 	       dpyTotalMax_1=abs(dpyTotal_1[iA,iB])
#-- 	    if dpzTotalMax_1 < abs(dpzTotal_1[iA,iB]):
#-- 	       dpzTotalMax_1=abs(dpzTotal_1[iA,iB])
#-- #            print 'Track %d: dpxTotalMax_1=%e, dpyTotalMax_1=%e, dpzTotalMax_1=%e' % \
#-- #	          (trackNumb_1,dpxTotalMax_1,dpyTotalMax_1,dpzTotalMax_1)
#--             timeEnd=os.times()
#-- 	    cpuTime_1[trackNumb_1]=1.e+6*(float(timeEnd[0])-float(timeStart[0]))  # CPU time , mks
#--             cpuTimeTotal += cpuTime_1[trackNumb_1] 
#-- #
#-- #------- End of approach_1 --------------
#-- #
#-- # for i in range(lastTrackNumber):
#-- #    print 'Track %d: larmor turns=%d, cpuTime(mks)=%e, time per turn(mks)=%6.1f' % \
#-- #          (i,larmorNumber[i],cpuTime[i],cpuTime[i]/larmorNumber[i])
#-- 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Approach_c: dragging with averaging over nLarmorAvrgng larmor rotation +
#             "Guiding center" method to calculate the transferred momenta
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#------- Start of calculations for approach_c --------------
#
         if runFlagApproach_c == 1:
            timeStart=os.times()
            if trackNumb_c == 0:
	       rhoFirstTurn=rhoCrrnt[iA,iB]
	       rhoLarmorFirstTurn=rho_larm[iA,iB]
#
# 6-vectors for ion and electron and distance 'b' between them for the first trajectory and 
# (T)rajectory with (M)aximal (T)ransferred (dpx) ( trajectory TMTdpx);
# (for checking only; indices 0-5 for electron, indices 6-11 for ion,
#  index=12 for 'b', index=13 for action and index=14 for b_gc):
#
               prtclCoorFirst_c=np.zeros((15,timePoints_c[trackNumb_c])) # first trajectory                 
            prtclCoorCrrnt_c=np.zeros((15,timePoints_c[trackNumb_c])) # current trajectory                
            prtclCoorMaxAbsDpx_c=np.zeros((15,timePoints_c[trackNumb_c]))      # trajectory TMTdpx
# Current distance from origin of the coordinate system to electron along the trajectory; cm
            bCrrnt_c=np.zeros(timePoints_c[trackNumb_c])           # cm
# Current log10 of two important ratios:
# First - ratio R_larmor/b; dimensionless
            larmR_bCrrnt_c=np.zeros(timePoints_c[trackNumb_c])       
# Second - ratio potential_energy/kinetic_energy; dimensionless
            uPot_enrgKinCrrnt_c=np.zeros(timePoints_c[trackNumb_c])  
# deltaPapprch_m=dpApprch_cCrrnt
            dpApprch_cCrrnt=np.zeros((3,timePoints_c[trackNumb_c]))
            for m in range(6): 
               z_ionCrrnt_c[m]=0.                     # Current initial zero-vector for ion
               z_elecCrrnt_c[m]=0.                    # Zeroing out of vector for electron
# Current initial vector for electron:
            z_elecCrrnt_c[Ix]=rhoCrrnt[iA,iB]+rho_larm[iA,iB]      # x, cm
            z_elecCrrnt_c[Iz]=-halfLintr[iA,iB]                    # z, cm
            z_elecCrrnt_c[Ipy]=m_elec*evTran[iA,iB]                # py, g*cm/sec
            z_elecCrrnt_c[Ipz]=m_elec*eVrmsLong                    # pz, g*cm/sec
	    z_elecCrrnt_gc=toGuidingCenter(z_elecCrrnt_c)          # transfer to system of guiding center
#	    if iA == 0 and iB == 0:
#	       print 'z_elecCrrnt_c: ', z_elecCrrnt_c 
#	       print 'z_elecCrrnt_gc: ', z_elecCrrnt_gc 
#-----------------------------------------------
# Main action - dragging of the current trajectories (for given i and j)
#
            for k in range(int(timePoints_c[trackNumb_c])):
#	       if k < 100: 
#	          print 'k=%d: x=%e, y=%e, z=%e' % \
#                      (k,1.e+4*z_elecCrrnt_gc[0],1.e+4*z_elecCrrnt_gc[2],1.e+4*z_elecCrrnt_gc[4]) 
#	       if k > timePoints_c[trackNumb_c]-100: 
#	          print 'k=%d: x=%e, y=%e, z=%e' % \
#                       (k,1.e+4*z_elecCrrnt_gc[0],1.e+4*z_elecCrrnt_gc[2],1.e+4*z_elecCrrnt_gc[4]) 
#
# dragging both paticles through first half of trajectory:
#
 	       z_elecCrrnt_gc=np.dot(matr_elec_c,z_elecCrrnt_gc)   # electron
 	       z_ionCrrnt_c=np.dot(matr_ion_c,z_ionCrrnt_c)        # ion
#
# dragging both paticles through interaction point:
#
	       dpIon,dpElec,action,b_gc=guidingCenterCollision(z_elecCrrnt_gc,z_ionCrrnt_c,timeStep_c) 
#	       if trackNumb_c == 0:
#	          print 'point %d: dpgcElec=%e, dpzElec=%e' % 
#                       (pointTrack_c[0],dpElec[1],dpElec[2])
	       for ic in range(3):
####
#### Taking into account transfer of momentum for both particles:
####
####	          z_ionCrrnt_c[2*ic+1] += dpIon[ic]   
####	          z_elecCrrnt_gc[2*ic+1] += dpElec[ic]
#   
# Current values to calculate deltaPapprch_c:  
#
	          dpApprch_cCrrnt[ic,k]=dpIon[ic]                  # g*cm/sec 
#
# dragging both paticles through second half of trajectory:
#
 	       z_elecCrrnt_gc=np.dot(matr_elec_c,z_elecCrrnt_gc)   # electron
 	       z_ionCrrnt_c=np.dot(matr_ion_c,z_ionCrrnt_c)        # ion
#	       crrntPoint=pointTrack_c[trackNumb_c]
#	       if iA == 0 and iB == 0 and crrntPoint < 10:
#                     print 'k, z_ionCrrnt_c', (k,z_ionCrrnt_c)
	       z_elecCrrnt_c=fromGuidingCenter(z_elecCrrnt_gc)     # transfer from system of guiding center 
# Current distance between ion and electron; cm:
 	       bCrrnt_c[k]=np.sqrt((z_ionCrrnt_c[0]-z_elecCrrnt_c[0])**2+ \
	                           (z_ionCrrnt_c[2]-z_elecCrrnt_c[2])**2+ \
			           (z_ionCrrnt_c[4]-z_elecCrrnt_c[4])**2)
# Current log10 of two important ratios:  
	       larmR_bCrrnt_c[k]=math.log10(rho_larm[iA,iB]/bCrrnt_c[k])      # dimensionless 
	       if maxLarmR_b_c < larmR_bCrrnt_c[k]:
	          maxLarmR_b_c=larmR_bCrrnt_c[k]
	       if minLarmR_b_c > larmR_bCrrnt_c[k]:
	          minLarmR_b_c=larmR_bCrrnt_c[k]
	       uPot_enrgKinCrrnt_c[k]=math.log10((q_elec**2/bCrrnt_c[k])/kinEnergy[iA,iB])
	       if maxUpot_enrgKin_c < uPot_enrgKinCrrnt_c[k]:
	          maxUpot_enrgKin_c=uPot_enrgKinCrrnt_c[k]
	       if minUpot_enrgKin_c > uPot_enrgKinCrrnt_c[k]:
	          minUpot_enrgKin_c=uPot_enrgKinCrrnt_c[k]
# To be prepared to draw future TMTdpx trajectory (for checking only):
               for ic in range(6):
                  prtclCoorCrrnt_c[ic,pointTrack_c[trackNumb_c]]=z_elecCrrnt_c[ic]  # 6-vector for electron
                  prtclCoorCrrnt_c[6+ic,pointTrack_c[trackNumb_c]]=z_ionCrrnt_c[ic] # 6-vector for ion 
	       prtclCoorCrrnt_c[12,pointTrack_c[trackNumb_c]]=bCrrnt_c[k]           # cm
	       prtclCoorCrrnt_c[13,pointTrack_c[trackNumb_c]]=action                # g*cm^2/sec
	       prtclCoorCrrnt_c[14,pointTrack_c[trackNumb_c]]=b_gc                  # cm
# To draw only first trajectory (for checking only):
               if trackNumb_c == 0:
                  for ic in range(6):
                     prtclCoorFirst_c[ic,pointTrack_c[trackNumb_c]]=z_elecCrrnt_c[ic]      # 6-vector for electron
                     prtclCoorFirst_c[6+ic,pointTrack_c[trackNumb_c]]=z_ionCrrnt_c[ic]     # 6-vector for ion 
	          prtclCoorFirst_c[12,pointTrack_c[trackNumb_c]]=bCrrnt_c[k]               # cm
	          prtclCoorFirst_c[13,pointTrack_c[trackNumb_c]]=action                    # g*cm^2/sec
	          prtclCoorFirst_c[14,pointTrack_c[trackNumb_c]]=b_gc                      # cm
	          if maxYcoorElec_c < abs(z_elecCrrnt_c[2]):
	             maxYcoorElec_c=abs(z_elecCrrnt_c[2])
	          if maxYcoorIon_c < abs(z_ionCrrnt_c[2]):
	             maxYcoorIon_c=abs(z_ionCrrnt_c[2])
               pointTrack_c[trackNumb_c] += 1
# End of dragging of the current trajectory	  
#-----------------------------------------------
# To draw the distribution of transferred dp (g*cm/sec):
               dpxTotal_c[iA,iB] +=-dpApprch_cCrrnt[0,k] 
               dpyTotal_c[iA,iB] +=-dpApprch_cCrrnt[1,k] 
               dpzTotal_c[iA,iB] +=-dpApprch_cCrrnt[2,k] 
#
# Total transferred dpx  for current track:
#
	       totAbsDpxApprch_c += abs(dpApprch_cCrrnt[0,k])
# 
# End of dragging of the current trajectory	  
#-----------------------------------------------
#
# Accumulate transferred momenta and other data: 
#
            if trackNumb_c == 0: 
# First definition of the total distance from origin of the
# coordinate system to electron along the trajectory; cm:
 	       b_c=bCrrnt_c                                        # cm
# First definition of the log10 of two important ratios; dimensionless:  
	       uPot_enrgKin_c=uPot_enrgKinCrrnt_c                  # ratio A 
	       if maxUpot_enrgKin_c < uPot_enrgKinCrrnt_c[k]:
	          maxUpot_enrgKin_c=uPot_enrgKinCrrnt_c[k]
	       if minUpot_enrgKin_c > uPot_enrgKinCrrnt_c[k]:
	          minUpot_enrgKin_c=uPot_enrgKinCrrnt_c[k]
	       larmR_b_c=larmR_bCrrnt_c                            # ratio B
	       if maxLarmR_b_c < larmR_bCrrnt_c[k]:
	          maxLarmR_b_c=larmR_bCrrnt_c[k]
	       if minLarmR_b_c > larmR_bCrrnt_c[k]:
	          minLarmR_b_c=larmR_bCrrnt_c[k]
#### First definition of the values to calculate deltaPapprch_c (g*cm/sec):
###               dpxApprch_c=dpApprch_cCrrnt[0,:] 
###               dpyApprch_c=dpApprch_cCrrnt[1,:]
###               dpzApprch_c=dpApprch_cCrrnt[2,:] 
            else:  
#### Total distance from origin of the coordinate system to electron along the trajectory:
 	       b_c=np.concatenate((b_c,bCrrnt_c),axis=0)           # cm
# Total log10 of two important ratios; dimensionless:  
	       uPot_enrgKin_c=np.concatenate((uPot_enrgKin_c,uPot_enrgKinCrrnt_c),axis=0)        
	       larmR_b_c=np.concatenate((larmR_b_c,larmR_bCrrnt_c),axis=0)                  
#### Total values to calculate deltaPapprch_c (g*cm/sec):
###               dpxApprch_c=np.concatenate((dpxApprch_c,dpApprch_cCrrnt[0,:]),axis=0)
###               dpyApprch_c=np.concatenate((dpyApprch_c,dpApprch_cCrrnt[1,:]),axis=0)
###               dpzApprch_c=np.concatenate((dpzApprch_c,dpApprch_cCrrnt[2,:]),axis=0)
#
# To draw TMTdpx trajectory (for checking only):
#
            if maxAbsDpxApprch_c < totAbsDpxApprch_c:
	       maxAbsDpxApprch_c=totAbsDpxApprch_c
	       indxAmaxAbsDpxApprch_c=iA
	       indxBmaxAbsDpxApprch_c=iB
	       trackNumbMaxAbsDpxApprch_c=trackNumb_c
	       prtclCoorMaxAbsDpx_c=prtclCoorCrrnt_c.copy()   
	       rhoMaxAbsDpxTurn_c=rhoCrrnt[iA,iB]
	       rhoLarmorMaxAbsDpxTurn_c=rho_larm[iA,iB]
###	       print 'iA=%d, iB=%d: TMT-track %d, points %d' % \
###	             (indxAmaxAbsDpxApprch_c,indxBmaxAbsDpxApprch_c, \
###		     trackNumbMaxAbsDpxApprch_c,pointTrack_c[trackNumbMaxAbsDpxApprch_c]) 
###	       print 'timePoints.shape: ', \
###	             (prtclCoorMaxAbsDpx_c.shape[0],prtclCoorMaxAbsDpx_c.shape[1])
# End of all calculations for approach_c
            sumPoints_c += pointTrack_c[trackNumb_c]
# 
# To draw the distribution of transferred dp:
#dpxTotalMax_c
	    if dpxTotalMax_c < abs(dpxTotal_c[iA,iB]):
	       dpxTotalMax_c=abs(dpxTotal_c[iA,iB])
	    if dpyTotalMax_c < abs(dpyTotal_c[iA,iB]):
	       dpyTotalMax_c=abs(dpyTotal_c[iA,iB])
	    if dpzTotalMax_c < abs(dpzTotal_c[iA,iB]):
	       dpzTotalMax_c=abs(dpzTotal_c[iA,iB])
#            print 'Track %d: dpxTotalMax_c=%e, dpyTotalMax_c=%e, dpzTotalMax_c=%e' % \
#	          (trackNumb_c,dpxTotalMax_c,dpyTotalMax_c,dpzTotalMax_c)
            timeEnd=os.times()
	    cpuTime_c[trackNumb_c]=1.e+6*(float(timeEnd[0])-float(timeStart[0]))  # CPU time , mks
            cpuTimeTotal += cpuTime_c[trackNumb_c] 
#
#------- End of approach_c --------------
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Approach_m: dragging with averaging over nLarmorAvrgng larmor rotation +
#             "Magnus expansion" method to calculate the transferred momenta
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#------- Start of calculations for approach_m --------------
         if runFlagApproach_m == 1:
            timeStart=os.times()
            if trackNumb_m == 0:
	       rhoFirstTurn=rhoCrrnt[iA,iB]
	       rhoLarmorFirstTurn=rho_larm[iA,iB]
#
# 6-vectors for ion and electron and distance 'b' between them for the first trajectory and 
# (T)rajectory with (M)aximal (T)ransferred (dpx) ( trajectory TMTdpx);
# (for checking only; indices 0-5 for electron, indices 6-11 for ion,
#  index=12 for 'b', index=13 for action and index=14 for dy_gc):
#
               prtclCoorFirst_m=np.zeros((15,timePoints_m[trackNumb_m]))       # first trajectory                 
            prtclCoorCrrnt_m=np.zeros((15,timePoints_m[trackNumb_m]))          # current trajectory                
            prtclCoorMaxAbsDpx_m=np.zeros((15,timePoints_m[trackNumb_m]))      # trajectory TMTdpx
# Current distance from origin of the coordinate system to electron along the trajectory; cm
            bCrrnt_m=np.zeros(timePoints_m[trackNumb_m])                       # cm
# Current log10 of two important ratios:
# First - ratio R_larmor/b; dimensionless
            larmR_bCrrnt_m=np.zeros(timePoints_m[trackNumb_m])       
# Second - ratio potential_energy/kinetic_energy; dimensionless
            uPot_enrgKinCrrnt_m=np.zeros(timePoints_m[trackNumb_m])  
# deltaPapprch_m=dpApprch_mCrrnt
            dpApprch_mCrrnt=np.zeros((3,timePoints_m[trackNumb_m]))
            for m in range(6): 
               z_ionCrrnt_m[m]=0.                     # Current initial zero-vector for ion
               z_elecCrrnt_m[m]=0.                    # Zeroing out of vector for electron
# Current initial vector for electron:
            z_elecCrrnt_m[Ix]=rhoCrrnt[iA,iB]+rho_larm[iA,iB]      # x, cm
            z_elecCrrnt_m[Iz]=-halfLintr[iA,iB]                    # z, cm
            z_elecCrrnt_m[Ipy]=m_elec*evTran[iA,iB]                # py, g*cm/sec
            z_elecCrrnt_m[Ipz]=m_elec*eVrmsLong                    # pz, g*cm/sec
	    z_elecCrrnt_gc=toGuidingCenter(z_elecCrrnt_m)          # transfer to system of guiding center
#	    if iA == 0 and iB == 0:
#	       print 'z_elecCrrnt_m: ', z_elecCrrnt_m 
#	       print 'z_elecCrrnt_gc: ', z_elecCrrnt_gc 
#-----------------------------------------------
# Main action - dragging of the current trajectories (for given i and j)
#
            for k in range(int(timePoints_m[trackNumb_m])):
#	       if k < 100: 
#	          print 'k=%d: x=%e, y=%e, z=%e' % \
#                      (k,1.e+4*z_elecCrrnt_gc[0],1.e+4*z_elecCrrnt_gc[2],1.e+4*z_elecCrrnt_gc[4]) 
#	       if k > timePoints_c[trackNumb_m]-100: 
#	          print 'k=%d: x=%e, y=%e, z=%e' % \
#                       (k,1.e+4*z_elecCrrnt_gc[0],1.e+4*z_elecCrrnt_gc[2],1.e+4*z_elecCrrnt_gc[4]) 
#
# dragging both paticles through first half of trajectory:
#
 	       z_elecCrrnt_gc=np.dot(matr_elec_m,z_elecCrrnt_gc)   # electron
 	       z_ionCrrnt_m=np.dot(matr_ion_m,z_ionCrrnt_m)        # ion
#
# dragging both paticles through interaction point:
#
	       dpIon,dpElec,action,dy_gc=MagnusExpansionCollision(z_elecCrrnt_gc,z_ionCrrnt_m,timeStep_m) 
#	       if trackNumb_m == 0:
#	          print 'point %d: dpgcElec=%e, dpzElec=%e' % 
#                       (pointTrack_m[0],dpElec[1],dpElec[2])
	       for ic in range(3):
####
#### Taking into account transfer of momentum for both particles:
####
####	          z_ionCrrnt_m[2*ic+1] += dpIon[ic]   
####	          z_elecCrrnt_gc[2*ic+1] += dpElec[ic]
#   
# Current values to calculate deltaPapprch_m:  
#
	          dpApprch_mCrrnt[ic,k]=dpIon[ic]                  # g*cm/sec
	       z_elecCrrnt_gc[2] += dy_gc                          # cm 
#
# dragging both paticles through second half of trajectory:
#
 	       z_elecCrrnt_gc=np.dot(matr_elec_m,z_elecCrrnt_gc)   # electron
 	       z_ionCrrnt_m=np.dot(matr_ion_m,z_ionCrrnt_m)        # ion
#	       crrntPoint=pointTrack_m[trackNumb_m]
#	       if iA == 0 and iB == 0 and crrntPoint < 10:
#                     print 'k, z_ionCrrnt_m', (k,z_ionCrrnt_m)
	       z_elecCrrnt_m=fromGuidingCenter(z_elecCrrnt_gc)     # transfer from system of guiding center 
# Current distance between ion and electron; cm:
 	       bCrrnt_m[k]=np.sqrt((z_ionCrrnt_m[0]-z_elecCrrnt_m[0])**2+ \
	                           (z_ionCrrnt_m[2]-z_elecCrrnt_m[2])**2+ \
			           (z_ionCrrnt_m[4]-z_elecCrrnt_m[4])**2)
# To be prepared to draw future TMTdpx trajectory (for checking only):
               for ic in range(6):
                  prtclCoorCrrnt_m[ic,pointTrack_m[trackNumb_m]]=z_elecCrrnt_m[ic]  # 6-vector for electron
                  prtclCoorCrrnt_m[6+ic,pointTrack_m[trackNumb_m]]=z_ionCrrnt_m[ic] # 6-vector for ion 
	       prtclCoorCrrnt_m[12,pointTrack_m[trackNumb_m]]=bCrrnt_m[k]           # cm
	       prtclCoorCrrnt_m[13,pointTrack_m[trackNumb_m]]=action                # g*cm^2/sec
	       prtclCoorCrrnt_m[14,pointTrack_m[trackNumb_m]]=dy_gc                 # cm
# To draw only first trajector3 (for checking only):
               if trackNumb_m == 0:
                  for ic in range(6):
                     prtclCoorFirst_m[ic,pointTrack_m[trackNumb_m]]=z_elecCrrnt_m[ic]      # 6-vector for electron
                     prtclCoorFirst_m[6+ic,pointTrack_m[trackNumb_m]]=z_ionCrrnt_m[ic]     # 6-vector for ion 
	          prtclCoorFirst_m[12,pointTrack_m[trackNumb_m]]=bCrrnt_m[k]               # cm
	          prtclCoorFirst_m[13,pointTrack_m[trackNumb_m]]=action                    # g*cm^2/sec
	          prtclCoorFirst_m[14,pointTrack_m[trackNumb_m]]=dy_gc                     # cm
	          if maxYcoorElec_m < abs(z_elecCrrnt_m[2]):
	             maxYcoorElec_m=abs(z_elecCrrnt_m[2])
	          if maxYcoorIon_m < abs(z_ionCrrnt_m[2]):
	             maxYcoorIon_m=abs(z_ionCrrnt_m[2])
               pointTrack_m[trackNumb_m] += 1
# End of dragging of the current trajectory	  
#-----------------------------------------------
# To draw the distribution of transferred dp (g*cm/sec):
               dpxTotal_m[iA,iB] +=-dpApprch_mCrrnt[0,k] 
               dpyTotal_m[iA,iB] +=-dpApprch_mCrrnt[1,k] 
               dpzTotal_m[iA,iB] +=-dpApprch_mCrrnt[2,k] 
#
# Total transferred dpx  for current track:
#
	       totAbsDpxApprch_m += abs(dpApprch_mCrrnt[0,k])
# 
# End of dragging of the current trajectory	  
#-----------------------------------------------
#
# Accumulate transferred momenta and other data: 
#
            if trackNumb_m == 0: 
# First definition of the total distance from origin of the
# coordinate system to electron along the trajectory; cm:
 	       b_m=bCrrnt_m                                        # cm
# First definition of the log10 of two important ratios; dimensionless:  
	       uPot_enrgKin_m=uPot_enrgKinCrrnt_m                  # ratio A 
	       if maxUpot_enrgKin_m < uPot_enrgKinCrrnt_m[k]:
	          maxUpot_enrgKin_m=uPot_enrgKinCrrnt_m[k]
	       if minUpot_enrgKin_m > uPot_enrgKinCrrnt_m[k]:
	          minUpot_enrgKin_m=uPot_enrgKinCrrnt_m[k]
	       larmR_b_m=larmR_bCrrnt_m                            # ratio B
	       if maxLarmR_b_m < larmR_bCrrnt_m[k]:
	          maxLarmR_b_m=larmR_bCrrnt_m[k]
	       if minLarmR_b_m > larmR_bCrrnt_m[k]:
	          minLarmR_b_m=larmR_bCrrnt_m[k]
#### First definition of the values to calculate deltaPapprch_m (g*cm/sec):
###               dpxApprch_m=dpApprch_mCrrnt[0,:] 
###               dpyApprch_m=dpApprch_mCrrnt[1,:]
###               dpzApprch_m=dpApprch_mCrrnt[2,:] 
            else:  
#### Total distance from origin of the coordinate system to electron along the trajectory:
 	       b_m=np.concatenate((b_m,bCrrnt_m),axis=0)           # cm
# Total log10 of two important ratios; dimensionless:  
	       uPot_enrgKin_m=np.concatenate((uPot_enrgKin_m,uPot_enrgKinCrrnt_m),axis=0)        
	       larmR_b_m=np.concatenate((larmR_b_m,larmR_bCrrnt_m),axis=0)                  
#### Total values to calculate deltaPapprch_m (g*cm/sec):
###               dpxApprch_m=np.concatenate((dpxApprch_m,dpApprch_mCrrnt[0,:]),axis=0)
###               dpyApprch_m=np.concatenate((dpyApprch_m,dpApprch_mCrrnt[1,:]),axis=0)
###               dpzApprch_m=np.concatenate((dpzApprch_m,dpApprch_mCrrnt[2,:]),axis=0)
#
# To draw TMTdpx trajectory (for checking only):
#
            if maxAbsDpxApprch_m < totAbsDpxApprch_m:
	       maxAbsDpxApprch_m=totAbsDpxApprch_m
	       indxAmaxAbsDpxApprch_m=iA
	       indxBmaxAbsDpxApprch_m=iB
	       trackNumbMaxAbsDpxApprch_m=trackNumb_m
	       prtclCoorMaxAbsDpx_m=prtclCoorCrrnt_m.copy()   
	       rhoMaxAbsDpxTurn_m=rhoCrrnt[iA,iB]
	       rhoLarmorMaxAbsDpxTurn_m=rho_larm[iA,iB]
###	       print 'iA=%d, iB=%d: TMT-track %d, points %d' % \
###	             (indxAmaxAbsDpxApprch_m,indxBmaxAbsDpxApprch_m, \
###		     trackNumbMaxAbsDpxApprch_m,pointTrack_m[trackNumbMaxAbsDpxApprch_m]) 
###	       print 'timePoints.shape: ', \
###	             (prtclCoorMaxAbsDpx_m.shape[0],prtclCoorMaxAbsDpx_m.shape[1])
# End of all calculations for approach_m
            sumPoints_m += pointTrack_m[trackNumb_m]
# 
# To draw the distribution of transferred dp:
#
	    if dpxTotalMax_m < abs(dpxTotal_m[iA,iB]):
	       dpxTotalMax_m=abs(dpxTotal_m[iA,iB])
	    if dpyTotalMax_m < abs(dpyTotal_m[iA,iB]):
	       dpyTotalMax_m=abs(dpyTotal_m[iA,iB])
	    if dpzTotalMax_m < abs(dpzTotal_m[iA,iB]):
	       dpzTotalMax_m=abs(dpzTotal_m[iA,iB])
            print 'Track %d: dpxTotalMax_m=%e, dpyTotalMax_m=%e, dpzTotalMax_m=%e' % \
	          (trackNumb_m,dpxTotalMax_m,dpyTotalMax_m,dpzTotalMax_m)
            timeEnd=os.times()
	    cpuTime_m[trackNumb_m]=1.e+6*(float(timeEnd[0])-float(timeStart[0]))  # CPU time , mks
            cpuTimeTotal += cpuTime_m[trackNumb_m] 
#
#------- End of approach_m --------------
#

#-- if runFlagApproach_1 == 1:
#--    print 'Approach_1: maxYcoorElec=%e mkm, maxYcoorIon=%e nm' % \
#--          (1.e+4*maxYcoorElec_1,1.e+7*maxYcoorIon_1)
#--    print 'Approach_1: for %d tracks number of points is %d' % (lastTrackNumber,sumPoints_1)


if runFlagApproach_c == 1:
   print 'Approach_c: for %d tracks number of points is %d' % (lastTrackNumber,sumPoints_c)
#   for i in range(b_cLenFirstTrack):
#      print 'action(%d)=%e' % (i,prtclCoor_c[13,i])

if runFlagApproach_m == 1:
   print 'Approach_m: for %d tracks number of points is %d' % (lastTrackNumber,sumPoints_m)

# for i in range(lastTrackNumber):
#    print 'Track %d: larmor turns=%d, cpuTime(mks)=%e, time per turn(mks)=%6.1f' % \
#          (i,larmorNumber[i],cpuTime[i],cpuTime[i]/larmorNumber[i])

print 'cpuTimeTotal(mksec) = %e' % cpuTimeTotal

#
# First track: to compare distances between particles for approach_1 and approach_c,3 (figure 325):
#

#-- if runFlagApproach_1 == 1:
#-- #   b_1Array=np.asarray(b_1)
#--    b_1LenFirstTrack=int(timePoints_1[0])
#--    print 'First track length for b: %d' % b_1LenFirstTrack


if runFlagApproach_c == 1:
#   b_cArray=np.asarray(b_c)
   b_cLenFirstTrack=int(timePoints_c[0])
   print 'First track length for b_c: %d' % b_cLenFirstTrack

if runFlagApproach_m == 1:
#   b_mArray=np.asarray(b_m)
   b_mLenFirstTrack=int(timePoints_m[0])
   print 'First track length for b_m: %d' % b_mLenFirstTrack

#
# Calculation of the difference for distance beteween electron and ion for approach_1 and approach_c:
#
'''
if runFlagApproach_1 == 1 and runFlagApproach_c == 1:

   diff_b2=np.zeros(b_cLenFirstTrack)
   print 'b_cLenFirstTrack=%d' % b_cLenFirstTrack

   k=0
   nStart=0
   for m in range(b_cLenFirstTrack): 
      for n in range(nStart,b_1LenFirstTrack):
         if prtclCoorFirst_1[4,n] == prtclCoorFirst_c[4,m]:
	       diff_b2[m]=prtclCoorFirst_1[12,n]-prtclCoorFirst_c[12,m]
               nStart=n-1
               break	 
         if prtclCoorFirst_1[4,n] > prtclCoorFirst_c[4,m]:
	    if n == 0:
	       diff_b2[m]=prtclCoorFirst_1[12,n]-prtclCoorFirst_c[12,m]
	    else:
	       bCurrent=prtclCoorFirst_1[12,n-1]+ \
	                (prtclCoorFirst_1[12,n]-prtclCoorFirst_1[12,n-1])* \
			(prtclCoorFirst_c[4,m]-prtclCoorFirst_1[4,n-1])/ \
			(prtclCoorFirst_1[4,n]-prtclCoorFirst_1[4,n-1])
               diff_b2[m]=bCurrent-prtclCoorFirst_c[14,m]
            nStart=n-1
            break	 
'''

###############################################
#
#                  Plotting
#
###############################################
#
# Checking of the first trajectory:
#
specFctr=1.e22                                                   # Use it in titles!
powerSpecFctr=22                                                 # Use it in titles!
X,Y=np.meshgrid(crrntA,crrntB) 

'''
if runFlagApproach_1 == 1 and plotFlagTracks == 1:
   pointsEndLarmor=int(larmorNumber[0])*stepsNumberOnGyro          # Number of points for drawing
#
# First electron's trajectory:
#
   turns=larmorNumber[0]                                           # Number of larmor turns for drawing 
   pointsTurns=turns*stepsNumberOnGyro                             # Number of points for drawing
   pointsStartLarmor=turns*stepsNumberOnGyro                       # Number of points for drawing
   lengthArrowElc=4
   pBegArrw=pointsTurns/2
   pEndArrw=pointsTurns/2+50

   fig40=plt.figure(40)
   ax40=fig40.gca(projection='3d')
   ax40.plot(1.e+4*prtclCoorFirst_1[0,0:pointsStartLarmor], \
             1.e+4*prtclCoorFirst_1[2,0:pointsStartLarmor], \
             1.e+4*prtclCoorFirst_1[4,0:pointsStartLarmor],'-r',linewidth=2)
   ax40.plot(1.e+4*prtclCoorFirst_1[0,0:lengthArrowElc],1.e+4*prtclCoorFirst_1[2,0:lengthArrowElc], \
             1.e+4*prtclCoorFirst_1[4,0:lengthArrowElc],'-b',linewidth=2)
   ax40.plot(1.e+4*prtclCoorFirst_1[0,pointsStartLarmor-lengthArrowElc:pointsStartLarmor], \
             1.e+4*prtclCoorFirst_1[2,pointsStartLarmor-lengthArrowElc:pointsStartLarmor], \
             1.e+4*prtclCoorFirst_1[4,pointsStartLarmor-lengthArrowElc:pointsStartLarmor], \
	     '-b',linewidth=2)
   plt.xlabel('x, $\mu m$',color='m',fontsize=16)
   plt.ylabel('y, $\mu m$',color='m',fontsize=16)
   ax40.set_zlabel('z, $\mu m$',color='m',fontsize=16)
   titleHeader="Approach-1. First Electron's Trajectory:"
   titleHeader += '\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m, $N_{Larm}=$%d'
   plt.title((titleHeader % (1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn,turns)), \
             color='m',fontsize=16)
#
# Trajectory of electron with maximal transferred px:
#
   turns=larmorNumber[trackNumbMaxAbsDpxApprch_1]                  # Number of larmor turns for drawing 
   pointsTurns=turns*stepsNumberOnGyro                             # Number of points for drawing
   pointsStartLarmor=turns*stepsNumberOnGyro                       # Number of points for drawing
   lengthArrowElc=4
   pBegArrw=pointsTurns/2
   pEndArrw=pointsTurns/2+50

   fig45=plt.figure(45)
   ax45=fig45.gca(projection='3d')
   ax45.plot(1.e+4*prtclCoorMaxAbsDpx_1[0,0:pointsStartLarmor], \
             1.e+4*prtclCoorMaxAbsDpx_1[2,0:pointsStartLarmor], \
             1.e+4*prtclCoorMaxAbsDpx_1[4,0:pointsStartLarmor],'-r',linewidth=2)
   ax45.plot(1.e+4*prtclCoorMaxAbsDpx_1[0,0:lengthArrowElc], \
             1.e+4*prtclCoorMaxAbsDpx_1[2,0:lengthArrowElc], \
             1.e+4*prtclCoorMaxAbsDpx_1[4,0:lengthArrowElc],'-b',linewidth=2)
   ax45.plot(1.e+4*prtclCoorMaxAbsDpx_1[0,pointsStartLarmor-lengthArrowElc:pointsStartLarmor], \
             1.e+4*prtclCoorMaxAbsDpx_1[2,pointsStartLarmor-lengthArrowElc:pointsStartLarmor], \
             1.e+4*prtclCoorMaxAbsDpx_1[4,pointsStartLarmor-lengthArrowElc:pointsStartLarmor], \
	     '-b',linewidth=2)
   plt.xlabel('x, $\mu m$',color='m',fontsize=16)
   plt.ylabel('y, $\mu m$',color='m',fontsize=16)
   ax45.set_zlabel('z, $\mu m$',color='m',fontsize=16)
   titleHeader='Approach-1. Trajectory of Electron with Maximal Transferred $p_x$:'
   titleHeader += '\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m, $N_{Larm}=$%d'
   plt.title((titleHeader % (1.e+4*rhoMaxAbsDpxTurn_1,1.e+4*rhoLarmorMaxAbsDpxTurn_1,turns)), \
             color='m',fontsize=16)
#
# First electron's trajectory (distance to ion):
#
   pointsTurns=pointTrack_1[0]
   plt.figure (50)
   plt.plot(1.e+4*prtclCoorFirst_1[4,0:pointsTurns],1.e+4*prtclCoorFirst_1[12,0:pointsTurns],'-xr')
   plt.xlabel('z, $\mu$m',color='m',fontsize=16)
   plt.ylabel('Distance to ion, $\mu$m',color='m',fontsize=16)
   plt.title("Approach-1: First Electron's Track", color='m',fontsize=16)
   plt.grid(True)
#
# First electron's trajectory (Start):
#
   turns=10                                                        # Number of larmor turns for drawing 
   pointsTurns=turns*stepsNumberOnGyro                             # Number of points for drawing
   pointsStartLarmor=turns*stepsNumberOnGyro                       # Number of points for drawing
   lengthArrowElc=4
   pBegArrw=pointsTurns/2
   pEndArrw=pointsTurns/2+50

   fig50=plt.figure(50)
   ax50=fig50.gca(projection='3d')
   ax50.plot(1.e+4*prtclCoorFirst_1[0,0:pointsStartLarmor], \
             1.e+4*prtclCoorFirst_1[2,0:pointsStartLarmor], \
             1.e+4*prtclCoorFirst_1[4,0:pointsStartLarmor],'-r',linewidth=2)
   ax50.plot(1.e+4*prtclCoorFirst_1[0,0:lengthArrowElc],1.e+4*prtclCoorFirst_1[2,0:lengthArrowElc], \
             1.e+4*prtclCoorFirst_1[4,0:lengthArrowElc],'-b',linewidth=2)
   ax50.plot(1.e+4*prtclCoorFirst_1[0,pointsStartLarmor-lengthArrowElc:pointsStartLarmor], \
             1.e+4*prtclCoorFirst_1[2,pointsStartLarmor-lengthArrowElc:pointsStartLarmor], \
             1.e+4*prtclCoorFirst_1[4,pointsStartLarmor-lengthArrowElc:pointsStartLarmor], \
	     '-b',linewidth=2)
   plt.xlabel('x, $\mu m$',color='m',fontsize=16)
   plt.ylabel('y, $\mu m$',color='m',fontsize=16)
   ax50.set_zlabel('z, $\mu m$',color='m',fontsize=16)
   titleHeader="Approach-1. First Electron's Trajectory (Start; $N_{Larm}=$%d):"
   titleHeader += '\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m'
   plt.title((titleHeader % (larmorNumber[0],1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)), \
             color='m',fontsize=16)
#
# First electron's trajectory (End):
#
   pBeg=pointsEndLarmor-pointsTurns

   fig60=plt.figure(60)
   ax60=fig60.gca(projection='3d')
   ax60.plot(1.e+4*prtclCoorFirst_1[0,pBeg:pointsEndLarmor], \
             1.e+4*prtclCoorFirst_1[2,pBeg:pointsEndLarmor], \
             1.e+4*prtclCoorFirst_1[4,pBeg:pointsEndLarmor],'-r',linewidth=2)
   ax60.plot(1.e+4*prtclCoorFirst_1[0,pointsEndLarmor-lengthArrowElc:pointsEndLarmor], \
             1.e+4*prtclCoorFirst_1[2,pointsEndLarmor-lengthArrowElc:pointsEndLarmor], \
             1.e+4*prtclCoorFirst_1[4,pointsEndLarmor-lengthArrowElc:pointsEndLarmor], \
	     '-b',linewidth=2)
   ax60.plot(1.e+4*prtclCoorFirst_1[0,pBeg:pBeg+lengthArrowElc], \
             1.e+4*prtclCoorFirst_1[2,pBeg:pBeg+lengthArrowElc], \
             1.e+4*prtclCoorFirst_1[4,pBeg:pBeg+lengthArrowElc],'-b',linewidth=2)
   plt.xlabel('x, $\mu m$',color='m',fontsize=16)
   plt.ylabel('y, $\mu m$',color='m',fontsize=16)
   ax60.set_zlabel('z, $\mu m$',color='m',fontsize=16)
   titleHeader="Approach-1. First Electron's Trajectory (End; $N_{Larm}=$%d):"
   titleHeader += '\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m'
   plt.title((titleHeader % (larmorNumber[0],1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)), \
             color='m',fontsize=16)
#
# First ion trajectory (Start):
#
   fig70=plt.figure(70)
   ax70=fig70.gca(projection='3d')
   ax70.plot(1.e+7*prtclCoorFirst_1[6,0:pointsStartLarmor],1.e+7*prtclCoorFirst_1[8,0:pointsStartLarmor], \
             1.e+7*prtclCoorFirst_1[10,0:pointsStartLarmor],'-b',linewidth=2)
   ax70.plot(1.e+7*prtclCoorFirst_1[6,pBegArrw:pEndArrw],1.e+7*prtclCoorFirst_1[8,pBegArrw:pEndArrw], \
             1.e+7*prtclCoorFirst_1[10,pBegArrw:pEndArrw],'-b',linewidth=4)
#     ax70.plot([arrwBegxIon,arrwEndxIon],[arrwBegyIon,arrwEndyIon], \
#               [arrwBegzIon,arrwEndzIon],color='r',alpha=0.8,lw=2)
   plt.xlabel('x, $nm$',color='m',fontsize=16)
   plt.ylabel('y, $nm$',color='m',fontsize=16)
   ax70.set_zlabel('z, $nm$',color='m',fontsize=16)
   plt.title('Approach-1. First Ion Trajectory (Start)',color='m',fontsize=16)

#   arrayRatioA=np.asarray(uPot_enrgKin_1)
#   print ('Length(arrayRatioA)=%d' % len(arrayRatioA))

if runFlagApproach_1 == 1 and plotFlagWorkingArrays ==1:
#
# Area of parameters A and B:
#
   plt.figure (110)
   plt.plot(uPot_enrgKin_1[0:sumPoints_1],larmR_b_1[0:sumPoints_1],'.r')
#   plt.plot(xAboundary,xBboundary,'-xb')
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   plt.title(('Approach-1: Map of Current Values of Parameters $A$ and $B$ (%d Tracks)' \
              % lastTrackNumber), color='m',fontsize=15)
   plt.text(-4.5,-3.75,'$b=|r_i-r_e|^{1/2}$',color='m',fontsize=30)
# plt.xlim([minA,maxA])
# plt.ylim([minB,maxB])
   plt.grid(True)

if runFlagApproach_c == 1 and plotFlagWorkingArrays ==1:
   plt.figure (115)
   plt.plot(uPot_enrgKin_c[0:sumPoints_c],larmR_b_c[0:sumPoints_c],'.r')
#   plt.plot(xAboundary,xBboundary,'-xb')
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   plt.title(('Approach-c: Map of Current Values of Parameters $A$ and $B$ (%d Tracks)' \
              % lastTrackNumber), color='m',fontsize=15)
   plt.text(-4.5,-3.75,'$b=|r_i-r_e|^{1/2}$',color='m',fontsize=30)
# plt.xlim([minA,maxA])
# plt.ylim([minB,maxB])
   plt.grid(True)

if runFlagApproach_1 == 1 and plotFlagDpTransf == 1:              
   specFctr=1.e22                                                   # Use it in titles!
   powerSpecFctr=22                                                 # Use it in titles!
   X,Y=np.meshgrid(crrntA,crrntB) 
#
# Transfered momentum px (surface):
#
   fig240=plt.figure(240)
   ax240=fig240.gca(projection='3d')
#    surf=ax240.plot_surface(X,Y,specFctr*dpxTotal_1,cmap=cm.coolwarm, \
#                            linewidth=0,antialiased=False)
   surf=ax240.plot_surface(X,Y,specFctr*dpxTotal_1,cmap=cm.jet,linewidth=0,antialiased=False)
#   titleHeader='Approach-1: Transfered Momentum $dP_x$ $(\cdot 10^{22}$; $g \cdot cm/sec)$'
#   titleHeader += '\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-22}$ $g \cdot cm/sec$)'
#   plt.title((titleHeader % (lastTrackNumber,specFctr*dpxTotalMax_1)), color='m',fontsize=16)
   titleHeader='Approach-1: Transfered Momentum $dP_x$ $(\cdot 10^{%2d}$; $g \cdot cm/sec)$'
   titleHeader +='\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-%2d}$ $g \cdot cm/sec)$)'
   plt.title((titleHeader % \
              (powerSpecFctr,lastTrackNumber,specFctr*dpxTotalMax_1,powerSpecFctr)), \
	     color='m',fontsize=14)
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   ax240.set_zlabel('$dP_x \cdot 10^{22}$; $g \cdot cm/sec$',color='m',fontsize=16)
   cb = fig240.colorbar(surf)
   plt.grid(True)
#
# Transfered momentum px (map):
#
   fig245=plt.figure(245)
   ax=fig245.add_subplot(111)                                       # for contours plotting
   mapDpx=ax.contourf(X,Y,specFctr*dpxTotal_1)   
   plt.plot(xAboundary[0:81],xBboundary[0:81],'-w')
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
#   titleHeader='Approach-1: Transfered Momentum $dP_x$ $(\cdot 10^{22}$; $g \cdot cm/sec)$'
#   titleHeader +='\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-22}$ $g \cdot cm/sec)$)'
#   plt.title((titleHeader % (lastTrackNumber,specFctr*dpxTotalMax_1)), color='m',fontsize=14)
   titleHeader='Approach-1: Transfered Momentum $dP_x$ $(\cdot 10^{%2d}$; $g \cdot cm/sec)$'
   titleHeader +='\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-%2d}$ $g \cdot cm/sec)$)'
   plt.title((titleHeader % \
              (powerSpecFctr,lastTrackNumber,specFctr*dpxTotalMax_1,powerSpecFctr)), \
	     color='m',fontsize=14)
#   plt.text(-4.9,-2.95,'Unallowable Tracks: Initial\nImpact Parameter > $R_{shield}$',color='w',fontsize=20)
#   plt.text(-4.65,-2.25,'Unallowable\nTracks: Initial\nImpact Parameter\nLarger than $R_{shield}$',color='w',fontsize=20)
   plt.text(-4.9,-2.5, \
            '        Area with\nUnallowable Tracks:\n    Initial Impact\nParameter Larger\n     than $R_{shield}$',\
            color='w',fontsize=20)
   plt.ylim([-3.,-.5])
   fig245.colorbar(mapDpx)
#
# Transfered momentum py (surface):
#
   fig250=plt.figure(250)
   ax250=fig250.gca(projection='3d')
#    surf=ax250.plot_surface(X,Y,specFctr*dpyTotal_1,cmap=cm.coolwarm, \
#                            linewidth=0,antialiased=False)
   surf=ax250.plot_surface(X,Y,specFctr*dpyTotal_1,cmap=cm.jet,linewidth=0,antialiased=False)
#   titleHeader='Approach-1: Transfered Momentum $dP_y$ $(\cdot 10^{22}$; $g \cdot cm/sec)$'
#   titleHeader += '\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-22}$ $g \cdot cm/sec$)'
#   plt.title((titleHeader % (lastTrackNumber,specFctr*dpyTotalMax_1)), color='m',fontsize=16)
   titleHeader='Approach-1: Transfered Momentum $dP_y$ $(\cdot 10^{%2d}$; $g \cdot cm/sec)$'
   titleHeader +='\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-%2d}$ $g \cdot cm/sec)$)'
   plt.title((titleHeader % \
              (powerSpecFctr,lastTrackNumber,specFctr*dpyTotalMax_1,powerSpecFctr)), \
	     color='m',fontsize=14)
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   ax250.set_zlabel('$dP_y \cdot 10^{22}$; $g \cdot cm/sec$',color='m',fontsize=16)
   cb = fig250.colorbar(surf)
   plt.grid(True)
#
# Transfered momentum py (map):
#
   fig255=plt.figure(255)
   ax=fig255.add_subplot(111)                                       # for contours plotting
   mapDpx=ax.contourf(X,Y,specFctr*dpyTotal_1)   
   plt.plot(xAboundary[0:81],xBboundary[0:81],'-w')
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
#   titleHeader='Approach-1: Transfered Momentum $dP_y$ $(\cdot 10^{22}$; $g \cdot cm/sec)$'
#   titleHeader +='\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-22}$ $g \cdot cm/sec)$)'
#   plt.title((titleHeader % (lastTrackNumber,specFctr*dpyTotalMax_1)), color='m',fontsize=14)
   titleHeader='Approach-1: Transfered Momentum $dP_y$ $(\cdot 10^{%2d}$; $g \cdot cm/sec)$'
   titleHeader +='\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-%2d}$ $g \cdot cm/sec)$)'
   plt.title((titleHeader % \
              (powerSpecFctr,lastTrackNumber,specFctr*dpyTotalMax_1,powerSpecFctr)), \
	     color='m',fontsize=14)
   plt.ylim([-3.,-.5])
#   plt.text(-4.9,-2.95,'Unallowable Tracks: Initial\nImpact Parameter > $R_{shield}$',color='w',fontsize=20)
#   plt.text(-4.65,-2.25,'Unallowable\nTracks: Initial\nImpact Parameter\nLarger than $R_{shield}$',color='w',fontsize=20)
   plt.text(-4.9,-2.5, \
            '        Area with\nUnallowable Tracks:\n    Initial Impact\nParameter Larger\n     than $R_{shield}$',\
            color='w',fontsize=20)
   fig255.colorbar(mapDpx)
#
# Transfered momentum pz (surface):
#
   fig260=plt.figure(260)
   ax260=fig260.gca(projection='3d')
#    surf=ax260.plot_surface(X,Y,specFctr*dpzTotal_1,cmap=cm.coolwarm, \
#                            linewidth=0,antialiased=False)
   surf=ax260.plot_surface(X,Y,specFctr*dpzTotal_1,cmap=cm.jet,linewidth=0,antialiased=False)
#   titleHeader='Approach-1: Transfered Momentum $dP_z$ $(\cdot 10^{22}$; $g \cdot cm/sec)$'
#   titleHeader += '\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-22}$ $g \cdot cm/sec$)'
#   plt.title((titleHeader % (lastTrackNumber,specFctr*dpzTotalMax_1)), color='m',fontsize=16)
   titleHeader='Approach-1: Transfered Momentum $dP_z$ $(\cdot 10^{%2d}$; $g \cdot cm/sec)$'
   titleHeader +='\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-%2d}$ $g \cdot cm/sec)$)'
   plt.title((titleHeader % \
              (powerSpecFctr,lastTrackNumber,specFctr*dpzTotalMax_1,powerSpecFctr)), \
	     color='m',fontsize=14)
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   ax260.set_zlabel('$dP_z \cdot 10^{22}$; $g \cdot cm/sec$',color='m',fontsize=16)
   cb = fig260.colorbar(surf)
   plt.grid(True)
#
# Transfered momentum pz (map):
#
   fig265=plt.figure(265)
   ax=fig265.add_subplot(111)                                       # for contours plotting
   mapDpx=ax.contourf(X,Y,specFctr*dpzTotal_1)   
   plt.plot(xAboundary[0:81],xBboundary[0:81],'-w')
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
#   titleHeader='Approach-1: Transfered Momentum $dP_z$ $(\cdot 10^{22}$; $g \cdot cm/sec)$'
#   titleHeader +='\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-22}$ $g \cdot cm/sec)$)'
#   plt.title((titleHeader % (lastTrackNumber,specFctr*dpzTotalMax_1)), color='m',fontsize=14)
   titleHeader='Approach-1: Transfered Momentum $dP_z$ $(\cdot 10^{%2d}$; $g \cdot cm/sec)$'
   titleHeader +='\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-%2d}$ $g \cdot cm/sec)$)'
   plt.title((titleHeader % \
              (powerSpecFctr,lastTrackNumber,specFctr*dpzTotalMax_1,powerSpecFctr)), \
	     color='m',fontsize=14)
   plt.ylim([-3.,-.5])
#   plt.text(-4.9,-2.95,'Unallowable Tracks: Initial\nImpact Parameter > $R_{shield}$',color='w',fontsize=20)
#   plt.text(-4.65,-2.25,'Unallowable\nTracks: Initial\nImpact Parameter\nLarger than $R_{shield}$',color='w',fontsize=20)
   plt.text(-4.9,-2.5, \
            '        Area with\nUnallowable Tracks:\n    Initial Impact\nParameter Larger\n     than $R_{shield}$',\
            color='w',fontsize=20)
   fig265.colorbar(mapDpx)
'''

if runFlagApproach_c == 1 and plotFlagTracks == 1:

#
# First electron's trajectory:
#
   turns=larmorNumber[0]                                           # Number of larmor turns for drawing 
   pointsTurns=turns                                               # Number of points for drawing
   lengthArrowElc=5
   pBegArrw=turns/2-lengthArrowElc/2
   pEndArrw=turns/2+lengthArrowElc/2

   fig140=plt.figure(140)
   ax140=fig140.gca(projection='3d')
   ax140.plot(1.e+4*prtclCoorFirst_c[0,0:pointsTurns], \
              1.e+4*prtclCoorFirst_c[2,0:pointsTurns], \
              1.e+4*prtclCoorFirst_c[4,0:pointsTurns],'-r',linewidth=2)
   ax140.plot(1.e+4*prtclCoorFirst_c[0,pBegArrw:pEndArrw], \
              1.e+4*prtclCoorFirst_c[2,pBegArrw:pEndArrw], \
              1.e+4*prtclCoorFirst_c[4,pBegArrw:pEndArrw],'-b',linewidth=4)
   plt.xlabel('x, $\mu m$',color='m',fontsize=16)
   plt.ylabel('y, $\mu m$',color='m',fontsize=16)
   ax140.set_zlabel('z, $\mu m$',color='m',fontsize=16)
   titleHeader="Approach-c. First Electron's Trajectory:"
   titleHeader += '\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m, $N_{Larm}=$%d'
   plt.title((titleHeader % (1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn,turns)), \
             color='m',fontsize=16)
#
# Trajectory of electron with maximal transferred px:
#
   turns=larmorNumber[trackNumbMaxAbsDpxApprch_c]                  # Number of larmor turns for drawing 
   pointsTurns=turns                                               # Number of points for drawing
   lengthArrowElc=8
   pBegArrw=turns/2-lengthArrowElc/2
   pEndArrw=turns/2+lengthArrowElc/2

   fig145=plt.figure(145)
   ax145=fig145.gca(projection='3d')
   ax145.plot(1.e+4*prtclCoorMaxAbsDpx_c[0,0:pointsTurns], \
              1.e+4*prtclCoorMaxAbsDpx_c[2,0:pointsTurns], \
              1.e+4*prtclCoorMaxAbsDpx_c[4,0:pointsTurns],'-r',linewidth=2)
   ax145.plot(1.e+4*prtclCoorMaxAbsDpx_c[0,pBegArrw:pEndArrw], \
              1.e+4*prtclCoorMaxAbsDpx_c[2,pBegArrw:pEndArrw], \
              1.e+4*prtclCoorMaxAbsDpx_c[4,pBegArrw:pEndArrw],'-b',linewidth=4)
   plt.xlabel('x, $\mu m$',color='m',fontsize=16)
   plt.ylabel('y, $\mu m$',color='m',fontsize=16)
   ax145.set_zlabel('z, $\mu m$',color='m',fontsize=16)
   titleHeader='Approach-c. Trajectory of Electron with Maximal Transferred $p_x$:'
   titleHeader += '\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m, $N_{Larm}=$%d'
   plt.title((titleHeader % (1.e+4*rhoMaxAbsDpxTurn_c,1.e+4*rhoLarmorMaxAbsDpxTurn_c,turns)), \
             color='m',fontsize=16)
#
# First electron's trajectory (action):
#
   specFctr=1.e31                                                  #       Use it in titles!
   specShift=5.86748911+1.478799548                                            # x e8! Use it in titles!

   plt.figure (160)
   plt.plot(1.e+4*prtclCoorFirst_c[4,0:pointsTurns], \
            specFctr*prtclCoorFirst_c[13,0:pointsTurns]-1.e8*specShift,'-xr',linewidth=2)
   plt.xlabel('z, $\mu$m',color='m',fontsize=16)
   plt.ylabel(('$10^{31}\cdot J$ $-$ %5.3f$\cdot10^8$, g$\cdot$cm$^2$/sec' % specShift), \
              color='m',fontsize=16)
   plt.title("Approach-c: Action $J$ (First Electron's Track)", color='m',fontsize=16)
   plt.grid(True)

'''
if runFlagApproach_c == 1 and plotFlagTracks == 1:
   fig325=plt.figure(325)
   plt.plot(1.e+4*prtclCoorFirst_c[4,0:b_cLenFirstTrack],1.e+4*diff_b2[0:b_cLenFirstTrack],'-xr',linewidth=2)
   plt.xlabel('z, $\mu m$',color='m',fontsize=16)
   plt.ylabel('$\Delta b=b_{Apprch-1}-b_{Apprch-c}$, $\mu m$',color='m',fontsize=16)
   plt.title('Approach-c: First Trajectory: $\Delta b$ (Difference for Distances $b$ between Particles)',color='m',fontsize=16)
#   plt.xlim([1.e+4*prtclCoorFirst_c[4,0]-50,1.e+4*prtclCoorFirst_c[4,b_cLenFirstTrack-1]+50])
   plt.ylim([-57.,-52.])
   plt.grid(True)
'''

if runFlagApproach_c == 1 and plotFlagDpTransf == 1:              

   specFctr=1.e22                                                  # Use it in titles!
   powerpecFctr=22                                                 # Use it in titles!
   X,Y=np.meshgrid(crrntA,crrntB) 
#
# Transfered momentum px (surface):
#
   fig340=plt.figure(340)
   ax340=fig340.gca(projection='3d')
#    surf=ax340.plot_surface(X,Y,specFctr*dpxTotal_c,cmap=cm.coolwarm, \
#                            linewidth=0,antialiased=False)
   surf=ax340.plot_surface(X,Y,specFctr*dpxTotal_c,cmap=cm.jet,linewidth=0,antialiased=False)
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   ax340.set_zlabel('$dP_x \cdot 10^{22}$; $g \cdot cm/sec$',color='m',fontsize=16)
   titleHeader='Approach-c: Transfered Momentum $dP_x$ $(\cdot 10^{22}$; $g \cdot cm/sec)$'
   titleHeader += '\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-22}$ $g \cdot cm/sec$)'
   plt.title((titleHeader % (lastTrackNumber,specFctr*dpxTotalMax_c)), color='m',fontsize=16)
   cb = fig340.colorbar(surf)
   plt.grid(True)
#
# Transfered momentum px (map):
#
   fig345=plt.figure(345)
   ax=fig345.add_subplot(111)                                       # for contours plotting
   mapDpx=ax.contourf(X,Y,specFctr*dpxTotal_c)   
   plt.plot(xAboundary[0:81],xBboundary[0:81],'-w')
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
#   titleHeader='Approach-c: Transfered Momentum $dP_x$ $(\cdot 10^{22}$; $g \cdot cm/sec)$'
#   titleHeader +='\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-22}$ $g \cdot cm/sec)$)'
#   plt.title((titleHeader % (lastTrackNumber,specFctr*dpxTotalMax_c)), color='m',fontsize=14)
   titleHeader='Approach-c: Transfered Momentum $dP_x$ $(\cdot 10^{%2d}$; $g \cdot cm/sec)$'
   titleHeader +='\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-%2d}$ $g \cdot cm/sec)$)'
   plt.title((titleHeader % \
              (powerSpecFctr,lastTrackNumber,specFctr*dpxTotalMax_c,powerSpecFctr)), \
	     color='m',fontsize=14)
#   plt.text(-4.9,-2.95,'Unallowable Tracks: Initial\nImpact Parameter > $R_{shield}$',color='w',fontsize=20)
#   plt.text(-4.65,-2.25,'Unallowable\nTracks: Initial\nImpact Parameter\nLarger than $R_{shield}$',color='w',fontsize=20)
   plt.text(-4.9,-2.5, \
            '        Area with\nUnallowable Tracks:\n    Initial Impact\nParameter Larger\n     than $R_{shield}$',\
            color='w',fontsize=20)
   plt.ylim([-3.,-.5])
   fig345.colorbar(mapDpx)

#
# Transfered momentum py (surface):
#
###   fig350=plt.figure(350)
###   ax350=fig350.gca(projection='3d')
####    surf=ax350.plot_surface(X,Y,specFctr*dpyTotal_c,cmap=cm.coolwarm, \
####                            linewidth=0,antialiased=False)
###   surf=ax350.plot_surface(X,Y,specFctr*dpyTotal_c,cmap=cm.jet,linewidth=0,antialiased=False)
###   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
###   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
###   ax350.set_zlabel('$dP_y \cdot 10^{22}$; $g \cdot cm/sec$',color='m',fontsize=16)
###   titleHeader='Approach-c: Transfered Momentum $dP_y$ $(\cdot 10^{22}$; $g \cdot cm/sec)$'
###   titleHeader += '\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-22}$ $g \cdot cm/sec$)'
###   plt.title((titleHeader % (lastTrackNumber,specFctr*dpyTotalMax_c)), color='m',fontsize=16)
###   cb = fig350.colorbar(surf)
###   plt.grid(True)
####
#### Transfered momentum py (map):
####
###   fig355=plt.figure(355)
###   ax=fig355.add_subplot(111)                                       # for contours plotting
###   mapDpy=ax.contourf(X,Y,specFctr*dpyTotal_c)   
###   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
###   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
###   titleHeader='Approach-c: Transfered Momentum $dP_y$ $(\cdot 10^{22}$; $g \cdot cm/sec)$'
###   titleHeader +='\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-22}$ $g \cdot cm/sec)$)'
###   plt.title((titleHeader % (lastTrackNumber,specFctr*dpyTotalMax_c)), color='m',fontsize=14)
###   fig355.colorbar(mapDpx)

#
# Transfered momentum pz (surface):
#
   fig360=plt.figure(360)
   ax360=fig360.gca(projection='3d')
#    surf=ax360.plot_surface(X,Y,specFctr*dpzTotal_c,cmap=cm.coolwarm, \
#                            linewidth=0,antialiased=False)
   surf=ax360.plot_surface(X,Y,specFctr*dpzTotal_c,cmap=cm.jet,linewidth=0,antialiased=False)
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   ax360.set_zlabel('$dP_z \cdot 10^{22}$; $g \cdot cm/sec$',color='m',fontsize=16)
#   titleHeader='Approach-c: Transfered Momentum $dP_z$ $(\cdot 10^{22}$; $g \cdot cm/sec)$'
#   titleHeader += '\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-22}$ $g \cdot cm/sec$)'
#   plt.title((titleHeader % (lastTrackNumber,specFctr*dpzTotalMax_c)), color='m',fontsize=16)
   titleHeader='Approach-c: Transfered Momentum $dP_z$ $(\cdot 10^{%2d}$; $g \cdot cm/sec)$'
   titleHeader +='\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-%2d}$ $g \cdot cm/sec)$)'
   plt.title((titleHeader % \
              (powerSpecFctr,lastTrackNumber,specFctr*dpzTotalMax_c,powerSpecFctr)), \
	     color='m',fontsize=14)
   cb = fig360.colorbar(surf)
   plt.grid(True)
#
# Transfered momentum pz (map):
#
   fig365=plt.figure(365)
   ax=fig365.add_subplot(111)                                       # for contours plotting
   mapDpz=ax.contourf(X,Y,specFctr*dpzTotal_c)   
   plt.plot(xAboundary[0:81],xBboundary[0:81],'-w')
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
#   titleHeader='Approach-c: Transfered Momentum $dP_z$ $(\cdot 10^{22}$; $g \cdot cm/sec)$'
#   titleHeader +='\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-22}$ $g \cdot cm/sec)$)'
#   plt.title((titleHeader % (lastTrackNumber,specFctr*dpzTotalMax_c)), color='m',fontsize=14)
   titleHeader='Approach-c: Transfered Momentum $dP_z$ $(\cdot 10^{%2d}$; $g \cdot cm/sec)$'
   titleHeader +='\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-%2d}$ $g \cdot cm/sec)$)'
   plt.title((titleHeader % \
              (powerSpecFctr,lastTrackNumber,specFctr*dpxTotalMax_c,powerSpecFctr)), \
	     color='m',fontsize=14)
#   plt.text(-4.9,-2.95,'Unallowable Tracks: Initial\nImpact Parameter > $R_{shield}$',color='w',fontsize=20)
#   plt.text(-4.65,-2.25,'Unallowable\nTracks: Initial\nImpact Parameter\nLarger than $R_{shield}$',color='w',fontsize=20)
   plt.text(-4.9,-2.5, \
            '        Area with\nUnallowable Tracks:\n    Initial Impact\nParameter Larger\n     than $R_{shield}$',\
            color='w',fontsize=20)
   plt.ylim([-3.,-.5])
   fig365.colorbar(mapDpx)

'''
if runFlagApproach_1 == 1 and runFlagApproach_c == 1 and plotFlagTracks == 1:

#
# First electron's trajectory:
#
   turns=larmorNumber[0]                                           # Number of larmor turns for drawing 
   pointsTurns=turns*stepsNumberOnGyro                             # Number of points for drawing
   pointsStartLarmor=turns*stepsNumberOnGyro                       # Number of points for drawing
   lengthArrowElc=4
   pBegArrw=turns/2-lengthArrowElc/2
   pEndArrw=turns/2+lengthArrowElc/2

   fig440=plt.figure(440)
   ax440=fig440.gca(projection='3d')
   ax440.plot(1.e+4*prtclCoorFirst_1[0,0:pointsStartLarmor], \
              1.e+4*prtclCoorFirst_1[2,0:pointsStartLarmor], \
              1.e+4*prtclCoorFirst_1[4,0:pointsStartLarmor],'-r',linewidth=2)
   ax440.plot(1.e+4*prtclCoorFirst_1[0,0:lengthArrowElc], \
              1.e+4*prtclCoorFirst_1[2,0:lengthArrowElc], \
              1.e+4*prtclCoorFirst_1[4,0:lengthArrowElc],'-r',linewidth=4)
   ax440.plot(1.e+4*prtclCoorFirst_1[0,pointsStartLarmor-lengthArrowElc:pointsStartLarmor], \
              1.e+4*prtclCoorFirst_1[2,pointsStartLarmor-lengthArrowElc:pointsStartLarmor], \
              1.e+4*prtclCoorFirst_1[4,pointsStartLarmor-lengthArrowElc:pointsStartLarmor], \
	      '-r',linewidth=4)
   ax440.plot(1.e+4*prtclCoorFirst_c[0,0:turns], \
              1.e+4*prtclCoorFirst_c[2,0:turns], \
              1.e+4*prtclCoorFirst_c[4,0:turns],'-b',linewidth=2)
   ax440.plot(1.e+4*prtclCoorFirst_c[0,pBegArrw:pEndArrw], \
              1.e+4*prtclCoorFirst_c[2,pBegArrw:pEndArrw], \
              1.e+4*prtclCoorFirst_c[4,pBegArrw:pEndArrw],'-b',linewidth=4)
   plt.xlabel('x, $\mu m$',color='m',fontsize=16)
   plt.ylabel('y, $\mu m$',color='m',fontsize=16)
   ax440.set_zlabel('z, $\mu m$',color='m',fontsize=16)
   titleHeader="First Electron's Trajectory:"
   titleHeader += '\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m, $N_{Larm}=$%d'
   titleHeader += '\n(Red $-$ Approach-1, Blue $-$ Approach-c)'
   plt.title((titleHeader % (1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn,turns)), \
             color='m',fontsize=16)
#
# Trajectory of electron with maximal transferred px:
#
   turns_1=larmorNumber[trackNumbMaxAbsDpxApprch_1]                  # Number of larmor turns for drawing 
   pointsTurns=turns_1*stepsNumberOnGyro                             # Number of points for drawing
   pointsStartLarmor=turns_1*stepsNumberOnGyro                       # Number of points for drawing
   lengthArrowElc=4
   turns_c=larmorNumber[trackNumbMaxAbsDpxApprch_c]                  # Number of larmor turns for drawing 
   pBegArrw=turns_c/2-lengthArrowElc
   pEndArrw=turns_c/2+lengthArrowElc
   print 'First: %d, max: %d' % (trackNumbMaxAbsDpxApprch_1,trackNumbMaxAbsDpxApprch_c)

   fig445=plt.figure(445)
   ax445=fig445.gca(projection='3d')
   ax445.plot(1.e+4*prtclCoorMaxAbsDpx_1[0,0:pointsStartLarmor], \
              1.e+4*prtclCoorMaxAbsDpx_1[2,0:pointsStartLarmor], \
              1.e+4*prtclCoorMaxAbsDpx_1[4,0:pointsStartLarmor],'-r',linewidth=2)
   ax445.plot(1.e+4*prtclCoorMaxAbsDpx_1[0,0:lengthArrowElc], \
              1.e+4*prtclCoorMaxAbsDpx_1[2,0:lengthArrowElc], \
              1.e+4*prtclCoorMaxAbsDpx_1[4,0:lengthArrowElc],'-r',linewidth=4)
   ax445.plot(1.e+4*prtclCoorMaxAbsDpx_1[0,pointsStartLarmor-lengthArrowElc:pointsStartLarmor], \
              1.e+4*prtclCoorMaxAbsDpx_1[2,pointsStartLarmor-lengthArrowElc:pointsStartLarmor], \
              1.e+4*prtclCoorMaxAbsDpx_1[4,pointsStartLarmor-lengthArrowElc:pointsStartLarmor], \
	      '-r',linewidth=4)
   ax445.plot(1.e+4*prtclCoorMaxAbsDpx_c[0,0:turns_c], \
              1.e+4*prtclCoorMaxAbsDpx_c[2,0:turns_c], \
              1.e+4*prtclCoorMaxAbsDpx_c[4,0:turns_c],'-b',linewidth=2)
   ax445.plot(1.e+4*prtclCoorMaxAbsDpx_c[0,pBegArrw:pEndArrw], \
              1.e+4*prtclCoorMaxAbsDpx_c[2,pBegArrw:pEndArrw], \
              1.e+4*prtclCoorMaxAbsDpx_c[4,pBegArrw:pEndArrw],'-b',linewidth=4)
   plt.xlabel('x, $\mu m$',color='m',fontsize=16)
   plt.ylabel('y, $\mu m$',color='m',fontsize=16)
   ax445.set_zlabel('z, $\mu m$',color='m',fontsize=16)
   titleHeader='Trajectory of Electron with Maximal Transferred $p_x$:'
   titleHeader += '\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m, $N_{Larm}=$%d'
   titleHeader += '\n(Red $-$ Approach-1, Blue $-$ Approach-c)'
   plt.title((titleHeader % (1.e+4*rhoMaxAbsDpxTurn_1,1.e+4*rhoLarmorMaxAbsDpxTurn_1,turns_1)), \
             color='m',fontsize=16)

   pointsTurns_1=pointTrack_1[0]
   pointsTurns_c=pointTrack_c[0]

   plt.figure (450)
   plt.plot(1.e+4*prtclCoorFirst_1[4,0:pointsTurns_1],1.e+4*prtclCoorFirst_1[12,0:pointsTurns_1], \
            '-xr',linewidth=2)
   plt.plot(1.e+4*prtclCoorFirst_c[4,0:pointsTurns_c],1.e+4*prtclCoorFirst_c[12,0:pointsTurns_c], \
            '-ob',linewidth=2)
   plt.xlabel('z, $\mu$m',color='m',fontsize=16)
   plt.ylabel('Distance to ion, $\mu$m',color='m',fontsize=16)
   plt.title("Electron's Track", color='m',fontsize=16)
   plt.legend(['Approach-1: With Larmor Rotation','Approach-c: Averaging of Rotation'], \
              loc='lower center',fontsize=16)
   plt.grid(True)

if runFlagApproach_1 == 1 and runFlagApproach_c == 1 and plotFlagDpTransf == 1:
#
# Difference of transfered momentum px (map):
#
   diffDpx_A2minusA1=np.zeros((nA,nB))
   for iA in range(nA):
      for iB in range(nB):
         if dpxTotal_1[iA,iB] !=0.:
            diffDpx_A2minusA1[iA,iB]=100.*(1.-dpxTotal_c[iA,iB]/dpxTotal_1[iA,iB])

   fig545=plt.figure(545)
   ax=fig545.add_subplot(111)                                       # for contours plotting
   mapDpx=ax.contourf(X,Y,diffDpx_A2minusA1)   
   plt.plot(xAboundary[0:81],xBboundary[0:81],'-w')
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   titleHeader='Difference of Transfered Momentum $dP_x$: 1-$Apprch_c$/$Apprch_1$ (%%)'
   titleHeader +='\nTracks: %d'
   plt.title((titleHeader % lastTrackNumber), color='m',fontsize=14)
#   plt.text(-4.9,-2.95,'Unallowable Tracks: Initial\nImpact Parameter > $R_{shield}$',color='w',fontsize=20)
#   plt.text(-4.65,-2.25,'Unallowable\nTracks: Initial\nImpact Parameter\nLarger than $R_{shield}$',color='w',fontsize=20)
   plt.text(-4.9,-2.5, \
            '        Area with\nUnallowable Tracks:\n    Initial Impact\nParameter Larger\n     than $R_{shield}$',\
            color='w',fontsize=20)
   plt.ylim([-3.,-.5])
   fig545.colorbar(mapDpx)
#
# Difference of transfered momentum pz (map):
#
   diffDpz_A2minusA1=np.zeros((nA,nB))
   for iA in range(nA):
      for iB in range(nB):
         if dpzTotal_1[iA,iB] !=0.:
            diffDpz_A2minusA1[iA,iB]=100.*(1.-dpzTotal_c[iA,iB]/dpzTotal_1[iA,iB])

   fig555=plt.figure(555)
   ax=fig555.add_subplot(111)                                       # for contours plotting
   mapDpx=ax.contourf(X,Y,diffDpz_A2minusA1)   
   plt.plot(xAboundary[0:81],xBboundary[0:81],'-w')
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   titleHeader='Difference of Transfered Momentum $dP_z$: 1-$Apprch_c$/$Apprch_1$ (%%)'
   titleHeader +='\nTracks: %d'
   plt.title((titleHeader % lastTrackNumber), color='m',fontsize=14)
#   plt.text(-4.9,-2.95,'Unallowable Tracks: Initial\nImpact Parameter > $R_{shield}$',color='w',fontsize=20)
#   plt.text(-4.65,-2.25,'Unallowable\nTracks: Initial\nImpact Parameter\nLarger than $R_{shield}$',color='w',fontsize=20)
   plt.text(-4.9,-2.5, \
            '        Area with\nUnallowable Tracks:\n    Initial Impact\nParameter Larger\n     than $R_{shield}$',\
            color='w',fontsize=20)
   plt.ylim([-3.,-.5])
   fig555.colorbar(mapDpx)
'''

if runFlagApproach_m == 1 and plotFlagTracks == 1:

#
# First electron's trajectory:
#
   turns=larmorNumber[0]                                           # Number of larmor turns for drawing 
   pointsTurns=turns                                               # Number of points for drawing
   lengthArrowElc=5
   pBegArrw=turns/2-lengthArrowElc/2
   pEndArrw=turns/2+lengthArrowElc/2

   fig640=plt.figure(640)
   ax640=fig640.gca(projection='3d')
   ax640.plot(1.e+4*prtclCoorFirst_m[0,0:pointsTurns], \
              1.e+4*prtclCoorFirst_m[2,0:pointsTurns], \
              1.e+4*prtclCoorFirst_m[4,0:pointsTurns],'-r',linewidth=2)
   ax640.plot(1.e+4*prtclCoorFirst_m[0,pBegArrw:pEndArrw], \
              1.e+4*prtclCoorFirst_m[2,pBegArrw:pEndArrw], \
              1.e+4*prtclCoorFirst_m[4,pBegArrw:pEndArrw],'-b',linewidth=4)
   plt.xlabel('x, $\mu m$',color='m',fontsize=16)
   plt.ylabel('y, $\mu m$',color='m',fontsize=16)
   ax640.set_zlabel('z, $\mu m$',color='m',fontsize=16)
   titleHeader="Approach-m. First Electron's Trajectory:"
   titleHeader += '\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m, $N_{Larm}=$%d'
   plt.title((titleHeader % (1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn,turns)), \
             color='m',fontsize=16)
#
# Trajectory of electron with maximal transferred px:
#
   turns=larmorNumber[trackNumbMaxAbsDpxApprch_m]                  # Number of larmor turns for drawing 
   pointsTurns=turns                                               # Number of points for drawing
   lengthArrowElc=8
   pBegArrw=turns/2-lengthArrowElc/2
   pEndArrw=turns/2+lengthArrowElc/2

   fig645=plt.figure(645)
   ax645=fig645.gca(projection='3d')
   ax645.plot(1.e+4*prtclCoorMaxAbsDpx_m[0,0:pointsTurns], \
              1.e+4*prtclCoorMaxAbsDpx_m[2,0:pointsTurns], \
              1.e+4*prtclCoorMaxAbsDpx_m[4,0:pointsTurns],'-r',linewidth=2)
   ax645.plot(1.e+4*prtclCoorMaxAbsDpx_m[0,pBegArrw:pEndArrw], \
              1.e+4*prtclCoorMaxAbsDpx_m[2,pBegArrw:pEndArrw], \
              1.e+4*prtclCoorMaxAbsDpx_m[4,pBegArrw:pEndArrw],'-b',linewidth=4)
   plt.xlabel('x, $\mu m$',color='m',fontsize=16)
   plt.ylabel('y, $\mu m$',color='m',fontsize=16)
   ax645.set_zlabel('z, $\mu m$',color='m',fontsize=16)
   titleHeader='Approach-m. Trajectory of Electron with Maximal Transferred $p_x$:'
   titleHeader += '\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m, $N_{Larm}=$%d'
   plt.title((titleHeader % (1.e+4*rhoMaxAbsDpxTurn_m,1.e+4*rhoLarmorMaxAbsDpxTurn_m,turns)), \
             color='m',fontsize=16)

#
# First electron's trajectory (action):
#
   specFctrJ=1.e31                                                  #       Use it in titles!
   specShiftJ=5.86748911-.000017                                    # x e8! Use it in titles!

   plt.figure (460)
   plt.plot(1.e+4*prtclCoorFirst_m[4,0:pointsTurns], \
            specFctrJ*prtclCoorFirst_m[13,0:pointsTurns]-1.e8*specShiftJ,'-xr',linewidth=2)
   plt.xlabel('z, $\mu$m',color='m',fontsize=16)
   plt.ylabel(('Action $J$: $10^{31}\cdot J$ $-$ %10.8f$\cdot10^8$, g$\cdot$cm$^2$/sec' % specShiftJ), \
              color='m',fontsize=16)
   plt.title("Approach-m: First Electron's Track", color='m',fontsize=16)
   plt.grid(True)

if runFlagApproach_m == 1 and plotFlagDpTransf == 1:              
 
   specFctr3=1.e26                                                   # Use it in titles!
   powerSpecFctr3=26                                                 # Use it in titles!
   X,Y=np.meshgrid(crrntA,crrntB) 
#
# Transfered momentum px (surface):
#
   fig540=plt.figure(540)
   ax540=fig540.gca(projection='3d')
#    surf=ax540.plot_surface(X,Y,specFctr3*dpxTotal_m,cmap=cm.coolwarm, \
#                            linewidth=0,antialiased=False)
   surf=ax540.plot_surface(X,Y,specFctr3*dpxTotal_m,cmap=cm.jet,linewidth=0,antialiased=False)
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   ax540.set_zlabel('$dP_x \cdot 10^{22}$; $g \cdot cm/sec$',color='m',fontsize=16)
#   titleHeader='Approach-m: Transfered Momentum $dP_x$ $(\cdot 10^{22}$; $g \cdot cm/sec)$'
#   titleHeader += '\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-22}$ $g \cdot cm/sec$)'
#   plt.title((titleHeader % (lastTrackNumber,specFctr3*dpxTotalMax_m)), color='m',fontsize=16)
   titleHeader='Approach-m: Transfered Momentum $dP_x$ $(\cdot 10^{%2d}$; $g \cdot cm/sec)$'
   titleHeader +='\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-%2d}$ $g \cdot cm/sec)$)'
   plt.title((titleHeader % \
              (powerSpecFctr3,lastTrackNumber,specFctr3*dpxTotalMax_m,powerSpecFctr3)), \
	     color='m',fontsize=14)
   cb = fig540.colorbar(surf)
   plt.grid(True)
#
# Transfered momentum px (map):
#
   fig845=plt.figure(845)
   ax=fig845.add_subplot(111)                                       # for contours plotting
   mapDpx=ax.contourf(X,Y,specFctr3*dpxTotal_m)   
   plt.plot(xAboundary[0:81],xBboundary[0:81],'-w')
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   titleHeader='Approach-m: Transfered Momentum $dP_x$ $(\cdot 10^{%2d}$; $g \cdot cm/sec)$'
   titleHeader +='\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-%2d}$ $g \cdot cm/sec)$)'
   plt.title((titleHeader % \
              (powerSpecFctr3,lastTrackNumber,specFctr3*dpxTotalMax_m,powerSpecFctr3)), \
	     color='m',fontsize=14)
   plt.ylim([-3.,-.5])
   plt.text(-4.9,-2.5, \
            '        Area with\nUnallowable Tracks:\n    Initial Impact\nParameter Larger\n     than $R_{shield}$',\
            color='w',fontsize=20)
   fig845.colorbar(mapDpx)

#
# Transfered momentum py (surface):
#
###   fig850=plt.figure(850)
###   ax850=fig850.gca(projection='3d')
####    surf=ax850.plot_surface(X,Y,specFctr3*dpyTotal_m,cmap=cm.coolwarm, \
####                            linewidth=0,antialiased=False)
###   surf=ax550.plot_surface(X,Y,specFctr3*dpyTotal_m,cmap=cm.jet,linewidth=0,antialiased=False)
###   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
###   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
###   ax850.set_zlabel('$dP_y \cdot 10^{22}$; $g \cdot cm/sec$',color='m',fontsize=16)
###   titleHeader='Approach-c: Transfered Momentum $dP_y$ $(\cdot 10^{22}$; $g \cdot cm/sec)$'
###   titleHeader += '\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-22}$ $g \cdot cm/sec$)'
###   plt.title((titleHeader % (lastTrackNumber,specFctr3*dpyTotalMax_m)), color='m',fontsize=16)
###   cb = fig850.colorbar(surf)
###   plt.grid(True)
####
#### Transfered momentum py (map):
####
###   fig855=plt.figure(855)
###   ax=fig855.add_subplot(111)                                       # for contours plotting
###   mapDpy=ax.contourf(X,Y,specFctr3*dpyTotal_m)   
###   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
###   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
###   titleHeader='Approach-m: Transfered Momentum $dP_y$ $(\cdot 10^{22}$; $g \cdot cm/sec)$'
###   titleHeader +='\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-22}$ $g \cdot cm/sec)$)'
###   plt.title((titleHeader % (lastTrackNumber,specFctr3*dpyTotalMax_m)), color='m',fontsize=14)
###   fig855.colorbar(mapDpx)

#
# Transfered momentum pz (surface):
#
   fig860=plt.figure(860)
   ax860=fig860.gca(projection='3d')
#    surf=ax850.plot_surface(X,Y,specFctr3*dpzTotal_m,cmap=cm.coolwarm, \
#                            linewidth=0,antialiased=False)
   surf=ax860.plot_surface(X,Y,specFctr3*dpzTotal_m,cmap=cm.jet,linewidth=0,antialiased=False)
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   ax860.set_zlabel('$dP_z \cdot 10^{22}$; $g \cdot cm/sec$',color='m',fontsize=16)
#   titleHeader='Approach-m: Transfered Momentum $dP_z$ $(\cdot 10^{22}$; $g \cdot cm/sec)$'
#   titleHeader += '\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-22}$ $g \cdot cm/sec$)'
#   plt.title((titleHeader % (lastTrackNumber,specFctr3*dpzTotalMax_m)), color='m',fontsize=16)
   titleHeader='Approach-m: Transfered Momentum $dP_z$ $(\cdot 10^{%2d}$; $g \cdot cm/sec)$'
   titleHeader +='\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-%2d}$ $g \cdot cm/sec)$)'
   plt.title((titleHeader % \
              (powerSpecFctr3,lastTrackNumber,specFctr3*dpzTotalMax_m,powerSpecFctr3)), \
	     color='m',fontsize=14)
   cb = fig860.colorbar(surf)
   plt.grid(True)
#
# Transfered momentum pz (map):
#
   fig865=plt.figure(865)
   ax=fig865.add_subplot(111)                                       # for contours plotting
   mapDpz=ax.contourf(X,Y,specFctr3*dpzTotal_m)   
   plt.plot(xAboundary[0:81],xBboundary[0:81],'-w')
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
#   titleHeader='Approach-m: Transfered Momentum $dP_z$ $(\cdot 10^{22}$; $g \cdot cm/sec)$'
#   titleHeader +='\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-22}$ $g \cdot cm/sec)$)'
#   plt.title((titleHeader % (lastTrackNumber,specFctr3*dpzTotalMax_m)), color='m',fontsize=14)
   titleHeader='Approach-m: Transfered Momentum $dP_z$ $(\cdot 10^{%2d}$; $g \cdot cm/sec)$'
   titleHeader +='\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-%2d}$ $g \cdot cm/sec)$)'
   plt.title((titleHeader % \
              (powerSpecFctr3,lastTrackNumber,specFctr3*dpzTotalMax_m,powerSpecFctr3)), \
	     color='m',fontsize=14)
   plt.ylim([-3.,-.5])
   plt.text(-4.9,-2.5, \
            '        Area with\nUnallowable Tracks:\n    Initial Impact\nParameter Larger\n     than $R_{shield}$',\
            color='w',fontsize=20)
   fig865.colorbar(mapDpx)

'''
if runFlagApproach_1 == 1 and runFlagApproach_m == 1 and plotFlagTracks == 1:
#
# First electron's trajectory (distances to ion - in lab.system and from "guidingCenterCollision"):
#
   pointsTurns_1=pointTrack_1[0]
   pointsTurns_m=pointTrack_m[0]

   plt.figure (650)
   plt.plot(1.e+4*prtclCoorFirst_1[4,0:pointsTurns_1],1.e+4*prtclCoorFirst_1[12,0:pointsTurns_1],'-xr')
   plt.plot(1.e+4*prtclCoorFirst_m[4,0:pointsTurns_m],1.e+4*prtclCoorFirst_m[12,0:pointsTurns_m],'-ob')
   plt.xlabel('z, $\mu$m',color='m',fontsize=16)
   plt.ylabel('Distance to ion, $\mu$m',color='m',fontsize=16)
   plt.title("First Electron's Track", color='m',fontsize=16)
   plt.legend(['Approach_1: With Larmor Rotation','Approach_m: "GuidingCenter" System'], \
              loc='best',fontsize=16)
   plt.grid(True)

if runFlagApproach_c == 1 and runFlagApproach_m == 1 and plotFlagTracks == 1:
#
# First electron's trajectory (distances to ion - in "guidingCenterCollision" system):
#
   pointsTurns_c=pointTrack_c[0]
   pointsTurns_m=pointTrack_m[0]

   plt.figure (660)
   plt.plot(1.e+4*prtclCoorFirst_c[4,0:pointsTurns_c], \
            1.e+4*prtclCoorFirst_c[12,0:pointsTurns_c],'-xr')
   plt.plot(1.e+4*prtclCoorFirst_m[4,0:pointsTurns_m], \
            1.e+4*prtclCoorFirst_m[12,0:pointsTurns_m],'ob')
   plt.xlabel('z, $\mu$m',color='m',fontsize=16)
   plt.ylabel('Distance to ion, $\mu$m',color='m',fontsize=16)
   plt.title('First Electron''s Track ("GuidingCenter" System)', color='m',fontsize=16)
   plt.legend(['Approach_c','Approach_m'],loc='best',fontsize=16)
   plt.grid(True)

if runFlagApproach_1 == 1 and runFlagApproach_m == 1 and plotFlagDpTransf == 1:
#
# Difference of transfered momentum px (map):
#
   diffDpx_A3minusA1=np.zeros((nA,nB))
   for iA in range(nA):
      for iB in range(nB):
         if dpxTotal_1[iA,iB] !=0.:
            diffDpx_A3minusA1[iA,iB]=100.*(1.-dpxTotal_m[iA,iB]/dpxTotal_1[iA,iB])

   fig745=plt.figure(745)
   ax=fig745.add_subplot(111)                                       # for contours plotting
   mapDpx=ax.contourf(X,Y,diffDpx_A3minusA1)   
   plt.plot(xAboundary[0:81],xBboundary[0:81],'-w')
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   titleHeader='Difference of Transfered Momentum $dP_x$: 1-$Apprch_m$/$Apprch_1$ (%%)'
   titleHeader +='\nTracks: %d'
   plt.title((titleHeader % lastTrackNumber), color='m',fontsize=14)
#   plt.text(-4.9,-2.95,'Unallowable Tracks: Initial\nImpact Parameter > $R_{shield}$',color='w',fontsize=20)
#   plt.text(-4.65,-2.25,'Unallowable\nTracks: Initial\nImpact Parameter\nLarger than $R_{shield}$',color='w',fontsize=20)
   plt.text(-4.9,-2.5, \
            '        Area with\nUnallowable Tracks:\n    Initial Impact\nParameter Larger\n     than $R_{shield}$',\
            color='w',fontsize=20)
   plt.ylim([-3.,-.5])
   fig745.colorbar(mapDpx)
#
# Difference of transfered momentum pz (map):
#
   diffDpz_A3minusA1=np.zeros((nA,nB))
   for iA in range(nA):
      for iB in range(nB):
         if dpzTotal_1[iA,iB] !=0.:
            diffDpz_A3minusA1[iA,iB]=100.*(1.-dpzTotal_m[iA,iB]/dpzTotal_1[iA,iB])

   fig755=plt.figure(755)
   ax=fig755.add_subplot(111)                                       # for contours plotting
   mapDpx=ax.contourf(X,Y,diffDpz_A3minusA1)   
   plt.plot(xAboundary[0:81],xBboundary[0:81],'-w')
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   titleHeader='Difference of Transfered Momentum $dP_z$: 1-$Apprch_m$/$Apprch_1$ (%%)'
   titleHeader +='\nTracks: %d'
   plt.title((titleHeader % lastTrackNumber), color='m',fontsize=14)
#   plt.text(-4.9,-2.95,'Unallowable Tracks: Initial\nImpact Parameter > $R_{shield}$',color='w',fontsize=20)
#   plt.text(-4.65,-2.25,'Unallowable\nTracks: Initial\nImpact Parameter\nLarger than $R_{shield}$',color='w',fontsize=20)
   plt.text(-4.9,-2.5, \
            '        Area with\nUnallowable Tracks:\n    Initial Impact\nParameter Larger\n     than $R_{shield}$',\
            color='w',fontsize=20)
   plt.ylim([-3.,-.5])
   fig755.colorbar(mapDpx)
'''

plt.show()

sys.exit()

