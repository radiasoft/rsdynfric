# from __future__ import division

#-------------------------------------
#
#        Started at 06/08/2018 (YuE)
#
# This script based on the previous script
# threeApproachesComparison_v6.py
#
# 07/13/2018: IT IS NOT FINISHED ("classical" approach works appropriate;
#                                 goog agreement with Magnus Expansion approach)!
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
alpha=0.                                # angle between V_i and magnetic field
trnsvT=0.5                              # transversal temperature, eV
longT=2.0e-4                         # longitudinal temperature, eV (was 2.0e-4)
eTempTran=trnsvT                        # to keep variable from previous script
eTempLong=longT                         # to keep variable from previous script
nField=1                             # number ov values  of the magnetic field
fieldB=np.zeros(nField)                 
fieldB[0]=3.e3                          # magnetic field, Gs
B_mag=fieldB[0]                         # to keep variable from previous script
omega_p=1.0e9                           # plasma frequency, 1/sec
n_e=omega_p**2*m_e/(4.*pi*q_e**2)       # plasma density, 3.1421e+08 cm-3

n_e1=8.e7                               # plasma density, cm-3
omega_p1=np.sqrt(4.*pi*n_e1*q_e**2/m_e) # plasma frequency, 5.0459e+08 1/s  
#
# Cooling system parameter:
#
coolLength=150.0        # typical length of the coolong section, cm

# 
# Calculated parameters of the electron beam:
#
V0=np.sqrt(2.*Ekin*eVtoErg/m_e)           # longitudinal velocity, cm/s
rmsTrnsvVe=np.sqrt(2.*trnsvT*eVtoErg/m_e) # RMS transversal velocity, cm/s
rmsLongVe=np.sqrt(2.*longT*eVtoErg/m_e)   # RMS longitudinal velocity, cm/s
# dens=curBeam*CtoPart/V0                 # density, 1/cm^3
# omega=np.sqrt(4.*pi*dens*q_e**2/m_e)    # plasma frequency, 1/s
cyclFreq=q_e*fieldB/(m_e*cLight)          # cyclotron frequency, 1/s
rmsRoLarm=rmsTrnsvVe*cyclFreq**(-1)       # RMS Larmor radius, cm
dens=omega_p**2*m_e/(4.*pi*q_e**2)        # density, 1/cm^3
likeDebyeR=(3./dens)**(1./3.)             # "Debye" sphere with 3 electrons, cm

coolPassTime=coolLength/V0             # time pass through cooling section, cm

powV0=round(np.log10(V0)) 
mantV0=V0/(10**powV0) 

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

print ('V0=%e cm/s, rmsTrnsvVe=%e cm/s, rmsLongVe=%e cm/s' % (V0,rmsTrnsvVe,rmsLongVe))

ergToEV = 1./1.60218e-12

rhoCrit=(m_e*(cLight/fieldB[0])**2)**(1./3)
print ('rhoCrit=%e cm' % rhoCrit)

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
print 'omega_Larmor= %e rad/sec, T_larm = %e sec, timeStep = %e sec' % \
      (omega_L,T_larm,timeStep)

nLarmorAvrgng=10               # number of averaged Larmor rotations 
#
# Data to integrate transferred momemta over the track:
#
timeStep_c=nLarmorAvrgng*stepsNumberOnGyro*timeStep                # sec      
print 'timeStep_c = %e s' % timeStep_c

eVrmsTran = np.sqrt(2.*eTempTran*eVtoErg/m_elec)                   # cm/sec
eVrmsLong = np.sqrt(2.*eTempLong*eVtoErg/m_elec)                   # cm/sec
kinEnergy = m_elec*(eVrmsTran**2+eVrmsLong**2)/2.     # kinetic energy; erg
print 'eVrmsTran = %e cm/sec, eVrmsLong = %e cm/sec, kinEnergy = %e eV' % \
      (eVrmsTran,eVrmsLong,ergToEV*kinEnergy)

ro_larmRMS = eVrmsTran/omega_L                                     # cm
print 'ro_larmRMS =%e mkm = ', 1.e4*ro_larmRMS
#
# Electrons are magnetized for impact parameter >> rhoCrit:
#
rhoCrit=math.pow(q_elec**2/(m_elec*omega_L**2),1./3)               # cm
print 'rhoCrit (mkm) = ' , 1.e+4*rhoCrit

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
# Factor to calculate transferred momenta for ion;
# it will be calculed also inside the routines "guidingCenterCollision"
# and "MagnusExpansionCollision":
#
dpFactor = q_elec**2*timeStep_c                              # g*cm^3/sec
print 'dpFactor = %e g*cm^3/s' % dpFactor
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
# Dimensions of dpIon, deElec are g*cm/sec:   
   dpIon[0]=-dpFactor_gc*deltaT*(vectrIon[0]-x_gc)/b_gc**3               
   dpIon[1]=-dpFactor_gc*deltaT*(vectrIon[2]-vectrElec_gc[2])/b_gc**3  
   dpIon[2]=-dpFactor_gc*deltaT*(vectrIon[4]-vectrElec_gc[4])/b_gc**3
   dpElec[1]=-dpIon[1] 
   dpElec[2]=-dpIon[2]
#    print 'dpIon[0]=%e, dpIon[1]=%e, dpIon[2]=%e' % \
#           (dpIon[0],dpIon[1],dpIon[2])
   return dpIon,dpElec,action,b_gc                                      

#
# "Magnus expansion" description of the collision during time interval 'deltaT'
# in the system coordinates of "guiding center" of electron
# input - 6-vectors for electron and ion before collision and time step deltaT; 
# output - transfered momenta to ion and electron and electron y_gc coordinate
# as well calculated parameters C1,C2,C3,b,D1,D2,q for testing: 
#
def MagnusExpansionCollision(vectrElec_gc,vectrIon,deltaT):

#    print 'Ion: x=%e, y=%e, z=%e' % (vectrIon[0],vectrIon[2],vectrIon[4])
#    print 'Electron: x=%e, y=%e, z=%e' % 
#          (vectrElec_gc[0],vectrElec_gc[4],vectrElec_gc[4])
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
   dpIon[1]=-2.*dpFactor_gc/q* \
           ((vectrIon[2]-vectrElec_gc[2])*D1-vectrIon[3]/M_ion*D2)
   dpIon[2]=-2.*dpFactor_gc/q*((vectrIon[4]-vectrElec_gc[4])*D1- \
                               (vectrIon[5]/M_ion-vectrElec_gc[5]/m_elec)*D2) 
   dpElec[1]=-dpIon[1] 
   dpElec[2]=-dpIon[2] 
   dy_gc=dpIon[0]/mOmegaLarm                                        # cm
#    print 'dpIon[0]=%e, dpIon[1]=%e, dpIon[2]=%e' % \
#           (dpIon[0],dpIon[1],dpIon[2])
   return dpIon,dpElec,action,dy_gc,C1,C2,C3,b,D1,D2,q                                      

sphereNe=3.
R_e=math.pow(sphereNe/n_e,1./3)                       # cm
print 'R_e (cm)=%e' % R_e

ro_Larm = eVrmsTran/omega_L                           # cm
print 'ro_Larm (cm)=%e' % ro_Larm

impctPrmtrMin=2.*ro_Larm

# rhoDependenceFlag = 1  # skip calculation of rho dependence if = 0!

# Taking into account the transfer of momenta for both particles
# (for "classical" only):
dpTransferFlag = 1        # no taking into account if = 0!

nVion=50
Vion=np.zeros(nVion)
VionRel=np.zeros(nVion)

vIonMin=4.e-3*eVrmsTran
vIonMax=10.*eVrmsTran
print 'VionMin=%e, vIonMax=%e' % (vIonMin,vIonMax)   

vIonMinRel=vIonMin/V0
vIonMaxRel=vIonMax/V0

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

'''
fig10=plt.figure(10)
plt.semilogx(impctPrmtrMax,impctPrmtrMaxCrrctdRel[:,0],'-r', \
           impctPrmtrMax,impctPrmtrMaxCrrctdRel[:,1],'-b', \
           impctPrmtrMax,impctPrmtrMaxCrrctdRel[:,2],'-g', \
           impctPrmtrMax,impctPrmtrMaxCrrctdRel[:,3],'-m',linewidth=2)
plt.grid(True)
hold=True
plt.xlabel('Maximal Impact parameter $R_{max}$, cm',color='m',fontsize=16)
plt.ylabel('$R_{max}^{Crrctd}/R_{Max}$',color='m',fontsize=16)
# plt.xlim([.9*min(impctPrmtrMax),1.1*max(impctPrmtrMax)])
plt.xlim([1.e-2,1.1*max(impctPrmtrMax)])
plt.ylim([.986,1.001])
titleHeader='$R_{max}^{Crrctd}=R_{Max} \cdot [1-(\pi\cdot N_{Larm} \cdot' 
titleHeader += '\Delta_{e||}/(\omega_{Larm} \cdot R_{max})]^{1/2}$'
plt.title(titleHeader,color='m',fontsize=16)
plt.legend([('$N_{Larm}=$%2d' % larmorTurnsMin[0]), \
           ('$N_{Larm}=$%2d' % larmorTurnsMin[1]), \
           ('$N_{Larm}=$%2d' % larmorTurnsMin[2]), \
           ('$N_{Larm}=$%2d' % larmorTurnsMin[3])],loc='lower center',fontsize=16)
#-----------------------------------------------------------------
'''


xLimit=[.9*VionRel[0],1.1*VionRel[nVion-1]]


fig209=plt.figure(209)
plt.loglog(VionRel,R_debye,'-r',VionRel,R_pass,'-b', \
                                VionRel,R_pass_1,'--b',linewidth=2)
plt.grid(True)
hold=True
plt.plot([VionRel[0],VionRel[nVion-1]],[R_e,R_e],color='m',linewidth=2)   
plt.xlabel('Relative Ion Velocity, $V_i/V_{e0}$',color='m',fontsize=16)
plt.ylabel('$R_{Debye}$, $R_{Pass}$, $R_e$, cm',color='m',fontsize=16)
titleHeader='Magnetized Collision: $R_{Debye}$, $R_{Pass}$, $R_e$: $V_{e0}=%5.3f\cdot10^{%2d}$cm/s'
plt.title(titleHeader % (mantV0,powV0),color='m',fontsize=16)
plt.xlim(xLimit)
yLimit=[1.e-3,10.]
plt.ylim(yLimit)
plt.plot([relVeTrnsv,relVeTrnsv],yLimit,'--m',linewidth=1)
plt.text(2.e-3,6.5e-4,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
plt.plot([relVeLong,relVeLong],yLimit,'--m',linewidth=1)
plt.text(4.4e-5,0.001175,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
plt.text(3.e-5,2.45e-3,'$R_e$',color='k',fontsize=16)
plt.text(3.e-5,5.e-2,'$R_{Debye}$',color='k',fontsize=16)
plt.text(3.e-5,1.8e-2,'$R_{Pass}$',color='k',fontsize=16)
plt.text(4.5e-5,4.8e-3,'$R_{Pass}$ $for$ $T_{e||}=0$',color='k',fontsize=16)

fig3151=plt.figure(3151)
plt.loglog(VionRel,impctPrmtrMax,'-r', VionRel,impctPrmtrMax_1,'--r', \
   [VionRel[0],VionRel[nVion-1]],[impctPrmtrMin,impctPrmtrMin],'-r',linewidth=2)
plt.grid(True)
hold=True
plt.xlabel('Relative Ion Velocity, $V_i/V_{e0}$',color='m',fontsize=16)
plt.ylabel('Impact Parameter, cm',color='m',fontsize=16)
titleHeader= \
          'Types of Collisions: $V_{e0}=%5.3f\cdot10^{%2d}$cm/s, $B=%6.1f$ Gs'
plt.title(titleHeader % (mantV0,powV0,fieldB[0]),color='m',fontsize=16)
plt.xlim(xLimit)
yLimit=[8.e-4,.6]
plt.ylim(yLimit)
plt.plot([relVeTrnsv,relVeTrnsv],yLimit,'--m',linewidth=1)
plt.text(2.e-3,5.9e-4,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
plt.plot([relVeLong,relVeLong],yLimit,'--m',linewidth=1)
plt.text(4.4e-5,.0018,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
plt.text(3.e-4,1.75e-3,'$R_{min}=2\cdot<rho_\perp>$',color='k',fontsize=16)
plt.text(7.e-4,5.e-2,'$R_{max}$',color='k',fontsize=16)
plt.text(2.85e-5,3.3e-3,'$R_{max}$ $for$ $T_{e||}=0$',color='k',fontsize=16)
plt.plot([VionRel[0],VionRel[nVion-1]],[20.*rhoCrit,20.*rhoCrit],color='k')   
plt.text(1.5e-4,7.e-3,'Magnetized Collisions',color='r',fontsize=25)
plt.text(1.e-4,10.e-4,'Adiabatic or Fast Collisions',color='r',fontsize=25)
plt.text(2.8e-5,.2,'Collisions are Screened',color='r',fontsize=25)
plt.text(1.6e-5,1.e-3,'$ \cong 20\cdot R_{Crit}$',color='k',fontsize=16)


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
#   print 'rhoMax(%d) = %e' % (i,rhoMax[i])
   log10rhoMax = math.log10(rhoMax[i])
   log10rhoStep = (log10rhoMax-log10rhoMin)/(nImpctPrmtr)
#   print 'Vion(%d) = %e, rhoMax = %e' % (i,Vion[i],rhoMax[i])
   for n in range(nImpctPrmtr):
      log10rhoCrrnt = log10rhoMin+(n+0.5)*log10rhoStep 
      rhoCrrnt = math.pow(10.,log10rhoCrrnt)
#      print '    rhoCrrnt(%d) = %e' % (n,rhoCrrnt)
      halfLintr[n,i] = np.sqrt(rhoMax[i]**2-rhoCrrnt**2)   # half length of interaction; cm
      timeHalfPath = halfLintr[n,i]/eVrmsLong     # 0.5 time of interaction; sec
      numbLarmor = int(2.*timeHalfPath/T_larm)             
      pointAlongTrack[n,i] = int(2.*timeHalfPath/timeStep_c)
      totalPoints += pointAlongTrack[n,i]
#      print '     %d: rhoCrrnt = %e, numbLarmor = %d, pointAlongTrack = %d' % \
#            (n,rhoCrrnt,numbLarmor,pointAlongTrack[n,i])
print 'totalPoints = %d' % totalPoints

arrayA=np.zeros(2*totalPoints)      
arrayB=np.zeros(2*totalPoints)     
bCrrnt_c = np.zeros(2*totalPoints)
#
# Variables for resting:
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
deltaEnrgIon_c = np.zeros((nImpctPrmtr,nVion))
#   
# "Magnus Expand" approach:
#
deltaPx_m = np.zeros((nImpctPrmtr,nVion))  
deltaPy_m = np.zeros((nImpctPrmtr,nVion))  
deltaPz_m = np.zeros((nImpctPrmtr,nVion))  
deltaEnrgIon_m = np.zeros((nImpctPrmtr,nVion))   
#   
# Comparison of approaches (ratio deltaEnrgIon_c/deltaEnrgIon_m):
#
dEion_c_m = np.zeros((nImpctPrmtr,nVion)) 
#
# Factor to calculate transferred energy to ion
# (the friction force is defined by this transfered energy): 
#
deFactor = 0.5/M_ion                               # 1/g

frctnForce_cSM = np.zeros(nVion)      # integration, using Simpson method
frctnForce_mSM = np.zeros(nVion)      # integration, using Simpson method

indx = 0
for i in range(nVion):
# Taking into account the corection of the maximal impact parameter
# on depence of preset number of minimal Larmor turns:
   rhoMax[i] = impctPrmtrMax[i]* \
   np.sqrt(1.- (pi*larmorTurns*eVrmsLong/omega_L/impctPrmtrMax[i])**2)
# Without taking into account the corection of the maximal impact parameter
# on depence of preset number of minimal Larmor turns:
   rhoMax[i] = impctPrmtrMax[i]
   log10rhoMax = math.log10(rhoMax[i])
   log10rhoStep = (log10rhoMax-log10rhoMin)/(nImpctPrmtr)
   print 'Vion(%d) = %e, rhoMax = %e' % (i,Vion[i],rhoMax[i])
   for n in range(nImpctPrmtr):
      log10rhoCrrnt = log10rhoMin+(n+0.5)*log10rhoStep 
      rhoCrrnt = math.pow(10.,log10rhoCrrnt)
#      rhoInit[i*nImpctPrmtr+n] = rhoCrrnt
      rhoInit[n,i] = rhoCrrnt
      halfLintr[n,i] = np.sqrt(rhoMax[i]**2-rhoCrrnt**2)   # half length of interaction; cm
      z_ionCrrnt = np.zeros(6)      # Zeroing out of vector for ion
      z_elecCrrnt = np.zeros(6)     # Zeroing out of vector for electron
# Zeroing out of "guiding center" vector for electron:
      z_elecCrrnt_gc = np.zeros(6)  
# Current values to transfered momemta 
# (second index numerates "classical", (if 0) and 
#                         "Magnus expantion" (if 1) approaches: 
      dpCrrnt = np.zeros((3,2))
# Intermediate arrays:
      dpIon_c = np.zeros(3) 
      dpIon_m = np.zeros(3) 
      dpElec_c = np.zeros(3) 
      dpElec_m = np.zeros(3) 
# Current initial vector for electron:
      z_elecCrrnt[Ix] = rhoCrrnt                     # x, cm
      z_elecCrrnt[Iz] = -halfLintr[n,i]              # z, cm
      z_elecCrrnt[Ipy] = m_elec*eVrmsTran            # py, g*cm/sec
      z_elecCrrnt[Ipz] = m_elec*eVrmsLong            # pz, g*cm/sec
# transfer to system of guiding center:
      z_elecCrrnt_gc=toGuidingCenter(z_elecCrrnt)  
#
# Main loop along the each track:
#
      for k in range(int(pointAlongTrack[n,i])):
#
# Dragging both particles through first half of the step of the track:
#
  	 z_elecCrrnt_gc = np.dot(matr_elec_c,z_elecCrrnt_gc) # electron
  	 z_ionCrrnt = np.dot(matr_ion_c,z_ionCrrnt)          # ion
# transfer from system of guiding center: 
	 z_elecCrrnt=fromGuidingCenter(z_elecCrrnt_gc)     
# Current distance between ion and electron; cm:
 	 bCrrnt_c[indx]=np.sqrt((z_ionCrrnt[0]-z_elecCrrnt[0])**2+ \
	                  (z_ionCrrnt[2]-z_elecCrrnt[2])**2+ \
			  (z_ionCrrnt[4]-z_elecCrrnt[4])**2)
# Current values of parameters A,B:  
	 arrayA[indx] = math.log10(ro_Larm/bCrrnt_c[indx])     
	 arrayB[indx] = math.log10((q_elec**2/bCrrnt_c[indx])/kinEnergy)
         indx += 1
#
# Dragging both particles through interaction during this step of track
# (for both approaches):
#
#    "Classical":
	 dpIon_c,dpElec_c,action,b_gc_c = \
	         guidingCenterCollision(z_elecCrrnt_gc,z_ionCrrnt,timeStep_c) 
#    "Magnus Expantion":
	 dpIon_m,dpElec_m,actionME,dy_gc_m,C1,C2,C3,b,D1,D2,q = \
	         MagnusExpansionCollision(z_elecCrrnt_gc,z_ionCrrnt,timeStep_c) 
# Save data for testing:
         b_gc[indxTest] = b_gc_c           # "classical" approach
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
	       z_ionCrrnt[2*ic+1] += dpIon_c[ic]   
	       z_elecCrrnt[2*ic+1] += dpElec_c[ic]
# transfer to system of guiding center:
         z_elecCrrnt_gc=toGuidingCenter(z_elecCrrnt)  
# Accumulation of the transfered momenta to ion along the track for both approaches:  
         for ic in range(3):
#	    if i == 0:
#	       print 'dpIon_c[%2d] = %20.14e, dpIon_m[%2d] = %20.14e' % \
#	             (ic,dpIon_c[ic],ic,dpIon_m[ic])
 	    dpCrrnt[ic,0] += dpIon_c[ic]       # "classical", g*cm/sec  
 	    dpCrrnt[ic,1] += dpIon_m[ic]       # "Magnus expansion", g*cm/sec  
#
# Dragging both particles through second half of the step of the track:
#
 	 z_elecCrrnt_gc = np.dot(matr_elec_c,z_elecCrrnt_gc)     # electron
 	 z_ionCrrnt = np.dot(matr_ion_c,z_ionCrrnt)              # ion
# transfer from system of guiding center: 
	 z_elecCrrnt=fromGuidingCenter(z_elecCrrnt_gc)     
# Current distance between ion and electron; cm:
 	 bCrrnt_c[indx]=np.sqrt((z_ionCrrnt[0]-z_elecCrrnt[0])**2+ \
	                  (z_ionCrrnt[2]-z_elecCrrnt[2])**2+ \
			  (z_ionCrrnt[4]-z_elecCrrnt[4])**2)
# Current values of parameters A,B:  
	 arrayA[indx] = math.log10(ro_Larm/bCrrnt_c[indx])     
	 arrayB[indx] = math.log10((q_elec**2/bCrrnt_c[indx])/kinEnergy)
         indx += 1
#
# Transferred momenta and energy along the track - "classical" approach:
#
      deltaPx_c[n,i] = abs(dpCrrnt[0,0])                         # dpx, g*cm/sec
#      if deltaPx_c[n,i] <= 0.: 
#         print 'deltaPx_c[%2d,%2d] = %e, dpCrrnt[%2d,%2d] = %e' % \
#	       (n,i,deltaPx_c[n,i],n,i,dpCrrnt[0,0])
      deltaPy_c[n,i] = abs(dpCrrnt[1,0])                         # dpy, g*cm/sec 
#      if deltaPy_c[n,i] <= 0.: 
#         print 'deltaPy_c[%2d,%2d] = %e' % (n,i,deltaPy_c[n,i])
      deltaPz_c[n,i] = abs(dpCrrnt[2,0])                         # dpz, g*cm/sec 
#      if deltaPz_c[n,i] <= 0.: 
#         print 'deltaPz_c[%2d,%2d] = %e' % (n,i,deltaPz_c[n,i])
      deltaEnrgIon_c[n,i] = (dpCrrnt[0,0]**2+dpCrrnt[1,0]**2+dpCrrnt[2,0]**2)* \
                            deFactor/eVtoErg                     # eV
#
# Transferred momenta and energy along the track - "Magnus expansion" approach:
#
      deltaPx_m[n,i] = abs(dpCrrnt[0,1]*dpFactor)                # dpx, g*cm/sec
#      if deltaPx_m[n,i] <= 0.: 
#         print 'deltaPx_m[%2d,%2d] = %e' % (n,i,deltaPx_m[n,i])
      deltaPy_m[n,i] = abs(dpCrrnt[1,1]*dpFactor) 
#      if deltaPy_m[n,i] <= 0.: 
#         print 'deltaPy_m[%2d,%2d] = %e' % (n,i,deltaPy_m[n,i])
      deltaPz_m[n,i] = abs(dpCrrnt[2,1]*dpFactor) 
#      if deltaPz_m[n,i] <= 0.: 
#         print 'deltaPz_m[%2d,%2d] = %e' % (n,i,deltaPz_m[n,i])
      deltaEnrgIon_m[n,i] = (dpCrrnt[0,1]**2+dpCrrnt[1,1]**2+dpCrrnt[2,1]**2)* \
                            deFactor/eVtoErg                     # eV
# Comparison of the approaches:
      dEion_c_m[n,i] = deltaEnrgIon_c[n,i]/deltaEnrgIon_m[n,i]-1.
#
# Integration using Simpson method:
#
      if (n > 0):
         frctnForce_cSM[i] +=  .5*(deltaEnrgIon_c[n,i]+deltaEnrgIon_c[n-1,i])* \
                               .5*(rhoInit[n,i]+rhoInit[n-1,i])* \
                               (rhoInit[n,i]-rhoInit[n-1,i])*n_e/100  # eV/m 
         frctnForce_mSM[i] +=  .5*(deltaEnrgIon_m[n,i]+deltaEnrgIon_m[n-1,i])* \
                               .5*(rhoInit[n,i]+rhoInit[n-1,i])* \
                               (rhoInit[n,i]-rhoInit[n-1,i])*n_e/100  # eV/m 

#
# Plotting of the tests:
#
print 'indxTestMax = %d' % indxTestMax
nn=np.arange(0,indxTestMax-1,1)

fig2020=plt.figure (2020)
plt.plot(nn,C1test[0:indxTestMax-1],'.r')
plt.xlabel('Points of Tracks',color='m',fontsize=16)
plt.ylabel('$C1$, $cm^2$',color='m',fontsize=16)
plt.title('$C1=[x_{gc}^2+y_{gc}^2+z_e^2+2J/(m_e \cdot \Omega_e)]^{0.5}$', \
          color='m',fontsize=16)
plt.xlim([-5000,indxTestMax+5000])
plt.grid(True)

fig2030=plt.figure (2030)
plt.plot(nn,1.e-5*C2test[0:indxTestMax-1],'.r')
plt.xlabel('Points of Tracks',color='m',fontsize=16)
plt.ylabel('$C2$, $\cdot 10^5$ $cm^2/s$',color='m',fontsize=16)
plt.title('$C2=2\cdot[V_{ix}\cdot(x_i-x_{gc})+V_{iy}\cdot(y_i-y_{gc})+(V_{iz}-V_{ez})\cdot(z_i-z_e)]$', \
          color='m',fontsize=16)
plt.xlim([-5000,indxTestMax+5000])
plt.grid(True)

fig2040=plt.figure (2040)
plt.plot(nn,1e-11*C3test[0:indxTestMax-1],'.r')
plt.xlabel('Points of Tracks',color='m',fontsize=16)
plt.ylabel('$C3$, $\cdot 10^{11}$ $cm^2/s^2$',color='m',fontsize=16)
plt.title('$C3=V_{ix}^2+V_{iy}^2+(V_{iz}-V_{ez})^2$',color='m',fontsize=16)
plt.xlim([-5000,indxTestMax+5000])
plt.grid(True)

fig2050=plt.figure (2050)
plt.plot(nn,b_ME[0:indxTestMax-1],'.r')
plt.xlabel('Points of Tracks',color='m',fontsize=16)
plt.ylabel('$b_{ME}$, $cm$',color='m',fontsize=16)
plt.title('Distance $b_{ME}$ between Particles for "Magnus Expansion" Approach', color='m',fontsize=16)
plt.text(3500,.425,('$b_{ME}=[C1+C2\cdot \delta t +C3 \cdot \delta t^2]^{0.5} (\delta t=%e$ $s)$' % timeStep_c), \
         color='m',fontsize=16)
plt.xlim([-5000,indxTestMax+5000])
plt.grid(True)

fig2055=plt.figure (2055)
plt.plot(nn,b_gc[0:indxTestMax-1],'.r')
plt.xlabel('Points of Tracks',color='m',fontsize=16)
plt.ylabel('$b_{GC}$, $cm$',color='m',fontsize=16)
plt.title('Distance $b_{GC}$ between Particles for "Guiding Center" Approach', color='m',fontsize=16)
plt.text(0,.425,'$b_{GC}=[(x_i-x_{gc})^2+(y_i-y_{gc})^2+(z_i-z_e)^2+2J/(m_e \cdot \Omega_e)]^{0.5}$', \
         color='m',fontsize=16)
plt.xlim([-5000,indxTestMax+5000])
plt.grid(True)

#
# Comparison of bCrrnt_c from "classical" with bTest from "Magnus expansion" approaches:
#
bCrrnt_cTest = np.zeros(indxTestMax)
bCrrnt_cTestRel = np.zeros(indxTestMax)
b_gc_ME_rel = np.zeros(indxTestMax)
print '  k      bCrrnt[2k]  bCrrnt[2k+1]   bCrrntTest       b_gc       b_ME      b_gc/b_ME  actn_gc_ME_rel\n' 
for k in range(indxTestMax):
   bCrrnt_cTest[k] = .5*(bCrrnt_c[2*k]+bCrrnt_c[2*k+1])
#   bCrrnt_cTestRel[k] = bCrrnt_cTest[k]/b_ME[k]
   b_gc_ME_rel[k] = b_gc[k]/b_ME[k]
   actn_gc_ME_rel[k] = 1.e6*(action_gc[k]/action_ME[k]-1)
   if (k < 50):
      print '%6d: %e %e %e %e %e %e %e' % \
            (k,bCrrnt_c[k],bCrrnt_c[k+1],bCrrnt_cTest[k],b_gc[k],b_ME[k],b_gc_ME_rel[k],actn_gc_ME_rel[k])
   
fig2060=plt.figure (2060)
# plt.semilogy(nn,bCrrnt_cTest[0:indxTestMax-1],'.r')
plt.plot(nn,bCrrnt_cTest[0:indxTestMax-1],'.r')
plt.xlabel('Points of Tracks',color='m',fontsize=16)
plt.ylabel('Test $b_{crrntTest}$, $cm$',color='m',fontsize=16)
plt.title('Test $b_{crrntTest} = .5 \cdot [b_{crrnt}(k)+b_{crrnt}(k+1)]$',color='m',fontsize=16)
plt.xlim([-5000,indxTestMax+5000])
# plt.ylim([.9*min(bCrrnt_cTest),1.1*max(bCrrnt_cTest)])
plt.grid(True)

fig2070=plt.figure (2070)
# plt.semilogy(nn,b_gc_ME_rel[0:indxTestMax-1],'.r')
plt.plot(nn,b_gc_ME_rel[0:indxTestMax-1],'.r')
plt.xlabel('Points of Tracks',color='m',fontsize=16)
plt.ylabel('$b_{GC}/b_{ME}$',color='m',fontsize=16)
plt.title('Comparison of Distances $b_{GC}$ and $b_{ME}$ between Particles',color='m',fontsize=16)
plt.xlim([-5000,indxTestMax+5000])
# plt.ylim([.9*min(b_gc_ME_rel),1.1*max(b_gc_ME_rel)])
plt.grid(True)

fig2080=plt.figure (2080)
# plt.semilogy(nn,actn_gc_ME_rel[0:indxTestMax-1],'.r')
plt.plot(nn,actn_gc_ME_rel[0:indxTestMax-1],'.r')
plt.xlabel('Points of Tracks',color='m',fontsize=16)
plt.ylabel('$10^6\cdot (J_{GC}/J_{ME}$ $-$ $1)$',color='m',fontsize=16)
plt.title('Comparison of Actions $J_{GC}$ and $J_{ME}$',color='m',fontsize=16)
plt.xlim([-5000,indxTestMax+5000])
plt.ylim([.9*min(actn_gc_ME_rel),1.1*max(actn_gc_ME_rel)])
plt.grid(True)


nn=np.arange(0,nVion*nImpctPrmtr,1)
halfLintrTest = np.zeros(nVion*nImpctPrmtr)
for i in range(nVion):
   for n in range(nImpctPrmtr):
     halfLintrTest[nVion*i+n] = halfLintr[i,n] 

fig2090=plt.figure (2090)
plt.semilogy(nn,halfLintrTest,'.r')
plt.xlabel('Points of Tracks',color='m',fontsize=16)
plt.ylabel('$0.5 \cdot L_{Intrctn}$, $cm$',color='m',fontsize=16)
plt.title('Total Length of Interaction: $L_{Intrctn}=2 \cdot [R_{max}^2-rho_{Init}^2)]^{0.5}$',color='m',fontsize=16)
plt.xlim([-100,nVion*nImpctPrmtr+100])
plt.ylim([.9*min(halfLintrTest),1.1*max(halfLintrTest)])
plt.grid(True)

#
# Fitting for figures with deltaEnrgIon_c (my own Least Squares Method - LSM;
# Python has own routine for LSM - see site  
#     http://scipy-cookbook.readthedocs.io/items/FittingData.html):
#
#
# Fittied function: 
#   
#  deltaEnrgIon = 10^fitA * rho^fitB,
#  so that
#
# log10(deltaEnrgIon) = fitB*log10(rho) + fitA
#

maxIndx = np.zeros(nVion)
fitA1 = np.zeros(nVion)
fitB1 = np.zeros(nVion)
fitA2 = np.zeros(nVion)
fitB2 = np.zeros(nVion)

#
# Preparing of the initial data:
# 
log10rhoInit = np.zeros((nImpctPrmtr,nVion))
log10deltaEnrgIon_c = np.zeros((nImpctPrmtr,nVion))

for i in range(nVion):
   indx = 0
   for n in range(nImpctPrmtr):
      if ((rhoInit[n,i] > 0.) and (deltaEnrgIon_c[n,i] > 0.)):
         log10rhoInit[indx,i] = np.log10(rhoInit[n,i])
         log10deltaEnrgIon_c[indx,i] = np.log10(deltaEnrgIon_c[n,i])
         indx += 1
   maxIndx[i] = indx
#    print 'maxIndx(%d) = %d' % (i,maxIndx[i]-1)

#
# Minimized functional:
#
# Func1 = {log10(deltaEnrgIon_c) - [fitB*log10(rho) + fitA]}^2
#
sumRho = np.zeros(nVion)
sumRho2 = np.zeros(nVion)
sumEnrg = np.zeros(nVion) 
sumRhoEnrg = np.zeros(nVion)
for i in range(nVion):
   for n in range(int(maxIndx[i])):
      sumRho[i] += log10rhoInit[n,i]
      sumRho2[i] += log10rhoInit[n,i]**2
      sumEnrg[i] += log10deltaEnrgIon_c[n,i]
      sumRhoEnrg[i] += log10rhoInit[n,i]*log10deltaEnrgIon_c[n,i]

   delta = maxIndx[i]*sumRho2[i]-sumRho[i]**2
   fitA1[i] = (sumRho2[i]*sumEnrg[i]-sumRho[i]*sumRhoEnrg[i])/delta
   fitB1[i] = (maxIndx[i]*sumRhoEnrg[i]-sumRho[i]*sumEnrg[i])/delta
#    print 'fitA1(%d) = %e, fitB1(%d) = %e' % (i,fitA1[i],i,fitB1[i])

rhoInitFit1 = np.zeros((maxIndx[0],nVion))
deltaEnrgIon_c_Fit1 = np.zeros((maxIndx[0],nVion))
funcHi2_1 = np.zeros(nVion)
for i in range(nVion):
   factorA = math.pow(10.,fitA1[i])
   for n in range(int(maxIndx[i])):
      rhoInitFit1[n,i] = math.pow(10.,log10rhoInit[n,i])
      deltaEnrgIon_c_Fit1[n,i] = factorA*math.pow(rhoInitFit1[n,i],fitB1[i])
      funcHi2_1[i] += (np.log10(deltaEnrgIon_c[n,i]-np.log10(deltaEnrgIon_c_Fit1[n,i])))**2  
#    print 'i=%2d: fitA1 = %e, fitB1 = %e, hi2_1 = %e' % \
#          (i,fitA1[i],fitB1[i],funcHi2_1[i])
#
# Minimized functional:
#
# Func2 = {1 - [fitB*log10(rho) + fitA]/log10(deltaEnrgIon_c)}^2
#
sumRhoEnrg = np.zeros(nVion)
sumRhoEnrg2 = np.zeros(nVion)
sumRho2Enrg2 = np.zeros(nVion)
sumEnrg = np.zeros(nVion) 
sumEnrg2 = np.zeros(nVion) 
for i in range(nVion):
   for n in range(int(maxIndx[i])):
      sumRhoEnrg[i] += log10rhoInit[n,i]/log10deltaEnrgIon_c[n,i]
      sumRhoEnrg2[i] += log10rhoInit[n,i]/log10deltaEnrgIon_c[n,i]**2
      sumRho2Enrg2[i] += (log10rhoInit[n,i]/log10deltaEnrgIon_c[n,i])**2
      sumEnrg[i] += 1./log10deltaEnrgIon_c[n,i]
      sumEnrg2[i] += 1./log10deltaEnrgIon_c[n,i]**2

   delta = sumEnrg2[i]*sumRho2Enrg2[i]-sumRhoEnrg2[i]**2
   fitA2[i] = (sumRho2Enrg2[i]*sumEnrg[i]-sumRhoEnrg[i]*sumRhoEnrg2[i])/delta
   fitB2[i] = (sumEnrg2[i]*sumRhoEnrg[i]-sumRhoEnrg2[i]*sumEnrg[i])/delta
#    print 'fitA2(%d) = %e, fitB2(%d) = %e' % (i,fitA2[i],i,fitB2[i])

rhoInitFit2 = np.zeros((maxIndx[0],nVion))
deltaEnrgIon_c_Fit2 = np.zeros((maxIndx[0],nVion))
funcHi2_2 = np.zeros(nVion)
for i in range(nVion):
   factorA = math.pow(10.,fitA2[i])
   for n in range(int(maxIndx[i])):
      rhoInitFit2[n,i] = math.pow(10.,log10rhoInit[n,i])
      deltaEnrgIon_c_Fit2[n,i] = factorA*math.pow(rhoInitFit2[n,i],fitB2[i])
      funcHi2_2[i] += (1.-np.log10(deltaEnrgIon_c_Fit2[n,i])/np.log10(deltaEnrgIon_c[n,i]))**2  
#    print 'i=%2d: fitA2 = %e, fitB2 = %e, hi2_2 = %e' % \
#          (i,fitA2[i],fitB2[i],funcHi2_2[i])

#      
# Main plotting:      
#      
fig110=plt.figure (110)
plt.plot(arrayA,arrayB,'.r')
plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
plt.title('Map of Parameters A,B', color='m',fontsize=16)
# plt.xlim([minA,maxA])
# plt.ylim([minB,maxB])
plt.grid(True)

nn=np.arange(0,2*totalPoints-1,1)
fig20=plt.figure (20)
plt.plot(nn,bCrrnt_c[0:2*totalPoints-1],'.r')
# plt.semilogy(nn,bCrrnt_c[0:2*totalPoints-1],'.r')
plt.xlabel('Points of Tracks',color='m',fontsize=16)
plt.ylabel('$b_{Lab.Sys}$, $cm$',color='m',fontsize=16)
plt.title('Distance $b_{Lab.Sys}$ between Particles in Lab.System', color='m',fontsize=16)
plt.xlim([-5000,2*totalPoints+5000])
# plt.xlim([0,2000])
plt.grid(True)

nn=np.arange(0,2*totalPoints-1,1)
fig30=plt.figure (30)
plt.plot(nn,arrayA[0:2*totalPoints-1],'.r',nn,arrayB[0:2*totalPoints-1],'.b')
plt.xlabel('Points of Tracks',color='m',fontsize=16)
plt.ylabel('$A$, $B$',color='m',fontsize=16)
plt.title('$A=log_{10}(q_e^2/b/E_{kin})$, $B=log_{10}(R_{Larm}/b)$', \
          color='m',fontsize=16)
plt.xlim([-5000,2*totalPoints+5000])
# plt.ylim([minB,maxB])
plt.grid(True)

xVionRel = np.zeros((nImpctPrmtr,nVion))
for i in range(nVion):
   for n in range(nImpctPrmtr):
       xVionRel[n,i] = VionRel[i]

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
plt.text(2.e-3,-.021,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
plt.plot([relVeLong,relVeLong],yLimit,'--m',linewidth=1)
plt.text(3.9e-5,.05,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)


'''
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
#      print'n_c=%2d: xRhoInitPx_c = %e, yDeltaPx_c = %e' % \
#            (indx_c,xRhoInitPx_c[indx_c],yDeltaPx_c[indx_c]) 
      indx_c += 1
   if deltaPx_m[n,0] > 0.:
      xRhoInitPx_m[indx_c] = rhoInit[n,0]
      yDeltaPx_m[indx_c] = deltaPx_m[n,0]
#      print'n_m=%2d: xRhoInitPx_m = %e, yDeltaPx_m = %e' % \
#            (indx_m,xRhoInitPx_m[indx_m],yDeltaPx_m[indx_m]) 
      indx_m += 1
maxIndx_c = indx_c-1
maxIndx_m = indx_m-1
# print 'maxIndx_c = %d, maxIndx_m = %d' % (maxIndx_c,maxIndx_m)

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
#      print'n_c=%2d: xRhoInitPz_c = %e, yDeltaPz_c = %e' % \
#            (indx_c,xRhoInitPz_c[indx_c],yDeltaPz_c[indx_c]) 
      indx_c += 1
   if deltaPz_m[n,0] > 0.:
      xRhoInitPz_m[indx_c] = rhoInit[n,0]
      yDeltaPz_m[indx_c] = deltaPz_m[n,0]
#      print'n_m=%2d: xRhoInitPz_m = %e, yDeltaPz_m = %e' % \
#            (indx_m,xRhoInitPz_m[indx_m],yDeltaPz_m[indx_m]) 
      indx_m += 1
maxIndx_c = indx_c-1
maxIndx_m = indx_m-1
# print 'maxIndx_c = %d, maxIndx_m = %d' % (maxIndx_c,maxIndx_m)

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


indxFigures = [0,9,12,18,19,23,27,29,31,34,39,49]
numbrFigures = [500,600,630,660,700,730,760,800,830,860,900,1000]
xPos = [.0024,.0024,.0024,.0027,.0028,.0028,.0028,.0028,.0029,.0028,.0029,.0048]
yPos = [6.4e-9,6.7e-9,6.4e-9,5.9e-9,6.2e-9,5.6e-9,5.8e-9,6.3e-9,5.8e-9,5.9e-9,5.8e-9,4.7e-9]

for i in range(12):
   VionCrrnt = V0*VionRel[indxFigures[i]]
   powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
   mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
   plt.figure(numbrFigures[i])
   plt.loglog(rhoInit[0:nImpctPrmtr,indxFigures[i]], \
              deltaEnrgIon_c[0:nImpctPrmtr,indxFigures[i]],'-xr', \
              rhoInitFit1[0:maxIndx[indxFigures[i]]+1,indxFigures[i]], \
	      deltaEnrgIon_c_Fit1[0:maxIndx[indxFigures[i]]+1,indxFigures[i]],'ob', \
              rhoInitFit2[0:maxIndx[indxFigures[i]]+1,indxFigures[i]], \
	      deltaEnrgIon_c_Fit2[0:maxIndx[indxFigures[i]]+1,indxFigures[i]],'or',linewidth=2)
   plt.xlabel('Track Initial Impact Parameter $rho_{Init}$, $cm$',color='m',fontsize=16)
   plt.ylabel('$\Delta E_{ion}$, $eV$', color='m',fontsize=16)
   titleHeader = 'Transferred Energy $\Delta E_{ion}$ to Single Ion for '
   titleHeader += ' $V_{ion}=%5.3f\cdot10^{%2d}$ $cm/s$'
   plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=16)
   plt.xlim([.95*rhoInit[0,indxFigures[i]],1.05*rhoInit[nImpctPrmtr-1,indxFigures[i]]])
   plt.ylim([.9*deltaEnrgIon_c[nImpctPrmtr-1,indxFigures[i]],1.1*deltaEnrgIon_c_Fit2[0,indxFigures[i]]])
   plt.legend(['Calculated Data',('Fitted Data (Func1): B = %5.3f' % abs(fitB1[indxFigures[i]])), \
              ('Fitted Data (Func2): B = B = %5.3f'% abs(fitB2[indxFigures[i]]))], \
             loc='lower center',fontsize=16)
   plt.text(xPos[i],yPos[i],'Fitted $\Delta E_{ion}$ are proportional to $1/rho_{Init}^B$', \
	    color='m',fontsize=16)
   plt.grid(True)



# indxFigures = [0,9,12,18,19,23,27,29,31,34,39,49]
# numbrFigures = [500,600,630,660,700,730,760,800,830,860,900,1000]
yPosText = [-2.115,-2.115,-2.115,-2.20,-2.115,-2.115,-2.115,-2.20,-2.115,-2.115,-2.115,-2.115]
fig1100=plt.figure (1100)
plt.semilogx(VionRel,fitB1,'-xb',VionRel,fitB2,'-xr',linewidth=2)
plt.xlabel('Relative Ion Velocity  $V_{ion}/V_0$',color='m',fontsize=16)
plt.ylabel('$B$', color='m',fontsize=16)
titleHeader = 'Dependence of Transferred Energy $\Delta E_{ion}$ to Single Ion on $rho^B$'
plt.title(titleHeader,color='m',fontsize=16)
plt.text(2.e-4,-1.93,('$V_{e0}=%5.3f\cdot10^{%2d}$cm/s' % (mantV0,powV0)),color='m',fontsize=16)
plt.xlim(xLimit)
yLimit=[min(fitB1[0],fitB2[0])-.025,max(fitB1[nImpctPrmtr-1],fitB2[nImpctPrmtr-1])+.025]
# print 'ylim[0] = %e, ylim[1] = %e' % (yLimit[0],yLimit[1])
plt.ylim(yLimit)
plt.plot([relVeTrnsv,relVeTrnsv],yLimit,'--m',linewidth=1)
plt.text(2.e-3,yLimit[0]-.021,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
plt.plot([relVeLong,relVeLong],yLimit,'--m',linewidth=1)
# plt.text(4.4e-5,yLimit[0]-.02,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)
plt.text(4.4e-5,-2.35,'$ \Delta V_{e||}/ sV_{e0}$',color='m',fontsize=14)
plt.legend(['Func1','Func2'],loc='lower right',fontsize=16)
plt.text(5.4e-4,-2.22,'Figures',color='k',fontsize=14)
for k in range(12):
   xPos = VionRel[indxFigures[k]]
   plt.plot([xPos,xPos],[-2.18,-2.125],'-k',linewidth=1)
   plt.text(xPos,yPosText[k],('%3d' % int(numbrFigures[k])),color='k',fontsize=10)
plt.grid(True)


# plt.show()

# sys.exit()

#----------------------------------------------------
#
# Integration of the transferred energy to ion from electron beam
#----------------------------------------------------
#
#  "Gauss-Kronrod" method of integration (GK)
#
#
# Points (psi_i) and weigths (w_i) to integrate for interval from -1 to 1;
# These data are from William H. Beyer. "Handbook of Mathematical Science".
# 5th Edition, CRC Press, Inc, 1978.
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
#----------------------------------------------------
#
# Data for GK:
#
#----------------------------------------------------

frctnForce_cGK = np.zeros(nVion)   # integration using "Gauss-Kronrod" method 
frctnForce_mGK = np.zeros(nVion)   # integration using "Gauss-Kronrod" method 

psi16=np.array([-0.9894009, -0.9445750, -0.8656312, -0.7554044, -0.6178762, \
                -0.4580168, -0.2816036, -0.0950125,  0.0950125,  0.2816036, \
		 0.4580168,  0.6178762,  0.7554044,  0.8656312,  0.9445750, \
		 0.9894009])
w16  =np.array([ 0.0271525,  0.0622535,  0.0951585,  0.1246290,  0.1495960, \
                 0.1691565,  0.1826034,  0.1894506,  0.1894506,  0.1826034, \
		 0.1691565,  0.1495960,  0.1246290,  0.0951585,  0.0622535, \
		 0.0271525])

nPointsGK = 16
rhoCrrntGK = np.zeros((nPointsGK,nVion))
deltaEnrgIon_cGK = np.zeros((nPointsGK,nVion))
deltaPx_cGK = np.zeros((nPointsGK,nVion))  
deltaPy_cGK = np.zeros((nPointsGK,nVion))  
deltaPz_cGK = np.zeros((nPointsGK,nVion))  
deltaEnrgIon_mGK = np.zeros((nPointsGK,nVion))
deltaPx_mGK = np.zeros((nPointsGK,nVion))  
deltaPy_mGK = np.zeros((nPointsGK,nVion))  
deltaPz_mGK = np.zeros((nPointsGK,nVion))  

indx = 0
totalPointsIntgrtn = 0
for i in range(nVion):
   rhoMaxCrrnt = impctPrmtrMax[i]* \
   np.sqrt(1.- (pi*larmorTurns*eVrmsLong/omega_L/impctPrmtrMax[i])**2)
   rhoMaxCrrnt = impctPrmtrMax[i]
#   log10rhoMax = math.log10(rhoMaxCrrnt)
#   log10rhoStep = (log10rhoMax-log10rhoMin)/(nPointsGK)
   print 'Vion(%d) = %e, rhoMin = %e, rhoMaxCrrnt = %e' % \
         (i,Vion[i],rhoMin,rhoMaxCrrnt)
   for n in range(nPointsGK):
#      log10rhoCrrntGK = log10rhoMin+(n+0.5)*log10rhoStep 
#      rhoCrrntGK = math.pow(10.,log10rhoCrrntGK)
      rhoCrrntGK[n,i] = psi16[n]*(rhoMaxCrrnt-rhoMin)/2 + \
                       (rhoMaxCrrnt+rhoMin)/2
#      print '    rhoCrrntGK[%2d,%2d] = %e' % (n,i,rhoCrrntGK[n,i])
      halfLintrCrrnt = np.sqrt(rhoMaxCrrnt**2-rhoCrrntGK[n,i]**2)   # half length of interaction; cm
      timeHalfPath = halfLintrCrrnt/eVrmsLong     # 0.5 time of interaction; sec
      numbLarmor = int(2.*timeHalfPath/T_larm)             
      pointAlongTrackCrrnt = int(2.*timeHalfPath/timeStep_c)
      totalPointsIntgrtn += pointAlongTrackCrrnt
      z_ionCrrnt = np.zeros(6)      # Zeroing out of vector for ion
      z_elecCrrnt = np.zeros(6)     # Zeroing out of vector for electron
# Zeroing out of "guiding center" vector for electron:
      z_elecCrrnt_gc = np.zeros(6)  
# Current values to transfered momemta 
# (second index numerates "classical", if 0, and 
#                         "Magnus expand" (if 1) approaches: 
      dpCrrnt = np.zeros((3,2))
# Intermediate arrays:
      dpIon_c = np.zeros(3) 
      dpIon_m = np.zeros(3) 
      dpElec_c = np.zeros(3) 
      dpElec_m = np.zeros(3) 
# Current initial vector for electron:
      z_elecCrrnt[Ix] = rhoCrrntGK[n,i]              # x, cm
      z_elecCrrnt[Iz] = -halfLintrCrrnt              # z, cm
      z_elecCrrnt[Ipy] = m_elec*eVrmsTran            # py, g*cm/sec
      z_elecCrrnt[Ipz] = m_elec*eVrmsLong            # pz, g*cm/sec
# transfer to system of guiding center:
      z_elecCrrnt_gc=toGuidingCenter(z_elecCrrnt)  
#
# Main loop along the each track:
#
      for k in range(int(pointAlongTrackCrrnt)):
#
# Dragging both particles through first half of the step of the track:
#
  	 z_elecCrrnt_gc = np.dot(matr_elec_c,z_elecCrrnt_gc) # electron
  	 z_ionCrrnt = np.dot(matr_ion_c,z_ionCrrnt)          # ion
# transfer from system of guiding center: 
	 z_elecCrrnt=fromGuidingCenter(z_elecCrrnt_gc)     
# Current distance between ion and electron; cm:
 	 bCrrnt_c[indx]=np.sqrt((z_ionCrrnt[0]-z_elecCrrnt[0])**2+ \
	                  (z_ionCrrnt[2]-z_elecCrrnt[2])**2+ \
			  (z_ionCrrnt[4]-z_elecCrrnt[4])**2)
         indx += 1
#
# Dragging both particles through interaction during this step of track
# (for both approaches):
#
#    "Classical":
	 dpIon_c,dpElec_c,action,b_gc_c = \
	         guidingCenterCollision(z_elecCrrnt_gc,z_ionCrrnt,timeStep_c) 
#    "Magnus Expand":
	 dpIon_m,dpElec_m,action,dy_gc_m,C1,C2,C3,b,D1,D2,q = \
	         MagnusExpansionCollision(z_elecCrrnt_gc,z_ionCrrnt,timeStep_c) 
#
# Taking into account transfer of momentum for both particles
# (for "classical" only):
#
         if (dpTransferFlag == 1):
	    for ic in range(3):
	       z_ionCrrnt[2*ic+1] += dpIon_c[ic]   
	       z_elecCrrnt[2*ic+1] += dpElec_c[ic]
# transfer to system of guiding center:
         z_elecCrrnt_gc=toGuidingCenter(z_elecCrrnt)  
# Accumulation of the transfered momenta along the track:  
         for ic in range(3):
#	    if i == 0:
#	       print 'dpIon_c[%2d] = %20.14e, dpIon_m[%2d] = %20.14e' % \
#	             (ic,dpIon_c[ic],ic,dpIon_m[ic])
 	    dpCrrnt[ic,0] += dpIon_c[ic]                      # g*cm/csec  
 	    dpCrrnt[ic,1] += dpIon_m[ic]                      # g*cm/csec  
#
# Dragging both particles through second half of the step of the track:
#
 	 z_elecCrrnt_gc = np.dot(matr_elec_c,z_elecCrrnt_gc)     # electron
 	 z_ionCrrnt = np.dot(matr_ion_c,z_ionCrrnt)              # ion
# transfer from system of guiding center: 
	 z_elecCrrnt=fromGuidingCenter(z_elecCrrnt_gc)     
# Current distance between ion and electron; cm:
 	 bCrrnt_c[indx]=np.sqrt((z_ionCrrnt[0]-z_elecCrrnt[0])**2+ \
	                  (z_ionCrrnt[2]-z_elecCrrnt[2])**2+ \
			  (z_ionCrrnt[4]-z_elecCrrnt[4])**2)
         indx += 1
#	 
# Transferred momenta and energy along the entire length of each track
# for both approaches:
#
      deltaPx_cGK[n,i] = abs(dpCrrnt[0,0]) 
      deltaPx_mGK[n,i] = abs(dpCrrnt[0,1]) 
#      if deltaPx_cGK[n,i] <= 0.: 
#         print 'deltaPx_cGK[%2d,%2d] = %e, dpCrrnt[%2d,%2d] = %e' % \
#	       (n,i,deltaPx_cGK[n,i],n,i,dpCrrnt[0,0])
#      if deltaPx_mGK[n,i] <= 0.: 
#         print 'deltaPx_mGK[%2d,%2d] = %e, dpCrrnt[%2d,%2d] = %e' % \
#	       (n,i,deltaPx_mGK[n,i],n,i,dpCrrnt[0,1])
      deltaPy_cGK[n,i] = abs(dpCrrnt[1,0]) 
      deltaPy_mGK[n,i] = abs(dpCrrnt[1,1]) 
#      if deltaPy_cGK[n,i] <= 0.: 
#         print 'deltaPy_cGK[%2d,%2d] = %e' % (n,i,deltaPy_cGK[n,i])
#      if deltaPy_mGK[n,i] <= 0.: 
#         print 'deltaPy_mGK[%2d,%2d] = %e' % (n,i,deltaPy_mGK[n,i])
      deltaPz_cGK[n,i] = abs(dpCrrnt[2,0]) 
      deltaPz_mGK[n,i] = abs(dpCrrnt[2,1]) 
#      if deltaPz_cGK[n,i] <= 0.: 
#         print 'deltaPz_cGK[%2d,%2d] = %e' % (n,i,deltaPz_cGK[n,i])
#      if deltaPz_mGK[n,i] <= 0.: 
#         print 'deltaPz_mGK[%2d,%2d] = %e' % (n,i,deltaPz_MGK[n,i])
      deltaEnrgIon_cGK[n,i] = (dpCrrnt[0,0]**2+dpCrrnt[1,0]**2+dpCrrnt[2,0]**2)* \
                            deFactor/eVtoErg                      # eV
      deltaPx_mGK[n,i] = abs(dpCrrnt[0,1]*dpFactor) 
#      if deltaPx_mGK[n,i] <= 0.: 
#         print 'deltaPx_mGK[%2d,%2d] = %e' % (n,i,deltaPx_mGK[n,i])
      deltaPy_mGK[n,i] = abs(dpCrrnt[1,1]*dpFactor) 
#      if deltaPy_mGK[n,i] <= 0.: 
#         print 'deltaPy_mGK[%2d,%2d] = %e' % (n,i,deltaPy_mGK[n,i])
      deltaPz_mGK[n,i] = abs(dpCrrnt[2,1]*dpFactor) 
#      if deltaPz_mGK[n,i] <= 0.: 
#         print 'deltaPz_mGK[%2d,%2d] = %e' % (n,i,deltaPz_mGK[n,i])
      deltaEnrgIon_mGK[n,i] = (dpCrrnt[0,1]**2+dpCrrnt[1,1]**2+dpCrrnt[2,1]**2)* \
                            deFactor/eVtoErg                      # eV
#
# Integration using "Gauss-Kronrod" method:
#
      frctnForce_cGK[i] += .5*(rhoMaxCrrnt-rhoMin)*w16[n]*deltaEnrgIon_cGK[n,i]* \
                          rhoCrrntGK[n,i]*n_e/100.                # eV/m
      frctnForce_mGK[i] += .5*(rhoMaxCrrnt-rhoMin)*w16[n]*deltaEnrgIon_mGK[n,i]* \
                          rhoCrrntGK[n,i]*n_e/100.                # eV/m

print '\ni          V_i           V_i/V0         ff_cSM          ff_mSM           ff_cGK           ff_mGK\n' 
for i in range(nVion):
   print '%2d    %10.6e    %10.6e    %10.6e    %10.6e    %10.6e    %10.6e' % \
         (i,Vion[i],Vion[i]/V0,frctnForce_cSM[i],frctnForce_mSM[i],frctnForce_cGK[i],frctnForce_mGK[i])

xVionRelIntgr = np.zeros((nPointsGK,nVion))
for i in range(nVion):
   for n in range(nPointsGK):
       xVionRelIntgr[n,i] = VionRel[i]

fig41=plt.figure (41)
for i in range(nVion):
    plt.semilogx(xVionRelIntgr[0:nPointsGK,i],rhoCrrntGK[0:nPointsGK,i],'.r')
plt.xlabel('Relative Ion Velocity, $V_i/V_{e0}$',color='m',fontsize=16)
plt.ylabel('$rho_{Init}$, cm',color='m',fontsize=16)
plt.title('Subdivisions for $rho_{Init}$ for Integration: Gauss-Kronrod Method', \
          color='m',fontsize=16)
plt.grid(True)
yLimit=[0.,max(rhoCrrntGK[0:nPointsGK,nVion-1])+.01]
plt.ylim(yLimit)
plt.plot([relVeTrnsv,relVeTrnsv],yLimit,'--m',linewidth=1)
plt.text(2.e-3,-.02,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
plt.plot([relVeLong,relVeLong],yLimit,'--m',linewidth=1)
plt.text(4.4e-5,.05,'$ \Delta V_{e||}/ V_{e0}$',color='m',fontsize=14)


indxFigures = [0,9,12,18,19,23,27,29,31,34,39,49]

'''
for i in range(2,10,1):
   plt.figure (2000+indxFigures[i])
   plt.loglog(rhoInit[:,indxFigures[i]],deltaEnrgIon_c[:,indxFigures[i]],'-xr', \
              rhoInit[:,indxFigures[i]],deltaEnrgIon_m[:,indxFigures[i]],'-xb',linewidth=2)
   plt.xlabel('Track Initial Impact Parameter $rho_{Init}$, $cm$',color='m',fontsize=16)
   plt.ylabel('$\Delta E_{ion}$, $eV$', color='m',fontsize=16)
   plt.title('Transferred Energy $\Delta E_{ion}$ to Single Ion', color='m',fontsize=16)
   plt.xlim([.95*rhoInit[0,indxFigures[i]],1.1*rhoInit[nImpctPrmtr-1,indxFigures[i]]])
   yLimit = [0.5*min(deltaEnrgIon_c[nImpctPrmtr-1,indxFigures[i]],deltaEnrgIon_m[nImpctPrmtr-1,indxFigures[i]]), \
             1.5*max(deltaEnrgIon_c[0,indxFigures[i]],deltaEnrgIon_m[0,indxFigures[i]])]
   plt.ylim(yLimit)
   plt.grid(True)
'''

for i in range(2,10,1):
   VionCrrnt = V0*VionRel[indxFigures[i]]
   powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
   mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
   plt.figure (3000+indxFigures[i])
   plt.semilogx(rhoInit[:,indxFigures[i]],dEion_c_m[:,indxFigures[i]],'-xr',linewidth=2)
   plt.xlabel('Track Initial Impact Parameter $rho_{Init}$, $cm$',color='m',fontsize=16)
   plt.ylabel('Difference of $\Delta E_{ion}$: "GC"/"ME" - 1', color='m',fontsize=16)
   titleHeader = 'Comparison of Approaches for $\Delta E_{ion}$:'
   titleHeader += ' $V_{ion}=%5.3f\cdot10^{%2d}$ $cm/s$'
   plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=16)
   plt.xlim([.95*rhoInit[0,indxFigures[i]],1.1*rhoInit[nImpctPrmtr-1,indxFigures[i]]])
   plt.ylim([min(dEion_c_m[:,indxFigures[i]])-.0005,max(dEion_c_m[:,indxFigures[i]])+.0005])
   plt.grid(True)

plt.show()

'''
# fig10.savefig('picturesCMA/correctedRmax_fig10cma.jpg')    
fig20.savefig('picturesCMA/particleDistance_ls_fig20cma.jpg')    
fig30.savefig('picturesCMA/parametersA-B_fig30cma.jpg')    
fig40.savefig('picturesCMA/initialImpactParameter_SM_fig40cma.jpg')    
fig41.savefig('picturesCMA/initialImpactParameter_GK_fig41cma.jpg')    
fig110.savefig('picturesCMA/mapA-B_fig110cma.jpg')    
fig209.savefig('picturesCMA/rDebye_rLikeDebye_rPass_fig209cma.jpg')    
fig3151.savefig('picturesCMA/impctPrmtr_fig3151cma.jpg')    
fig500.savefig('picturesCMA/deltaEtransf_indxPlot-0_fig500cma.jpg')    
fig600.savefig('picturesCMA/deltaEtransf_indxPlot-fig9_600cma.jpg')    
fig630.savefig('picturesCMA/deltaEtransf_indxPlot-12_fig630cma.jpg')    
fig660.savefig('picturesCMA/deltaEtransf_indxPlot-18_fig660cma.jpg')    
fig700.savefig('picturesCMA/deltaEtransf_indxPlot-19_fig700cma.jpg')    
fig730.savefig('picturesCMA/deltaEtransf_indxPlot-23_fig730cma.jpg')    
fig760.savefig('picturesCMA/deltaEtransf_indxPlot-27_fig760cma.jpg')    
fig800.savefig('picturesCMA/deltaEtransf_indxPlot-29_fig800cma.jpg')    
fig830.savefig('picturesCMA/deltaEtransf_indxPlot-31_fig830cma.jpg')    
fig860.savefig('picturesCMA/deltaEtransf_indxPlot-34_fig860cma.jpg')    
fig900.savefig('picturesCMA/deltaEtransf_indxPlot-39_fig900cma.jpg')    
fig1000.savefig('picturesCMA/deltaEtransf_indxPlot-49_fig1000cma.jpg')    
fig1100.savefig('picturesCMA/exponentB_on_ionVelocity_fig1100cma.jpg')    
fig2020.savefig('picturesCMA/magnusExpansion_C1_fig2020cma.jpg')    
fig2030.savefig('picturesCMA/magnusExpansion_C2_fig2030cma.jpg')    
fig2040.savefig('picturesCMA/magnusExpansion_C3_fig2040cma.jpg')    
fig2050.savefig('picturesCMA/particleDistance_me_fig2050cma.jpg')    
fig2055.savefig('picturesCMA/particleDistance_gc_fig2055cma.jpg')    
fig2070.savefig('picturesCMA/particleDistanceComprsn_gc_me_fig2070cma.jpg')    
fig2080.savefig('picturesCMA/actionComprsn_gc_me_fig2080cma.jpg')    
fig2090.savefig('picturesCMA/totalLengthIntrsctn_fig2090cma.jpg')    
fig3012.savefig('picturesCMA/deltaEcomprsn_indxPlot-12_fig3012cma.jpg')    
fig3018.savefig('picturesCMA/deltaEcomprsn_indxPlot-18_fig3018cma.jpg')    
fig3019.savefig('picturesCMA/deltaEcomprsn_indxPlot-19_fig3019cma.jpg')    
fig3023.savefig('picturesCMA/deltaEcomprsn_indxPlot-23_fig3023cma.jpg')    
fig3027.savefig('picturesCMA/deltaEcomprsn_indxPlot-27_fig3027cma.jpg')    
fig3029.savefig('picturesCMA/deltaEcomprsn_indxPlot-29_fig3029cma.jpg')    
fig3031.savefig('picturesCMA/deltaEcomprsn_indxPlot-31_fig3031cma.jpg')    
fig3034.savefig('picturesCMA/deltaEcomprsn_indxPlot-34_fig3034cma.jpg')    
'''

'''
#
# Opening the output file: 
#
apprchClsscl_file='resultsClassicalApproach.dat'
print 'Open output file "%s"...' % apprchClsscl_file
apprchClsscl_flag=0
try:
   outfile = open(apprchClsscl_file,'w')
   apprchClsscl_file_flag=1
except:
   print 'Problem to open output file "%s"' % apprchClsscl_file

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
   
sys.exit()

