# from __future__ import division

#-------------------------------------
#
#        Started at 06/08/2018 (YuE)
#
# This script based on the previous script
# threeApproachesComparison_v6.py
#
# 06/25/2018: IT IS NOT FINISHED ("classical" approach works appropriate)!
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

print 'V0=%e cm/s, rmsTrnsvVe=%e cm/s, rmsLongVe=%e cm/s' % \
       (V0,rmsTrnsvVe,rmsLongVe)

ergToEV = 1./1.60218e-12

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
# output - transfered momenta to ion and electron and electron y_gc coordinate: 
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
   C1=np.sqrt((vectrIon[0]-x_gc)**2+ \
              (vectrIon[2]-vectrElec_gc[2])**2+ \
              (vectrIon[4]-vectrElec_gc[4])**2+2.*action/mOmegaLarm)   # cm^2
   C2=2.*((vectrIon[0]-x_gc)*vectrIon[1]/M_ion+ \
          (vectrIon[2]-vectrElec_gc[2])*vectrIon[3]/M_ion+ \
          (vectrIon[4]-vectrElec_gc[4])* \
	  (vectrIon[5]/M_ion-vectrElec_gc[5]/m_elec))  # cm^2/sec
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
   return dpIon,dpElec,action,dy_gc                                      

sphereNe=3.
R_e=math.pow(sphereNe/n_e,1./3)                       # cm
print 'R_e (cm)=%e' % R_e

ro_Larm = eVrmsTran/omega_L                           # cm
print 'ro_Larm (cm)=%e' % ro_Larm

impctPrmtrMin=2.*ro_Larm

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
plt.title(titleHeader % (mantV0,powV0,fieldB[0]),color='m',fontsize=16)
plt.xlim(xLimit)
yLimit=[8.e-4,.6]
plt.ylim(yLimit)
plt.plot([relVeTrnsv,relVeTrnsv],yLimit,'--m',linewidth=1)
plt.text(2.e-3,5.9e-4,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
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
# Factor to calculate transferred energy to ion
# (the friction force is defined by this transfered energy): 
#
deFactor = 0.5/M_ion                               # 1/g

indx = 0
for i in range(nVion):
   rhoMax[i] = impctPrmtrMax[i]* \
   np.sqrt(1.- (pi*larmorTurns*eVrmsLong/omega_L/impctPrmtrMax[i])**2)
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
# (second index numerates "classical", if 0, and 
#                         "Magnus expand" (if 1) approaches: 
      dpCrrnt = np.zeros((3,2))
# Intermediate arrays:
      dpIon_c = np.zeros(3) 
      dpIon_m = np.zeros(3) 
      dpElec_c = np.zeros(3) 
      dpElec_m = np.zeros(3) 
# Current initial vector for electron:
      z_elecCrrnt[Ix] = rhoCrrnt                     # x, cm
      z_elecCrrnt[Iz] = -halfLintr[n,i]                   # z, cm
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
#    "Magnus Expand":
	 dpIon_m,dpElec_m,action,dy_gc_m = \
	         MagnusExpansionCollision(z_elecCrrnt_gc,z_ionCrrnt,timeStep_c) 
####
#### Taking into account transfer of momentum for both particles:
####
####	     for ic in range(3):
####	        z_ionCrrnt[2*ic+1] += dpIon[ic]   
####	        z_elecCrrnt[2*ic+1] += dpElec[ic]
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
# Current values of parameters A,B:  
	 arrayA[indx] = math.log10(ro_Larm/bCrrnt_c[indx])     
	 arrayB[indx] = math.log10((q_elec**2/bCrrnt_c[indx])/kinEnergy)
         indx += 1
# Transferred momenta and energy along the track for both approaches:
      deltaPx_c[n,i] = abs(dpCrrnt[0,0]) 
#      if deltaPx_c[n,i] <= 0.: 
#         print 'deltaPx_c[%2d,%2d] = %e, dpCrrnt[%2d,%2d] = %e' % \
#	       (n,i,deltaPx_c[n,i],n,i,dpCrrnt[0,0])
      deltaPy_c[n,i] = abs(dpCrrnt[1,0]) 
#      if deltaPy_c[n,i] <= 0.: 
#         print 'deltaPy_c[%2d,%2d] = %e' % (n,i,deltaPy_c[n,i])
      deltaPz_c[n,i] = abs(dpCrrnt[2,0]) 
#      if deltaPz_c[n,i] <= 0.: 
#         print 'deltaPz_c[%2d,%2d] = %e' % (n,i,deltaPz_c[n,i])
      deltaEnrgIon_c[n,i] = (dpCrrnt[0,0]**2+dpCrrnt[1,0]**2+dpCrrnt[2,0]**2)* \
                            deFactor/eVtoErg                      # eV
      deltaPx_m[n,i] = abs(dpCrrnt[0,1]*dpFactor) 
#      if deltaPx_m[n,i] <= 0.: 
#         print 'deltaPx_m[%2d,%2d] = %e' % (n,i,deltaPx_m[n,i])
      deltaPy_m[n,i] = abs(dpCrrnt[1,1]*dpFactor) 
#      if deltaPy_m[n,i] <= 0.: 
#         print 'deltaPy_m[%2d,%2d] = %e' % (n,i,deltaPy_m[n,i])
      deltaPz_m[n,i] = abs(dpCrrnt[2,1]*dpFactor) 
#      if deltaPz_m[n,i] <= 0.: 
#         print 'deltaPz_m[%2d,%2d] = %e' % (n,i,deltaPz_m[n,i])
      deltaEnrgIon_m[n,i] = (dpCrrnt[0,1]**2+dpCrrnt[1,1]**2+dpCrrnt[2,1]**2)* \
                            deFactor/eVtoErg                      # eV

#
# Fitting for figures with deltaEnrgIon_c (my own Least Squares Method - LSM;
# Python has own routine for LSM - see site  
#     http://scipy-cookbook.readthedocs.io/items/FittingData.html):
#
log10rhoInit = np.zeros((nImpctPrmtr,nVion))
log10deltaEnrgIon_c = np.zeros((nImpctPrmtr,nVion))

maxIndx = np.zeros(nVion)
fitA = np.zeros(nVion)
fitB = np.zeros(nVion)

for i in range(nVion):
   indx = 0
   for n in range(nImpctPrmtr):
      if ((rhoInit[n,i] > 0.) and (deltaEnrgIon_c[n,i] > 0.)):
         log10rhoInit[indx,i] = np.log10(rhoInit[n,i])
         log10deltaEnrgIon_c[indx,i] = np.log10(deltaEnrgIon_c[n,i])
         indx += 1
   maxIndx[i] = indx-1
   print 'maxIndx(%d) = %d' % (i,maxIndx[i])


for i in range(nVion):
   sumRho = np.zeros(nVion)
   sumRho2 = np.zeros(nVion)
   sumEnrg = np.zeros(nVion) 
   sumRhoEnrg = np.zeros(nVion)
   for n in range(int(maxIndx[i])):
      sumRho[0] += log10rhoInit[n,i]
      sumRho2[0] += log10rhoInit[n,i]**2
      sumEnrg[0] += log10deltaEnrgIon_c[n,i]
      sumRhoEnrg[0] += log10rhoInit[n,i]*log10deltaEnrgIon_c[n,i]

   delta = maxIndx[i]*sumRho2[0]-sumRho[0]**2
   fitA[i] = (sumRho2[0]*sumEnrg[0]-sumRho[0]*sumRhoEnrg[0])/delta
   fitB[i] = (maxIndx[i]*sumRhoEnrg[0]-sumRho[0]*sumEnrg[0])/delta
   print 'fitA(%d) = %e, fitB(%d) = %e' % (i,fitA[i],i,fitB[i])

rhoInitFit = np.zeros((maxIndx[0],nVion))
deltaEnrgIon_c_Fit = np.zeros((maxIndx[0],nVion))
for i in range(nVion):
   for n in range(int(maxIndx[i])):
      rhoInitFit[n,i] = math.pow(10.,log10rhoInit[n,i])
      deltaEnrgIon_c_Fit[n,i] = math.pow(10.,fitA[i])* \
                             math.pow(rhoInitFit[n,i],fitB[i])

#      
# Plotting:      
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
#   plt.plot(xAboundary,xBboundary,'-xb')
plt.xlabel('Points of Tracks',color='m',fontsize=16)
plt.ylabel('$b_{crrnt}$, cm',color='m',fontsize=16)
plt.title('Distance between Particles $b_{crrnt}$', color='m',fontsize=16)
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

# nn=np.arange(0,nImpctPrmtr*nVion-1,1)
# fig40=plt.figure (40)
# plt.plot(nn,rhoInit[0:nImpctPrmtr*nVion-1],'.r')
# plt.xlabel('Point of Tracks',color='m',fontsize=16)
# plt.ylabel('$rho_{Init}$, cm',color='m',fontsize=16)
# plt.xlim([-100,nImpctPrmtr*nVion+100])
# plt.title('Track Initial Impact Parameter $rho_{Init}$', color='m',fontsize=16)
# plt.grid(True)

xVionRel = np.zeros((nImpctPrmtr,nVion))
for i in range(nVion):
   for n in range(nImpctPrmtr):
       xVionRel[n,i] = VionRel[i]

fig40=plt.figure (40)
for i in range(nVion):
    plt.semilogx(xVionRel[0:nImpctPrmtr,i],rhoInit[0:nImpctPrmtr,i],'.r')
plt.xlabel('Relative Ion Velocity, $V_i/V_{e0}$',color='m',fontsize=16)
plt.ylabel('$rho_{Init}$, cm',color='m',fontsize=16)
plt.title('Track Initial Impact Parameter $rho_{Init}$', color='m',fontsize=16)
plt.grid(True)
yLimit=[0.,.405]
plt.ylim(yLimit)
plt.plot([relVeTrnsv,relVeTrnsv],yLimit,'--m',linewidth=1)
plt.text(2.e-3,-.021,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)


# for i in range(nVion):
#    plt.figure (50+i)
#    plt.loglog(rhoInit[:,i],deltaEnrgIon[:,i],'-r',linewidth=2)
#    plt.xlabel('Track Initial Impact Parameter $rho_{Init}$', \
#               color='m',fontsize=16)
#    plt.ylabel('$\Delta E$, erg', color='m',fontsize=16)
#    plt.title('Transferred Energy $\Delta E$', color='m',fontsize=16)
#    plt.grid(True)
#    time.sleep(5)


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
indxPlot=0
VionCrrnt = V0*VionRel[indxPlot]
powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
fig500=plt.figure (500)
plt.loglog(rhoInit[0:nImpctPrmtr-1,indxPlot], \
           deltaEnrgIon_c[0:nImpctPrmtr-1,indxPlot],'-xr', \
           rhoInitFit[0:maxIndx[indxPlot],indxPlot], \
	   deltaEnrgIon_c_Fit[0:maxIndx[indxPlot],indxPlot],'ob',linewidth=2)
plt.xlabel('Track Initial Impact Parameter $rho_{Init}$', \
           color='m',fontsize=16)
plt.ylabel('$\Delta E_{ion}$, $eV$', color='m',fontsize=16)
titleHeader = 'Transferred Energy $\Delta E_{ion}$ to Single Ion for '
titleHeader += ' $V_{ion}=%5.3f\cdot10^{%2d}$ $cm/s$'
plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=16)
plt.xlim([.95*rhoInit[0,0],1.05*rhoInit[nImpctPrmtr-1,0]])
plt.ylim([3.e-11,2.e-8])
plt.legend(['Calculated Data','Fitted Data'],loc='lower center',fontsize=16)
plt.text(2.5e-3,1.e-8, \
         ('Fitted $\Delta E_{ion}$ are proportional to $1/rho_{Init}^{%5.3f}$' % \
	  abs(fitB[indxPlot])),color='m',fontsize=16)
plt.grid(True)

indxPlot=9
VionCrrnt = V0*VionRel[indxPlot]
powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
fig600=plt.figure (600)
plt.loglog(rhoInit[0:nImpctPrmtr-1,indxPlot], \
           deltaEnrgIon_c[0:nImpctPrmtr-1,indxPlot],'-xr', \
           rhoInitFit[0:maxIndx[indxPlot],indxPlot], \
	   deltaEnrgIon_c_Fit[0:maxIndx[indxPlot],indxPlot],'ob',linewidth=2)
plt.xlabel('Track Initial Impact Parameter $rho_{Init}$', \
           color='m',fontsize=16)
plt.ylabel('$\Delta E_{ion}$, $eV$', color='m',fontsize=16)
titleHeader = 'Transferred Energy $\Delta E_{ion}$ to Single Ion for '
titleHeader += ' $V_{ion}=%5.3f\cdot10^{%2d}$ $cm/s$'
plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=16)
plt.xlim([.95*rhoInit[0,0],1.05*rhoInit[nImpctPrmtr-1,0]])
plt.ylim([8.e-11,2.e-8])
plt.legend(['Calculated Data','Fitted Data'],loc='lower center',fontsize=16)
plt.text(2.5e-3,1.e-8, \
         ('Fitted $\Delta E_{ion}$ are proportional to $1/rho_{Init}^{%5.3f}$' % \
	  abs(fitB[indxPlot])),color='m',fontsize=16)
plt.grid(True)

yLimit=[1.5e-10,2.e-8]

indxPlot=19
VionCrrnt = V0*VionRel[indxPlot]
powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
fig700=plt.figure (700)
plt.loglog(rhoInit[0:nImpctPrmtr-1,indxPlot], \
           deltaEnrgIon_c[0:nImpctPrmtr-1,indxPlot],'-xr', \
           rhoInitFit[0:maxIndx[indxPlot],indxPlot], \
	   deltaEnrgIon_c_Fit[0:maxIndx[indxPlot],indxPlot],'ob',linewidth=2)
plt.xlabel('Track Initial Impact Parameter $rho_{Init}$', \
           color='m',fontsize=16)
plt.ylabel('$\Delta E_{ion}$, $eV$', color='m',fontsize=16)
titleHeader = 'Transferred Energy $\Delta E_{ion}$ to Single Ion for '
titleHeader += ' $V_{ion}=%5.3f\cdot10^{%2d}$ $cm/s$'
plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=16)
plt.xlim([.95*rhoInit[0,0],1.05*rhoInit[nImpctPrmtr-1,0]])
plt.ylim(yLimit)  # for "indxPlot=19"
plt.legend(['Calculated Data','Fitted Data'],loc='lower center',fontsize=16)
plt.text(2.5e-3,1.e-8, \
         ('Fitted $\Delta E_{ion}$ are proportional to $1/rho_{Init}^{%5.3f}$' % \
	  abs(fitB[indxPlot])),color='m',fontsize=16)
plt.grid(True)

indxPlot=29
VionCrrnt = V0*VionRel[indxPlot]
powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
fig800=plt.figure (800)
plt.loglog(rhoInit[0:nImpctPrmtr-1,indxPlot], \
           deltaEnrgIon_c[0:nImpctPrmtr-1,indxPlot],'-xr', \
           rhoInitFit[0:maxIndx[indxPlot],indxPlot], \
	   deltaEnrgIon_c_Fit[0:maxIndx[indxPlot],indxPlot],'ob',linewidth=2)
plt.xlabel('Track Initial Impact Parameter $rho_{Init}$', \
           color='m',fontsize=16)
plt.ylabel('$\Delta E_{ion}$, $eV$', color='m',fontsize=16)
titleHeader = 'Transferred Energy $\Delta E_{ion}$ to Single Ion for '
titleHeader += ' $V_{ion}=%5.3f\cdot10^{%2d}$ $cm/s$'
plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=16)
plt.xlim([.95*rhoInit[0,0],1.05*rhoInit[nImpctPrmtr-1,0]])
plt.ylim(yLimit)  # for "indxPlot=29"
plt.legend(['Calculated Data','Fitted Data'],loc='lower center',fontsize=16)
plt.text(2.5e-3,1.e-8, \
         ('Fitted $\Delta E_{ion}$ are proportional to $1/rho_{Init}^{%5.3f}$' % \
	  abs(fitB[indxPlot])),color='m',fontsize=16)
plt.grid(True)

indxPlot=39
VionCrrnt = V0*VionRel[indxPlot]
powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
fig900=plt.figure (900)
plt.loglog(rhoInit[0:nImpctPrmtr-1,indxPlot], \
           deltaEnrgIon_c[0:nImpctPrmtr-1,indxPlot],'-xr', \
           rhoInitFit[0:maxIndx[indxPlot],indxPlot], \
	   deltaEnrgIon_c_Fit[0:maxIndx[indxPlot],indxPlot],'ob',linewidth=2)
plt.xlabel('Track Initial Impact Parameter $rho_{Init}$', \
           color='m',fontsize=16)
plt.ylabel('$\Delta E_{ion}$, $eV$', color='m',fontsize=16)
titleHeader = 'Transferred Energy $\Delta E_{ion}$ to Single Ion for '
titleHeader += ' $V_{ion}=%5.3f\cdot10^{%2d}$ $cm/s$'
plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=16)
plt.xlim([.95*rhoInit[0,0],1.05*rhoInit[nImpctPrmtr-1,0]])
plt.ylim(yLimit)  # for "indxPlot=39"
plt.legend(['Calculated Data','Fitted Data'],loc='lower center',fontsize=16)
plt.text(2.5e-3,1.e-8, \
         ('Fitted $\Delta E_{ion}$ are proportional to $1/rho_{Init}^{%5.3f}$' % \
	  abs(fitB[indxPlot])),color='m',fontsize=16)
plt.grid(True)

indxPlot=49
VionCrrnt = V0*VionRel[indxPlot]
powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
fig1000=plt.figure (1000)
plt.loglog(rhoInit[0:nImpctPrmtr-1,indxPlot], \
           deltaEnrgIon_c[0:nImpctPrmtr-1,indxPlot],'-xr', \
           rhoInitFit[0:maxIndx[indxPlot],indxPlot], \
	   deltaEnrgIon_c_Fit[0:maxIndx[indxPlot],indxPlot],'ob',linewidth=2)
plt.xlabel('Track Initial Impact Parameter $rho_{Init}$', \
           color='m',fontsize=16)
plt.ylabel('$\Delta E_{ion}$, $eV$', color='m',fontsize=16)
titleHeader = 'Transferred Energy $\Delta E_{ion}$ to Single Ion for '
titleHeader += ' $V_{ion}=%5.3f\cdot10^{%2d}$ $cm/s$'
plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=16)
plt.xlim([.95*rhoInit[0,0],1.05*rhoInit[nImpctPrmtr-1,0]])
plt.ylim(yLimit)  # for "indxPlot=49"
plt.legend(['Calculated Data','Fitted Data'],loc='lower center',fontsize=16)
plt.text(2.5e-3,1.e-8, \
         ('Fitted $\Delta E_{ion}$ are proportional to $1/rho_{Init}^{%5.3f}$' % \
	  abs(fitB[indxPlot])),color='m',fontsize=16)
plt.grid(True)

fig1100=plt.figure (1100)
plt.semilogx(VionRel,fitB,'-xr',linewidth=2)
plt.xlabel('Relative Ion Velocity  $V_{ion}/V_0$',color='m',fontsize=16)
plt.ylabel('$B$', color='m',fontsize=16)
titleHeader = \
    'Dependence of $\Delta E$ on $rho^B$ ($V_{e0}=%5.3f\cdot10^{%2d}$cm/s)'
plt.title(titleHeader % (mantV0,powV0),color='m',fontsize=16)
plt.xlim(xLimit)
yLimit=[-2.3,-1.95]
plt.ylim(yLimit)
plt.plot([relVeTrnsv,relVeTrnsv],yLimit,'--m',linewidth=1)
plt.text(2.e-3,-2.3175,'$ \Delta V_{e\perp}/ V_{e0}$',color='m',fontsize=14)
plt.grid(True)

plt.show()


'''
fig10.savefig('picturesCMA/correctedRmax_fig10cma.jpg')    
fig20.savefig('picturesCMA/particles_distance_fig20cma.jpg')    
fig30.savefig('picturesCMA/parametersA-B_fig30cma.jpg')    
fig40.savefig('picturesCMA/initialImpactParameter_fig40cma.jpg')    
fig110.savefig('picturesCMA/mapA-B_fig110cma.jpg')    
fig209.savefig('picturesCMA/rDebye_rLikeDebye_rPass_fig209cma.jpg')    
fig3151.savefig('picturesCMA/impctPrmtr_3151cma.jpg')    
fig500.savefig('picturesCMA/deltaEtransf_indxPlot-0_500cma.jpg')    
fig600.savefig('picturesCMA/deltaEtransf_indxPlot-9_600cma.jpg')    
fig700.savefig('picturesCMA/deltaEtransf_indxPlot-19_700cma.jpg')    
fig800.savefig('picturesCMA/deltaEtransf_indxPlot-29_800cma.jpg')    
fig900.savefig('picturesCMA/deltaEtransf_indxPlot-39_900cma.jpg')    
fig1000.savefig('picturesCMA/deltaEtransf_indxPlot-49_1000cma.jpg')    
fig1100.savefig('picturesCMA/exponentB_on_ionVelocity_1100cma.jpg')    
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

'''
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
   
print 'bEdges[0]=%e, bEdges[nTotal]=%e (all sizes in mkm)' % \
      (1.e+4*bEdges[0],1.e+4*bEdges[nTotal])

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
rhoCrit=math.pow(q_elec**2/(m_elec*omega_L**2),1./3)    # cm
alpha=2*q_elec**2*omega_L/(m_elec*eVrmsLong**3)         # dimensionless
beta=2*q_elec**2/(m_elec*eVrmsLong**2)                  # cm
print 'Rshield(mkm)=%e, rhoCrit(mkm)=%f, alpha=%f' % \
      (1.e+4*Rshield[0],1.e+4*rhoCrit,alpha)

#
# Array A=log10(Upot/Ekin)=
#        =log10([q_e^2/rho]/[m_e((V_transverse^2+V_longitudinal^2)/2]): 
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

evTran=np.zeros((nA,nB))                                # cm/sec
rho_larm=np.zeros((nA,nB))                              # cm
kinEnergy=np.zeros((nA,nB))                             # erg
rhoCrrnt=np.zeros((nA,nB))                              # cm
#-----
# To check the solution of the  equation for V_transversal, using comparison 
# two formulas to calculate the impact parameter:
rhoCheck=np.zeros((nA,nB))                              # cm
#-----
halfLintr=np.zeros((nA,nB))                             # cm
timePath=np.zeros((nA,nB))                              # sec
numbLarmor=np.zeros((nA,nB))                            # dimensionless

# To plot area of the impact parameter b; 
# First index =0 for b < Rshield[0], 
#             =1 for Rshield[0] < b < Rshield[1], 
#             =2 for Rshield[1] < b < Rshield[2]: 

mapAimpctParmtr=np.zeros((3,nA*nB))
mapBimpctParmtr=np.zeros((3,nA*nB))
rhoCrrntMax=0.
rhoCrrntMax_A=0
rhoCrrntMax_B=0
trackNumb=-1                                    # Tracks are enumerated from 0!
track_1=-1                                      # Tracks are enumerated from 0!
track_c=-1                                      # Tracks are enumerated from 0!

for iA in range(nA):
   crrntA[iA]=minA+stepA*iA                     # log10(Upot/Ekin)
   for iB in range(nB):
      crrntB[iB]=minB+stepB*iB                  # log10(Rlarm/rho)
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

trackNumb_1=-1                                                     # Tracks will be enumerated from 0!

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

maxAbsDpxApprch_c=0.
maxAbsDpxApprch_m=0.

totAbsDpxApprch_c=0.
totAbsDpxApprch_m=0.

indxAmaxAbsDpxApprch_c=0
indxAmaxAbsDpxApprch_m=0

indxBmaxAbsDpxApprch_c=0
indxBmaxAbsDpxApprch_m=0

trackNumbMaxAbsDpxApprch_c=0
trackNumbMaxAbsDpxApprch_m=0

#
# Working arrays:
#
larmorNumber=np.zeros(nA*nB)
larmorRadius=np.zeros(nA*nB)

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

maxYcoorElec_c=0.
maxYcoorIon_c=0.
maxYcoorElec_m=0.
maxYcoorIon_m=0.

cpuTimeTotal=0

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
# First track: to compare distances between particles for approach_c,m (figure 325):
#

if runFlagApproach_c == 1:
#   b_cArray=np.asarray(b_c)
   b_cLenFirstTrack=int(timePoints_c[0])
   print 'First track length for b_c: %d' % b_cLenFirstTrack

if runFlagApproach_m == 1:
#   b_mArray=np.asarray(b_m)
   b_mLenFirstTrack=int(timePoints_m[0])
   print 'First track length for b_m: %d' % b_mLenFirstTrack

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

plt.show()

sys.exit()


'''
