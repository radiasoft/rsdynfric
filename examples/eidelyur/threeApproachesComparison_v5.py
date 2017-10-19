# from __future__ import division

#-------------------------------------
#
#        Started at 09/12/2017 (YuE)
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


eVtoErg=1.602e-12                                                  # energy from eV to erg (from CI to CGS)
# Indices:
(Ix, Ipx, Iy, Ipy, Iz, Ipz) = range(6)

#
# Initial parameters:
#
Z_ion = qe*2.997e+9                                                # charge of ion (proton), CGSE units of the charge
M_ion = mp*1.e+3                                                   # mass of ion (proton), g
q_elec = qe*2.997e+9                                               # charge of electron, CGSE units of the charge (without sign!)
m_elec = me*1.e+3                                                  # mass of electron, g

B_mag=1.e+3                                                        # Gs
eTempTran = 0.5                                                    # transversal temperature of electrons, eV
eTempLong = 2.e-4                                                  # longitudinal temperature of electrons, eV

stepsNumberOnGyro = 25                                             # number of the steps on each Larmour period

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
# Convertion from electron's "coordinates" to guiding-center coordinates:
# For electron z_e=(x_e,px_e,y_e,py_e,z_e,pz_e) --> zgc_e=(phi,p_phi,y_gc,py_gc,z_e,pz_e);
# z_e and zgc_e are 6-vectors 
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
# Convertion from guiding-center coordinates to electron's "coordinates":
# For each electron zgc_e=(phi,p_phi,y_gc,py_gc,z_e,pz_e) --> z_e=(x_e,px_e,y_e,py_e,z_e,pz_e);
# zgc_c and z_e are 6-vectors
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
   return dpIon,dpElec,action                                      

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
          (vectrIon[4]-vectrElec_gc[4])**2+rhoLarm_gc**2)**(3/2)                # cm^3
   action=vectrElec_gc[1]+dpFactor_gc*numer*rhoLarm_gc/(omega_L*denom)          # g*cm^2/sec
   C1=np.sqrt((vectrIon[0]-x_gc)**2+ \
              (vectrIon[2]-vectrElec_gc[2])**2+ \
              (vectrIon[4]-vectrElec_gc[4])**2+2.*action/mOmegaLarm)            # cm^2
   C2=2.*((vectrIon[0]-x_gc)*vectrIon[1]/M_ion+(vectrIon[2]-vectrElec_gc[2])*vectrIon[3]/M_ion+ \
          (vectrIon[4]-vectrElec_gc[4])*(vectrIon[5]/M_ion-vectrElec_gc[5]/m_elec))               # cm^2/sec
   C3=(vectrIon[1]/M_ion)**2+(vectrIon[3]/M_ion)**2+ \
      (vectrIon[5]/M_ion-vectrElec_gc[5]/m_elec)**2                # cm^2/sec^2
   b=np.sqrt(C1+C2*deltaT+C3*deltaT**2)                            # cm
   D1=(2.*C3*deltaT+C2)/b-C2/np.sqrt(C1)                           # cm/sec
   D2=(C2*deltaT+2.*C1)/b-np.sqrt(C1)                              # cm
   q=4.*C1*C3-C2**2                                                # cm^4/sec^2
   dpIon[0]=-2.*dpFactor_gc/q*((vectrIon[0]-x_gc)*D1-vectrIon[1]/M_ion*D2)                # g*cm/sec
   dpIon[1]=-2.*dpFactor_gc/q*((vectrIon[2]-vectrElec_gc[2])*D1-vectrIon[3]/M_ion*D2)     # g*cm/sec
   dpIon[2]=-2.*dpFactor_gc/q*((vectrIon[4]-vectrElec_gc[4])*D1- \
                                (vectrIon[5]/M_ion-vectrElec_gc[5]/m_elec)*D2)            # g*cm/sec
   dpElec[1]=-dpIon[1]                                             # g*cm/sec
   dpElec[2]=-dpIon[2]                                             # g*cm/sec
   dy_gc=dpIon[0]/mOmegaLarm                                       # cm
#    print 'dpIon[0]=%e, dpIon[1]=%e, dpIon[2]=%e' % (dpIon[0],dpIon[1],dpIon[2])
   return dpIon,dpElec,action,dy_gc                                      

z_elecCrrnt=np.zeros(6)                                            # 6-vector for electron (for Approach_1)
z_ionCrrnt=np.zeros(6)                                             # 6-vector for ion (for Approach_1)
z_elecCrrnt_2=np.zeros(6)                                          # 6-vector for electron (for Approach_2)
z_ionCrrnt_2=np.zeros(6)                                           # 6-vector for ion (for Approach_2)
z_elecCrrnt_gc=np.zeros(6)                     # 6-vector for electron in #guiding center" system (for Approaches_2,3)
z_elecCrrnt_3=np.zeros(6)                                          # 6-vector for electron (for Approach_3)
z_ionCrrnt_3=np.zeros(6)                                           # 6-vector for ion (for Approach_3)

#===============================================================================
#
# Case: the same longidudinal (eVrmsLong).
# The selection of main parameters is sescribed in the
# document "notesToChoiceCodeParameters.docx"
#

rhoMin=1.e-4                                                       # R_crit=1 mkm; cm
rShield=100.e-4                                                    # max value of the R_shield=50 mkm; cm
rhoMax=rShield                                                     # cm

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
# All data will be placed in nTotal bins (sorted by the range of impact parameters):
#
bTot=np.zeros(nTotal)                                              # cm
larmR_bTot=np.zeros(nTotal)                                        # ratio R_larmor/b; dimensionless
uPot_enrgKinTot=np.zeros(nTotal)                                   # ratio eVca/Ekin; dimensionless
population=np.zeros(nTotal)                                        # number of data in each bin
#
# Transfered momenta for different approaches:
#      deltaP=q_e^2*timeStep*abs(dpApprch_NTot|).
# First index numerates the x,y and z-component od deltaP:
#
dpApprch_1Tot=np.zeros((3,nTotal))                                 # 1/cm^2
dpApprch_2Tot=np.zeros((3,nTotal))                                 # 1/cm^2
dpApprch_3Tot=np.zeros((3,nTotal))                                 # 1/cm^2

sumPoints=0
sumPoints_2=0
sumPoints_3=0

#
# Electrons are magnetized for impact parameter >> rhoCrit:
#
rhoCrit=math.pow(q_elec**2/(m_elec*omega_L**2),1./3)               # cm
alpha=2*q_elec**2*omega_L/(m_elec*eVrmsLong**3)                    # dimensionless
print 'rhoCrit(mkm)=%f, alpha=%f' % (1.e+4*rhoCrit, alpha)

#
# Array A=log10(Upot/Ekin): 
#
nA=5
crrntA=np.zeros(nA)
minA=-5.
maxA=0.
stepA=(maxA-minA)/(nA-1)

#
# Array B=log10(R_larm/b): 
#
nB=5
crrntB=np.zeros(nB)
minB=-3.
maxB=-.5
stepB=(maxB-minB)/(nB-1)

#
# Working arrays:
#
larmorNumber=np.zeros(nA*nB)
larmorRadius=np.zeros(nA*nB)

pointTrack=np.zeros(nA*nB)  
timePoints=np.zeros(nA*nB)
cpuTime=np.zeros(nA*nB)

pointTrack_2=np.zeros(nA*nB)  
timePoints_2=np.zeros(nA*nB)
cpuTime_2=np.zeros(nA*nB)

pointTrack_3=np.zeros(nA*nB)  
timePoints_3=np.zeros(nA*nB)
cpuTime_3=np.zeros(nA*nB)

evTran=np.zeros((nA,nB))                                         # cm/sec
rho_larm=np.zeros((nA,nB))                                       # cm
kinEnergy=np.zeros((nA,nB))                                      # erg
rhoCrrnt=np.zeros((nA,nB))                                       # cm
halfLintr=np.zeros((nA,nB))                                      # cm
timePath=np.zeros((nA,nB))                                       # sec
numbLarmor=np.zeros((nA,nB))                                     # dimensionless

# 
# For approach_1 (without averaging over larmor rotation):
#
matr_elec=solenoid_eMatrix(B_mag,timeStep)                         # matrix for electron for timeStep in magnetic field
matr_ion=drift_Matrix(M_ion,timeStep)                              # matrix for ion (with mass M_ion) for timeStep
trackNumb=-1                                                       # Tracks will be enumerated from 0!
dpFactor=q_elec**2*timeStep                                        # g*cm^3/sec
minLarmR_b_1=1.e8                                                  # dimensionless     
maxLarmR_b_1=-1.e8                                                 # dimensionless
minUpot_enrgKin_1=1.e8                                             # dimensionless      
maxUpot_enrgKin_1=-1.e8                                            # dimensionless      

# 
# For approach_2 (with averaging over nLarmorAvrgng larmor rotation) +
# "Magnus ex[ansion" method to calculate the transferred momenmta):
#
nLarmorAvrgng=1                                                    # number of averaged Larmor rotations 
timeStep_2=nLarmorAvrgng*stepsNumberOnGyro*timeStep                # time step for approach_2
matr_elec_2=guidingCenter_Matrix(.5*timeStep_2)                    # matrix for electron for timeStep_2/2 
matr_ion_2=drift_Matrix(M_ion,.5*timeStep_2)                       # matrix for ion (with mass M_ion) for timeStep_2/2
trackNumb_2=-1                                                     # Tracks will be enumerated from 0!
dpFactor_2=q_elec**2*timeStep_2                                    # g*cm^3/sec
minLarmR_b_2=1.e8                                                  # dimensionless     
maxLarmR_b_2=-1.e8                                                 # dimensionless
minUpot_enrgKin_2=1.e8                                             # dimensionless      
maxUpot_enrgKin_2=-1.e8                                            # dimensionless      

# 
# For approach_3 (with averaging over nLarmorAvrgng larmor rotation +
# "Guiding center" method to calculate the transferred momenmta):
#
nLarmorAvrgng=1                                                    # number of averaged Larmor rotations 
timeStep_3=nLarmorAvrgng*stepsNumberOnGyro*timeStep                # time step for approach_3
matr_elec_3=guidingCenter_Matrix(.5*timeStep_3)                    # matrix for electron for timeStep_3/2 
matr_ion_3=drift_Matrix(M_ion,.5*timeStep_3)                       # matrix for ion (with mass M_ion) for timeStep_3/2
trackNumb_3=-1                                                     # Tracks will be enumerated from 0!
dpFactor_3=q_elec**2*timeStep_3                                    # g*cm^3/sec
minLarmR_b_3=1.e8                                                  # dimensionless     
maxLarmR_b_3=-1.e8                                                 # dimensionless
minUpot_enrgKin_3=1.e8                                             # dimensionless      
maxUpot_enrgKin_3=-1.e8                                            # dimensionless      

print 'timeStep=%e, timeStep_2=%e, timeStep_3=%e' % (timeStep,timeStep_2,timeStep_3)

# 
# To draw track with maximal transfer of px:
#
maxAbsDpxApprch=0.
maxAbsDpxApprch_2=0.
maxAbsDpxApprch_3=0.

totAbsDpxApprch=0.
totAbsDpxApprch_2=0.
totAbsDpxApprch_3=0.

indxAmaxAbsDpxApprch=0
indxAmaxAbsDpxApprch_2=0
indxAmaxAbsDpxApprch_3=0

indxBmaxAbsDpxApprch=0
indxBmaxAbsDpxApprch_2=0
indxBmaxAbsDpxApprch_3=0

trackNumbMaxAbsDpxApprch=0
trackNumbMaxAbsDpxApprch_2=0
trackNumbMaxAbsDpxApprch_3=0

runFlagApproach_1=1                                                # run approach1 if = 1, otherwise don't run
runFlagApproach_2=1                                                # run approach2 if = 1, otherwise don't run
runFlagApproach_3=1                                                # run approach3 if = 1, otherwise don't run

plotFlagWorkingArrays=1                                            # Plot working arrays if = 1, orherwise don't plot
plotFlagTracks=1                                                   # Plot everything about tracks if = 1, orherwise don't plot
plotFlagDpTransf=0                                                 # Plot everything about dp transf. if = 1, orherwise don't plot

maxYcoorElec=0.
maxYcoorIon=0.
cpuTimeTotal=0
for iA in range(nA):
   crrntA[iA]=minA+stepA*iA                                        # log10(Upot/Ekin)
   for iB in range(nB):
      crrntB[iB]=minB+stepB*iB                                     # log10(Rlarm/b)
      q=-alpha*math.pow(10.,crrntB[iB]-crrntA[iA])
#
# Equation: y^3+y+q=0:
#
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
      evTran[iA,iB]=eVrmsLong*root                                 # cm/sec
      rho_larm[iA,iB]=evTran[iA,iB]/omega_L                        # cm
      kinEnergy[iA,iB]=m_elec*(evTran[iA,iB]**2+eVrmsLong**2)/2.   # erg
      rhoCrrnt[iA,iB]=rhoCrit+rho_larm[iA,iB]                      # cm
      halfLintr[iA,iB]=rho_larm[iA,iB]/math.pow(10.,crrntB[iB])    # cm
      timePath[iA,iB]=halfLintr[iA,iB]/eVrmsLong                   # sec
      numbLarmor[iA,iB]=int(timePath[iA,iB]/T_larm)                # dimensionless
      print 'iA=%d, iB=%d, trackNumber=%d, numbLarmor[iA,iB]=%d: Lintr(mkm)=%e, rhoLarm(mkm)=%e' % \
	    (iA,iB,trackNumb,numbLarmor[iA,iB],1.e4*halfLintr[iA,iB],1.e4*rho_larm[iA,iB])   
      if numbLarmor[iA,iB] >= 40:
         trackNumb += 1                                            # Tracks are enumerated from 0!  
         trackNumb_2 += 1                                          # Tracks are enumerated from 0!  
         trackNumb_3 += 1                                          # Tracks are enumerated from 0!  
         larmorNumber[trackNumb]=numbLarmor[iA,iB]
         larmorRadius[trackNumb]=rho_larm[iA,iB]                          # cm
         timePoints[trackNumb]=int(numbLarmor[iA,iB]*stepsNumberOnGyro)   # for approach_1; dimensionless
         timePoints_2[trackNumb_2]=int(numbLarmor[iA,iB]/nLarmorAvrgng)   # for approach_2; dimensionless
         timePoints_3[trackNumb_3]=int(numbLarmor[iA,iB]/nLarmorAvrgng)   # for approach_3; dimensionless
	 print 'timePoints=%d, timePoints_2=%d, timePoints_3=%d, numbLarmor[iA,iB]=%d' % \
	       (timePoints[trackNumb],timePoints_2[trackNumb_2],timePoints_2[trackNumb_3],larmorNumber[trackNumb])   
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Approach_1: dragging without averaging over larmor rotation
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# To draw only first trajectories (for checking only):
#
#------- Start of calculations for approach_1 --------------
#
         if runFlagApproach_1 == 1:
            timeStart=os.times()
            if trackNumb == 0:
	       rhoFirstTurn=rhoCrrnt[iA,iB]
	       rhoLarmorFirstTurn=rho_larm[iA,iB]
# 6-vectors for ion and electron and distance 'b' between them for the first trajectory and 
# trajectory with maximal transferred dpx (TMTdpx)
# (for checking only; indices 0-5 for electron, indices 6-11 for ion and index=12 for 'b'):
               prtclCoor=np.zeros((13,timePoints[trackNumb]))                   # first trajectory                 
               prtclCoorCrrnt=np.zeros((13,timePoints[trackNumb]))              # current trajectory                
               prtclCoorMaxAbsDpx=np.zeros((13,timePoints[trackNumb]))          # trajectory TMTdpx                    
# Current distance from origin of the coordinate system to electron along the trajectory; cm
            bCrrnt=np.zeros(timePoints[trackNumb])                 # cm
# Current log10 of two important ratios; dimensionless:
            larmR_bCrrnt=np.zeros(timePoints[trackNumb])           # ratio R_larmor/b; dimensionless
            uPot_enrgKinCrrnt=np.zeros(timePoints[trackNumb])      # ratio potential_energy/kinetic_energy; dimensionless
            dpApprch_1Crrnt=np.zeros((3,timePoints[trackNumb]))    # 1/cm^2; deltaPapprch_1=dpFactor*dpApprch_1Crrnt; dpFactor=q_e^2*timeStep
            for m in range(6): 
               z_ionCrrnt[m]=0.                                    # Current initial zero-vector for ion
               z_elecCrrnt[m]=0.                                   # Zeroing out of vector for electron
# Current initial vector for electron:
            z_elecCrrnt[Ix]=rhoCrrnt[iA,iB]                        # x, cm
            z_elecCrrnt[Iz]=-halfLintr[iA,iB]                      # z, cm
            z_elecCrrnt[Ipy]=m_elec*evTran[iA,iB]                  # py, g*cm/sec
            z_elecCrrnt[Ipz]=m_elec*eVrmsLong                      # pz, g*cm/sec
#-----------------------------------------------
# Main action - dragging of the current trajectories (for given i and j)
#
            for k in range(int(timePoints[trackNumb])):
 	       z_elecCrrnt=np.dot(matr_elec,z_elecCrrnt)           # electron's dragging
 	       z_ionCrrnt=np.dot(matr_ion,z_ionCrrnt)              # ion's dragging
# Current distance between ion and electron; cm:   
 	       bCrrnt[k]=np.sqrt((z_ionCrrnt[0]-z_elecCrrnt[0])**2+ \
	                         (z_ionCrrnt[2]-z_elecCrrnt[2])**2+ \
			         (z_ionCrrnt[4]-z_elecCrrnt[4])**2)
# Current log10 of two important ratios:  
	       larmR_bCrrnt[k]=math.log10(rho_larm[iA,iB]/bCrrnt[k])      # dimensionless 
	       if maxLarmR_b_1 < larmR_bCrrnt[k]:
	          maxLarmR_b_1=larmR_bCrrnt[k]
	       if minLarmR_b_1 > larmR_bCrrnt[k]:
	          minLarmR_b_1=larmR_bCrrnt[k]
	       uPot_enrgKinCrrnt[k]=math.log10((q_elec**2/bCrrnt[k])/kinEnergy[iA,iB])      # dimensionless 
	       if maxUpot_enrgKin_1 < uPot_enrgKinCrrnt[k]:
	          maxUpot_enrgKin_1=uPot_enrgKinCrrnt[k]
	       if minUpot_enrgKin_1 > uPot_enrgKinCrrnt[k]:
	          minUpot_enrgKin_1=uPot_enrgKinCrrnt[k]
# Current values to calculate deltaPapprch_1=dpFactor*dpApprch_1Crrnt; dpFactor=q_e^2*timeStep:  
               for ic in range(3):
	          dpApprch_1Crrnt[ic,k]= abs(z_ionCrrnt[2*ic]-z_elecCrrnt[2*ic])/ bCrrnt[k]**3      # 1/cm^2;  
# To draw future TMTdpx trajectory (for checking only):
               for ic in range(6):
                  prtclCoorCrrnt[ic,pointTrack[trackNumb]]=z_elecCrrnt[ic]           # 6-vector for electron
                  prtclCoorCrrnt[6+ic,pointTrack[trackNumb]]=z_ionCrrnt[ic]          # 6-vector for ion 
	       prtclCoorCrrnt[12,pointTrack[trackNumb]]=bCrrnt[k]                    # cm
# To draw only first trajectory (for checking only):
               if trackNumb == 0:
                  for ic in range(6):
                     prtclCoor[ic,pointTrack[trackNumb]]=z_elecCrrnt[ic]             # 6-vector for electron
                     prtclCoor[6+ic,pointTrack[trackNumb]]=z_ionCrrnt[ic]            # 6-vector for ion 
	          prtclCoor[12,pointTrack[trackNumb]]=bCrrnt[k]                      # cm
	          if maxYcoorElec < abs(z_elecCrrnt[2]):
	             maxYcoorElec=abs(z_elecCrrnt[2])
	          if maxYcoorIon < abs(z_ionCrrnt[2]):
	             maxYcoorIon=abs(z_ionCrrnt[2])
#	          crrntPoint=pointTrack[trackNumb]
#	          if crrntPoint < 100 or crrntPoint > 45400:
#                     print 'Point %d: x=%e, y=%e, z=%e' % \
#	                   (crrntPoint,1.e+7*prtclCoor[6,crrntPoint],1.e+7*prtclCoor[8,crrntPoint], \
#		                       1.e+7*prtclCoor[10,crrntPoint])
#
# Taking into account transfer of momentum for both particles:
#
               for ic in range(3):
                  z_ionCrrnt[2*ic+1] += -dpFactor*dpApprch_1Crrnt[ic,k]
                  z_elecCrrnt[2*ic+1] += dpFactor*dpApprch_1Crrnt[ic,k]
               pointTrack[trackNumb] += 1
#
# Total transferred dpx  for current track:
#
	       totAbsDpxApprch += abs(dpFactor*dpApprch_1Crrnt[0,k])
# 
# End of dragging of the current trajectory	  
#-----------------------------------------------
#
# Accumulate transferred momenta and other data: 
#
            if trackNumb == 0: 
# First definition of the total distance from origin of the coordinate system to electron along the trajectory; cm:
 	       b=bCrrnt                                            # cm
# First definition of the log10 of two important ratios; dimensionless:  
	       larmR_b=larmR_bCrrnt                                # dimensionless 
	       uPot_enrgKin=uPot_enrgKinCrrnt                      # dimensionless 
# First definition of the values to calculate deltaPapprch_1:
               dpxApprch_1=dpFactor*dpApprch_1Crrnt[0,:]           # g*cm/sec  
               dpyApprch_1=dpFactor*dpApprch_1Crrnt[1,:]           # g*cm/sec  
               dpzApprch_1=dpFactor*dpApprch_1Crrnt[2,:]           # g*cm/sec  
            else:  
# Total distance from origin of the coordinate system to electron along the trajectory:
 	       b=np.concatenate((b,bCrrnt),axis=0)                 # cm
# Total log10 of two important ratios; dimensionless:  
	       larmR_b=np.concatenate((larmR_b,larmR_bCrrnt),axis=0)                  
	       uPot_enrgKin=np.concatenate((uPot_enrgKin,uPot_enrgKinCrrnt),axis=0)        
# Total values to calculate deltaPapprch_1:
               dpxApprch_1=np.concatenate((dpxApprch_1,dpFactor*dpApprch_1Crrnt[0,:]),axis=0)       # g*cm/sec  
               dpyApprch_1=np.concatenate((dpyApprch_1,dpFactor*dpApprch_1Crrnt[1,:]),axis=0)       # g*cm/sec  
               dpzApprch_1=np.concatenate((dpzApprch_1,dpFactor*dpApprch_1Crrnt[2,:]),axis=0)       # g*cm/sec  
# 
# To draw TMTdpx trajectory (for checking only):
#
            if maxAbsDpxApprch < totAbsDpxApprch:
	       maxAbsDpxApprch=totAbsDpxApprch
	       indxAmaxAbsDpxApprch=iA
	       indxBmaxAbsDpxApprch=iB
	       trackNumbMaxAbsDpxApprch=trackNumb
	       prtclCoorMaxAbsDpx=prtclCoorCrrnt.copy()   
	       rhoMaxAbsDpxTurn=rhoCrrnt[iA,iB]
	       rhoLarmorMaxAbsDpxTurn=rho_larm[iA,iB]
	       print 'iA=%d, iB=%d: track %d, points %d' % \
	             (indxAmaxAbsDpxApprch,indxBmaxAbsDpxApprch,trackNumbMaxAbsDpxApprch,pointTrack[trackNumbMaxAbsDpxApprch]) 
	       print 'timePoints.shape: ', (prtclCoorMaxAbsDpx.shape[0],prtclCoorMaxAbsDpx.shape[1])                                            
# End of all calculations for approach_1
            lastTrackNumber=trackNumb+1                            # quantity of tracks = trackNumb + 1!     
            sumPoints += pointTrack[trackNumb]
            timeEnd=os.times()
	    cpuTime[trackNumb]=1.e+6*(float(timeEnd[0])-float(timeStart[0]))         # CPU time , mks
            cpuTimeTotal += cpuTime[trackNumb] 
#
#------- End of approach_1 --------------
#
# for i in range(lastTrackNumber):
#    print 'Track %d: larmor turns=%d, cpuTime(mks)=%e, time per turn(mks)=%6.1f' % \
#          (i,larmorNumber[i],cpuTime[i],cpuTime[i]/larmorNumber[i])

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Approach_2: dragging with averaging over nLarmorAvrgng larmor rotation +
#             "Guiding center" method to calculate the transferred momenta
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#
#------- Start of calculations for approach_2 --------------
#
         if runFlagApproach_2 == 1:
            timeStart=os.times()
            if trackNumb_2 == 0:
	       rhoFirstTurn=rhoCrrnt[iA,iB]
	       rhoLarmorFirstTurn=rho_larm[iA,iB]
# 6-vectors for ion and electron and distance 'b' between them for the first trajectories 
# (for checking only; indices 0-5 for electron, indices 6-11 for ion, index=12 for 'b' and index=13 for action):
               prtclCoor_2=np.zeros((14,timePoints_2[trackNumb_2]))                    
               prtclCoorCrrnt_2=np.zeros((14,timePoints_2[trackNumb_2]))             # current trajectory                
               prtclCoorMaxAbsDp_2=np.zeros((14,timePoints_2[trackNumb_2]))          # trajectory TMTdpx                    
# Current distance from origin of the coordinate system to electron along the trajectory; cm
            bCrrnt_2=np.zeros(timePoints_2[trackNumb_2])           # cm
# Current log10 of two important ratios; dimensionless:
            larmR_bCrrnt_2=np.zeros(timePoints_2[trackNumb_2])     # ratio R_larmor/b; dimensionless
            uPot_enrgKinCrrnt_2=np.zeros(timePoints_2[trackNumb_2])        # ratio potential_energy/kinetic_energy; dimensionless
            dpApprch_2Crrnt=np.zeros((3,timePoints_2[trackNumb_2]))        # g*cm/sec
            for m in range(6): 
               z_ionCrrnt_2[m]=0.                                  # Current initial zero-vector for ion
               z_elecCrrnt_2[m]=0.                                 # Zeroing out of vector for electron
# Current initial vector for electron:
            z_elecCrrnt_2[Ix]=rhoCrrnt[iA,iB]+larmorRadius[trackNumb_2]   # x, cm
            z_elecCrrnt_2[Iz]=-halfLintr[iA,iB]                           # z, cm
            z_elecCrrnt_2[Ipy]=m_elec*evTran[iA,iB]                       # py, g*cm/sec
            z_elecCrrnt_2[Ipz]=m_elec*eVrmsLong                    # pz, g*cm/sec
	    z_elecCrrnt_gc=toGuidingCenter(z_elecCrrnt_2)          # transfer to system of guiding center
#	    if iA == 0 and iB == 0:
#	       print 'z_elecCrrnt_2: ', z_elecCrrnt_2 
#	       print 'z_elecCrrnt_gc: ', z_elecCrrnt_gc 
#-----------------------------------------------
# Main action - dragging of the current trajectories (for given i and j)
#
            for k in range(int(timePoints_2[trackNumb_2])):
#	       if k < 100: 
#	          print 'k=%d: x=%e, y=%e, z=%e' % (k,1.e+4*z_elecCrrnt_gc[0],1.e+4*z_elecCrrnt_gc[2],1.e+4*z_elecCrrnt_gc[4]) 
#	       if k > timePoints_2[trackNumb_2]-100: 
#	          print 'k=%d: x=%e, y=%e, z=%e' % (k,1.e+4*z_elecCrrnt_gc[0],1.e+4*z_elecCrrnt_gc[2],1.e+4*z_elecCrrnt_gc[4]) 
# 	       z_elecCrrnt_gc=matr_elec_2.dot(z_elecCrrnt_gc)      # electron's dragging for first nhalf timeStep
# 	       z_ionCrrnt_2=matr_ion_2.dot(z_ionCrrnt_2)           # ion's dragging for first half timeStep
 	       z_elecCrrnt_gc=np.dot(matr_elec_2,z_elecCrrnt_gc)   # electron's dragging for first nhalf timeStep
 	       z_ionCrrnt_2=np.dot(matr_ion_2,z_ionCrrnt_2)        # ion's dragging for first half timeStep
# dragging both paticles through interaction point:
	       dpIon,dpElec,action=guidingCenterCollision(z_elecCrrnt_gc,z_ionCrrnt_2,timeStep_2) 
#	       if trackNumb_2 == 0:
#	          print 'point %d: dpgcElec=%e, dpzElec=%e' % (pointTrack_2[0],dpElec[1],dpElec[2])
	       for ic in range(3):
	          z_ionCrrnt_2[2*ic+1] += dpIon[ic]   
	          z_elecCrrnt_gc[2*ic+1] += dpElec[ic]   
# Current values to calculate deltaPapprch_2:  
	          dpApprch_2Crrnt[ic,k]=dpIon[ic]                  # g*cm/sec 
# 	       z_elecCrrnt_gc=matr_elec_2.dot(z_elecCrrnt_gc)      # electron's dragging for second half timeStep
# 	       z_ionCrrnt_2=matr_ion_2.dot(z_ionCrrnt_2)           # ion's dragging for second half timeStep
 	       z_elecCrrnt_gc=np.dot(matr_elec_2,z_elecCrrnt_gc)   # electron's dragging for first nhalf timeStep
 	       z_ionCrrnt_2=np.dot(matr_ion_2,z_ionCrrnt_2)        # ion's dragging for first half timeStep
#	       crrntPoint=pointTrack_2[trackNumb_2]
#	       if iA == 0 and iB == 0 and crrntPoint < 10:
#                     print 'k, z_ionCrrnt_2', (k,z_ionCrrnt_2)
	       z_elecCrrnt_2=fromGuidingCenter(z_elecCrrnt_gc)     # transfer from system of guiding center 
# Current distance between ion and electron; cm:
 	       bCrrnt_2[k]=np.sqrt((z_ionCrrnt_2[0]-z_elecCrrnt_2[0])**2+ \
	                           (z_ionCrrnt_2[2]-z_elecCrrnt_2[2])**2+ \
			           (z_ionCrrnt_2[4]-z_elecCrrnt_2[4])**2)
# Current log10 of two important ratios:  
	       larmR_bCrrnt_2[k]=math.log10(rho_larm[iA,iB]/bCrrnt_2[k])  # dimensionless 
	       if maxLarmR_b_2 < larmR_bCrrnt_2[k]:
	          maxLarmR_b_2=larmR_bCrrnt_2[k]
	       if minLarmR_b_2 > larmR_bCrrnt_2[k]:
	          minLarmR_b_2=larmR_bCrrnt_2[k]
	       uPot_enrgKinCrrnt_2[k]=math.log10((q_elec**2/bCrrnt_2[k])/kinEnergy[iA,iB])  # dimensionless 
	       if maxUpot_enrgKin_2 < uPot_enrgKinCrrnt_2[k]:
	          maxUpot_enrgKin_2=uPot_enrgKinCrrnt_2[k]
	       if minUpot_enrgKin_2 > uPot_enrgKinCrrnt_2[k]:
	          minUpot_enrgKin_2=uPot_enrgKinCrrnt_2[k]
# To draw future TMTdpx trajectory (for checking only):
               for ic in range(6):
                  prtclCoorCrrnt_2[ic,pointTrack_2[trackNumb_2]]=z_elecCrrnt_2[ic]   # 6-vector for electron
                  prtclCoorCrrnt_2[6+ic,pointTrack_2[trackNumb_2]]=z_ionCrrnt_2[ic]  # 6-vector for ion 
	       prtclCoorCrrnt_2[12,pointTrack_2[trackNumb_2]]=bCrrnt_2[k]            # cm
	       prtclCoorCrrnt_2[13,pointTrack_2[trackNumb_2]]=action                 # g*cm^2/sec
# To draw only first trajectory (for checking only):
               if trackNumb_2 == 0:
                  for ic in range(6):
                     prtclCoor_2[ic,pointTrack_2[trackNumb_2]]=z_elecCrrnt_2[ic]     # 6-vector for electron
                     prtclCoor_2[6+ic,pointTrack_2[trackNumb_2]]=z_ionCrrnt_2[ic]    # 6-vector for ion 
	          prtclCoor_2[12,pointTrack_2[trackNumb_2]]=bCrrnt_2[k]              # cm
	          prtclCoor_2[13,pointTrack_2[trackNumb_2]]=action                   # g*cm^2/sec
#	          crrntPoint=pointTrack_2[trackNumb_2]
#	          if crrntPoint < 100 or crrntPoint > 1700:
#                     print 'Point %d: x=%e, y=%e, z=%e' % \
#	                   (crrntPoint,1.e+7*prtclCoor_2[6,crrntPoint],1.e+7*prtclCoor_2[8,crrntPoint], \
#		                       1.e+7*prtclCoor_2[10,crrntPoint])
	       
               pointTrack_2[trackNumb_2] += 1
#
# Total transferred dpx  for current track:
#
	       totAbsDpxApprch_2 += abs(dpFactor*dpApprch_2Crrnt[0,k])
# End of dragging of the current trajectory	  
#-----------------------------------------------
#
# Accumulate transferred momenta and other data: 
#
            if trackNumb_2 == 0: 
# First definition of the total distance from origin of the coordinate system to electron along the trajectory:
 	       b_2=bCrrnt_2                                        # cm 
# First definition of the log10 of two important ratios; dimensionless:  
	       larmR_b_2=larmR_bCrrnt_2                            
	       uPot_enrgKin_2=uPot_enrgKinCrrnt_2                 
# First definition of the values deltaPapprch_2:
               dpxApprch_2=dpApprch_2Crrnt[0,:]                    # g*cm/sec  
               dpyApprch_2=dpApprch_2Crrnt[1,:]                    # g*cm/sec  
               dpzApprch_2=dpApprch_2Crrnt[2,:]                    # g*cm/sec  
            else:  
# Total distance from origin of the coordinate system to electron along the trajectory:
 	       b_2=np.concatenate((b_2,bCrrnt_2),axis=0)           # cm
# Total log10 of two important ratios; dimensionless :  
	       larmR_b_2=np.concatenate((larmR_b_2,larmR_bCrrnt_2),axis=0)                  
	       uPot_enrgKin_2=np.concatenate((uPot_enrgKin_2,uPot_enrgKinCrrnt_2),axis=0)        
# Total values deltaPapprch_2:
               dpxApprch_2=np.concatenate((dpxApprch_2,dpApprch_2Crrnt[0,:]),axis=0)                   # g*cm/sec  
               dpyApprch_2=np.concatenate((dpyApprch_2,dpApprch_2Crrnt[1,:]),axis=0)                   # g*cm/sec
               dpzApprch_2=np.concatenate((dpzApprch_2,dpApprch_2Crrnt[2,:]),axis=0)                   # g*cm/sec 
#             print 'trackNumb_2:%d: shapes: b=%d, larmR_b_2=%d, uPot=%d, dpx=%d, dpy=%d, dpz=%d' % \
# 	            (trackNumb_2,b.shape[0],larmR_b_2.shape[0],uPot_enrgKin_2.shape[0], \
# 		    dpxApprch_2.shape[0],dpyApprch_2.shape[0],dpzApprch_2.shape[0])  
# 
# To draw TMTdpx trajectory (for checking only):
#
            if maxAbsDpxApprch_2 < totAbsDpxApprch_2:
	       maxAbsDpxApprch_2=totAbsDpxApprch_2
	       indxAmaxAbsDpxApprch_2=iA
	       indxBmaxAbsDpxApprch_2=iB
	       trackNumbMaxAbsDpxApprch_2=trackNumb
	       prtclCoorMaxAbsDpx_2=prtclCoorCrrnt_2.copy()   
	       rhoMaxAbsDpxTurn_2=rhoCrrnt[iA,iB]
	       rhoLarmorMaxAbsDpxTurn_2=rho_larm[iA,iB]
	       print 'iA=%d, iB=%d: track %d, points %d' % (indxAmaxAbsDpxApprch_2,\
	             indxBmaxAbsDpxApprch_2,trackNumbMaxAbsDpxApprch_2,pointTrack[trackNumbMaxAbsDpxApprch_2]) 
	       print 'timePoints.shape: ', (prtclCoorMaxAbsDpx_2.shape[0],prtclCoorMaxAbsDpx_2.shape[1])                                            
# End of all calculations for approach_2
            lastTrackNumber_2=trackNumb_2+1                        # quantity of tracks = trackNumb_2 + 1!     
            sumPoints_2 += pointTrack_2[trackNumb_2]
            timeEnd=os.times()
	    cpuTime_2[trackNumb_2]=1.e+6*(float(timeEnd[0])-float(timeStart[0]))        # CPU time , mks
            cpuTimeTotal += cpuTime_2[trackNumb_2] 
#
#------- End of approach_2 --------------
#

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Approach_3: dragging with averaging over nLarmorAvrgng larmor rotation +
#             "Magnus expansion" method to calculate the transferred momenta
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#
#------- Start of calculations for approach_3 --------------
#
         if runFlagApproach_3 == 1:
            timeStart=os.times()
            if trackNumb_3 == 0:
	       rhoFirstTurn=rhoCrrnt[iA,iB]
	       rhoLarmorFirstTurn=rho_larm[iA,iB]
# 6-vectors for ion and electron and distance 'b' between them for the first trajectories 
# (for checking only; indices 0-5 for electron, indices 6-11 for ion, index=12 for 'b', 
#                     index=13 for action and index=14 for dy_gc):
               prtclCoor_3=np.zeros((15,timePoints_3[trackNumb_3]))                    
               prtclCoorCrrnt_3=np.zeros((15,timePoints_3[trackNumb_3]))             # current trajectory                
               prtclCoorMaxAbsDp_3=np.zeros((15,timePoints_3[trackNumb_3]))          # trajectory TMTdpx                    
# Current distance from origin of the coordinate system to electron along the trajectory; cm
            bCrrnt_3=np.zeros(timePoints_3[trackNumb_3])           # cm
# Current log10 of two important ratios; dimensionless:
            larmR_bCrrnt_3=np.zeros(timePoints_3[trackNumb_3])     # ratio R_larmor/b; dimensionless
            uPot_enrgKinCrrnt_3=np.zeros(timePoints_3[trackNumb_3])       # ratio potential_energy/kinetic_energy; dimensionless
            dpApprch_3Crrnt=np.zeros((3,timePoints_3[trackNumb_3]))       # g*cm/sec
            for m in range(6): 
               z_ionCrrnt_3[m]=0.                                  # Current initial zero-vector for ion
               z_elecCrrnt_3[m]=0.                                 # Zeroing out of vector for electron
# Current initial vector for electron:
            z_elecCrrnt_3[Ix]=rhoCrrnt[iA,iB]+larmorRadius[trackNumb_3]   # x, cm
            z_elecCrrnt_3[Iz]=-halfLintr[iA,iB]                           # z, cm
            z_elecCrrnt_3[Ipy]=m_elec*evTran[iA,iB]                       # py, g*cm/sec
            z_elecCrrnt_3[Ipz]=m_elec*eVrmsLong                    # pz, g*cm/sec
	    z_elecCrrnt_gc=toGuidingCenter(z_elecCrrnt_3)          # transfer to system of guiding center
#	    if iA == 0 and iB == 0:
#	       print 'z_elecCrrnt_3: ', z_elecCrrnt_3 
#	       print 'z_elecCrrnt_gc: ', z_elecCrrnt_gc 
#-----------------------------------------------
# Main action - dragging of the current trajectories (for given i and j)
#
            for k in range(int(timePoints_3[trackNumb_3])):
#	       if k < 100: 
#	          print 'k=%d: x=%e, y=%e, z=%e' % (k,1.e+4*z_elecCrrnt_gc[0],1.e+4*z_elecCrrnt_gc[2],1.e+4*z_elecCrrnt_gc[4]) 
#	       if k > timePoints_2[trackNumb_3]-100: 
#	          print 'k=%d: x=%e, y=%e, z=%e' % (k,1.e+4*z_elecCrrnt_gc[0],1.e+4*z_elecCrrnt_gc[2],1.e+4*z_elecCrrnt_gc[4]) 
# 	       z_elecCrrnt_gc=matr_elec_3.dot(z_elecCrrnt_gc)      # electron's dragging for first half timeStep
# 	       z_ionCrrnt_3=matr_ion_3.dot(z_ionCrrnt_3)           # ion's dragging for first half timeStep
 	       z_elecCrrnt_gc=np.dot(matr_elec_3,z_elecCrrnt_gc)   # electron's dragging for first half timeStep
 	       z_ionCrrnt_3=np.dot(matr_ion_3,z_ionCrrnt_3)        # ion's dragging for first half timeStep
# dragging both paticles through interaction point:
	       dpIon,dpElec,action,dy_gc=MagnusExpansionCollision(z_elecCrrnt_gc,z_ionCrrnt_3,timeStep_3) 
#	       if trackNumb_3 == 0:
#	          print 'point %d: dpxIon=%e, dpyIon=%e, dpzIon=%e, dy_gc=%e' % \
#		        (pointTrack_3[0],dpIon[0],dpIon[1],dpElec[2],dy_gc)
	       for ic in range(3):
	          z_ionCrrnt_3[2*ic+1] += dpIon[ic]   
	          z_elecCrrnt_gc[2*ic+1] += dpElec[ic]   
# Current values to calculate deltaPapprch_3:  
	          dpApprch_3Crrnt[ic,k]=dpIon[ic]                  # g*cm/sec 
	       z_elecCrrnt_gc[2] += dy_gc                          # cm   
# 	       z_elecCrrnt_gc=matr_elec_2.dot(z_elecCrrnt_gc)      # electron's dragging for second half timeStep
# 	       z_ionCrrnt_2=matr_ion_2.dot(z_ionCrrnt_2)           # ion's dragging for second half timeStep
 	       z_elecCrrnt_gc=np.dot(matr_elec_3,z_elecCrrnt_gc)   # electron's dragging for second half timeStep
 	       z_ionCrrnt_2=np.dot(matr_ion_3,z_ionCrrnt_3)        # ion's dragging for second half timeStep
#	       crrntPoint=pointTrack_3[trackNumb_3]
#	       if iA == 0 and iB == 0 and crrntPoint < 10:
#                     print 'k, z_ionCrrnt_3', (k,z_ionCrrnt_3)
	       z_elecCrrnt_3=fromGuidingCenter(z_elecCrrnt_gc)     # transfer from system of guiding center 
# Current distance between ion and electron; cm:
 	       bCrrnt_3[k]=np.sqrt((z_ionCrrnt_3[0]-z_elecCrrnt_3[0])**2+ \
	                           (z_ionCrrnt_3[2]-z_elecCrrnt_3[2])**2+ \
			           (z_ionCrrnt_3[4]-z_elecCrrnt_3[4])**2)
# Current log10 of two important ratios:  
	       larmR_bCrrnt_3[k]=math.log10(rho_larm[iA,iB]/bCrrnt_3[k])  # dimensionless 
	       if maxLarmR_b_3 < larmR_bCrrnt_3[k]:
	          maxLarmR_b_3=larmR_bCrrnt_3[k]
	       if minLarmR_b_3 > larmR_bCrrnt_3[k]:
	          minLarmR_b_3=larmR_bCrrnt_3[k]
	       uPot_enrgKinCrrnt_3[k]=math.log10((q_elec**2/bCrrnt_3[k])/kinEnergy[iA,iB])  # dimensionless 
	       if maxUpot_enrgKin_3 < uPot_enrgKinCrrnt_3[k]:
	          maxUpot_enrgKin_3=uPot_enrgKinCrrnt_3[k]
	       if minUpot_enrgKin_3 > uPot_enrgKinCrrnt_3[k]:
	          minUpot_enrgKin_3=uPot_enrgKinCrrnt_3[k]
# To draw future TMTdpx trajectory (for checking only):
               for ic in range(6):
                  prtclCoorCrrnt_3[ic,pointTrack_3[trackNumb_3]]=z_elecCrrnt_3[ic]   # 6-vector for electron
                  prtclCoorCrrnt_3[6+ic,pointTrack_3[trackNumb_3]]=z_ionCrrnt_3[ic]  # 6-vector for ion 
	          prtclCoorCrrnt_3[12,pointTrack_3[trackNumb_3]]=bCrrnt_3[k]         # cm
		  prtclCoorCrrnt_3[13,pointTrack_3[trackNumb_3]]=action              # g*cm^2/sec
		  prtclCoorCrrnt_3[14,pointTrack_3[trackNumb_3]]=dy_gc               # cm
# To draw only first trajectory (for checking only):
               if trackNumb_3 == 0:
                  for ic in range(6):
                     prtclCoor_3[ic,pointTrack_3[trackNumb_3]]=z_elecCrrnt_3[ic]     # 6-vector for electron
                     prtclCoor_3[6+ic,pointTrack_3[trackNumb_3]]=z_ionCrrnt_3[ic]    # 6-vector for ion 
	          prtclCoor_3[12,pointTrack_3[trackNumb_3]]=bCrrnt_3[k]              # cm
		  prtclCoor_3[13,pointTrack_3[trackNumb_3]]=action                   # g*cm^2/sec
		  prtclCoor_3[14,pointTrack_3[trackNumb_3]]=dy_gc                    # cm
#	          crrntPoint=pointTrack_3[trackNumb_3]
#	          if crrntPoint < 100 or crrntPoint > 1700:
#                     print 'Point %d: x=%e, y=%e, z=%e' % \
#	                   (crrntPoint,1.e+7*prtclCoor_3[6,crrntPoint],1.e+7*prtclCoor_3[8,crrntPoint], \
#		                       1.e+7*prtclCoor_3[10,crrntPoint])
	       
               pointTrack_3[trackNumb_3] += 1
#
# Total transferred dpx  for current track:
#
	       totAbsDpxApprch_3 += abs(dpFactor*dpApprch_3Crrnt[0,k])
# End of dragging of the current trajectory	  
#-----------------------------------------------
#
# Accumulate transferred momenta and other data: 
#
            if trackNumb_3 == 0: 
# First definition of the total distance from origin of the coordinate system to electron along the trajectory:
 	       b_3=bCrrnt_3                                        # cm 
# First definition of the log10 of two important ratios; dimensionless:  
	       larmR_b_3=larmR_bCrrnt_3                             
	       uPot_enrgKin_3=uPot_enrgKinCrrnt_3                   
# First definition of the values deltaPapprch_3:
               dpxApprch_3=dpApprch_3Crrnt[0,:]                    # g*cm/sec  
               dpyApprch_3=dpApprch_3Crrnt[1,:]                    # g*cm/sec  
               dpzApprch_3=dpApprch_3Crrnt[2,:]                    # g*cm/sec  
            else:  
# Total distance from origin of the coordinate system to electron along the trajectory:
 	       b_3=np.concatenate((b_3,bCrrnt_3),axis=0)           # cm
# Total log10 of two important ratios; dimensionless :  
	       larmR_b_3=np.concatenate((larmR_b_3,larmR_bCrrnt_3),axis=0)                  
	       uPot_enrgKin_3=np.concatenate((uPot_enrgKin_3,uPot_enrgKinCrrnt_3),axis=0)        
# Total values deltaPapprch_3:
               dpxApprch_3=np.concatenate((dpxApprch_3,dpApprch_3Crrnt[0,:]),axis=0)                   # g*cm/sec  
               dpyApprch_3=np.concatenate((dpyApprch_3,dpApprch_3Crrnt[1,:]),axis=0)                   # g*cm/sec
               dpzApprch_3=np.concatenate((dpzApprch_3,dpApprch_3Crrnt[2,:]),axis=0)                   # g*cm/sec 
#             print 'trackNumb_3:%d: shapes: b=%d, larmR_b_3=%d, uPot=%d, dpx=%d, dpy=%d, dpz=%d' % \
# 	            (trackNumb_3,b.shape[0],larmR_b_3.shape[0],uPot_enrgKin_3.shape[0], \
# 		     dpxApprch_3.shape[0],dpyApprch_3.shape[0],dpzApprch_3.shape[0])  
# 
# To draw TMTdpx trajectory (for checking only):
#
            if maxAbsDpxApprch_3 < totAbsDpxApprch_3:
	       maxAbsDpxApprch_3=totAbsDpxApprch_3
	       indxAmaxAbsDpxApprch_3=iA
	       indxBmaxAbsDpxApprch_3=iB
	       trackNumbMaxAbsDpxApprch_3=trackNumb
	       prtclCoorMaxAbsDpx_3=prtclCoorCrrnt_3.copy()   
	       rhoMaxAbsDpxTurn_3=rhoCrrnt[iA,iB]
	       rhoLarmorMaxAbsDpxTurn_3=rho_larm[iA,iB]
	       print 'iA=%d, iB=%d: track %d, points %d' % (indxAmaxAbsDpxApprch_3,\
	             indxBmaxAbsDpxApprch_3,trackNumbMaxAbsDpxApprch_3,pointTrack[trackNumbMaxAbsDpxApprch_3]) 
	       print 'timePoints.shape: ', (prtclCoorMaxAbsDpx_3.shape[0],prtclCoorMaxAbsDpx_3.shape[1])                                            
# End of all calculations for approach_3
            lastTrackNumber_3=trackNumb_3+1                           # quantity of tracks = trackNumb_2 + 1!     
            sumPoints_3 += pointTrack_3[trackNumb_3]
            timeEnd=os.times()
	    cpuTime_3[trackNumb_3]=1.e+6*(float(timeEnd[0])-float(timeStart[0]))     # CPU time , mks
            cpuTimeTotal += cpuTime_3[trackNumb_3] 
# for k in range (int(pointTrack_3[0])):
#   print 'Point %d: dpxIon=%e, dpyIon=%e, dpzIon=%e' % \
#		     (k,dpxApprch_3[k],dpyApprch_3[k],dpzApprch_3[k])
#
#------- End of approach_3 --------------
#

if runFlagApproach_1 == 1:
   print 'Approach_1: maxYcoorElec=%e mkm, maxYcoorIon=%e nm' % (1.e+4*maxYcoorElec,1.e+7*maxYcoorIon)
   print 'Approach_1: for %d tracks number of points is %d' % (lastTrackNumber,sumPoints)
if runFlagApproach_2 == 1:
   print 'Approach_2: for %d tracks number of points is %d' % (lastTrackNumber_2,sumPoints_2)
#   for i in range(b_2LenFirstTrack):
#      print 'action(%d)=%e' % (i,prtclCoor_2[13,i])
if runFlagApproach_3 == 1:
   print 'Approach_3: for %d tracks number of points is %d' % (lastTrackNumber_3,sumPoints_3)

# for i in range(lastTrackNumber):
#    print 'Track %d: larmor turns=%d, cpuTime(mks)=%e, time per turn(mks)=%6.1f' % \
#          (i,larmorNumber[i],cpuTime[i],cpuTime[i]/larmorNumber[i])

print 'cpuTimeTotal(mksec) = %e' % cpuTimeTotal

#
# First track: to compare distances between particles for approach_1 and approach_2,3 (figure 325):
#
if runFlagApproach_1 == 1:
   bArray=np.asarray(b)
   bLenFirstTrack=int(timePoints[0])
   print 'First track length for b: %d' % bLenFirstTrack

if runFlagApproach_2 == 1:
   b_2Array=np.asarray(b_2)
   b_2LenFirstTrack=int(timePoints_2[0])
   print 'First track length for b_2: %d' % b_2LenFirstTrack

if runFlagApproach_3 == 1:
   b_3Array=np.asarray(b_3)
   b_3LenFirstTrack=int(timePoints_3[0])
   print 'First track length for b_3: %d' % b_3LenFirstTrack

#
# Calculation of the difference for distance beteween electron and ion for approach_1 and approach_2:
#
if runFlagApproach_1 == 1 and runFlagApproach_2 == 1:

   diff_b2=np.zeros(b_2LenFirstTrack)
   print 'b_2LenFirstTrack=%d' % b_2LenFirstTrack

   k=0
   nStart=0
   for m in range(b_2LenFirstTrack): 
      for n in range(nStart,bLenFirstTrack):
         if prtclCoor[4,n] == prtclCoor_2[4,m]:
	       diff_b2[m]=bArray[n]-b_2Array[m]
               nStart=n-1
               break	 
         if prtclCoor[4,n] > prtclCoor_2[4,m]:
	    if n == 0:
	       diff_b2[m]=bArray[n]-b_2Array[m]
	    else:
	       bCurrent=bArray[n-1]+(bArray[n]-bArray[n-1])*(prtclCoor_2[4,m]-prtclCoor[4,n-1])/(prtclCoor[4,n]-prtclCoor[4,n-1])
               diff_b2[m]=bCurrent-b_2Array[m]
            nStart=n-1
            break	 
#
# First track: to compare distances between particles for approach_1 and approach_3 (figure 335):
#
if runFlagApproach_1 == 1 and runFlagApproach_3 == 1:

   diff_b3=np.zeros(b_3LenFirstTrack)
   print 'b_3LenFirstTrack=%d' % b_3LenFirstTrack
#
# Calculation of the difference for distance beteween electron and ion for approach_1 and approach_3:
#
   k=0
   nStart=0
   for m in range(b_3LenFirstTrack): 
      for n in range(nStart,bLenFirstTrack):
         if prtclCoor[4,n] == prtclCoor_3[4,m]:
	       diff_b3[m]=bArray[n]-b_3Array[m]
               nStart=n-1
               break	 
         if prtclCoor[4,n] > prtclCoor_3[4,m]:
	    if n == 0:
	       diff_b3[m]=bArray[n]-b_3Array[m]
	    else:
	       bCurrent=bArray[n-1]+(bArray[n]-bArray[n-1])*(prtclCoor_3[4,m]-prtclCoor[4,n-1])/(prtclCoor[4,n]-prtclCoor[4,n-1])
               diff_b3[m]=bCurrent-b_3Array[m]
            nStart=n-1
            break	 

nBins=80

timeStart=os.times()
###################################################
#
#       Data processing of the first approach:
#: layout of arrays xA=uPot_enrgKin, yB=larmR_b to nBins channels
# and arrays zApprch1dpx, zApprch1dpy, zApprch1dpz to (nBins x nBins) channels
#
###################################################   

if runFlagApproach_1 == 1:
   xA=np.zeros(nBins)
   xAedges=np.zeros(nBins+1)
   xAnumb=np.zeros(nBins)
   xAstep=(maxUpot_enrgKin_1-minUpot_enrgKin_1)/nBins
   for i in range(nBins+1):
      xAedges[i]=minUpot_enrgKin_1+xAstep*i

   yB=np.zeros(nBins)
   yBedges=np.zeros(nBins+1)
   yBnumb=np.zeros(nBins)
   yBstep=(maxLarmR_b_1-minLarmR_b_1)/nBins
   for i in range(nBins+1):
      yBedges[i]=minLarmR_b_1+yBstep*i

   zApprch1dpNumb=np.zeros((nBins,nBins))
   zApprch1dpx=np.zeros((nBins,nBins))
   zApprch1dpy=np.zeros((nBins,nBins))
   zApprch1dpz=np.zeros((nBins,nBins))

   for nPoint in range(int(sumPoints)):
      for iA in range(nBins):
         searchAflag=0
         if (xAedges[iA] <= uPot_enrgKin[nPoint] < xAedges[iA+1]):
            if xAnumb[iA] == 0:
	       xA[iA]=uPot_enrgKin[nPoint]                         # log10(Upot/Ekin)
	    else:
	       xA[iA]=(xA[iA]*xAnumb[iA]+uPot_enrgKin[nPoint])/(xAnumb[iA]+1)        # averaging inside bin iA
            xAnumb[iA] += 1
	    searchAflag=1
	    break
      if searchAflag == 0:
         xA[nBins-1]=(xA[nBins-1]*xAnumb[nBins-1]+uPot_enrgKin[nPoint])/(xAnumb[nBins-1]+1)            # averaging inside bin iA 
         xAnumb[nBins-1] += 1
      for iB in range(nBins):
         searchBflag=0
         if (yBedges[iB] <= larmR_b[nPoint] < yBedges[iB+1]):
            if yBnumb[iB] == 0:
	       yB[iB]=larmR_b[nPoint]                              # log10(Rlarm/b)
	    else:
	       yB[iB]=(yB[iB]*yBnumb[iB]+larmR_b[nPoint])/(yBnumb[iB]+1)             # averaging inside bin iB  
            yBnumb[iB] += 1
	    searchBflag=1
	    break
      if searchBflag == 0:
         yB[nBins-1]=(yB[nBins-1]*yBnumb[nBins-1]+larmR_b[nPoint])/(yBnumb[nBins-1]+1)                 # averaging inside bin iB 
         yBnumb[nBins-1] += 1
      if zApprch1dpNumb[iA,iB] == 0:
         zApprch1dpx[iA,iB]=dpxApprch_1[nPoint]
         zApprch1dpy[iA,iB]=dpyApprch_1[nPoint]
         zApprch1dpz[iA,iB]=dpzApprch_1[nPoint]
      else:
         zApprch1dpx[iA,iB]= \
         (zApprch1dpx[iA,iB]*zApprch1dpNumb[iA,iB]+dpxApprch_1[nPoint])/(zApprch1dpNumb[iA,iB]+1)  # averaging inside rectangle iA,iB
         zApprch1dpy[iA,iB]= \
         (zApprch1dpy[iA,iB]*zApprch1dpNumb[iA,iB]+dpyApprch_1[nPoint])/(zApprch1dpNumb[iA,iB]+1)  # averaging inside rectangle iA,iB
         zApprch1dpz[iA,iB]= \
         (zApprch1dpz[iA,iB]*zApprch1dpNumb[iA,iB]+dpzApprch_1[nPoint])/(zApprch1dpNumb[iA,iB]+1)  # averaging inside rectangle iA,iB
      zApprch1dpNumb[iA,iB] += 1

# Some checkings:
   sumAnumb=0
   sumBnumb=0
   sumCnumb=0
   for iA in range(nBins):
      sumAnumb += xAnumb[iA]
      sumBnumb += yBnumb[iA]
      for iB in range(nBins):
         sumCnumb += zApprch1dpNumb[iA,iB]
#         print 'iA=%d, iB=%d: xA=%f (%d), xB=%f (%d), zApprch1dpx=%e' % \
#               (iA,iB,xA[iA],xAnumb[iA],yB[iB],yBnumb[iB],zApprch1dpx[iA,iB])
   
   print 'sumA=%d, sumB=%d, sumC=%d; sumPoints=%d' % (sumAnumb,sumBnumb,sumCnumb,int(sumPoints)) 


###################################################
#
#       Data processing of the second approach:
#: layout of arrays xA_2=uPot_enrgKin_2, yB_2=larmR_b_2 to nBins channels
# and arrays zApprch2dpx, zApprch2dpy, zApprch2dpz to (nBins x nBins) channels
#
###################################################   

if runFlagApproach_2 == 1:
   xA_2=np.zeros(nBins)
   xAedges_2=np.zeros(nBins+1)
   xAnumb_2=np.zeros(nBins)
   xAstep=(maxUpot_enrgKin_2-minUpot_enrgKin_2)/nBins
   for i in range(nBins+1):
      xAedges_2[i]=minUpot_enrgKin_2+xAstep*i

   yB_2=np.zeros(nBins)
   yBedges_2=np.zeros(nBins+1)
   yBnumb_2=np.zeros(nBins)
   yBstep=(maxLarmR_b_2-minLarmR_b_2)/nBins
   for i in range(nBins+1):
      yBedges_2[i]=minLarmR_b_2+yBstep*i

   zApprch2dpNumb=np.zeros((nBins,nBins))
   zApprch2dpx=np.zeros((nBins,nBins))
   zApprch2dpy=np.zeros((nBins,nBins))
   zApprch2dpz=np.zeros((nBins,nBins))

   for nPoint in range(int(sumPoints_2)):
      for iA in range(nBins):
         searchAflag=0
         if (xAedges_2[iA] <= uPot_enrgKin_2[nPoint] < xAedges_2[iA+1]):
            if xAnumb_2[iA] == 0:
	       xA_2[iA]=uPot_enrgKin_2[nPoint]                     # log10(Upot/Ekin)
	    else:
	       xA_2[iA]=(xA_2[iA]*xAnumb_2[iA]+uPot_enrgKin_2[nPoint])/(xAnumb_2[iA]+1)                # averaging inside bin iA_2
            xAnumb_2[iA] += 1
	    searchAflag=1
	    break
      if searchAflag == 0:
         xA_2[nBins-1]=(xA_2[nBins-1]*xAnumb_2[nBins-1]+uPot_enrgKin_2[nPoint])/(xAnumb_2[nBins-1]+1)  # averaging inside bin iA_2 
         xAnumb_2[nBins-1] += 1
      for iB in range(nBins):
         searchBflag=0
         if (yBedges_2[iB] <= larmR_b_2[nPoint] < yBedges_2[iB+1]):
            if yBnumb_2[iB] == 0:
	       yB_2[iB]=larmR_b_2[nPoint]                          # log10(Rlarm/b)
	    else:
	       yB_2[iB]=(yB_2[iB]*yBnumb_2[iB]+larmR_b_2[nPoint])/(yBnumb_2[iB]+1)                     # averaging inside bin iB_2  
            yBnumb_2[iB] += 1
	    searchBflag=1
	    break
      if searchBflag == 0:
         yB_2[nBins-1]=(yB_2[nBins-1]*yBnumb_2[nBins-1]+larmR_b_2[nPoint])/(yBnumb_2[nBins-1]+1)       # averaging inside bin iB_2 
         yBnumb_2[nBins-1] += 1
      if zApprch2dpNumb[iA,iB] == 0:
         zApprch2dpx[iA,iB]=dpxApprch_2[nPoint]
         zApprch2dpy[iA,iB]=dpyApprch_2[nPoint]
         zApprch2dpz[iA,iB]=dpzApprch_2[nPoint]
      else:
         zApprch2dpx[iA,iB]= \
         (zApprch2dpx[iA,iB]*zApprch2dpNumb[iA,iB]+dpxApprch_2[nPoint])/(zApprch2dpNumb[iA,iB]+1)  # averaging inside rectangle iA,iB
         zApprch2dpy[iA,iB]= \
         (zApprch2dpy[iA,iB]*zApprch2dpNumb[iA,iB]+dpyApprch_2[nPoint])/(zApprch2dpNumb[iA,iB]+1)  # averaging inside rectangle iA,iB
         zApprch2dpz[iA,iB]= \
         (zApprch2dpz[iA,iB]*zApprch2dpNumb[iA,iB]+dpzApprch_2[nPoint])/(zApprch2dpNumb[iA,iB]+1)  # averaging inside rectangle iA,iB
      zApprch2dpNumb[iA,iB] += 1

# Some checkings:
   sumAnumb_2=0
   sumBnumb_2=0
   sumCnumb_2=0
   for iA in range(nBins):
      sumAnumb_2 += xAnumb_2[iA]
      sumBnumb_2 += yBnumb_2[iA]
      for iB in range(nBins):
         sumCnumb_2 += zApprch2dpNumb[iA,iB]
#         print 'iA=%d, iB=%d: xA=%f (%d), xB=%f (%d), zApprch1dpx=%e' % \
#               (iA,iB,xA_2[iA],xAnumb_2[iA],yB_2[iB],yBnumb_2[iB],zApprch2dpx[iA,iB])
   
   print 'sumA_2=%d, sumB_2=%d, sumC_2=%d; sumPoints_2=%d' % (sumAnumb_2,sumBnumb_2,sumCnumb_2,int(sumPoints_2)) 

###################################################
#
#       Data processing of the third approach:
#: layout of arrays xA_3=uPot_enrgKin_3, yB_3=larmR_b_3 to nBins channels
# and arrays zApprch3dpx, zApprch3dpy, zApprch3dpz to (nBins x nBins) channels
#
###################################################   

if runFlagApproach_3 == 1:
   xA_3=np.zeros(nBins)
   xAedges_3=np.zeros(nBins+1)
   xAnumb_3=np.zeros(nBins)
   xAstep=(maxUpot_enrgKin_3-minUpot_enrgKin_3)/nBins
   for i in range(nBins+1):
      xAedges_3[i]=minUpot_enrgKin_3+xAstep*i

   yB_3=np.zeros(nBins)
   yBedges_3=np.zeros(nBins+1)
   yBnumb_3=np.zeros(nBins)
   yBstep=(maxLarmR_b_3-minLarmR_b_3)/nBins
   for i in range(nBins+1):
      yBedges_3[i]=minLarmR_b_3+yBstep*i

   zApprch3dpNumb=np.zeros((nBins,nBins))
   zApprch3dpx=np.zeros((nBins,nBins))
   zApprch3dpy=np.zeros((nBins,nBins))
   zApprch3dpz=np.zeros((nBins,nBins))

   for nPoint in range(int(sumPoints_3)):
      for iA in range(nBins):
         searchAflag=0
         if (xAedges_3[iA] <= uPot_enrgKin_3[nPoint] < xAedges_3[iA+1]):
            if xAnumb_3[iA] == 0:
	       xA_3[iA]=uPot_enrgKin_3[nPoint]                     # log10(Upot/Ekin)
	    else:
	       xA_3[iA]=(xA_3[iA]*xAnumb_3[iA]+uPot_enrgKin_3[nPoint])/(xAnumb_3[iA]+1)                # averaging inside bin iA_3
            xAnumb_3[iA] += 1
	    searchAflag=1
	    break
      if searchAflag == 0:
         xA_3[nBins-1]=(xA_3[nBins-1]*xAnumb_3[nBins-1]+uPot_enrgKin_3[nPoint])/(xAnumb_3[nBins-1]+1)  # averaging inside bin iA_3 
         xAnumb_3[nBins-1] += 1
      for iB in range(nBins):
         searchBflag=0
         if (yBedges_3[iB] <= larmR_b_3[nPoint] < yBedges_3[iB+1]):
            if yBnumb_3[iB] == 0:
	       yB_3[iB]=larmR_b_3[nPoint]                          # log10(Rlarm/b)
	    else:
	       yB_3[iB]=(yB_3[iB]*yBnumb_3[iB]+larmR_b_3[nPoint])/(yBnumb_3[iB]+1)                     # averaging inside bin iB_3  
            yBnumb_3[iB] += 1
	    searchBflag=1
	    break
      if searchBflag == 0:
         yB_3[nBins-1]=(yB_3[nBins-1]*yBnumb_3[nBins-1]+larmR_b_3[nPoint])/(yBnumb_3[nBins-1]+1)       # averaging inside bin iB_3 
         yBnumb_3[nBins-1] += 1
      if zApprch3dpNumb[iA,iB] == 0:
         zApprch3dpx[iA,iB]=dpxApprch_3[nPoint]
         zApprch3dpy[iA,iB]=dpyApprch_3[nPoint]
         zApprch3dpz[iA,iB]=dpzApprch_3[nPoint]
      else:
         zApprch3dpx[iA,iB]= \
         (zApprch3dpx[iA,iB]*zApprch3dpNumb[iA,iB]+dpxApprch_3[nPoint])/(zApprch3dpNumb[iA,iB]+1)  # averaging inside rectangle iA,iB
         zApprch3dpy[iA,iB]= \
         (zApprch3dpy[iA,iB]*zApprch3dpNumb[iA,iB]+dpyApprch_3[nPoint])/(zApprch3dpNumb[iA,iB]+1)  # averaging inside rectangle iA,iB
         zApprch3dpz[iA,iB]= \
         (zApprch3dpz[iA,iB]*zApprch3dpNumb[iA,iB]+dpzApprch_3[nPoint])/(zApprch3dpNumb[iA,iB]+1)  # averaging inside rectangle iA,iB
#	 print 'nPoint=%d: dpzApprch_3(%d,%d)=%e --> zApprch3dpz=%e' % (nPoint,iA,iB,dpzApprch_3[nPoint],zApprch3dpz[iA,iB])
      zApprch3dpNumb[iA,iB] += 1

# Some checkings:
   sumAnumb_3=0
   sumBnumb_3=0
   sumCnumb_3=0
   for iA in range(nBins):
      sumAnumb_3 += xAnumb_3[iA]
      sumBnumb_3 += yBnumb_3[iA]
      for iB in range(nBins):
         sumCnumb_3 += zApprch3dpNumb[iA,iB]
#         print 'iA=%d, iB=%d: xA=%f (%d), xB=%f (%d), zApprch1dpx=%e' % \
#               (iA,iB,xA_3[iA],xAnumb_3[iA],yB_3[iB],yBnumb_3[iB],zApprch3dpx[iA,iB])
   
   print 'sumA_3=%d, sumB_3=%d, sumC_3=%d; sumPoints_3=%d' % (sumAnumb_3,sumBnumb_3,sumCnumb_3,int(sumPoints_3)) 

timeEnd=os.times()
runTime=1.e+6*(float(timeEnd[0])-float(timeStart[0]))              # CPU time , mks
print 'runTime(mksec) = %e' % runTime

#
# Checking of the first trajectories:
#
turns=10                                                           # Number of larmor turns for drawing 
pointsTurns=turns*stepsNumberOnGyro                                # Number of points for drawing
pointsStartLarmor=turns*stepsNumberOnGyro                          # Number of points for drawing
lengthArrowElc=4
pBegArrw=pointsTurns/2
pEndArrw=pointsTurns/2+50

                                                 # Plot everything about dp transf. if = 1, orherwise don't plot

if runFlagApproach_1 == 1 and plotFlagTracks == 1:
   pointsEndLarmor=int(larmorNumber[0])*stepsNumberOnGyro             # Number of points for drawing

   fig10=plt.figure(10)
   ax10=fig10.gca(projection='3d')
   ax10.plot(1.e+4*prtclCoor[0,0:pointsStartLarmor],1.e+4*prtclCoor[2,0:pointsStartLarmor], \
             1.e+4*prtclCoor[4,0:pointsStartLarmor],'-r',linewidth=2)
   ax10.plot(1.e+4*prtclCoor[0,0:lengthArrowElc],1.e+4*prtclCoor[2,0:lengthArrowElc], \
             1.e+4*prtclCoor[4,0:lengthArrowElc],'-b',linewidth=2)
   ax10.plot(1.e+4*prtclCoor[0,pointsStartLarmor-lengthArrowElc:pointsStartLarmor], \
             1.e+4*prtclCoor[2,pointsStartLarmor-lengthArrowElc:pointsStartLarmor], \
             1.e+4*prtclCoor[4,pointsStartLarmor-lengthArrowElc:pointsStartLarmor],'-b',linewidth=2)
   plt.xlabel('x, $\mu m$',color='m',fontsize=16)
   plt.ylabel('y, $\mu m$',color='m',fontsize=16)
   ax10.set_zlabel('z, $\mu m$',color='m',fontsize=16)
   plt.title( \
      ('Approach-1: First Electron Trajectory (Start; $N_{Larm}=$%d):\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m' \
       % (larmorNumber[0],1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)),color='m',fontsize=16)

   fig15=plt.figure(15)
   ax15=fig15.gca(projection='3d')
   ax15.plot(1.e+7*prtclCoor[6,0:pointsStartLarmor],1.e+7*prtclCoor[8,0:pointsStartLarmor], \
             1.e+7*prtclCoor[10,0:pointsStartLarmor],'-b',linewidth=2)
   ax15.plot(1.e+7*prtclCoor[6,pBegArrw:pEndArrw],1.e+7*prtclCoor[8,pBegArrw:pEndArrw], \
             1.e+7*prtclCoor[10,pBegArrw:pEndArrw],'-b',linewidth=4)
#     ax15.plot([arrwBegxIon,arrwEndxIon],[arrwBegyIon,arrwEndyIon],[arrwBegzIon,arrwEndzIon],color='r',alpha=0.8,lw=2)
   plt.xlabel('x, $nm$',color='m',fontsize=16)
   plt.ylabel('y, $nm$',color='m',fontsize=16)
   ax15.set_zlabel('z, $nm$',color='m',fontsize=16)
   plt.title('Approach-1: First Ion Trajectory (Start)',color='m',fontsize=16)


   pBeg=pointsEndLarmor-pointsTurns

   fig20=plt.figure(20)
   ax20=fig20.gca(projection='3d')
   ax20.plot(1.e+4*prtclCoor[0,pBeg:pointsEndLarmor], \
             1.e+4*prtclCoor[2,pBeg:pointsEndLarmor], \
             1.e+4*prtclCoor[4,pBeg:pointsEndLarmor],'-r',linewidth=2)
   ax20.plot(1.e+4*prtclCoor[0,pointsEndLarmor-lengthArrowElc:pointsEndLarmor], \
             1.e+4*prtclCoor[2,pointsEndLarmor-lengthArrowElc:pointsEndLarmor], \
             1.e+4*prtclCoor[4,pointsEndLarmor-lengthArrowElc:pointsEndLarmor],'-b',linewidth=2)
   ax20.plot(1.e+4*prtclCoor[0,pBeg:pBeg+lengthArrowElc], 1.e+4*prtclCoor[2,pBeg:pBeg+lengthArrowElc], \
             1.e+4*prtclCoor[4,pBeg:pBeg+lengthArrowElc],'-b',linewidth=2)
   plt.xlabel('x, $\mu m$',color='m',fontsize=16)
   plt.ylabel('y, $\mu m$',color='m',fontsize=16)
   ax20.set_zlabel('z, $\mu m$',color='m',fontsize=16)
   plt.title( \
      ('Approach-1: First Electron Trajectory (End; $N_{Larm}=$%d):\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m' \
       % (larmorNumber[0],1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)),color='m',fontsize=16)

   pEndArrw=pointsEndLarmor-pointsTurns/2
   pBegArrw=pEndArrw-50

   fig25=plt.figure(25)
   ax25=fig25.gca(projection='3d')
   ax25.plot(1.e+7*prtclCoor[6,pBeg:pointsEndLarmor],1.e+7*prtclCoor[8,pBeg:pointsEndLarmor], \
             1.e+7*prtclCoor[10,pBeg:pointsEndLarmor],'-b',linewidth=2)
   ax25.plot(1.e+7*prtclCoor[6,pBegArrw:pEndArrw],1.e+7*prtclCoor[8,pBegArrw:pEndArrw], \
             1.e+7*prtclCoor[10,pBegArrw:pEndArrw],'-b',linewidth=4)
   plt.xlabel('x, $nm$',color='m',fontsize=16)
   plt.ylabel('y, $nm$',color='m',fontsize=16)
   ax25.set_zlabel('z, $nm$',color='m',fontsize=16)
   plt.title('Approach-1: First Ion Trajectory (End)',color='m',fontsize=16)

   plt.figure(30)
   plt.plot(range(int(pointsEndLarmor)),1.e4*prtclCoor[4,0:pointsEndLarmor],'-r',linewidth=2)
   plt.xlabel('Points',color='m',fontsize=16)
   plt.ylabel('z, $\mu$m',color='m',fontsize=16)
   plt.title(('Approach-1: First Electron Trajectory ($N_{Larm}=$%d):\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m' % \
              (larmorNumber[0],1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)),color='m',fontsize=16)
   plt.grid(True)

   plt.figure(40)
   plt.plot(1.e4*prtclCoor[4,0:pointsEndLarmor],1.e4*prtclCoor[12,0:pointsEndLarmor],'-r',linewidth=2)
   plt.xlabel('z, $\mu$m',color='m',fontsize=16)
   plt.ylabel('Distance Between Particles, $\mu$m',color='m',fontsize=16)
   plt.title(('Approach-1: First Electron Trajectory ($N_{Larm}=$%d):\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m' % \
              (larmorNumber[0],1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)),color='m',fontsize=16)
   plt.grid(True)

   pointsDraw=int(pointTrack[trackNumbMaxAbsDpxApprch])             # Number of points for drawing

   plt.figure(1030)
   plt.plot(range(int(pointsDraw)),1.e4*prtclCoorMaxAbsDpx[4,0:pointsDraw],'-r',linewidth=2)
   plt.xlabel('Points',color='m',fontsize=16)
   plt.ylabel('z, $\mu$m',color='m',fontsize=16)
   plt.title(('Approach-1: First Electron Track %d ($N_{Larm}=$%d):\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m' % \
              (trackNumbMaxAbsDpxApprch,larmorNumber[trackNumbMaxAbsDpxApprch],1.e+4*rhoMaxAbsDpxTurn, \
	       1.e+4*rhoLarmorMaxAbsDpxTurn)),color='m',fontsize=16)
   plt.grid(True)

   plt.figure(1040)
   plt.plot(1.e4*prtclCoorMaxAbsDpx[4,0:pointsDraw],1.e4*prtclCoorMaxAbsDpx[12,0:pointsDraw],'.r',linewidth=2)
   plt.xlabel('z, $\mu$m',color='m',fontsize=16)
   plt.ylabel('Distance Between Particles, $\mu$m',color='m',fontsize=16)
   plt.title(('Approach-1: Electron Track %d ($N_{Larm}=$%d):\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m' % \
              (trackNumbMaxAbsDpxApprch,larmorNumber[trackNumbMaxAbsDpxApprch],1.e+4*rhoMaxAbsDpxTurn, \
	       1.e+4*rhoLarmorMaxAbsDpxTurn)), color='m',fontsize=16)
   plt.grid(True)

   plt.figure(55)
   plt.plot(uPot_enrgKin,larmR_b,'.r')
   plt.xlim([minUpot_enrgKin_1-.1,maxUpot_enrgKin_1+.1])
   plt.ylim([minLarmR_b_1-.1,maxLarmR_b_1+.1])
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   plt.title(('Approach-1: Map for Transfered Momenta Calculations\nTracks: %d' % lastTrackNumber), color='m',fontsize=16)
   plt.grid(True)

if runFlagApproach_2 == 1 and plotFlagTracks == 1:
   pointsGCtot_2=int(pointTrack_2[0])
   print 'pointsGCtot_2=%d' % pointsGCtot_2

   turns_2=10*turns
   pBegArrw_2=pointsGCtot_2-turns_2/2
   pEndArrw_2=pBegArrw_2+25

   fig225=plt.figure(225)
   ax225=fig225.gca(projection='3d')
   ax225.plot(1.e+7*prtclCoor_2[6,pointsGCtot_2-turns_2:pointsGCtot_2],1.e+7*prtclCoor_2[8,pointsGCtot_2-turns_2:pointsGCtot_2], \
              1.e+7*prtclCoor_2[10,pointsGCtot_2-turns_2:pointsGCtot_2],'-b',linewidth=2)
   ax225.plot(1.e+7*prtclCoor_2[6,pBegArrw_2:pEndArrw_2],1.e+7*prtclCoor_2[8,pBegArrw_2:pEndArrw_2], \
              1.e+7*prtclCoor_2[10,pBegArrw_2:pEndArrw_2],'-b',linewidth=4)
   plt.xlabel('x, $nm$',color='m',fontsize=16)
   plt.ylabel('y, $nm$',color='m',fontsize=16)
   ax225.set_zlabel('z, $nm$',color='m',fontsize=16)
   plt.title('Approach-2: First Ion Trajectory (End)',color='m',fontsize=16)

if runFlagApproach_1 == 1 and runFlagApproach_2 == 1 and plotFlagTracks == 1:
   fig325=plt.figure(325)
   plt.plot(1.e+4*prtclCoor_2[4,0:b_2LenFirstTrack],1.e+4*diff_b2[0:b_2LenFirstTrack],'.r',linewidth=2)
   plt.xlabel('z, $\mu m$',color='m',fontsize=16)
   plt.ylabel('$\Delta b=b_{Apprch-1}-b_{Apprch-2}$, $\mu m$',color='m',fontsize=16)
   plt.title('First Trajectory: $\Delta b$ (Difference for Distances $b$ between Particles)',color='m',fontsize=16)
   plt.xlim([1.e+4*prtclCoor_2[4,0]-50,1.e+4*prtclCoor_2[4,b_2LenFirstTrack-1]+50])
   plt.ylim([-.5,.5])
   plt.grid(True)

if plotFlagWorkingArrays == 1:

   fig1000=plt.figure(1000)
   plt.plot(range(nB),1.e+4*halfLintr[0,:],'-xr',linewidth=2)
   plt.plot(range(nB),1.e+4*halfLintr[1,:],'-xb',linewidth=2)
   plt.plot(range(nB),1.e+4*halfLintr[2,:],'-xm',linewidth=2)
   plt.plot(range(nB),1.e+4*halfLintr[3,:],'-xg',linewidth=2)
   plt.plot(range(nB),1.e+4*halfLintr[4,:],'-xk',linewidth=2)
#   plt.title(('Approach-1: TransferedMomentum $dP_x$\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{24}$ $g \cdot cm/sec$)' % \
#              (lastTrackNumber,dpxDataMax)), color='m',fontsize=16)
   plt.xlabel('$iB$',color='m',fontsize=16)
   plt.ylabel('Half Length of Trajectory, $\mu$m',color='m',fontsize=16)
   plt.legend(['$iA=0$','$iA=1$','$iA=2$','$iA=3$','$iA=4$'],loc='upper right',fontsize=16)
   plt.grid(True)
#
# Normalization of the arrays zApprch1dpx,zApprch1dpy,zApprch1dpz
# to improve visual presentation of the results:
#
if runFlagApproach_1 == 1 and plotFlagDpTransf == 1:
   visPrezApprch1dpx=np.zeros((nBins,nBins))
   visPrezApprch1dpy=np.zeros((nBins,nBins))
   visPrezApprch1dpz=np.zeros((nBins,nBins))

   dpxDataMax=-1.e24                                               # maximum for array of zApprch1dpx
   indxAdpxDataMax=0                                               # index A of the maximum for array of zApprch1dpx
   indxBdpxDataMax=0                                               # index B of the maximum for array of zApprch1dpx
   dpyDataMax=-1.e24                                               # maximum for array of zApprch1dpx
   indxAdpyDataMax=0                                               # index A of the maximum for array of zApprch1dpx
   indxBdpyDataMax=0                                               # index B of the maximum for array of zApprch1dpx
   dpzDataMax=-1.e24                                               # maximum for array of zApprch1dpx
   indxAdpzDataMax=0                                               # index A of the maximum for array of zApprch1dpx
   indxBdpzDataMax=0                                               # index B of the maximum for array of zApprch1dpx

   normFactor=1.e24
   for iA in range(nBins):
      for iB in range(nBins):
         visPrezApprch1dpx[iA,iB]=normFactor*zApprch1dpx[iA,iB]
         visPrezApprch1dpy[iA,iB]=normFactor*zApprch1dpy[iA,iB]
         visPrezApprch1dpz[iA,iB]=normFactor*zApprch1dpz[iA,iB]
         if dpxDataMax < abs(normFactor*zApprch1dpx[iA,iB]):
            dpxDataMax=abs(normFactor*zApprch1dpx[iA,iB])
            indxAdpxDataMax=iA
            indxBdpxDataMax=iB
         if dpyDataMax < abs(normFactor*zApprch1dpy[iA,iB]):
            dpyDataMax=abs(normFactor*zApprch1dpy[iA,iB])
            indxAdpyDataMax=iA
            indxBdpyDataMax=iB
         if dpzDataMax < abs(normFactor*zApprch1dpz[iA,iB]):
            dpzDataMax=abs(normFactor*zApprch1dpz[iA,iB])
            indxAdpzDataMax=iA
            indxBdpzDataMax=iB

   print '|Maximum| (approach-1): dpx=%e,  dpy=%e,  dpz=%e' % (dpxDataMax,dpyDataMax,dpzDataMax)

   X,Y=np.meshgrid(xA,yB) 

   fig70=plt.figure(70)
   ax70=fig70.gca(projection='3d')
#    surf=ax70.plot_surface(X,Y,visPrezApprch1dpx,cmap=cm.coolwarm,linewidth=0,antialiased=False)
   surf=ax70.plot_surface(X,Y,visPrezApprch1dpx,cmap=cm.jet,linewidth=0,antialiased=False)
   plt.title(('Approach-1: Transfered Momentum $dP_x$\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{24}$ $g \cdot cm/sec$)' % \
              (lastTrackNumber,dpxDataMax)), color='m',fontsize=16)
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   ax70.set_zlabel('$dP_x \cdot 10^{-24}$; $g \cdot cm/sec$',color='m',fontsize=16)
   cb = fig70.colorbar(surf)
   plt.grid(True)

   fig75=plt.figure(75)
   ax=fig75.add_subplot(111)                                       # for contours plotting
#    mapDpx=ax.contourf(X,Y,zApprch1dpx)   
   mapDpx=ax.contourf(X,Y,visPrezApprch1dpx)   
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   plt.title( \
      ('Approach-1: Transfered Momentum $dP_x$ $(\cdot 10^{-24}$; $g \cdot cm/sec)$\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-24}$)' % \
          (lastTrackNumber,dpxDataMax)), color='m',fontsize=14)
   fig75.colorbar(mapDpx)

   fig80=plt.figure(80)
   ax80=fig80.gca(projection='3d')
#    surf=ax80.plot_surface(X,Y,visPrezApprch1dpx,cmap=cm.coolwarm,linewidth=0,antialiased=False)
   surf=ax80.plot_surface(X,Y,visPrezApprch1dpy,cmap=cm.jet,linewidth=0,antialiased=False)
   plt.title(('Approach-1: Transfered Momentum $dP_y$\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{24}$ $g \cdot cm/sec$)' % \
              (lastTrackNumber,dpyDataMax)), color='m',fontsize=16)
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   ax80.set_zlabel('$dP_y \cdot 10^{-24}$; $g \cdot cm/sec$',color='m',fontsize=16)
   cb = fig80.colorbar(surf)
   plt.grid(True)

   fig85=plt.figure(85)
   ax=fig85.add_subplot(111)                                       # for contours plotting
#    mapDpx=ax.contourf(X,Y,zApprch1dpx)   
   mapDpx=ax.contourf(X,Y,visPrezApprch1dpy)   
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   plt.title( \
      ('Approach-1: Transfered Momentum $dP_y$ $(\cdot 10^{-24}$; $g \cdot cm/sec)$\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-24}$)' % \
        (lastTrackNumber,dpyDataMax)), color='m',fontsize=14)
   fig85.colorbar(mapDpx)

   fig90=plt.figure(90)
   ax90=fig90.gca(projection='3d')
#    surf=ax90.plot_surface(X,Y,visPrezApprch1dpx,cmap=cm.coolwarm,linewidth=0,antialiased=False)
   surf=ax90.plot_surface(X,Y,visPrezApprch1dpz,cmap=cm.jet,linewidth=0,antialiased=False)
   plt.title( \
      ('Approach-1: Transfered Momentum $dP_z$\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{24}$ $g \cdot cm/sec$)' % \
       (lastTrackNumber,dpzDataMax)), color='m',fontsize=16)
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   ax90.set_zlabel('$dP_z \cdot 10^{-24}$; $g \cdot cm/sec$',color='m',fontsize=16)
   cb = fig90.colorbar(surf)
   plt.grid(True)

   fig95=plt.figure(95)
   ax=fig95.add_subplot(111)                                       # for contours plotting
#    mapDpx=ax.contourf(X,Y,zApprch1dpx)   
   mapDpx=ax.contourf(X,Y,visPrezApprch1dpz)   
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   plt.title( \
      ('Approach-1: Transfered Momentum $dP_z$ $(\cdot 10^{-24}$; $g \cdot cm/sec)$\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-24}$)' % \
       (lastTrackNumber,dpzDataMax)), color='m',fontsize=14)
   fig95.colorbar(mapDpx)

if runFlagApproach_2 == 1 and plotFlagTracks == 1:
   pointsGCturns=turns                                             # Number of points for drawing ("GC")

   pointsDraw=int(pointTrack_2[trackNumbMaxAbsDpxApprch_2])           # Number of points for drawing

   plt.figure(1130)
   plt.plot(range(int(pointsDraw)),1.e4*prtclCoorMaxAbsDpx_2[4,0:pointsDraw],'-xr',linewidth=2)
   plt.xlabel('Points',color='m',fontsize=16)
   plt.ylabel('z, $\mu$m',color='m',fontsize=16)
   plt.title(('Approach-2: First Electron Track %d ($N_{Larm}=$%d):\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m' % \
              (trackNumbMaxAbsDpxApprch_2,larmorNumber[trackNumbMaxAbsDpxApprch_2],1.e+4*rhoMaxAbsDpxTurn_2, \
	       1.e+4*rhoLarmorMaxAbsDpxTurn_2)),color='m',fontsize=16)
   plt.grid(True)

   plt.figure(1140)
   plt.plot(1.e4*prtclCoorMaxAbsDpx_2[4,0:pointsDraw],1.e4*prtclCoorMaxAbsDpx_2[12,0:pointsDraw],'-xr',linewidth=2)
   plt.xlabel('z, $\mu$m',color='m',fontsize=16)
   plt.ylabel('Distance Between Particles, $\mu$m',color='m',fontsize=16)
   plt.title(('Approach-2: Track %d ($N_{Larm}=$%d):\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m' % \
              (trackNumbMaxAbsDpxApprch_2,larmorNumber[trackNumbMaxAbsDpxApprch_2],\
	       1.e+4*rhoMaxAbsDpxTurn_2,1.e+4*rhoLarmorMaxAbsDpxTurn_2)), color='m',fontsize=16)
   plt.grid(True)

   pBegArrw=pointsGCturns/2-1

if runFlagApproach_1 == 1 and runFlagApproach_2 == 1 and plotFlagTracks == 1:
   fig110=plt.figure(110)
   ax110=fig110.gca(projection='3d')
   ax110.plot(1.e+4*prtclCoor[0,0:pointsStartLarmor],1.e+4*prtclCoor[2,0:pointsStartLarmor], \
              1.e+4*prtclCoor[4,0:pointsStartLarmor],'-r',linewidth=2)
   ax110.plot(1.e+4*prtclCoor[0,0:lengthArrowElc],1.e+4*prtclCoor[2,0:lengthArrowElc], \
              1.e+4*prtclCoor[4,0:lengthArrowElc],'-b',linewidth=2)
   ax110.plot(1.e+4*prtclCoor[0,pointsStartLarmor-lengthArrowElc:pointsStartLarmor], \
              1.e+4*prtclCoor[2,pointsStartLarmor-lengthArrowElc:pointsStartLarmor], \
              1.e+4*prtclCoor[4,pointsStartLarmor-lengthArrowElc:pointsStartLarmor],'-b',linewidth=2)
   ax110.plot(1.e+4*prtclCoor_2[0,0:pointsGCturns],1.e+7*prtclCoor_2[2,0:pointsGCturns], \
              1.e+4*prtclCoor_2[4,0:pointsGCturns],'-b',linewidth=2)
   ax110.plot(1.e+4*prtclCoor_2[0,pBegArrw:pBegArrw+3],1.e+4*prtclCoor_2[2,pBegArrw:pBegArrw+3], \
              1.e+4*prtclCoor_2[4,pBegArrw:pBegArrw+3],'-b',linewidth=4)
   plt.xlabel('x, $\mu m$',color='m',fontsize=16)
   plt.ylabel('y, ($\mu m$',color='m',fontsize=16)
   ax110.set_zlabel('z, $\mu m$',color='m',fontsize=16)
   plt.title(('First Electron Trajectory (Start; $N_{Larm}=$%d):\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m' \
           % (larmorNumber[0],1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)),color='m',fontsize=16)
   ax110.text(0.,0.,-5425.,'Blue is Trajectory in the System of "Guiding Center"',color='b',fontsize=16)


#    for m in range(pointsEndLarmor-50,pointsEndLarmor):
#       print 'Electron: z(%d)=%e' % (m,1.e+4*prtclCoor[4,m])
#    for m in range(pointsGCtot_2-turns,pointsGCtot_2):
#       print 'Electron_gc: z(%d)=%e' % (m,1.e+4*prtclCoor_2[4,m])

   print 'Last point of the first track: z(%d)=%e mkm, z_gc(%d)=%e mkm' % \
         (pointsEndLarmor,1.e+4*prtclCoor[4,pointsEndLarmor-1],pointsGCtot_2,1.e+4*prtclCoor_2[4,pointsGCtot_2-1])

   fig120=plt.figure(120) 
   ax120=fig120.gca(projection='3d')
   ax120.plot(1.e+4*prtclCoor[0,pBeg:pointsEndLarmor], \
              1.e+4*prtclCoor[2,pBeg:pointsEndLarmor], \
              1.e+4*prtclCoor[4,pBeg:pointsEndLarmor],'-r',linewidth=2)
   ax120.plot(1.e+4*prtclCoor[0,pointsEndLarmor-lengthArrowElc:pointsEndLarmor], \
              1.e+4*prtclCoor[2,pointsEndLarmor-lengthArrowElc:pointsEndLarmor], \
              1.e+4*prtclCoor[4,pointsEndLarmor-lengthArrowElc:pointsEndLarmor],'-r',linewidth=4)
   ax120.plot(1.e+4*prtclCoor[0,pBeg:pBeg+lengthArrowElc], 1.e+4*prtclCoor[2,pBeg:pBeg+lengthArrowElc], \
              1.e+4*prtclCoor[4,pBeg:pBeg+lengthArrowElc],'-r',linewidth=4)
   ax120.plot(1.e+4*prtclCoor_2[0,pointsGCtot_2-turns:pointsGCtot_2],1.e+4*prtclCoor_2[2,pointsGCtot_2-turns:pointsGCtot_2], \
              1.e+4*prtclCoor_2[4,pointsGCtot_2-turns:pointsGCtot_2],'-b',linewidth=2)
   ax120.plot(1.e+4*prtclCoor_2[0,pointsGCtot_2-6:pointsGCtot_2-3],1.e+4*prtclCoor_2[2,pointsGCtot_2-6:pointsGCtot_2-3], \
              1.e+4*prtclCoor_2[4,pointsGCtot_2-6:pointsGCtot_2-3],'-b',linewidth=4)
   plt.xlabel('x, $\mu m$',color='m',fontsize=16)
   plt.ylabel('y, $\mu m$',color='m',fontsize=16)
   ax120.set_zlabel('z, $\mu m$',color='m',fontsize=16)
   plt.title(('First Electron Trajectory (End; $N_{Larm}=$%d):\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m' % \
              (larmorNumber[0],1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)),color='m',fontsize=16)
   ax120.text(0.,0.,25.,'Blue is Trajectory in the System of "Guiding Center"',color='b',fontsize=16)


   pEndArrw=pointsEndLarmor-pointsTurns/2
   pBegArrw=pEndArrw-50

#    print 'pointsGCtot_2=%d, turns=%d' % (pointsGCtot_2,turns)
#    for m in range(pointsGCtot_2-turns,pointsGCtot_2):
#       print 'Point %d: x=%e, y=%e, z=%e' % (m,1.e+7*prtclCoor_2[6,m],1.e+7*prtclCoor_2[8,m],1.e+7*prtclCoor_2[10,m])


   fig125=plt.figure(125)
   ax125=fig125.gca(projection='3d')
   ax125.plot(1.e+7*prtclCoor[6,pBeg:pointsEndLarmor],1.e+7*prtclCoor[8,pBeg:pointsEndLarmor], \
              1.e+7*prtclCoor[10,pBeg:pointsEndLarmor],'-r',linewidth=2)
   ax125.plot(1.e+7*prtclCoor[6,pBegArrw:pEndArrw],1.e+7*prtclCoor[8,pBegArrw:pEndArrw], \
              1.e+7*prtclCoor[10,pBegArrw:pEndArrw],'-r',linewidth=4)
   ax125.plot(1.e+7*prtclCoor_2[6,pointsGCtot_2-turns:pointsGCtot_2],1.e+7*prtclCoor_2[8,pointsGCtot_2-turns:pointsGCtot_2], \
              1.e+7*prtclCoor_2[10,pointsGCtot_2-turns:pointsGCtot_2],'-b',linewidth=2)
   ax125.plot(1.e+7*prtclCoor_2[6,pointsGCtot_2-turns+3:pointsGCtot_2-3],1.e+7*prtclCoor_2[8,pointsGCtot_2+3-turns:pointsGCtot_2-3], \
              1.e+7*prtclCoor_2[10,pointsGCtot_2-turns+3:pointsGCtot_2-3],'-b',linewidth=4)
   plt.xlabel('x, $nm$',color='m',fontsize=16)
   plt.ylabel('y, $nm$',color='m',fontsize=16)
   ax125.set_zlabel('z, $nm$',color='m',fontsize=16)
   plt.title('Approach-2: First Ion Trajectory (End)',color='m',fontsize=16)

if runFlagApproach_2 == 1:
   pointCoor_2=np.zeros(pointsGCtot_2)
   xCoor_2=np.zeros(pointsGCtot_2)
   yCoor_2=np.zeros(pointsGCtot_2)
   zCoor_2=np.zeros(pointsGCtot_2)
   stepNzCoor=50

   k=0
   for m in range(0,pointsGCtot_2,stepNzCoor):
      pointCoor_2[k]=k*stepNzCoor*stepsNumberOnGyro
      xCoor_2[k]=1.e4*prtclCoor_2[0,m]                             # mkm
      yCoor_2[k]=1.e4*prtclCoor_2[2,m]                             # mkm
      zCoor_2[k]=1.e4*prtclCoor_2[4,m]                             # mkm
      k += 1
      NzCoor_2=k

if runFlagApproach_3 == 1:
   pointsGCtot_3=int(pointTrack_3[0])
   pointCoor_3=np.zeros(pointsGCtot_3)
   xCoor_3=np.zeros(pointsGCtot_3)                             
   yCoor_3=np.zeros(pointsGCtot_3)
   zCoor_3=np.zeros(pointsGCtot_3)
   stepNzCoor=50

   k=0
   for m in range(0,pointsGCtot_3,stepNzCoor):
      pointCoor_3[k]=k*stepNzCoor*stepsNumberOnGyro
      xCoor_3[k]=1.e4*prtclCoor_3[0,m]                             # mkm
      yCoor_3[k]=1.e4*prtclCoor_3[2,m]                             # mkm
      zCoor_3[k]=1.e4*prtclCoor_3[4,m]                             # mkm
      k += 1
      NzCoor_3=k

if runFlagApproach_3 == 1 and plotFlagTracks == 1:
   pointsGCturns=turns                                             # Number of points for drawing ("GC")
   pointsDraw=int(pointTrack_3[0])                                 # Number of points for drawing

   plt.figure(230)
   plt.plot(range(int(pointsDraw)),1.e4*prtclCoor_3[4,0:pointsDraw],'-xr',linewidth=2)
   plt.xlabel('Points',color='m',fontsize=16)
   plt.ylabel('z, $\mu$m',color='m',fontsize=16)
   plt.title(('Approach-3: First Electron Trajectory ($N_{Larm}=$%d):\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m' % \
              (larmorNumber[0],1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)),color='m',fontsize=16)
   plt.grid(True)

   plt.figure(240)
   plt.plot(1.e4*prtclCoor_3[4,0:pointsDraw],1.e4*prtclCoor_3[12,0:pointsDraw],'-xr',linewidth=2)
   plt.xlabel('z, $\mu$m',color='m',fontsize=16)
   plt.ylabel('Distance Between Particles, $\mu$m',color='m',fontsize=16)
   plt.title(('Approach-3: First Electron Trajectory ($N_{Larm}=$%d):\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m' % \
              (larmorNumber[0],1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)),color='m',fontsize=16)
   plt.grid(True)

   pointsDraw=int(pointTrack_3[trackNumbMaxAbsDpxApprch_3])        # Number of points for drawing

   plt.figure(1230)
   plt.plot(range(int(pointsDraw)),1.e4*prtclCoorMaxAbsDpx_3[4,0:pointsDraw],'-xr',linewidth=2)
   plt.xlabel('Points',color='m',fontsize=16)
   plt.ylabel('z, $\mu$m',color='m',fontsize=16)
   plt.title(('Approach-3: Electron Track %d ($N_{Larm}=$%d):\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m' % \
              (trackNumbMaxAbsDpxApprch_3,larmorNumber[trackNumbMaxAbsDpxApprch_3],1.e+4*rhoMaxAbsDpxTurn_3, \
	       1.e+4*rhoLarmorMaxAbsDpxTurn_3)),color='m',fontsize=16)
   plt.grid(True)

   plt.figure(1240)
   plt.plot(1.e4*prtclCoorMaxAbsDpx_3[4,0:pointsDraw],1.e4*prtclCoorMaxAbsDpx_3[12,0:pointsDraw],'-xr',linewidth=2)
   plt.xlabel('z, $\mu$m',color='m',fontsize=16)
   plt.ylabel('Distance Between Particles, $\mu$m',color='m',fontsize=16)
   plt.title(('Approach-3: Track %d ($N_{Larm}=$%d):\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m' % \
              (trackNumbMaxAbsDpxApprch_3,larmorNumber[trackNumbMaxAbsDpxApprch_3],\
	       1.e+4*rhoMaxAbsDpxTurn_3,1.e+4*rhoLarmorMaxAbsDpxTurn_3)), color='m',fontsize=16)
   plt.grid(True)

if runFlagApproach_1 == 1 and runFlagApproach_2 == 1 and plotFlagTracks == 1:
   plt.figure(135)
   plt.plot(range(int(pointsEndLarmor)),1.e4*prtclCoor[4,0:pointsEndLarmor],'-r',linewidth=2)
   plt.plot(pointCoor_2[0:NzCoor_2],zCoor_2[0:NzCoor_2],'xb',markersize=10)
   plt.xlim([-100,pointsEndLarmor+100])
   plt.ylim([-5600,250])
   plt.xlabel('Points',color='m',fontsize=16)
   plt.ylabel('z, $\mu$m',color='m',fontsize=16)
   plt.title(('First Electron Trajectory ($N_{Larm}=$%d):\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m' % \
           (larmorNumber[0],1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)),color='m',fontsize=16)
   plt.text(1000,-300,'Blue is Trajectory in the System of "Guiding Center"',color='b',fontsize=16)
   plt.grid(True)

if runFlagApproach_1 == 1 and runFlagApproach_3 == 1 and plotFlagTracks == 1:
   plt.figure(235)
   plt.plot(range(int(pointsEndLarmor)),1.e4*prtclCoor[4,0:pointsEndLarmor],'-r',linewidth=2)
   plt.plot(pointCoor_3[0:NzCoor_3],zCoor_3[0:NzCoor_3],'xb',markersize=10)
#   plt.xlim([-100,pointsEndLarmor+100])
#   plt.ylim([-5600,250])
   plt.xlabel('Points',color='m',fontsize=16)
   plt.ylabel('z, $\mu$m',color='m',fontsize=16)
   plt.title(('First Electron Trajectory ($N_{Larm}=$%d):\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m' % \
           (larmorNumber[0],1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)),color='m',fontsize=16)
   plt.text(1000,-300,'Blue is Trajectory in the System of "Guiding Center"',color='b',fontsize=16)
   plt.grid(True)

if runFlagApproach_1 == 1 and runFlagApproach_2 == 1 and runFlagApproach_3 == 1 and plotFlagTracks == 1:

   pointsDraw_1=int(pointTrack[0])                                 # Number of points for drawing (approach-1)
   pointsDraw_2=int(pointTrack_2[0])                               # Number of points for drawing (approach-2)
   pointsDraw_3=int(pointTrack_3[0])                               # Number of points for drawing (approach-3)

# First trajectory:

   plt.figure(2000)
   plt.plot(range(pointsDraw_1),1.e4*prtclCoor[4,0:pointsDraw_1],'-r',linewidth=2)
   plt.plot(pointCoor_2,zCoor_2,'xk',markersize=10)
   plt.plot(pointCoor_3,zCoor_3,'xm',markersize=10)
   plt.xlabel('Points',color='m',fontsize=16)
   plt.ylabel('z, $\mu$m',color='m',fontsize=16)
   plt.title(('First Electron Trajectory ($N_{Larm}=$%d):\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m' % \
              (larmorNumber[0],1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)),color='m',fontsize=16)
   plt.legend(['Approach-1','Approach-2','Approach-3'],loc='upper left',fontsize=16)
   plt.grid(True)

   pointsDraw_1max=int(pointTrack[trackNumbMaxAbsDpxApprch])       # Number of points for drawing
   pointsDraw_2max=int(pointTrack_2[trackNumbMaxAbsDpxApprch_2])   # Number of points for drawing
   pointsDraw_3max=int(pointTrack_3[trackNumbMaxAbsDpxApprch_3])   # Number of points for drawing
 
   print 'pointsDraw_1max=%d, pointsDraw_2max=%d, pointsDraw_3max=%d' % (pointsDraw_1max,pointsDraw_2max,pointsDraw_3max)  
   pointCoor_2max=np.zeros(pointsDraw_2max)
   for k in range(pointsDraw_2max):
      pointCoor_2max[k]=stepsNumberOnGyro*k
      
   pointCoor_3max=np.zeros(pointsDraw_3max)
   for k in range(pointsDraw_3max):
      pointCoor_3max[k]=stepsNumberOnGyro*k
      
   plt.figure(2100)
   plt.plot(range(pointsDraw_1max),1.e4*prtclCoorMaxAbsDpx[4,0:pointsDraw_1max],'-xr',linewidth=2)
   plt.plot(pointCoor_2max,1.e4*prtclCoorMaxAbsDpx_2[4,0:pointsDraw_2max],'xk',markersize=10)
   plt.plot(pointCoor_3max,1.e4*prtclCoorMaxAbsDpx_3[4,0:pointsDraw_3max],'xm',markersize=10)
   plt.xlabel('Points',color='m',fontsize=16)
   plt.ylabel('z, $\mu$m',color='m',fontsize=16)
   plt.title(('Electron Track %d ($N_{Larm}=$%d):\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m' % \
              (trackNumbMaxAbsDpxApprch_3,larmorNumber[trackNumbMaxAbsDpxApprch_3],1.e+4*rhoMaxAbsDpxTurn_3, \
	       1.e+4*rhoLarmorMaxAbsDpxTurn_3)),color='m',fontsize=16)
   plt.legend(['Approach-1','Approach-2','Approach-3'],loc='upper left',fontsize=16)
   plt.grid(True)

   fig3000=plt.figure(3000)
   plt.plot(1.e+4*prtclCoor_2[4,0:b_2LenFirstTrack],1.e+4*diff_b2[0:b_2LenFirstTrack],'.r',linewidth=2)
   plt.plot(1.e+4*prtclCoor_3[4,0:b_3LenFirstTrack],1.e+4*diff_b3[0:b_3LenFirstTrack],'.b',linewidth=2)
   plt.xlabel('z, $\mu m$',color='m',fontsize=16)
   plt.ylabel('$\Delta b$, $\mu m$',color='m',fontsize=16)
   plt.title('First Trajectory: $\Delta b$ (Difference for Distances $b$ between Particles)',color='m',fontsize=16)
   plt.legend(['$b_{Apprch-1}-b_{Apprch-2}$','$b_{Apprch-1}-b_{Apprch-3}$'],loc='upper left',fontsize=16)
#   plt.xlim([1.e+4*prtclCoor_3[4,0]-50,1.e+4*prtclCoor_3[4,b_3LenFirstTrack-1]+50])
   plt.ylim([-.5,.5])
   plt.grid(True)

#
# Normalization of the arrays zApprch2dpx,zApprch2dpy,zApprch2dpz
# to improve visual presentation of the results:
#
if runFlagApproach_2 == 1 and plotFlagDpTransf == 1:
   visPrezApprch2dpx=np.zeros((nBins,nBins))
   visPrezApprch2dpy=np.zeros((nBins,nBins))
   visPrezApprch2dpz=np.zeros((nBins,nBins))

   dpxDataMax_2=-1.e24                                             # maximum for array of zApprch2dpx
   indxAdpxDataMax_2=0                                             # index A_2 of the maximum for array of zApprch2dpx
   indxBdpxDataMax_2=0                                             # index B_2 of the maximum for array of zApprch2dpx
   dpyDataMax_2=-1.e24                                             # maximum for array of zApprch2dpx
   indxAdpyDataMax_2=0                                             # index A_2 of the maximum for array of zApprch2dpx
   indxBdpyDataMax_2=0                                             # index B_2 of the maximum for array of zApprch2dpx
   dpzDataMax_2=-1.e24                                             # maximum for array of zApprch2dpx
   indxAdpzDataMax_2=0                                             # index A_2 of the maximum for array of zApprch2dpx
   indxBdpzDataMax_2=0                                             # index B_2 of the maximum for array of zApprch2dpx

   normFactor=1.e24
   for iA in range(nBins):
      for iB in range(nBins):
         visPrezApprch2dpx[iA,iB]=normFactor*zApprch2dpx[iA,iB]
         visPrezApprch2dpy[iA,iB]=normFactor*zApprch2dpy[iA,iB]
         visPrezApprch2dpz[iA,iB]=normFactor*zApprch2dpz[iA,iB]
         if dpxDataMax_2 < abs(normFactor*zApprch2dpx[iA,iB]):
            dpxDataMax_2=abs(normFactor*zApprch2dpx[iA,iB])
            indxAdpxDataMax_2=iA
            indxBdpxDataMax_2=iB
         if dpyDataMax_2 < abs(normFactor*zApprch2dpy[iA,iB]):
            dpyDataMax_2=abs(normFactor*zApprch2dpy[iA,iB])
            indxAdpyDataMax_2=iA
            indxBdpyDataMax_2=iB
         if dpzDataMax_2 < abs(normFactor*zApprch2dpz[iA,iB]):
            dpzDataMax_2=abs(normFactor*zApprch2dpz[iA,iB])
            indxAdpzDataMax_2=iA
            indxBdpzDataMax_2=iB

   print '|Maximum| (approach-2): dpx=%e,  dpy=%e,  dpz=%e' % (dpxDataMax_2,dpyDataMax_2,dpzDataMax_2)

   X2,Y2=np.meshgrid(xA_2,yB_2) 

   fig170=plt.figure(170)
   ax170=fig170.gca(projection='3d')
#    surf=ax170.plot_surface(X2,Y2,zApprch2dpx,cmap=cm.coolwarm,linewidth=0,antialiased=False)
   surf=ax170.plot_surface(X2,Y2,visPrezApprch2dpx,cmap=cm.jet,linewidth=0,antialiased=False)
   plt.title(('Approach-2: Transfered Momentum $dP_x$\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{24}$ $g \cdot cm/sec$)' % \
              (lastTrackNumber_2,dpxDataMax_2)), color='m',fontsize=16)
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   ax170.set_zlabel('$dP_x \cdot 10^{-24}$; $g \cdot cm/sec$',color='m',fontsize=16)
   cb = fig170.colorbar(surf)
   plt.grid(True)

   fig175=plt.figure(175)
   ax=fig175.add_subplot(111)                                      # for contours plotting
#    mapDpx=ax.contourf(X2,Y2,zApprch2dpx)   
   mapDpx=ax.contourf(X2,Y2,visPrezApprch2dpx)   
#    mapDpx=ax.contourf(X2,Y2,dpxApprch_2,levels)                  # Next 4 commands not work  
#    Contourrange=[int(NlarmCutofDown+1)]
#    mapTurnCutoff=ax.contour(X2,Y2,Nlarm,Contourrange,format='%d',colors=('w'),linewidths=(2)) 
#    plt.clabel(mapTurnCutoff,inline=1,fontsize=14,manual=[(-3,-1.5)])  
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   plt.title(('Approach-2: Transfered Momentum $dP_x$ $(\cdot 10^{-24}$; $g \cdot cm/sec)$\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-24}$)' % \
              (lastTrackNumber_2,dpxDataMax_2)), color='m',fontsize=14)
   fig175.colorbar(mapDpx)


   fig180=plt.figure(180)
   ax180=fig180.gca(projection='3d')
#    surf=ax180.plot_surface(X2,Y2,visPrezApprch2dpx,cmap=cm.coolwarm,linewidth=0,antialiased=False)
   surf=ax180.plot_surface(X2,Y2,visPrezApprch2dpy,cmap=cm.jet,linewidth=0,antialiased=False)
   plt.title(('Approach-2: Transfered Momentum $dP_y$\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{24}$ $g \cdot cm/sec$)' % \
              (lastTrackNumber_2,dpyDataMax_2)), color='m',fontsize=16)
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   ax180.set_zlabel('$dP_y \cdot 10^{-24}$; $g \cdot cm/sec$',color='m',fontsize=16)
   cb = fig180.colorbar(surf)
#    cbar=fig180.colorbar(surf,ticks=[0,1000,2000,3000,4000,5000,6000,7000])  # Next 2 commands not work
#    cbar.ax.set_yticklabels(['0','1000','2000','3000','4000','5000','6000','7000'])
#    labels=np.arange(0,8000,1000)               # Next 4 commands not work
#    location=labels
#    cb.set_ticks(location)
#    cb.set_ticklabels(labels)
#    tick_locator = ticker.MaxNLocator(nbins=10) # Next 3 commands not work
#    cb.locator = tick_locator
#    cb.update_ticks()
   plt.grid(True)

   fig185=plt.figure(185)
   ax=fig185.add_subplot(111)                                      # for contours plotting
#    mapDpx=ax.contourf(X2,Y2,zApprch2dpy)   
   mapDpy=ax.contourf(X2,Y2,visPrezApprch2dpy)   
#    mapDpy=ax.contourf(X2,Y2,dpxApprch_2,levels)                  # Next 4 commands not work  
#    Contourrange=[int(NlarmCutofDown+1)]
#    mapTurnCutoff=ax.contour(X2,Y2,Nlarm,Contourrange,format='%d',colors=('w'),linewidths=(2)) 
#    plt.clabel(mapTurnCutoff,inline=1,fontsize=14,manual=[(-3,-1.5)])  
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   plt.title(('Approach-2: Transfered Momentum $dP_y$ $(\cdot 10^{-24}$; $g \cdot cm/sec)$\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-24}$)' % \
              (lastTrackNumber_2,dpyDataMax_2)), color='m',fontsize=14)
   cb = fig185.colorbar(mapDpy)


   fig190=plt.figure(190)
   ax190=fig190.gca(projection='3d')
#    surf=ax190.plot_surface(X2,Y2,visPrezApprch2dpx,cmap=cm.coolwarm,linewidth=0,antialiased=False)
   surf=ax190.plot_surface(X2,Y2,visPrezApprch2dpz,cmap=cm.jet,linewidth=0,antialiased=False)
   plt.title(('Approach-2: Transfered Momentum $dP_z$\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{24}$ $g \cdot cm/sec$)' % \
              (lastTrackNumber_2,dpzDataMax_2)), color='m',fontsize=16)
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   ax190.set_zlabel('$dP_z \cdot 10^{-24}$; $g \cdot cm/sec$',color='m',fontsize=16)
   cb = fig190.colorbar(surf)
   plt.grid(True)

   fig195=plt.figure(195)
   ax=fig195.add_subplot(111)                                      # for contours plotting
#    mapDpx=ax.contourf(X2,Y2,zApprch2dpz)   
   mapDpz=ax.contourf(X2,Y2,visPrezApprch2dpz)   
#    mapDpx=ax.contourf(X2,Y2,dpxApprch_2,levels)                  # Next 4 commands not work  
#    Contourrange=[int(NlarmCutofDown+1)]
#    mapTurnCutoff=ax.contour(X2,Y2,Nlarm,Contourrange,format='%d',colors=('w'),linewidths=(2)) 
#    plt.clabel(mapTurnCutoff,inline=1,fontsize=14,manual=[(-3,-1.5)])  
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   plt.title(('Approach-2: Transfered Momentum $dP_z$ $(\cdot 10^{-24}$; $g \cdot cm/sec)$\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-24}$)' % \
              (lastTrackNumber_2,dpzDataMax_2)), color='m',fontsize=14)
   fig195.colorbar(mapDpz)


if runFlagApproach_3 == 1:
   pointsGCturns=turns                                             # Number of points for drawing ("GC")

   pBegArrw=pointsGCturns/2-1

if runFlagApproach_1 == 1 and runFlagApproach_3 == 1 and plotFlagTracks == 1:
   fig210=plt.figure(210)
   ax210=fig210.gca(projection='3d')
   ax210.plot(1.e+4*prtclCoor[0,0:pointsStartLarmor],1.e+4*prtclCoor[2,0:pointsStartLarmor], \
              1.e+4*prtclCoor[4,0:pointsStartLarmor],'-r',linewidth=2)
   ax210.plot(1.e+4*prtclCoor[0,0:lengthArrowElc],1.e+4*prtclCoor[2,0:lengthArrowElc], \
              1.e+4*prtclCoor[4,0:lengthArrowElc],'-b',linewidth=2)
   ax210.plot(1.e+4*prtclCoor[0,pointsStartLarmor-lengthArrowElc:pointsStartLarmor], \
              1.e+4*prtclCoor[2,pointsStartLarmor-lengthArrowElc:pointsStartLarmor], \
              1.e+4*prtclCoor[4,pointsStartLarmor-lengthArrowElc:pointsStartLarmor],'-b',linewidth=2)
   ax210.plot(1.e+4*prtclCoor_3[0,0:pointsGCturns],1.e+7*prtclCoor_3[2,0:pointsGCturns], \
              1.e+4*prtclCoor_3[4,0:pointsGCturns],'-b',linewidth=2)
   ax210.plot(1.e+4*prtclCoor_3[0,pBegArrw:pBegArrw+3],1.e+4*prtclCoor_3[2,pBegArrw:pBegArrw+3], \
              1.e+4*prtclCoor_3[4,pBegArrw:pBegArrw+3],'-b',linewidth=4)
   plt.xlabel('x, $\mu m$',color='m',fontsize=16)
   plt.ylabel('y, ($\mu m$',color='m',fontsize=16)
   ax210.set_zlabel('z, $\mu m$',color='m',fontsize=16)
   plt.title(('First Electron Trajectory (Start; $N_{Larm}=$%d):\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m' \
           % (larmorNumber[0],1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)),color='m',fontsize=16)
   ax210.text(0.,0.,-5425.,'Blue is Trajectory in the System of "Guiding Center"',color='b',fontsize=16)


#    for m in range(pointsEndLarmor-50,pointsEndLarmor):
#       print 'Electron: z(%d)=%e' % (m,1.e+4*prtclCoor[4,m])
#    for m in range(pointsGCtot_3-turns,pointsGCtot_3):
#       print 'Electron_gc: z(%d)=%e' % (m,1.e+4*prtclCoor_3[4,m])

   print 'Last point of the first track: z(%d)=%e mkm, z_gc(%d)=%e mkm' % \
         (pointsEndLarmor,1.e+4*prtclCoor[4,pointsEndLarmor-1],pointsGCtot_3,1.e+4*prtclCoor_3[4,pointsGCtot_3-1])

   fig220=plt.figure(220) 
   ax220=fig220.gca(projection='3d')
   ax220.plot(1.e+4*prtclCoor[0,pBeg:pointsEndLarmor], \
              1.e+4*prtclCoor[2,pBeg:pointsEndLarmor], \
              1.e+4*prtclCoor[4,pBeg:pointsEndLarmor],'-r',linewidth=2)
   ax220.plot(1.e+4*prtclCoor[0,pointsEndLarmor-lengthArrowElc:pointsEndLarmor], \
              1.e+4*prtclCoor[2,pointsEndLarmor-lengthArrowElc:pointsEndLarmor], \
              1.e+4*prtclCoor[4,pointsEndLarmor-lengthArrowElc:pointsEndLarmor],'-r',linewidth=4)
   ax220.plot(1.e+4*prtclCoor[0,pBeg:pBeg+lengthArrowElc], 1.e+4*prtclCoor[2,pBeg:pBeg+lengthArrowElc], \
              1.e+4*prtclCoor[4,pBeg:pBeg+lengthArrowElc],'-r',linewidth=4)
   ax220.plot(1.e+4*prtclCoor_3[0,pointsGCtot_3-turns:pointsGCtot_3],1.e+4*prtclCoor_3[2,pointsGCtot_3-turns:pointsGCtot_3], \
              1.e+4*prtclCoor_3[4,pointsGCtot_3-turns:pointsGCtot_3],'-b',linewidth=2)
   ax220.plot(1.e+4*prtclCoor_3[0,pointsGCtot_3-6:pointsGCtot_3-3],1.e+4*prtclCoor_3[2,pointsGCtot_3-6:pointsGCtot_3-3], \
              1.e+4*prtclCoor_3[4,pointsGCtot_3-6:pointsGCtot_3-3],'-b',linewidth=4)
   plt.xlabel('x, $\mu m$',color='m',fontsize=16)
   plt.ylabel('y, $\mu m$',color='m',fontsize=16)
   ax220.set_zlabel('z, $\mu m$',color='m',fontsize=16)
   plt.title(('First Electron Trajectory (End; $N_{Larm}=$%d):\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m' % \
              (larmorNumber[0],1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)),color='m',fontsize=16)
   ax220.text(0.,0.,25.,'Blue is Trajectory in the System of "Guiding Center"',color='b',fontsize=16)


   pEndArrw=pointsEndLarmor-pointsTurns/2
   pBegArrw=pEndArrw-50

#    print 'pointsGCtot_3=%d, turns=%d' % (pointsGCtot_3,turns)
#    for m in range(pointsGCtot_3-turns,pointsGCtot_3):
#       print 'Point %d: x=%e, y=%e, z=%e' % (m,1.e+7*prtclCoor_3[6,m],1.e+7*prtclCoor_3[8,m],1.e+7*prtclCoor_3[10,m])


   fig225=plt.figure(225)
   ax225=fig225.gca(projection='3d')
   ax225.plot(1.e+7*prtclCoor[6,pBeg:pointsEndLarmor],1.e+7*prtclCoor[8,pBeg:pointsEndLarmor], \
              1.e+7*prtclCoor[10,pBeg:pointsEndLarmor],'-r',linewidth=2)
   ax225.plot(1.e+7*prtclCoor[6,pBegArrw:pEndArrw],1.e+7*prtclCoor[8,pBegArrw:pEndArrw], \
              1.e+7*prtclCoor[10,pBegArrw:pEndArrw],'-r',linewidth=4)
   ax225.plot(1.e+7*prtclCoor_3[6,pointsGCtot_3-turns:pointsGCtot_3],1.e+7*prtclCoor_3[8,pointsGCtot_3-turns:pointsGCtot_3], \
              1.e+7*prtclCoor_3[10,pointsGCtot_3-turns:pointsGCtot_3],'xb',linewidth=2)
   plt.xlabel('x, $nm$',color='m',fontsize=16)
   plt.ylabel('y, $nm$',color='m',fontsize=16)
   ax225.set_zlabel('z, $nm$',color='m',fontsize=16)
   plt.title('First Ion Trajectory (End)',color='m',fontsize=16)

if runFlagApproach_3 == 1:
   pointCoor_3=np.zeros(pointsGCtot_3)
   zCoor_3=np.zeros(pointsGCtot_3)
   stepNzCoor=50

   k=0
   for m in range(0,pointsGCtot_3,stepNzCoor):
      pointCoor_3[k]=k*stepNzCoor*stepsNumberOnGyro
      zCoor_3[k]=1.e4*prtclCoor_3[4,m]
      k += 1
      NzCoor_3=k

if runFlagApproach_2 == 1 and plotFlagDpTransf == 1:
   plt.figure(155)
   plt.plot(uPot_enrgKin_2,larmR_b_2,'.r')
   plt.xlim([minUpot_enrgKin_2-.1,maxUpot_enrgKin_2+.1])
   plt.ylim([minLarmR_b_2-.1,maxLarmR_b_2+.1])
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   plt.title(('Approach-2: Map for Transfered Momenta Calculations\nTracks: %d' % lastTrackNumber_2), color='m',fontsize=16)
   plt.grid(True)

if runFlagApproach_3 == 1 and plotFlagDpTransf == 1:
   plt.figure(255)
   plt.plot(uPot_enrgKin_3,larmR_b_3,'.r')
   plt.xlim([minUpot_enrgKin_3-.1,maxUpot_enrgKin_3+.1])
   plt.ylim([minLarmR_b_3-.1,maxLarmR_b_3+.1])
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   plt.title(('Approach-3: Map for Transfered Momenta Calculations\nTracks: %d' % lastTrackNumber_3), color='m',fontsize=16)
   plt.grid(True)

if runFlagApproach_1 == 1 and runFlagApproach_3 == 1 and plotFlagTracks == 1:
   fig425=plt.figure(425)
   plt.plot(1.e+4*prtclCoor_3[4,0:b_3LenFirstTrack],1.e+4*diff_b3[0:b_3LenFirstTrack],'.r',linewidth=2)
   plt.xlabel('z, $\mu m$',color='m',fontsize=16)
   plt.ylabel('$\Delta b=b_{Apprch-1}-b_{Apprch-3}$, $\mu m$',color='m',fontsize=16)
   plt.title('First Trajectory: $\Delta b$ (Difference for Distances $b$ between Particles)',color='m',fontsize=16)
   plt.xlim([1.e+4*prtclCoor_3[4,0]-50,1.e+4*prtclCoor_3[4,b_3LenFirstTrack-1]+50])
   plt.ylim([-.5,.5])
   plt.grid(True)

if runFlagApproach_2 == 1 and runFlagApproach_3 == 1 and plotFlagTracks == 1:
   plt.figure(430)
   plt.plot(pointCoor_2[0:NzCoor_2],xCoor_2[0:NzCoor_2],'-r',markersize=10)
   plt.plot(pointCoor_3[0:NzCoor_3],xCoor_3[0:NzCoor_3],'xb',markersize=10)
#   plt.xlim([-100,pointsEndLarmor+100])
#   plt.ylim([-5600,250])
   plt.xlabel('Points',color='m',fontsize=16)
   plt.ylabel('x, $\mu$m',color='m',fontsize=16)
   plt.title(('First Electron Trajectory ($N_{Larm}=$%d):\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m' % \
           (larmorNumber[0],1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)),color='m',fontsize=16)
#   plt.text(1000,-300,'Blue is Trajectory in the System of "Guiding Center"',color='b',fontsize=16)
   plt.grid(True)

   plt.figure(440)
   plt.plot(pointCoor_2[0:NzCoor_2],yCoor_2[0:NzCoor_2],'-r',markersize=10)
   plt.plot(pointCoor_3[0:NzCoor_3],yCoor_3[0:NzCoor_3],'xb',markersize=10)
#   plt.xlim([-100,pointsEndLarmor+100])
   plt.ylim([-.000005,.000025])
   plt.xlabel('Points',color='m',fontsize=16)
   plt.ylabel('y, $\mu$m',color='m',fontsize=16)
   plt.title(('First Electron Trajectory ($N_{Larm}=$%d):\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m' % \
           (larmorNumber[0],1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)),color='m',fontsize=16)
#   plt.text(1000,-300,'Blue is Trajectory in the System of "Guiding Center"',color='b',fontsize=16)
   plt.grid(True)

   plt.figure(450)
   plt.plot(pointCoor_2[0:NzCoor_2],zCoor_2[0:NzCoor_2],'-r',markersize=10)
   plt.plot(pointCoor_3[0:NzCoor_3],zCoor_3[0:NzCoor_3],'xb',markersize=10)
#   plt.xlim([-100,pointsEndLarmor+100])
#   plt.ylim([-5600,250])
   plt.xlabel('Points',color='m',fontsize=16)
   plt.ylabel('z, $\mu$m',color='m',fontsize=16)
   plt.title(('First Electron Trajectory ($N_{Larm}=$%d):\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m' % \
           (larmorNumber[0],1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)),color='m',fontsize=16)
#   plt.text(1000,-300,'Blue is Trajectory in the System of "Guiding Center"',color='b',fontsize=16)
   plt.grid(True)

   fig460=plt.figure(460)
   plt.plot(1.e+4*prtclCoor_2[4,0:pointTrack_2[0]],1.e+30*prtclCoor_2[13,0:pointTrack_2[0]]-2.384782e+6,'.r',linewidth=2)
   plt.plot(1.e+4*prtclCoor_3[4,0:pointTrack_3[0]],1.e+30*prtclCoor_3[13,0:pointTrack_3[0]]-2.384782e+6,'xb',linewidth=2)
   plt.xlabel('z, $\mu m$',color='m',fontsize=16)
   plt.ylabel('$J$, $g \cdot cm^2/sec$',color='m',fontsize=16)
   plt.title('First Electron Trajectory: Action $J$',color='m',fontsize=16)
   plt.grid(True)

if runFlagApproach_3 == 1 and plotFlagTracks == 1:
   plt.figure(500)
   plt.plot(1.e4*prtclCoor_3[4,0:pointTrack_3[0]],1.e4*prtclCoor_3[13,0:pointTrack_3[0]],'.r',linewidth=2)
#   plt.xlim([-100,pointsEndLarmor+100])
#   plt.ylim([-5600,250])
   plt.ylabel('$dy_{gc}$, $\mu$m',color='m',fontsize=16)
   plt.xlabel('z, $\mu$m',color='m',fontsize=16)
#   plt.title(('First Electron Trajectory ($N_{Larm}=$%d):\nImpact Parameter=%5.2f $\mu$m, $R_{Larm}$=%5.2f $\mu$m' % \
#           (larmorNumber[0],1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)),color='m',fontsize=16)
#   plt.text(1000,-300,'Blue is Trajectory in the System of "Guiding Center"',color='b',fontsize=16)
   plt.grid(True)

#
# Normalization of the arrays zApprch3dpx,zApprch3dpy,zApprch3dpz
# to improve visual presentation of the results:
#
if runFlagApproach_3 == 1 and plotFlagDpTransf == 1:
   visPrezApprch3dpx=np.zeros((nBins,nBins))
   visPrezApprch3dpy=np.zeros((nBins,nBins))
   visPrezApprch3dpz=np.zeros((nBins,nBins))

   dpxDataMax_3=-1.e24                                             # maximum for array of zApprch3dpx
   indxAdpxDataMax_3=0                                             # index A_3 of the maximum for array of zApprch3dpx
   indxBdpxDataMax_3=0                                             # index B_3 of the maximum for array of zApprch3dpx
   dpyDataMax_3=-1.e24                                             # maximum for array of zApprch3dpx
   indxAdpyDataMax_3=0                                             # index A_3 of the maximum for array of zApprch3dpx
   indxBdpyDataMax_3=0                                             # index B_3 of the maximum for array of zApprch3dpx
   dpzDataMax_3=-1.e24                                             # maximum for array of zApprch3dpx
   indxAdpzDataMax_3=0                                             # index A_3 of the maximum for array of zApprch3dpx
   indxBdpzDataMax_3=0                                             # index B_3 of the maximum for array of zApprch3dpx

   normFactor=1.e24
   for iA in range(nBins):
      for iB in range(nBins):
         visPrezApprch3dpx[iA,iB]=normFactor*zApprch3dpx[iA,iB]
         visPrezApprch3dpy[iA,iB]=normFactor*zApprch3dpy[iA,iB]
         visPrezApprch3dpz[iA,iB]=normFactor*zApprch3dpz[iA,iB]
         if dpxDataMax_3 < abs(normFactor*zApprch3dpx[iA,iB]):
            dpxDataMax_3=abs(normFactor*zApprch3dpx[iA,iB])
            indxAdpxDataMax_3=iA
            indxBdpxDataMax_3=iB
         if dpyDataMax_3 < abs(normFactor*zApprch3dpy[iA,iB]):
            dpyDataMax_3=abs(normFactor*zApprch3dpy[iA,iB])
            indxAdpyDataMax_3=iA
            indxBdpyDataMax_3=iB
#	 print 'zApprch3dpz[%d,%d]=%e, max=%e' % (iA,iB,normFactor*zApprch3dpz[iA,iB],dpzDataMax_3)   
         if dpzDataMax_3 < abs(normFactor*zApprch3dpz[iA,iB]):
            dpzDataMax_3=abs(normFactor*zApprch3dpz[iA,iB])
            indxAdpzDataMax_3=iA
            indxBdpzDataMax_3=iB
   print '|Maximum| (approach-3): dpx=%e,  dpy=%e,  dpz=%e' % (dpxDataMax_3,dpyDataMax_3,dpzDataMax_3)


   X3,Y3=np.meshgrid(xA_3,yB_3) 

   fig270=plt.figure(270)
   ax270=fig270.gca(projection='3d')
#    surf=ax270.plot_surface(X3,Y3,zApprch2dpx,cmap=cm.coolwarm,linewidth=0,antialiased=False)
   surf=ax270.plot_surface(X3,Y3,visPrezApprch3dpx,cmap=cm.jet,linewidth=0,antialiased=False)
   plt.title(('Approach-3: Transfered Momentum $dP_x$\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{24}$ $g \cdot cm/sec$)' % \
              (lastTrackNumber_3,dpxDataMax_3)), color='m',fontsize=16)
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   ax270.set_zlabel('$dP_x \cdot 10^{-24}$; $g \cdot cm/sec$',color='m',fontsize=16)
   cb = fig270.colorbar(surf)
   plt.grid(True)

   fig275=plt.figure(275)
   ax=fig275.add_subplot(111)                                      # for contours plotting
#    mapDpx=ax.contourf(X3,Y3,zApprch3dpx)   
   mapDpx=ax.contourf(X3,Y3,visPrezApprch3dpx)   
#    mapDpx=ax.contourf(X3,Y3,dpxApprch_3,levels)                  # Next 4 commands not work  
#    Contourrange=[int(NlarmCutofDown+1)]
#    mapTurnCutoff=ax.contour(X3,Y3,Nlarm,Contourrange,format='%d',colors=('w'),linewidths=(2)) 
#    plt.clabel(mapTurnCutoff,inline=1,fontsize=14,manual=[(-3,-1.5)])  
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   plt.title(('Approach-3: Transfered Momentum $dP_x$ $(\cdot 10^{-24}$; $g \cdot cm/sec)$\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-24}$)' % \
              (lastTrackNumber_3,dpxDataMax_3)), color='m',fontsize=14)
   fig275.colorbar(mapDpx)


   fig280=plt.figure(280)
   ax280=fig280.gca(projection='3d')
#    surf=ax280.plot_surface(X3,Y3,visPrezApprch3dpx,cmap=cm.coolwarm,linewidth=0,antialiased=False)
   surf=ax280.plot_surface(X3,Y3,visPrezApprch3dpy,cmap=cm.jet,linewidth=0,antialiased=False)
   plt.title(('Approach-3: Transfered Momentum $dP_y$\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{24}$ $g \cdot cm/sec$)' % \
              (lastTrackNumber_3,dpyDataMax_3)), color='m',fontsize=16)
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   ax280.set_zlabel('$dP_y \cdot 10^{-24}$; $g \cdot cm/sec$',color='m',fontsize=16)
   cb = fig280.colorbar(surf)
#    cbar=fig280.colorbar(surf,ticks=[0,1000,2000,3000,4000,5000,6000,7000])  # Next 2 commands not work
#    cbar.ax.set_yticklabels(['0','1000','2000','3000','4000','5000','6000','7000'])
#    labels=np.arange(0,8000,1000)                                 # Next 4 commands not work
#    location=labels
#    cb.set_ticks(location)
#    cb.set_ticklabels(labels)
#    tick_locator = ticker.MaxNLocator(nbins=10) # Next 3 commands not work
#    cb.locator = tick_locator
#       cb.update_ticks()
   plt.grid(True)

   fig285=plt.figure(285)
   ax=fig285.add_subplot(111)                                      # for contours plotting
#    mapDpx=ax.contourf(X3,Y3,zApprch3dpy)   
   mapDpy=ax.contourf(X3,Y3,visPrezApprch3dpy)   
#    mapDpx=ax.contourf(X3,Y3,dpxApprch_3,levels)                  # Next 4 commands not work  
#    Contourrange=[int(NlarmCutofDown+1)]
#    mapTurnCutoff=ax.contour(X3,Y3,Nlarm,Contourrange,format='%d',colors=('w'),linewidths=(2)) 
#    plt.clabel(mapTurnCutoff,inline=1,fontsize=14,manual=[(-3,-1.5)])  
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   plt.title(('Approach-3: Transfered Momentum $dP_y$ $(\cdot 10^{-24}$; $g \cdot cm/sec)$\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-24}$)' % \
              (lastTrackNumber_3,dpyDataMax_3)), color='m',fontsize=14)
   fig285.colorbar(mapDpy)


   fig290=plt.figure(290)
   ax290=fig290.gca(projection='3d')
#    surf=ax290.plot_surface(X3,Y3,visPrezApprch3dpx,cmap=cm.coolwarm,linewidth=0,antialiased=False)
   surf=ax290.plot_surface(X3,Y3,visPrezApprch3dpz,cmap=cm.jet,linewidth=0,antialiased=False)
   plt.title(('Approach-3: Transfered Momentum $dP_z$\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{24}$ $g \cdot cm/sec$)' % \
              (lastTrackNumber_3,dpzDataMax_3)), color='m',fontsize=16)
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   ax290.set_zlabel('$dP_z \cdot 10^{-24}$; $g \cdot cm/sec$',color='m',fontsize=16)
   cb = fig290.colorbar(surf)
   plt.grid(True)

   fig295=plt.figure(295)
   ax=fig295.add_subplot(111)                                      # for contours plotting
#    mapDpx=ax.contourf(X3,Y3,zApprch3dpz)   
   mapDpz=ax.contourf(X3,Y3,visPrezApprch3dpz)   
#    mapDpx=ax.contourf(X3,Y3,dpxApprch_3,levels)                  # Next 4 commands not work  
#    Contourrange=[int(NlarmCutofDown+1)]
#    mapTurnCutoff=ax.contour(X3,Y3,Nlarm,Contourrange,format='%d',colors=('w'),linewidths=(2)) 
#    plt.clabel(mapTurnCutoff,inline=1,fontsize=14,manual=[(-3,-1.5)])  
   plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
   plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
   plt.title(('Approach-3: Transfered Momentum $dP_z$ $(\cdot 10^{-24}$; $g \cdot cm/sec)$\nTracks: %d (|Maximum| = %5.1f $\cdot 10^{-24}$)' % \
              (lastTrackNumber_3,dpzDataMax_3)), color='m',fontsize=14)
   fig295.colorbar(mapDpz)

plt.show()   

sys.exit()

###########################################################
#
# Writing 6-vector for electron along first trajectory:
#
###########################################################
outputFile='elec6vector.dat'
print 'Open output file "%s"...' % outputFile
outfileFlag=0
try:
   outfile = open(outputFile,'w')
   outfileFlag=1
except:
   print 'Problem to open output file "%s"' % outputFile

outfile.write ('\n                              First trajectory of electron (total points %d)' % pointsEndLarmor)
outfile.write ('\n      x, cm            y, cm            z, cm         px, g*cm/sec     py, g*cm/sec     pz, g*cm/sec\n\n')

for it in range(pointsEndLarmor):
   strLine=        '  {: 12.6e}'.format(prtclCoor[0,it])+',   '+'{: 12.6e}'.format(prtclCoor[2,it])+ \
                   ',   '+'{: 12.6e}'.format(prtclCoor[4,it])+',   '
   strLine=strLine+'{: 12.6e}'.format(prtclCoor[1,it])+',   '+'{: 12.6e}'.format(prtclCoor[3,it])+ \
                   ',   '+'{: 12.6e}'.format(prtclCoor[5,it])
   outfile.write ('%s\n' % strLine)

outfile.close()
print 'Close the written output file "%s"' % outputFile

###########################################################
#
# Writing 6-vector for ion along first trajectory:
#
###########################################################
outputFile='ion6vector.dat'
print 'Open output file "%s"...' % outputFile
outfileFlag=0
try:
   outfile = open(outputFile,'w')
   outfileFlag=1
except:
   print 'Problem to open output file "%s"' % outputFile

outfile.write ('\n                                First trajectory of ion (total points %d)' % pointsEndLarmor)
outfile.write ('\n      x, cm            y, cm            z, cm         px, g*cm/sec     py, g*cm/sec     pz, g*cm/sec\n\n')

for it in range(pointsEndLarmor):
   strLine=        '  {: 12.6e}'.format(prtclCoor[6,it])+',   '+'{: 12.6e}'.format(prtclCoor[8,it])+',   '+'{: 12.6e}'.format(prtclCoor[10,it])+',   '
   strLine=strLine+'{: 12.6e}'.format(prtclCoor[7,it])+',   '+'{: 12.6e}'.format(prtclCoor[9,it])+',   '+'{: 12.6e}'.format(prtclCoor[11,it])
   outfile.write ('%s\n' % strLine)

outfile.close()
print 'Close the written output file "%s"' % outputFile

# sys.exit()   

###################################################
#
# Writing the results to output file for later processing and vizualization: 
#
###################################################   
#
# Opening the output file for approach-1: 
#
outputFile='dpApprch1.dat'
print 'Open output file "%s"...' % outputFile
outfileFlag=0
try:
   outfile = open(outputFile,'w')
   outfileFlag=1
except:
   print 'Problem to open output file "%s"' % outputFile

nInLine=10

outfile.write ('\n       Tracks: %d , Normalization Factor: %e\n' % (lastTrackNumber,normFactor))

#
# Writing the array xA=log10(Upot/Ekin) to output file: 
#
outfile.write ('\n       xA=log10(eVca/Ekin) ( Entries: %d with %d per line)\n\n' % (nBins,nInLine))

k=0
for m in range(nBins): 
   strVal='{:f}'.format(xA[m])
   if k == 0:
      xA_line=strVal
   else:
      xA_line=xA_line+', '+strVal
   k += 1
   if k == nInLine:
      outfile.write ('%s\n' % xA_line)
      k=0
if k != 0:
   outfile.write ('%s\n' % xA_line)

#
# Writing the array yB=log10(Rlarm/b) to output file: 
#
outfile.write ('\n       yB=log10(Rlarm/b) ( Entries: %d with %d per line)\n\n' % (nBins,nInLine))

k=0
for m in range(nBins): 
   strVal='{:f}'.format(yB[m])
   if k == 0:
      yB_line=strVal
   else:
      yB_line=yB_line+', '+strVal
   k += 1
   if k == nInLine:
      outfile.write ('%s\n' % yB_line)
      k=0
if k != 0:
   outfile.write ('%s\n' % yB_line)

'''
#
# Writing the array zApprch1dpx to output file (without of the skiping of the repeated zero values): 
#
outfile.write ('\n    zApprch1dpx[iA,iB] (1/cm**2; Entries: %d x %d )' % (nBins,nBins)) 
outfile.write \
('\nFormat: for iB=0: iA=0 --> nBins-1; then for iB=1: iA=0 --> nBins-1 and so on till iB=nBins-1\n\n')

k=0
for mB in range(2): 
   for mA in range(nBins):
      strVal='{:e}'.format(zApprch1dpx[mA,mB])
      if k == 0:
         zRes_line=strVal
      else:
         zRes_line=zRes_line+', '+strVal
      k += 1
      if k == nInLine:
         outfile.write ('%s\n' % zRes_line)
         k=0
if k != 0:
   outfile.write ('%s\n' % zRes_line)
'''

#
# Writing the array zApprch1dpx to output file (skiping of the repeated zero values): 
#

nInLine=8

outfile.write ('\n    zApprch1dpx[iA,iB] (1/cm**2; Entries: %d x %d with %d per line)' % (nBins,nBins,nInLine)) 
outfile.write \
('\nFormat: for iB=0: iA=0 --> nBins-1, then for iB=1: iA=0 --> nBins-1 and so on till for iB=nBins-1: iA=0 --> nBins-1\n\n')

k=0
countsRes=0
zRes_line=''
nZero=0
valPrev=-1.
for mB in range(nBins): 
   for mA in range(nBins):
      valCrrnt=zApprch1dpx[mA,mB]
      strVal='{:e}'.format(valCrrnt)
      if valPrev == -1.:
         zRes_line=zRes_line+strVal
         k += 1
         countsRes += 1
	 if valCrrnt == 0.:
	    nZero=1
      else:	 
         if valPrev != 0. and valCrrnt != 0.:
            if countsRes == 1:
               zRes_line=zRes_line+', '+strVal+', '
	    else:
               zRes_line=zRes_line+strVal+', '
            k += 1
            countsRes += 1
         if valPrev != 0. and valCrrnt == 0.:
            if countsRes == 1:
               zRes_line=zRes_line+', '+strVal
	    else:
               zRes_line=zRes_line+strVal
            k += 1
            nZero += 1
         if valPrev == 0. and valCrrnt == 0.:
            nZero += 1
         if valPrev == 0. and valCrrnt != 0.:
	    if k < nInLine:
               zRes_line=zRes_line+'('+'{:d}'.format(nZero)+'), '+strVal+', '
	       k += 1
               countsRes += nZero+1
	       nZero=0
	    else:
               zRes_line=zRes_line+'('+'{:d}'.format(nZero)+'),'
               outfile.write ('%s\n' % zRes_line)
               zRes_line=''
	       k=0
               zRes_line=zRes_line+strVal+', '
	       k += 1
               countsRes += nZero+1
	       nZero=0
      if k >= nInLine:
         if nZero == 0:
            outfile.write ('%s\n' % zRes_line)
            zRes_line=''
	    k=0 
#      print 'mA=%d, mB=%d,valCrrnt=%e, valPrev=%e; nZero=%d, k=%d, countsRes=%d' % (mA,mB,valCrrnt,valPrev,nZero,k,countsRes)
      valPrev=valCrrnt        
if nZero !=0:   
   zRes_line=zRes_line+'('+'{:d}'.format(nZero)+')'
   if nZero != nBins*nBins:
      countsRes += nZero
   else:
      countsRes=nBins*nBins
   outfile.write ('%s\n' % zRes_line)

print 'Total counts for zApprch1dpx data =%d' % countsRes
if countsRes != nBins*nBins:
   print 'Something wrong in the writing of zApprch1dpx data to output file!'

#
# Writing the array zApprch1dpy to output file (skiping of the repeated zero values): 
#

nInLine=8

outfile.write ('\n    zApprch1dpy[iA,iB] (1/cm**2; Entries: %d x %d with %d per line)' % (nBins,nBins,nInLine)) 
outfile.write \
('\nFormat: for iB=0: iA=0 --> nBins-1, then for iB=1: iA=0 --> nBins-1 and so on till for iB=nBins-1: iA=0 --> nBins-1\n\n')

k=0
countsRes=0
zRes_line=''
nZero=0
valPrev=-1.
for mB in range(nBins): 
   for mA in range(nBins):
      valCrrnt=zApprch1dpy[mA,mB]
      strVal='{:e}'.format(valCrrnt)
      if valPrev == -1.:
         zRes_line=zRes_line+strVal
         k += 1
         countsRes += 1
	 if valCrrnt == 0.:
	    nZero=1
      else:	 
         if valPrev != 0. and valCrrnt != 0.:
            if countsRes == 1:
               zRes_line=zRes_line+', '+strVal+', '
	    else:
               zRes_line=zRes_line+strVal+', '
            k += 1
            countsRes += 1
         if valPrev != 0. and valCrrnt == 0.:
            if countsRes == 1:
               zRes_line=zRes_line+', '+strVal
	    else:
               zRes_line=zRes_line+strVal
            k += 1
            nZero += 1
         if valPrev == 0. and valCrrnt == 0.:
            nZero += 1
         if valPrev == 0. and valCrrnt != 0.:
	    if k < nInLine:
               zRes_line=zRes_line+'('+'{:d}'.format(nZero)+'), '+strVal+', '
	       k += 1
               countsRes += nZero+1
	       nZero=0
	    else:
               zRes_line=zRes_line+'('+'{:d}'.format(nZero)+'),'
               outfile.write ('%s\n' % zRes_line)
               zRes_line=''
	       k=0
               zRes_line=zRes_line+strVal+', '
	       k += 1
               countsRes += nZero+1
	       nZero=0
      if k >= nInLine:
         if nZero == 0:
            outfile.write ('%s\n' % zRes_line)
            zRes_line=''
	    k=0 
#      print 'mA=%d, mB=%d,valCrrnt=%e, valPrev=%e; nZero=%d, k=%d, countsRes=%d' % (mA,mB,valCrrnt,valPrev,nZero,k,countsRes)
      valPrev=valCrrnt        
if nZero !=0:   
   zRes_line=zRes_line+'('+'{:d}'.format(nZero)+')'
   if nZero != nBins*nBins:
      countsRes += nZero
   else:
      countsRes=nBins*nBins
   outfile.write ('%s\n' % zRes_line)

print 'Total counts for zApprch1dpy data =%d' % countsRes
if countsRes != nBins*nBins:
   print 'Something wrong in the writing of zApprch1dpy data to output file!'

#
# Writing the array zApprch1dpz to output file (skiping of the repeated zero values): 
#

nInLine=8

outfile.write ('\n    zApprch1dpz[iA,iB] (1/cm**2; Entries: %d x %d with %d per line)' % (nBins,nBins,nInLine)) 
outfile.write \
('\nFormat: for iB=0: iA=0 --> nBins-1, then for iB=1: iA=0 --> nBins-1 and so on till for iB=nBins-1: iA=0 --> nBins-1\n\n')

k=0
countsRes=0
zRes_line=''
nZero=0
valPrev=-1.
for mB in range(nBins): 
   for mA in range(nBins):
      valCrrnt=zApprch1dpz[mA,mB]
      strVal='{:e}'.format(valCrrnt)
      if valPrev == -1.:
         zRes_line=zRes_line+strVal
         k += 1
         countsRes += 1
	 if valCrrnt == 0.:
	    nZero=1
      else:	 
         if valPrev != 0. and valCrrnt != 0.:
            if countsRes == 1:
               zRes_line=zRes_line+', '+strVal+', '
	    else:
               zRes_line=zRes_line+strVal+', '
            k += 1
            countsRes += 1
         if valPrev != 0. and valCrrnt == 0.:
            if countsRes == 1:
               zRes_line=zRes_line+', '+strVal
	    else:
               zRes_line=zRes_line+strVal
            k += 1
            nZero += 1
         if valPrev == 0. and valCrrnt == 0.:
            nZero += 1
         if valPrev == 0. and valCrrnt != 0.:
	    if k < nInLine:
               zRes_line=zRes_line+'('+'{:d}'.format(nZero)+'), '+strVal+', '
	       k += 1
               countsRes += nZero+1
	       nZero=0
	    else:
               zRes_line=zRes_line+'('+'{:d}'.format(nZero)+'),'
               outfile.write ('%s\n' % zRes_line)
               zRes_line=''
	       k=0
               zRes_line=zRes_line+strVal+', '
	       k += 1
               countsRes += nZero+1
	       nZero=0
      if k >= nInLine:
         if nZero == 0:
            outfile.write ('%s\n' % zRes_line)
            zRes_line=''
	    k=0 
#      print 'mA=%d, mB=%d,valCrrnt=%e, valPrev=%e; nZero=%d, k=%d, countsRes=%d' % (mA,mB,valCrrnt,valPrev,nZero,k,countsRes)
      valPrev=valCrrnt        
if nZero !=0:   
   zRes_line=zRes_line+'('+'{:d}'.format(nZero)+')'
   if nZero != nBins*nBins:
      countsRes += nZero
   else:
      countsRes=nBins*nBins
   outfile.write ('%s\n' % zRes_line)

outfile.write ('\n')
            

print 'Total counts for zApprch1dpz data =%d' % countsRes
if countsRes != nBins*nBins:
   print 'Something wrong in the writing of zApprch1dpz data to output file!'

outfile.close()
print 'Close the written output file "%s"' % outputFile

#
# Opening the output file for approach-2: 
#
outputFile='dpApprch2.dat'
print 'Open output file "%s"...' % outputFile
outfileFlag=0
try:
   outfile = open(outputFile,'w')
   outfileFlag=1
except:
   print 'Problem to open output file "%s"' % outputFile

nInLine=10

outfile.write ('\n       Tracks: %d , Normalization Factor: %e\n' % (lastTrackNumber,normFactor))

#
# Writing the array xA_2=log10(Upot/Ekin) to output file: 
#
outfile.write ('\n       xA_2=log10(eVca/Ekin) ( Entries: %d with %d per line)\n\n' % (nBins,nInLine))

k=0
for m in range(nBins): 
   strVal='{:f}'.format(xA_2[m])
   if k == 0:
      xA_line=strVal
   else:
      xA_line=xA_line+', '+strVal
   k += 1
   if k == nInLine:
      outfile.write ('%s\n' % xA_line)
      k=0
if k != 0:
   outfile.write ('%s\n' % xA_line)

#
# Writing the array yB_2=log10(Rlarm/b) to output file: 
#
outfile.write ('\n       yB_2=log10(Rlarm/b) ( Entries: %d with %d per line)\n\n' % (nBins,nInLine))

k=0
for m in range(nBins): 
   strVal='{:f}'.format(yB_2[m])
   if k == 0:
      yB_line=strVal
   else:
      yB_line=yB_line+', '+strVal
   k += 1
   if k == nInLine:
      outfile.write ('%s\n' % yB_line)
      k=0
if k != 0:
   outfile.write ('%s\n' % yB_line)

'''
#
# Writing the array zApprch2dpx to output file (without of the skiping of the repeated zero values): 
#
outfile.write ('\n    zApprch2dpx[iA_2,iB_2] (1/cm**2; Entries: %d x %d )' % (nBins,nBins)) 
outfile.write \
('\nFormat: for iB_2=0: iA_2=0 --> nBins-1; then for iB_2=1: iA=0 --> nBins-1 and so on till iB_2=nBins-1\n\n')

k=0
for mB_2 in range(2): 
   for mA_2 in range(nBins):
      strVal='{:e}'.format(zApprch2dpx[mA_2,mB_2])
      if k == 0:
         zRes_line=strVal
      else:
         zRes_line=zRes_line+', '+strVal
      k += 1
      if k == nInLine:
         outfile.write ('%s\n' % zRes_line)
         k=0
if k != 0:
   outfile.write ('%s\n' % zRes_line)
'''

#
# Writing the array zApprch2dpx to output file (skiping of the repeated zero values): 
#

nInLine=8

outfile.write ('\n    zApprch2dpx[iA_2,iB_2] (1/cm**2; Entries: %d x %d with %d per line)' % (nBins,nBins,nInLine)) 
outfile.write \
('\nFormat: for iB_2=0: iA_2=0 --> nBins-1, then for iB_2=1: iA_2=0 --> nBins-1 and so on till for iB_2=nBins-1: iA_2=0 --> nBins-1\n\n')

k=0
countsRes=0
zRes_line=''
nZero=0
valPrev=-1.
for mB_2 in range(nBins): 
   for mA_2 in range(nBins):
      valCrrnt=zApprch2dpx[mA_2,mB_2]
      strVal='{:e}'.format(valCrrnt)
      if valPrev == -1.:
         zRes_line=zRes_line+strVal
         k += 1
         countsRes += 1
	 if valCrrnt == 0.:
	    nZero=1
      else:	 
         if valPrev != 0. and valCrrnt != 0.:
            if countsRes == 1:
               zRes_line=zRes_line+', '+strVal+', '
	    else:
               zRes_line=zRes_line+strVal+', '
            k += 1
            countsRes += 1
         if valPrev != 0. and valCrrnt == 0.:
            if countsRes == 1:
               zRes_line=zRes_line+', '+strVal
	    else:
               zRes_line=zRes_line+strVal
            k += 1
            nZero += 1
         if valPrev == 0. and valCrrnt == 0.:
            nZero += 1
         if valPrev == 0. and valCrrnt != 0.:
	    if k < nInLine:
               zRes_line=zRes_line+'('+'{:d}'.format(nZero)+'), '+strVal+', '
	       k += 1
               countsRes += nZero+1
	       nZero=0
	    else:
               zRes_line=zRes_line+'('+'{:d}'.format(nZero)+'),'
               outfile.write ('%s\n' % zRes_line)
               zRes_line=''
	       k=0
               zRes_line=zRes_line+strVal+', '
	       k += 1
               countsRes += nZero+1
	       nZero=0
      if k >= nInLine:
         if nZero == 0:
            outfile.write ('%s\n' % zRes_line)
            zRes_line=''
	    k=0 
#      print 'mA_2=%d, mB_2=%d,valCrrnt=%e, valPrev=%e; nZero=%d, k=%d, countsRes=%d' % \
#            (mA_2,mB_2,valCrrnt,valPrev,nZero,k,countsRes)
      valPrev=valCrrnt        
if nZero !=0:   
   zRes_line=zRes_line+'('+'{:d}'.format(nZero)+')'
   if nZero != nBins*nBins:
      countsRes += nZero
   else:
      countsRes=nBins*nBins
   outfile.write ('%s\n' % zRes_line)

print 'Total counts for zApprch2dpx data =%d' % countsRes
if countsRes != nBins*nBins:
   print 'Something wrong in the writing of zApprch2dpx data to output file!'

#
# Writing the array zApprch2dpy to output file (skiping of the repeated zero values): 
#

nInLine=8

outfile.write ('\n    zApprch2dpy[iA_2,iB_2] (1/cm**2; Entries: %d x %d with %d per line)' % (nBins,nBins,nInLine)) 
outfile.write \
('\nFormat: for iB_2=0: iA_2=0 --> nBins-1, then for iB_2=1: iA_2=0 --> nBins-1 and so on till for iB_2=nBins-1: iA_2=0 --> nBins-1\n\n')

k=0
countsRes=0
zRes_line=''
nZero=0
valPrev=-1.
for mB_2 in range(nBins): 
   for mA_2 in range(nBins):
      valCrrnt=zApprch2dpy[mA_2,mB_2]
      strVal='{:e}'.format(valCrrnt)
      if valPrev == -1.:
         zRes_line=zRes_line+strVal
         k += 1
         countsRes += 1
	 if valCrrnt == 0.:
	    nZero=1
      else:	 
         if valPrev != 0. and valCrrnt != 0.:
            if countsRes == 1:
               zRes_line=zRes_line+', '+strVal+', '
	    else:
               zRes_line=zRes_line+strVal+', '
            k += 1
            countsRes += 1
         if valPrev != 0. and valCrrnt == 0.:
            if countsRes == 1:
               zRes_line=zRes_line+', '+strVal
	    else:
               zRes_line=zRes_line+strVal
            k += 1
            nZero += 1
         if valPrev == 0. and valCrrnt == 0.:
            nZero += 1
         if valPrev == 0. and valCrrnt != 0.:
	    if k < nInLine:
               zRes_line=zRes_line+'('+'{:d}'.format(nZero)+'), '+strVal+', '
	       k += 1
               countsRes += nZero+1
	       nZero=0
	    else:
               zRes_line=zRes_line+'('+'{:d}'.format(nZero)+'),'
               outfile.write ('%s\n' % zRes_line)
               zRes_line=''
	       k=0
               zRes_line=zRes_line+strVal+', '
	       k += 1
               countsRes += nZero+1
	       nZero=0
      if k >= nInLine:
         if nZero == 0:
            outfile.write ('%s\n' % zRes_line)
            zRes_line=''
	    k=0 
#       print 'mA_2=%d, mB_2=%d,valCrrnt=%e, valPrev=%e; nZero=%d, k=%d, countsRes=%d' % \
#             (mA_2,mB_2,valCrrnt,valPrev,nZero,k,countsRes)
      valPrev=valCrrnt        
if nZero !=0:   
   zRes_line=zRes_line+'('+'{:d}'.format(nZero)+')'
   if nZero != nBins*nBins:
      countsRes += nZero
   else:
      countsRes=nBins*nBins
   outfile.write ('%s\n' % zRes_line)

print 'Total counts for zApprch2dpy data =%d' % countsRes
if countsRes != nBins*nBins:
   print 'Something wrong in the writing of zApprch2dpy data to output file!'

#
# Writing the array zApprch2dpz to output file (skiping of the repeated zero values): 
#

nInLine=8

outfile.write ('\n    zApprch2dpz[iA_2,iB_2] (1/cm**2; Entries: %d x %d with %d per line)' % (nBins,nBins,nInLine)) 
outfile.write \
('\nFormat: for iB_2=0: iA_2=0 --> nBins-1, then for iB_2=1: iA_2=0 --> nBins-1 and so on till for iB_2=nBins-1: iA_2=0 --> nBins-1\n\n')

k=0
countsRes=0
zRes_line=''
nZero=0
valPrev=-1.
for mB_2 in range(nBins): 
   for mA_2 in range(nBins):
      valCrrnt=zApprch2dpz[mA_2,mB_2]
      strVal='{:e}'.format(valCrrnt)
      if valPrev == -1.:
         zRes_line=zRes_line+strVal
         k += 1
         countsRes += 1
	 if valCrrnt == 0.:
	    nZero=1
      else:	 
         if valPrev != 0. and valCrrnt != 0.:
            if countsRes == 1:
               zRes_line=zRes_line+', '+strVal+', '
	    else:
               zRes_line=zRes_line+strVal+', '
            k += 1
            countsRes += 1
         if valPrev != 0. and valCrrnt == 0.:
            if countsRes == 1:
               zRes_line=zRes_line+', '+strVal
	    else:
               zRes_line=zRes_line+strVal
            k += 1
            nZero += 1
         if valPrev == 0. and valCrrnt == 0.:
            nZero += 1
         if valPrev == 0. and valCrrnt != 0.:
	    if k < nInLine:
               zRes_line=zRes_line+'('+'{:d}'.format(nZero)+'), '+strVal+', '
	       k += 1
               countsRes += nZero+1
	       nZero=0
	    else:
               zRes_line=zRes_line+'('+'{:d}'.format(nZero)+'),'
               outfile.write ('%s\n' % zRes_line)
               zRes_line=''
	       k=0
               zRes_line=zRes_line+strVal+', '
	       k += 1
               countsRes += nZero+1
	       nZero=0
      if k >= nInLine:
         if nZero == 0:
            outfile.write ('%s\n' % zRes_line)
            zRes_line=''
	    k=0 
#      print 'mA_2=%d, mB_2=%d,valCrrnt=%e, valPrev=%e; nZero=%d, k=%d, countsRes=%d' % \
#            (mA_2,mB_2,valCrrnt,valPrev,nZero,k,countsRes)
      valPrev=valCrrnt        
if nZero !=0:   
   zRes_line=zRes_line+'('+'{:d}'.format(nZero)+')'
   if nZero != nBins*nBins:
      countsRes += nZero
   else:
      countsRes=nBins*nBins
   outfile.write ('%s\n' % zRes_line)

outfile.write ('\n')
            

print 'Total counts for zApprch2dpz data =%d' % countsRes
if countsRes != nBins*nBins:
   print 'Something wrong in the writing of zApprch2dpz data to output file!'

outfile.close()
print 'Close the written output file "%s"' % outputFile

sys.exit()

###################################################
#
# Reading the data from output files for checking (these part of code  will be used in 
# the script threeApproachesVisualozation.py for later processing and vizualization): 
#
###################################################   

#
# Opening the input file 'dpApprch1.dat': 
#
inputFile='dpApprch1.dat'
print 'Open input file "%s"...' % inputFile
inpfileFlag=0
try:
   inpfile = open(inputFile,'r')
   inpfileFlag=1
except:
   print 'Problem to open input file "%s"' % inputFile

#
# Reading the results from input file: 
#
xAheaderLineNumber=0                                               # Serial number of line with header for xA-Data
yBheaderLineNumber=0                                               # Serial number of line with header for yB-Data
dataDpxHeaderLineNumber=0                                          # Serial number of line with header for dataDpx
dataDpyHeaderLineNumber=0                                          # Serial number of line with header for dataDpy
dataDpzHeaderLineNumber=0                                          # Serial number of line with header for dataDpz
lastLineNumber=0                                                   # Number of the last line 
xAdataFlag=0                                                       # =1 when xA-Data already read
yBdataFlag=0                                                       # =1 when yB-Data already read 
dataDpxFlag=0                                                      # =1 when dataDpx already read
dataDpyFlag=0                                                      # =1 when dataDpy already read
dataDpzFlag=0                                                      # =1 when dataDpz already read
xAdata=[]                                                          # Array of xA-Data
yBdata=[]                                                          # Array of yB-Data
dataDpx=[]                                                         # Array of dataDpx
dataDpy=[]                                                         # Array of dataDpy
dataDpz=[]                                                         # Array of dataDpz

lines=0                                                            # Number of current line from input file   
linesFull=0                                                        # Number of fully filled rows with each type of data
dataNumber=0                                                       # Number of current value of any types of Data
while True:
   lineData=inpfile.readline()
#    print 'line=%d: %s' % (lines,lineData)
   if not lineData:
      break
   lines += 1
   if lines == 2:
      words=lineData.split()
      indx1=words.index('Tracks:')
      lastTrackNumber=int(words[indx1+1])
      indx2=words.index('Factor:')
      normFactor=float(words[indx2+1])
      print 'lastTrackNumber=%d, normFactor=%e' % (lastTrackNumber,normFactor)   
# Header for xA-Data:
   if lines == 4: 
      words=lineData.split()
      indx1=words.index('Entries:')
      xAentries=int(words[indx1+1])
      indx2=words.index('with')
      perLine=int(words[indx2+1])
#       print 'xAdata-Header from %d: words =%s, index1=%d, entries=%d, index2=%d, perLine=%d' % (lines,words,indx1,xAentries,indx2,perLine)
      xAdata=np.zeros(xAentries)
      linesFull=xAentries//perLine
      entriesRem=xAentries-perLine*linesFull
      yBheaderLineNumber=linesFull+7
      if entriesRem > 0:
         yBheaderLineNumber += 1
      print 'yBheaderLineNumber=%d' % yBheaderLineNumber
   if xAdataFlag == 0:
# xA-Data:
      if lines > 5 and lines <= yBheaderLineNumber-2:
         words=lineData.split()
         nWords=len(words)
#          print 'xA-Data from %d: words=%s, number of entries = %d' % (lines,words,nWords)
         for m in range(nWords):
            wordCrrnt=words[m].split(",")
            xAdata[dataNumber]=float(wordCrrnt[0])
	    dataNumber += 1
      if lines == yBheaderLineNumber-2:
         xAdataFlag=1   
         print 'xA-Data(%d entries) already read' % xAentries 
# Header for yB-Data:
   if lines == yBheaderLineNumber:
      words=lineData.split()
      indx1=words.index('Entries:')
      yBentries=int(words[indx1+1])
      indx2=words.index('with')
      perLine=int(words[indx2+1])
#       print 'yBdata-Header from %d: words =%s, index1=%d, entries=%d, index2=%d, perLine=%d' % (lines,words,indx1,yBentries,indx2,perLine)
      yBdata=np.zeros(yBentries)
      linesFull=yBentries//perLine
      entriesRem=yBentries-perLine*linesFull
      dataDpxHeaderLineNumber=yBheaderLineNumber+linesFull+3
      if entriesRem > 0:
         dataDpxHeaderLineNumber += 1
#      print 'dataDpxHeaderLineNumber=%d' % dataDpxHeaderLineNumber
      dataNumber=0
   if xAdataFlag == 1 and yBdataFlag == 0:
# yB-Data:
      if lines >  yBheaderLineNumber+1 and lines <= dataDpxHeaderLineNumber-2:
         words=lineData.split()
         nWords=len(words)
#          print 'yB-Data from %d: words=%s, number of entries = %d' % (lines,words,nWords)
         for m in range(nWords):
            wordCrrnt=words[m].split(",")
            yBdata[dataNumber]=float(wordCrrnt[0])
	    dataNumber += 1
      if lines == dataDpxHeaderLineNumber-2:
         yBdataFlag=1   
         print 'yB-Data(%d entries) already read' % yBentries  
# Header for dataDpx:
   if lines == dataDpxHeaderLineNumber:
      words=lineData.split()
      indx1=words.index('Entries:')
      entriesA=int(words[indx1+1])
      indx2=words.index('x')
      entriesB=int(words[indx2+1])
      indx3=words.index('with')
      perLine=int(words[indx3+1])
#       print 'dataDpx-Header from %d: words =%s, index1=%d, entriesA=%d, index2=%d, entriesB=%d, index3=%d, perLine=%d' % \
#             (lines,words,indx1,entriesA,indx2,entriesB,indx3,perLine)
      dataDpx=np.zeros((entriesA,entriesB))
# If data are written without skiping of the repeated zero values:      
#       linesFull=entriesA*entriesB//perLine
#       entriesRem=entriesA*entriesB-perLine*linesFull
      mA=0
      mB=0
      total_mAmB=0
   if yBdataFlag == 1 and dataDpxFlag == 0:    
# dataDpx: (Format: at the begining, the index mA for "xA-direction" increases and only after that index mB for "yB-direction")
      if lines >  dataDpxHeaderLineNumber+2:
#            print 'line %d: "%s' % (lines,lineData)
         words=lineData.split()
         nWords=len(words)
#          print 'Data from %d: words=%s, number of entries = %d' % (lines,words,nWords)
	 if nWords == 0:
	    dataDpxFlag=1
	    dataDpyHeaderLineNumber=lines+1   
#             print 'dataDpyHeaderLineNumber=%d' % dataDpyHeaderLineNumber
	 else:
            for m in range(nWords):
               wordCrrnt=words[m].split(",")
#	        print 'wordCrrnt: ', wordCrrnt
 	       indxBrsktOpn=wordCrrnt[0].find('(')
 	       if indxBrsktOpn > 0:
 	          indxBrsktCls=wordCrrnt[0].find(')')
# Nonrepeated zero values:
 	       if indxBrsktOpn < 0:
# 	          print 'nonZero value=%e' % float(wordCrrnt[0])
                  dataDpx[mA,mB]=float(wordCrrnt[0])
                  total_mAmB += 1
                  mA += 1
	          if mA  == entriesA:
	             mA=0
		     mB += 1
 	       else:
		  wordHelp=wordCrrnt[0]
 	          nZero=int(wordHelp[indxBrsktOpn+1:indxBrsktCls])
# 	           print 'Number of zero=%d' % nZero
                  for nZ in range(nZero):
		     dataDpx[mA,mB]=0.
                     total_mAmB += 1
                     mA += 1
	             if mA  == entriesA:
	                mA=0
		        mB += 1
         if dataDpxFlag == 1:
            print 'dataDpx(%d x %d entries) already read (total %d)' % (entriesA,entriesB,total_mAmB) 
# Header for dataDpy:
   if lines == dataDpyHeaderLineNumber:
      words=lineData.split()
      indx1=words.index('Entries:')
      entriesA=int(words[indx1+1])
      indx2=words.index('x')
      entriesB=int(words[indx2+1])
      indx3=words.index('with')
      perLine=int(words[indx3+1])
#       print 'dataDpy-Header from %d: words =%s, index1=%d, entriesA=%d, index2=%d, entriesB=%d, index3=%d, perLine=%d' % \
#             (lines,words,indx1,entriesA,indx2,entriesB,indx3,perLine)
      dataDpy=np.zeros((entriesA,entriesB))
# If data are written without skiping of the repeated zero values:      
#       linesFull=entriesA*entriesB//perLine
#       entriesRem=entriesA*entriesB-perLine*linesFull
      mA=0
      mB=0
      total_mAmB=0
   if dataDpxFlag == 1 and dataDpyFlag == 0:    
# dataDpy: (Format: at the begining, the index mA for "xA-direction" increases and only after that index mB for "yB-direction")
      if lines >  dataDpyHeaderLineNumber+2:
         words=lineData.split()
         nWords=len(words)
#          print 'Data from %d: words=%s, number of entries = %d' % (lines,words,nWords)
	 if nWords == 0:
	    dataDpyFlag=1
	    dataDpzHeaderLineNumber=lines+1   
#             print 'dataDpzHeaderLineNumber=%d' % dataDpzHeaderLineNumber
	 else:
            for m in range(nWords):
               wordCrrnt=words[m].split(",")
#	        print 'wordCrrnt: ', wordCrrnt
 	       indxBrsktOpn=wordCrrnt[0].find('(')
 	       if indxBrsktOpn > 0:
 	          indxBrsktCls=wordCrrnt[0].find(')')
# Nonrepeated zero values:
 	       if indxBrsktOpn < 0:
# 	          print 'nonZero value=%e' % float(wordCrrnt[0])
                  dataDpy[mA,mB]=float(wordCrrnt[0])
                  total_mAmB += 1
                  mA += 1
	          if mA  == entriesA:
	             mA=0
		     mB += 1
 	       else:
		  wordHelp=wordCrrnt[0]
 	          nZero=int(wordHelp[indxBrsktOpn+1:indxBrsktCls])
# 	           print 'Number of zero=%d' % nZero
                  for nZ in range(nZero):
		     dataDpy[mA,mB]=0.
                     total_mAmB += 1
                     mA += 1
	             if mA  == entriesA:
	                mA=0
		        mB += 1
         if dataDpyFlag == 1:
            print 'dataDpy(%d x %d entries) already read (total %d)' % (entriesA,entriesB,total_mAmB) 
# Header for dataDpz:
   if lines == dataDpzHeaderLineNumber:
      words=lineData.split()
      indx1=words.index('Entries:')
      entriesA=int(words[indx1+1])
      indx2=words.index('x')
      entriesB=int(words[indx2+1])
      indx3=words.index('with')
      perLine=int(words[indx3+1])
#       print 'dataDpz-Header from %d: words =%s, index1=%d, entriesA=%d, index2=%d, entriesB=%d, index3=%d, perLine=%d' % \
#             (lines,words,indx1,entriesA,indx2,entriesB,indx3,perLine)
      dataDpz=np.zeros((entriesA,entriesB))
# If data are written without skiping of the repeated zero values:      
#       linesFull=entriesA*entriesB//perLine
#       entriesRem=entriesA*entriesB-perLine*linesFull
      mA=0
      mB=0
      total_mAmB=0
   if dataDpyFlag == 1 and dataDpzFlag == 0:    
# dataDpz: (Format: at the begining, the index mA for "xA-direction" increases and only after that index mB for "yB-direction")
      if lines >  dataDpzHeaderLineNumber+2:
         words=lineData.split()
         nWords=len(words)
#          print 'Data from %d: words=%s, number of entries = %d' % (lines,words,nWords)
	 if nWords == 0:
	    dataDpzFlag=1
# Not necessary:	    
# 	    dataDpzHeaderLineNumber=lines+1   
#             print 'dataDpzHeaderLineNumber=%d' % dataDpzHeaderLineNumber
	 else:
            for m in range(nWords):
               wordCrrnt=words[m].split(",")
#	        print 'wordCrrnt: ', wordCrrnt
 	       indxBrsktOpn=wordCrrnt[0].find('(')
 	       if indxBrsktOpn > 0:
 	          indxBrsktCls=wordCrrnt[0].find(')')
# Nonrepeated zero values:
 	       if indxBrsktOpn < 0:
# 	          print 'nonZero value=%e' % float(wordCrrnt[0])
                  dataDpz[mA,mB]=float(wordCrrnt[0])
                  total_mAmB += 1
                  mA += 1
	          if mA  == entriesA:
	             mA=0
		     mB += 1
 	       else:
		  wordHelp=wordCrrnt[0]
 	          nZero=int(wordHelp[indxBrsktOpn+1:indxBrsktCls])
# 	           print 'Number of zero=%d' % nZero
                  for nZ in range(nZero):
		     dataDpz[mA,mB]=0.
                     total_mAmB += 1
                     mA += 1
	             if mA  == entriesA:
	                mA=0
		        mB += 1
         if dataDpzFlag == 1:
            print 'dataDpz(%d x %d entries) already read (total %d)' % (entriesA,entriesB,total_mAmB) 
	    break 
   
inpfile.close()
print 'Close input file "%s"' % inputFile
#
# Comparison written and readed data for checking writing/reading:
#
print 'Comparison written and readed data for checking writing/reading:' 
nBad=0
for mA in range(xAentries):
   if abs(xA[mA]-xAdata[mA]) > 1.e-6:
      nBad += 1
      print 'Something wrong with xA data for %d: xA=%16.10e, xAdata=%16.10e' % (mA,xA[mA],xAdata[mA])
if nBad == 0:
   print 'Everything is OK for xA data' 
   
nBad=0
for mB in range(yBentries):
   if abs(yB[mA]-yBdata[mA]) > 1.e-6:
      nBad += 1
      print 'Something wrong with yB data for %d: yB=%16.10e, yBdata=%16.10e' % (mB,yB[mB],yBdata[mB])
if nBad == 0:
   print 'Everything is OK for yB data' 
   
nBad=0
for mB in range(entriesB):
   for mA in range(entriesA):
      if abs(zApprch1dpx[mA,mB]-dataDpx[mA,mB]) > 1.e-3:
         nBad += 1
         print 'Something wrong with zApprch1dpx data for mA=%d and mB=%d: zApprch1dpx=%16.10e, dataDpx=%16.10e' % \
	       (mA,mB,zApprch1dpx[mA,mB],dataDpx[mA,mB])
if nBad == 0:
   print 'Everything is OK for zApprch1dpx data' 
	 
nBad=0
for mB in range(entriesB):
   for mA in range(entriesA):
      if abs(zApprch1dpy[mA,mB]-dataDpy[mA,mB]) > 1.e-3:
         nBad += 1
         print 'Something wrong with zApprch1dpy data for mA=%d and mB=%d: zApprch1dpy=%16.10e, dataDpy=%16.10e' % \
	       (mA,mB,zApprch1dpy[mA,mB],dataDpy[mA,mB])
if nBad == 0:
   print 'Everything is OK for zApprch1dpy data' 
	 
nBad=0
for mB in range(entriesB):
   for mA in range(entriesA):
      if abs(zApprch1dpz[mA,mB]-dataDpz[mA,mB]) > 1.e-3:
         nBad += 1
         print 'Something wrong with zApprch1dpz data for mA=%d and mB=%d: zApprch1dpz=%16.10e, dataDpz=%16.10e' % \
	       (mA,mB,zApprch1dpz[mA,mB],dataDpz[mA,mB])
if nBad == 0:
   print 'Everything is OK for zApprch1dpz data' 
	 
#
# Comparison written and readed data for checking writing/reading:
#
print 'Comparison written and readed data for checking writing/reading:' 
nBad=0
for mA in range(xAentries):
   if abs(xA[mA]-xAdata[mA]) > 1.e-6:
      nBad += 1
      print 'Something wrong with xA data for %d: xA=%16.10e, xAdata=%16.10e' % (mA,xA[mA],xAdata[mA])
if nBad == 0:
   print 'Everything is OK for xA data' 
   
nBad=0
for mB in range(yBentries):
   if abs(yB[mA]-yBdata[mA]) > 1.e-6:
      nBad += 1
      print 'Something wrong with yB data for %d: yB=%16.10e, yBdata=%16.10e' % (mB,yB[mB],yBdata[mB])
if nBad == 0:
   print 'Everything is OK for yB data' 
   
nBad=0
for mB in range(entriesB):
   for mA in range(entriesA):
      if abs(zApprch1dpx[mA,mB]-dataDpx[mA,mB]) > 1.e-3:
         nBad += 1
         print 'Something wrong with zApprch1dpx data for mA=%d and mB=%d: zApprch1dpx=%16.10e, dataDpx=%16.10e' % \
	       (mA,mB,zApprch1dpx[mA,mB],dataDpx[mA,mB])
if nBad == 0:
   print 'Everything is OK for zApprch1dpx data' 
	 
nBad=0
for mB in range(entriesB):
   for mA in range(entriesA):
      if abs(zApprch1dpy[mA,mB]-dataDpy[mA,mB]) > 1.e-3:
         nBad += 1
         print 'Something wrong with zApprch1dpy data for mA=%d and mB=%d: zApprch1dpy=%16.10e, dataDpy=%16.10e' % \
	       (mA,mB,zApprch1dpy[mA,mB],dataDpy[mA,mB])
if nBad == 0:
   print 'Everything is OK for zApprch1dpy data' 
	 
nBad=0
for mB in range(entriesB):
   for mA in range(entriesA):
      if abs(zApprch1dpz[mA,mB]-dataDpz[mA,mB]) > 1.e-3:
         nBad += 1
         print 'Something wrong with zApprch1dpz data for mA=%d and mB=%d: zApprch1dpz=%16.10e, dataDpz=%16.10e' % \
	       (mA,mB,zApprch1dpz[mA,mB],dataDpz[mA,mB])
if nBad == 0:
   print 'Everything is OK for zApprch1dpz data' 

#
# Opening the input file 'dpApprch2.dat': 
#
inputFile='dpApprch2.dat'
print 'Open input file "%s"...' % inputFile
inpfileFlag=0
try:
   inpfile = open(inputFile,'r')
   inpfileFlag=1
except:
   print 'Problem to open input file "%s"' % inputFile

#
# Reading the results from input file: 
#
xA_2headerLineNumber=0                                             # Serial number of line with header for xA_2data
yB_2headerLineNumber=0                                             # Serial number of line with header for yB_2data
dataDpxHeaderLineNumber_2=0                                        # Serial number of line with header for dataDpx_2
dataDpyHeaderLineNumber_2=0                                        # Serial number of line with header for dataDpy_2
dataDpzHeaderLineNumber_2=0                                        # Serial number of line with header for dataDpz_2
lastLineNumber_2=0                                                 # Number of the last line 
xA_2dataFlag=0                                                     # =1 when xA_2data already read
yB_2dataFlag=0                                                     # =1 when yB_2data already read 
dataDpxFlag_2=0                                                    # =1 when dataDpx_2 already read
dataDpyFlag_2=0                                                    # =1 when dataDpy_2 already read
dataDpzFlag_2=0                                                    # =1 when dataDpz_2 already read
xA_2data=[]                                                        # Array of xA_2data
yB_2data=[]                                                        # Array of yB_2data
dataDpx_2=[]                                                       # Array of dataDpx_2
dataDpy_2=[]                                                       # Array of dataDpy_2
dataDpz_2=[]                                                       # Array of dataDpz_2

lines=0                                                            # Number of current line from input file   
linesFull=0                                                        # Number of fully filled rows with each type of data
dataNumber=0                                                       # Number of current value of any types of Data
while True:
   lineData=inpfile.readline()
#    print 'line=%d: %s' % (lines,lineData)
   if not lineData:
      break
   lines += 1
   if lines == 2:
      words=lineData.split()
      indx1=words.index('Tracks:')
      lastTrackNumber=int(words[indx1+1])
      indx2=words.index('Factor:')
      normFactor=float(words[indx2+1])
      print 'lastTrackNumber=%d, normFactor=%e' % (lastTrackNumber,normFactor)   
# Header for xA-Data:
   if lines == 4: 
      words=lineData.split()
      indx1=words.index('Entries:')
      xAentries=int(words[indx1+1])
      indx2=words.index('with')
      perLine=int(words[indx2+1])
#       print 'xAdata-Header from %d: words =%s, index1=%d, entries=%d, index2=%d, perLine=%d' % (lines,words,indx1,xAentries,indx2,perLine)
      xA_2data=np.zeros(xAentries)
      linesFull=xAentries//perLine
      entriesRem=xAentries-perLine*linesFull
      yB_2headerLineNumber=linesFull+7
      if entriesRem > 0:
         yB_2headerLineNumber += 1
      print 'yB_2headerLineNumber=%d' % yB_2headerLineNumber
   if xA_2dataFlag == 0:
# xA_2data:
      if lines > 5 and lines <= yB_2headerLineNumber-2:
         words=lineData.split()
         nWords=len(words)
#          print 'xA_2data from %d: words=%s, number of entries = %d' % (lines,words,nWords)
         for m in range(nWords):
            wordCrrnt=words[m].split(",")
            xA_2data[dataNumber]=float(wordCrrnt[0])
	    dataNumber += 1
      if lines == yB_2headerLineNumber-2:
         xA_2dataFlag=1   
         print 'xA_2data(%d entries) already read' % xAentries 
# Header for yB_2data:
   if lines == yBheaderLineNumber:
      words=lineData.split()
      indx1=words.index('Entries:')
      yBentries=int(words[indx1+1])
      indx2=words.index('with')
      perLine=int(words[indx2+1])
#       print 'yB_2data-Header from %d: words =%s, index1=%d, entries=%d, index2=%d, perLine=%d' % (lines,words,indx1,yBentries,indx2,perLine)
      yB_2data=np.zeros(yBentries)
      linesFull=yBentries//perLine
      entriesRem=yBentries-perLine*linesFull
      dataDpxHeaderLineNumber_2=yB_2headerLineNumber+linesFull+3
      if entriesRem > 0:
         dataDpxHeaderLineNumber_2 += 1
#      print 'dataDpxHeaderLineNumber_2=%d' % dataDpxHeaderLineNumber_2
      dataNumber=0
   if xA_2dataFlag == 1 and yB_2dataFlag == 0:
# yB-Data:
      if lines >  yB_2headerLineNumber+1 and lines <= dataDpxHeaderLineNumber_2-2:
         words=lineData.split()
         nWords=len(words)
#          print 'yB-Data from %d: words=%s, number of entries = %d' % (lines,words,nWords)
         for m in range(nWords):
            wordCrrnt=words[m].split(",")
            yB_2data[dataNumber]=float(wordCrrnt[0])
	    dataNumber += 1
      if lines == dataDpxHeaderLineNumber_2-2:
         yB_2dataFlag=1   
         print 'yB_2data(%d entries) already read' % yBentries  
# Header for dataDpx_2:
   if lines == dataDpxHeaderLineNumber_2:
      words=lineData.split()
      indx1=words.index('Entries:')
      entriesA=int(words[indx1+1])
      indx2=words.index('x')
      entriesB=int(words[indx2+1])
      indx3=words.index('with')
      perLine=int(words[indx3+1])
#       print 'dataDpx_2-Header from %d: words =%s, index1=%d, entriesA=%d, index2=%d, entriesB=%d, index3=%d, perLine=%d' % \
#             (lines,words,indx1,entriesA,indx2,entriesB,indx3,perLine)
      dataDpx_2=np.zeros((entriesA,entriesB))
# If data are written without skiping of the repeated zero values:      
#       linesFull=entriesA*entriesB//perLine
#       entriesRem=entriesA*entriesB-perLine*linesFull
      mA=0
      mB=0
      total_mAmB=0
   if yB_2dataFlag == 1 and dataDpxFlag_2 == 0:    
# dataDpx: (Format: at the begining, the index mA for "xA-direction" increases and only after that index mB for "yB-direction")
      if lines >  dataDpxHeaderLineNumber_2+2:
#            print 'line %d: "%s' % (lines,lineData)
         words=lineData.split()
         nWords=len(words)
#          print 'Data from %d: words=%s, number of entries = %d' % (lines,words,nWords)
	 if nWords == 0:
	    dataDpxFlag_2=1
	    dataDpyHeaderLineNumber_2=lines+1   
#             print 'dataDpyHeaderLineNumber_2=%d' % dataDpyHeaderLineNumber_2
	 else:
            for m in range(nWords):
               wordCrrnt=words[m].split(",")
#	        print 'wordCrrnt: ', wordCrrnt
 	       indxBrsktOpn=wordCrrnt[0].find('(')
 	       if indxBrsktOpn > 0:
 	          indxBrsktCls=wordCrrnt[0].find(')')
# Nonrepeated zero values:
 	       if indxBrsktOpn < 0:
# 	          print 'nonZero value=%e' % float(wordCrrnt[0])
                  dataDpx_2[mA,mB]=float(wordCrrnt[0])
                  total_mAmB += 1
                  mA += 1
	          if mA  == entriesA:
	             mA=0
		     mB += 1
 	       else:
		  wordHelp=wordCrrnt[0]
 	          nZero=int(wordHelp[indxBrsktOpn+1:indxBrsktCls])
# 	           print 'Number of zero=%d' % nZero
                  for nZ in range(nZero):
		     dataDpx_2[mA,mB]=0.
                     total_mAmB += 1
                     mA += 1
	             if mA  == entriesA:
	                mA=0
		        mB += 1
         if dataDpxFlag_2 == 1:
            print 'dataDpx_2(%d x %d entries) already read (total %d)' % (entriesA,entriesB,total_mAmB) 
# Header for dataDpy_2:
   if lines == dataDpyHeaderLineNumber_2:
      words=lineData.split()
      indx1=words.index('Entries:')
      entriesA=int(words[indx1+1])
      indx2=words.index('x')
      entriesB=int(words[indx2+1])
      indx3=words.index('with')
      perLine=int(words[indx3+1])
#       print 'dataDpy_2-Header from %d: words =%s, index1=%d, entriesA=%d, index2=%d, entriesB=%d, index3=%d, perLine=%d' % \
#             (lines,words,indx1,entriesA,indx2,entriesB,indx3,perLine)
      dataDpy_2=np.zeros((entriesA,entriesB))
# If data are written without skiping of the repeated zero values:      
#       linesFull=entriesA*entriesB//perLine
#       entriesRem=entriesA*entriesB-perLine*linesFull
      mA=0
      mB=0
      total_mAmB=0
   if dataDpxFlag_2 == 1 and dataDpyFlag_2 == 0:    
# dataDpy: (Format: at the begining, the index mA for "xA-direction" increases and only after that index mB for "yB-direction")
      if lines >  dataDpyHeaderLineNumber_2+2:
         words=lineData.split()
         nWords=len(words)
#          print 'Data from %d: words=%s, number of entries = %d' % (lines,words,nWords)
	 if nWords == 0:
	    dataDpyFlag_2=1
	    dataDpzHeaderLineNumber_2=lines+1   
#             print 'dataDpzHeaderLineNumber_2=%d' % dataDpzHeaderLineNumber_2
	 else:
            for m in range(nWords):
               wordCrrnt=words[m].split(",")
#	        print 'wordCrrnt: ', wordCrrnt
 	       indxBrsktOpn=wordCrrnt[0].find('(')
 	       if indxBrsktOpn > 0:
 	          indxBrsktCls=wordCrrnt[0].find(')')
# Nonrepeated zero values:
 	       if indxBrsktOpn < 0:
# 	          print 'nonZero value=%e' % float(wordCrrnt[0])
                  dataDpy_2[mA,mB]=float(wordCrrnt[0])
                  total_mAmB += 1
                  mA += 1
	          if mA  == entriesA:
	             mA=0
		     mB += 1
 	       else:
		  wordHelp=wordCrrnt[0]
 	          nZero=int(wordHelp[indxBrsktOpn+1:indxBrsktCls])
# 	           print 'Number of zero=%d' % nZero
                  for nZ in range(nZero):
		     dataDpy_2[mA,mB]=0.
                     total_mAmB += 1
                     mA += 1
	             if mA  == entriesA:
	                mA=0
		        mB += 1
         if dataDpyFlag_2 == 1:
            print 'dataDpy_2(%d x %d entries) already read (total %d)' % (entriesA,entriesB,total_mAmB) 
# Header for dataDpz_2:
   if lines == dataDpzHeaderLineNumber_2:
      words=lineData.split()
      indx1=words.index('Entries:')
      entriesA=int(words[indx1+1])
      indx2=words.index('x')
      entriesB=int(words[indx2+1])
      indx3=words.index('with')
      perLine=int(words[indx3+1])
#       print 'dataDpz-Header from %d: words =%s, index1=%d, entriesA=%d, index2=%d, entriesB=%d, index3=%d, perLine=%d' % \
#             (lines,words,indx1,entriesA,indx2,entriesB,indx3,perLine)
      dataDpz_2=np.zeros((entriesA,entriesB))
# If data are written without skiping of the repeated zero values:      
#       linesFull=entriesA*entriesB//perLine
#       entriesRem=entriesA*entriesB-perLine*linesFull
      mA=0
      mB=0
      total_mAmB=0
   if dataDpyFlag_2 == 1 and dataDpzFlag_2 == 0:    
# dataDpz: (Format: at the begining, the index mA for "xA-direction" increases and only after that index mB for "yB-direction")
      if lines >  dataDpzHeaderLineNumber_2+2:
         words=lineData.split()
         nWords=len(words)
#          print 'Data from %d: words=%s, number of entries = %d' % (lines,words,nWords)
	 if nWords == 0:
	    dataDpzFlag_2=1
# Not necessary:	    
# 	    dataDpzHeaderLineNumber_2=lines+1   
#             print 'dataDpzHeaderLineNumber_2=%d' % dataDpzHeaderLineNumber_2
	 else:
            for m in range(nWords):
               wordCrrnt=words[m].split(",")
#	        print 'wordCrrnt: ', wordCrrnt
 	       indxBrsktOpn=wordCrrnt[0].find('(')
 	       if indxBrsktOpn > 0:
 	          indxBrsktCls=wordCrrnt[0].find(')')
# Nonrepeated zero values:
 	       if indxBrsktOpn < 0:
# 	          print 'nonZero value=%e' % float(wordCrrnt[0])
                  dataDpz_2[mA,mB]=float(wordCrrnt[0])
                  total_mAmB += 1
                  mA += 1
	          if mA  == entriesA:
	             mA=0
		     mB += 1
 	       else:
		  wordHelp=wordCrrnt[0]
 	          nZero=int(wordHelp[indxBrsktOpn+1:indxBrsktCls])
# 	           print 'Number of zero=%d' % nZero
                  for nZ in range(nZero):
		     dataDpz_2[mA,mB]=0.
                     total_mAmB += 1
                     mA += 1
	             if mA  == entriesA:
	                mA=0
		        mB += 1
         if dataDpzFlag_2 == 1:
            print 'dataDpz_2(%d x %d entries) already read (total %d)' % (entriesA,entriesB,total_mAmB) 
	    break 
   
inpfile.close()
print 'Close input file "%s"' % inputFile

#
# Comparison written and readed data for checking writing/reading:
#
print 'Comparison written and readed data for checking writing/reading:' 
nBad=0
for mA in range(xAentries):
   if abs(xA_2[mA]-xA_2data[mA]) > 1.e-6:
      nBad += 1
      print 'Something wrong with xA data for %d: xA=%16.10e, xAdata=%16.10e' % (mA,xA[mA],xAdata[mA])
if nBad == 0:
   print 'Everything is OK for xA data' 
   
nBad=0
for mB in range(yBentries):
   if abs(yB_2[mA]-yB_2data[mA]) > 1.e-6:
      nBad += 1
      print 'Something wrong with yB data for %d: yB=%16.10e, yBdata=%16.10e' % (mB,yB[mB],yBdata[mB])
if nBad == 0:
   print 'Everything is OK for yB data' 
   
nBad=0
for mB in range(entriesB):
   for mA in range(entriesA):
      if abs(zApprch2dpx[mA,mB]-dataDpx_2[mA,mB]) > 1.e-3:
         nBad += 1
         print 'Something wrong with zApprch2dpx data for mA=%d and mB=%d: zApprch2dpx=%16.10e, dataDpx_2=%16.10e' % \
	       (mA,mB,zApprch2dpx[mA,mB],dataDpx_2[mA,mB])
if nBad == 0:
   print 'Everything is OK for zApprch2dpx data' 
	 
nBad=0
for mB in range(entriesB):
   for mA in range(entriesA):
      if abs(zApprch2dpy[mA,mB]-dataDpy_2[mA,mB]) > 1.e-3:
         nBad += 1
         print 'Something wrong with zApprch2dpy data for mA=%d and mB=%d: zApprch2dpy=%16.10e, dataDpy_2=%16.10e' % \
	       (mA,mB,zApprch2dpy[mA,mB],dataDpy_2[mA,mB])
if nBad == 0:
   print 'Everything is OK for zApprch2dpy data' 
	 
nBad=0
for mB in range(entriesB):
   for mA in range(entriesA):
      if abs(zApprch2dpz[mA,mB]-dataDpz_2[mA,mB]) > 1.e-3:
         nBad += 1
         print 'Something wrong with zApprch2dpz data for mA=%d and mB=%d: zApprch2dpz=%16.10e, dataDpz_2=%16.10e' % \
	       (mA,mB,zApprch2dpz[mA,mB],dataDpz_2[mA,mB])
if nBad == 0:
   print 'Everything is OK for zApprch2dpz data' 

sys.exit()

