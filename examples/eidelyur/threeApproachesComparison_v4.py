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
# Matrix describes the electron - collision during time interval 'deltaT'
# in the system coordinates of "guiding center" of electron
# input - 6-vectors for electron and ion before collision and time step deltaT; 
# output - both vectors after collision: 
#
def guidingCenterCollision(vectrElec,vectrIon,deltaT):
   mOmegaLarm=m_elec*omega_L                                       # g/sec
   denom=((vectrIon[0]-vectrElec[3]/mOmegaLarm)**2+ \
          (vectrIon[2]-vectrElec[2])**2+(vectrIon[4]-vectrElec[4])**2)**(3/2)        # cm^3
   omega_gc=Z_ion*q_elec**2/(mOmegaLarm*denom)                     # 1/sec
   phase=omega_gc*deltaT                                           # dimensionless
   sinOmega_gc=math.sin(phase)
   cosOmega_gc=1.-math.cos(phase)
   dx_gc=(vectrIon[0]-vectrElec[3]/mOmegaLarm)*cosOmega_gc+ \
         (vectrIon[2]-vectrElec[2])*sinOmega_gc                    # cm 
   dy_gc=-(vectrIon[0]-vectrElec[3]/mOmegaLarm)*sinOmega_gc+ \
          (vectrIon[2]-vectrElec[2])*cosOmega_gc                   # cm
   dpz=mOmegaLarm*phase*(vectrIon[4]-vectrElec[4])                 # g*cm/sec
   vectrIon[1] +=  mOmegaLarm*dy_gc                                # g*cm/sec 
   vectrIon[3] += -mOmegaLarm*dx_gc                                # g*cm/sec
   vectrIon[5] += -dpz                                             # g*cm/sec
   vectrElec[2] += dy_gc                                           # cm
   vectrElec[3] += mOmegaLarm*dx_gc                                # g*cm/sec 
   vectrIon[5] +=  dpz                                             # g*cm/sec
   return vectrElec,vectrIon                                       

z_elecCrrnt=np.zeros(6)                                            # 6-vector for electron (for Approach_1)
z_ionCrrnt=np.zeros(6)                                             # 6-vector for ion (for Approach_1)
z_elecCrrnt_2=np.zeros(6)                                          # 6-vector for electron (for Approach_2)
z_ionCrrnt_2=np.zeros(6)                                           # 6-vector for ion (for Approach_2)
z_elecCrrnt_prev=np.zeros(6)                                       # 6-vector for electron for previous time step (for Approach_2)
z_elecCrrnt_gc=np.zeros(6)                                         # 6-vector for electron in #guiding center" system (for Approach_2)

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

mPrev=0
bPrev=0.
sumPoints=0
sumPoints_2=0

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

pointTrack=np.zeros(nA*nB)  
larmorNumber=np.zeros(nA*nB)
cpuTime=np.zeros(nA*nB)
pointTrack_2=np.zeros(nA*nB)  
larmorNumber_2=np.zeros(nA*nB)
cpuTime_2=np.zeros(nA*nB)


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
# For approach_2 (with averaging over nLarmorAvrgng larmor rotation):
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
      eVtran=eVrmsLong*root                                        # cm/sec
      rho_larm=eVtran/omega_L                                      # cm
      kinEnergy=m_elec*(eVtran**2+eVrmsLong**2)/2.                 # erg
      rhoCrrnt=rhoCrit+rho_larm                                    # cm
      halfLintr=rho_larm/math.pow(10.,crrntB[iB])                  # cm
      timePath=halfLintr/eVrmsLong                                 # sec
      numbLarmor=int(timePath/T_larm)                              # dimensionless
      if numbLarmor >= 40:
         timeStart=os.times()
         trackNumb += 1                                            # Tracks are enumerated from 0!  
         trackNumb_2 += 1                                          # Tracks are enumerated from 0!  
         larmorNumber[trackNumb]=numbLarmor
	 print 'iA=%d, iB=%d, trackNumber_1=%d, numbLarmor=%d; trackNumber_2=%d' % (iA,iB,trackNumb,numbLarmor,trackNumb_2)   
         timePoints=int(numbLarmor*stepsNumberOnGyro)              # for approach_1; dimensionless
         timePoints_2=int(numbLarmor//nLarmorAvrgng)               # for approach_2; dimensionless
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Approach_1: dragging without averaging over larmor rotation
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# To draw only first trajectories (for checking only):
#
#------- Start of calculations for approach_1 --------------
#
         if trackNumb == 0:
	    rhoFirstTurn=rhoCrrnt
	    rhoLarmorFirstTurn=rho_larm
# 6-vectors for ion and electron and distance 'b' between them for the first trajectories 
# (for checking only; indices 0-5 for electron, indices 6-11 for ion and index=12 for 'b'):
            prtclCoor=np.zeros((13,timePoints))                    
# Current distance from origin of the coordinate system to electron along the trajectory; cm
         bCrrnt=np.zeros(timePoints)                               # cm
# Current log10 of two important ratios; dimensionless:
         larmR_bCrrnt=np.zeros(timePoints)                         # ratio R_larmor/b; dimensionless
         uPot_enrgKinCrrnt=np.zeros(timePoints)                    # ratio potential_energy/kinetic_energy; dimensionless
         dpApprch_1Crrnt=np.zeros((3,timePoints))                  # 1/cm^2; deltaPapprch_1=dpFactor*dpApprch_1Crrnt; dpFactor=q_e^2*timeStep
         for m in range(6): 
            z_ionCrrnt[m]=0.                                       # Current initial zero-vector for ion
            z_elecCrrnt[m]=0.                                      # Zeroing out of vector for electron
# Current initial vector for electron:
         z_elecCrrnt[Ix]=rhoCrrnt                                  # x, cm
         z_elecCrrnt[Iz]=-halfLintr                                # z, cm
         z_elecCrrnt[Ipy]=m_elec*eVtran                            # py, g*cm/sec
         z_elecCrrnt[Ipz]=m_elec*eVrmsLong                         # pz, g*cm/sec
#-----------------------------------------------
# Main action - dragging of the current trajectories (for given i and j)
#
         for k in range(timePoints):
 	    z_elecCrrnt=matr_elec.dot(z_elecCrrnt)                 # electron's dragging
 	    z_ionCrrnt=matr_ion.dot(z_ionCrrnt)                    # ion's dragging
# Current distance from origin of the coordinate system to electron along the it's trajectory 
# (It means that ion is fixed in the origin); cm:
#  	    bCrrnt[k]=np.sqrt(z_elecCrrnt[0]**2+z_elecCrrnt[2]**2+z_elecCrrnt[4]**2)
# Current distance between ion and electron; cm:
 	    bCrrnt[k]=np.sqrt((z_ionCrrnt[0]-z_elecCrrnt[0])**2+ \
	                      (z_ionCrrnt[2]-z_elecCrrnt[2])**2+ \
			      (z_ionCrrnt[4]-z_elecCrrnt[4])**2)
# Current log10 of two important ratios:  
	    larmR_bCrrnt[k]=math.log10(rho_larm/bCrrnt[k])         # dimensionless 
	    if maxLarmR_b_1 < larmR_bCrrnt[k]:
	       maxLarmR_b_1=larmR_bCrrnt[k]
	    if minLarmR_b_1 > larmR_bCrrnt[k]:
	       minLarmR_b_1=larmR_bCrrnt[k]
	    uPot_enrgKinCrrnt[k]=math.log10((q_elec**2/bCrrnt[k])/kinEnergy)         # dimensionless 
	    if maxUpot_enrgKin_1 < uPot_enrgKinCrrnt[k]:
	       maxUpot_enrgKin_1=uPot_enrgKinCrrnt[k]
	    if minUpot_enrgKin_1 > uPot_enrgKinCrrnt[k]:
	       minUpot_enrgKin_1=uPot_enrgKinCrrnt[k]
# Current values to calculate deltaPapprch_1=dpFactor*dpApprch_1Crrnt; dpFactor=q_e^2*timeStep:  
            for ic in range(3):
	       dpApprch_1Crrnt[ic,k]= abs(z_elecCrrnt[2*ic])/ bCrrnt[k]**2           # 1/cm^2;  
# To draw only first trajectory (for checking only):
            if trackNumb == 0:
               for ic in range(6):
                  prtclCoor[ic,pointTrack[trackNumb]]=z_elecCrrnt[ic]                # 6-vector for electron
                  prtclCoor[6+ic,pointTrack[trackNumb]]=z_ionCrrnt[ic]               # 6-vector for ion 
	       prtclCoor[12,pointTrack[trackNumb]]=bCrrnt[k]                         # cm
#  	 print 'i=%d, j=%d, k=%d: bPrev=%e (m=%d), bCrrnt=%e, x=%e, y=%e, z=%e' % \
#  	       (i,j,k,1.e+4*bPrev,mPrev,1.e+4*bCrrnt[k],1.e+4*z_elecCrrnt[0],1.e+4*z_elecCrrnt[2],1.e+4*z_elecCrrnt[4])       
#
# Taking into  account transfer of momentum for both particles:
#
            for ic in range(3):
               z_ionCrrnt[2*ic+1] += -dpFactor*dpApprch_1Crrnt[ic,k]
               z_elecCrrnt[2*ic+1] += dpFactor*dpApprch_1Crrnt[ic,k]
            pointTrack[trackNumb] += 1
# 
# End of gragging of the current trajectory	  
#-----------------------------------------------
# 
         if trackNumb == 0: 
# First definition of the total distance from origin of the coordinate system to electron along the trajectory; cm:
 	    b=bCrrnt                                               # cm
# First definition of the total log10 of two important ratios; dimensionless:  
	    larmR_b=larmR_bCrrnt                                   # dimensionless 
	    uPot_enrgKin=uPot_enrgKinCrrnt                         # dimensionless 
# First definition of the values to calculate deltaPapprch_1=dpFactor*dpApprch_1Crrnt; dpFactor=q_e^2*timeStep:
            dpxApprch_1=dpApprch_1Crrnt[0,:]                       # 1/cm^2;  
            dpyApprch_1=dpApprch_1Crrnt[1,:]                       # 1/cm^2;  
            dpzApprch_1=dpApprch_1Crrnt[2,:]                       # 1/cm^2;  
         else:  
# Total distance from origin of the coordinate system to electron along the trajectory:
 	    b=np.concatenate((b,bCrrnt),axis=0)                    # cm
# Total log10 of two important ratios; dimensionless :  
	    larmR_b=np.concatenate((larmR_b,larmR_bCrrnt),axis=0)                  
	    uPot_enrgKin=np.concatenate((uPot_enrgKin,uPot_enrgKinCrrnt),axis=0)        
# Total values to calculate deltaPapprch_1=dpFactor*dpApprch_1Crrnt; dpFactor=q_e^2*timeStep:
            dpxApprch_1=np.concatenate((dpxApprch_1,dpApprch_1Crrnt[0,:]),axis=0)                 # 1/cm^2;  
            dpyApprch_1=np.concatenate((dpyApprch_1,dpApprch_1Crrnt[1,:]),axis=0)                 # 1/cm^2;  
            dpzApprch_1=np.concatenate((dpzApprch_1,dpApprch_1Crrnt[2,:]),axis=0)                 # 1/cm^2;  
#          print 'trackNumb:%d: shapes: b=%d, larmR_b=%d, uPot=%d, dpx=%d, dpy=%d, dpz=%d' % \
# 	          (trackNumb,b.shape[0],larmR_b.shape[0],uPot_enrgKin.shape[0], \
# 		   dpxApprch_1.shape[0],dpyApprch_1.shape[0],dpzApprch_1.shape[0])  
#
         lastTrackNumber=trackNumb+1                               # quantity of tracks = trackNumb + 1!     
         sumPoints += pointTrack[trackNumb]
         timeEnd=os.times()
	 cpuTime[trackNumb]=1.e+6*(float(timeEnd[0])-float(timeStart[0]))            # CPU time , mks
         cpuTimeTotal += cpuTime[trackNumb] 
#
#------- End of calculations for approach_1 --------------
#
# for i in range(lastTrackNumber):
#    print 'Track %d: larmor turns=%d, cpuTime(mks)=%e, time per turn(mks)=%6.1f' % \
#          (i,larmorNumber[i],cpuTime[i],cpuTime[i]/larmorNumber[i])

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Approach_2: dragging with averaging over nLarmorAvrgng larmor rotation
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#
#------- Start of calculations for approach_2 --------------
#
         if trackNumb_2 == 0:
	    rhoFirstTurn=rhoCrrnt
	    rhoLarmorFirstTurn=rho_larm
# 6-vectors for ion and electron and distance 'b' between them for the first trajectories 
# (for checking only; indices 0-5 for electron, indices 6-11 for ion and index=12 for 'b'):
         prtclCoor_2=np.zeros((13,timePoints_2))                    
# Current distance from origin of the coordinate system to electron along the trajectory; cm
         bCrrnt_2=np.zeros(timePoints_2)                           # cm
# Current log10 of two important ratios; dimensionless:
         larmR_bCrrnt_2=np.zeros(timePoints_2)                     # ratio R_larmor/b; dimensionless
         uPot_enrgKinCrrnt_2=np.zeros(timePoints_2)                # ratio potential_energy/kinetic_energy; dimensionless
         dpApprch_2Crrnt=np.zeros((3,timePoints_2))                # g*cm/sec
         timeEnd=os.times()
         for m in range(6): 
            z_ionCrrnt_2[m]=0.                                     # Current initial zero-vector for ion
            z_elecCrrnt_2[m]=0.                                    # Zeroing out of vector for electron
# Current initial vector for electron:
         z_elecCrrnt_2[Ix]=rhoCrrnt                                # x, cm
         z_elecCrrnt_2[Iz]=-halfLintr                              # z, cm
         z_elecCrrnt_2[Ipy]=m_elec*eVtran                          # py, g*cm/sec
         z_elecCrrnt_2[Ipz]=m_elec*eVrmsLong                       # pz, g*cm/sec
	 z_elecCrrnt_gc=toGuidingCenter(z_elecCrrnt_2)             # transfer to system of guiding center
#-----------------------------------------------
# Main action - dragging of the current trajectories (for given i and j)
#
         for k in range(timePoints_2):
 	    z_elecCrrnt_gc=matr_elec_2.dot(z_elecCrrnt_gc)         # electron's dragging for first nhalf timeStep
 	    z_ionCrrnt_2=matr_ion_2.dot(z_ionCrrnt_2)              # ion's dragging for first half timeStep
 	    z_elecCrrnt_gc=matr_elec_2.dot(z_elecCrrnt_gc)         # electron's dragging for second half timeStep
# gragging both paticles through interaction point:
	    z_elecCrrnt_gc,z_ionCrrnt=guidingCenterCollision(z_elecCrrnt_gc,z_ionCrrnt_2,timeStep_2)    
 	    z_ionCrrnt_2=matr_ion_2.dot(z_ionCrrnt_2)              # ion's dragging for second half timeStep
	    z_elecCrrnt_2=fromGuidingCenter(z_elecCrrnt_gc)        # transfer from system of guiding center 
# Current distance between ion and electron; cm:
 	    bCrrnt_2[k]=np.sqrt((z_ionCrrnt_2[0]-z_elecCrrnt_2[0])**2+ \
	                        (z_ionCrrnt_2[2]-z_elecCrrnt_2[2])**2+ \
			        (z_ionCrrnt_2[4]-z_elecCrrnt_2[4])**2)
# Current log10 of two important ratios:  
	    larmR_bCrrnt_2[k]=math.log10(rho_larm/bCrrnt_2[k])     # dimensionless 
	    if maxLarmR_b_2 < larmR_bCrrnt_2[k]:
	       maxLarmR_b_2=larmR_bCrrnt_2[k]
	    if minLarmR_b_2 > larmR_bCrrnt_2[k]:
	       minLarmR_b_2=larmR_bCrrnt_2[k]
	    uPot_enrgKinCrrnt_2[k]=math.log10((q_elec**2/bCrrnt_2[k])/kinEnergy)     # dimensionless 
	    if maxUpot_enrgKin_2 < uPot_enrgKinCrrnt_2[k]:
	       maxUpot_enrgKin_2=uPot_enrgKinCrrnt_2[k]
	    if minUpot_enrgKin_2 > uPot_enrgKinCrrnt_2[k]:
	       minUpot_enrgKin_2=uPot_enrgKinCrrnt_2[k]
# Current values to calculate deltaPapprch_2:  
#            for ic in range(3):
#	       dpApprch_2Crrnt[ic,k]= z_elecCrrnt_prev[2*ic]-z_elecCrrnt_2[2*ic]     # g*scm/sec 
	    z_elecCrrnt_prev=z_elecCrrnt_2.copy() 
# To draw only first trajectory (for checking only):
            if trackNumb_2 == 0:
               for ic in range(6):
                  prtclCoor_2[ic,pointTrack_2[trackNumb_2]]=z_elecCrrnt_2[ic]        # 6-vector for electron
                  prtclCoor_2[6+ic,pointTrack_2[trackNumb_2]]=z_ionCrrnt_2[ic]       # 6-vector for ion 
	       prtclCoor_2[12,pointTrack_2[trackNumb_2]]=bCrrnt_2[k]                 # cm
            pointTrack_2[trackNumb_2] += 1
#
# End of gragging of the current trajectory	  
#-----------------------------------------------
#
         if trackNumb_2 == 0: 
# First definition of the total distance from origin of the coordinate system to electron along the trajectory; cm:
 	    b_2=bCrrnt_2                                           # cm 
### # First definition of the total log10 of two important ratios; dimensionless:  
	    larmR_b_2=larmR_bCrrnt_2                               # dimensionless 
	    uPot_enrgKin_2=uPot_enrgKinCrrnt_2                     # dimensionless 
# First definition of the values deltaPapprch_2:
            dpxApprch_2=dpApprch_2Crrnt[0,:]                       # g*cm/sec  
            dpyApprch_2=dpApprch_2Crrnt[1,:]                       # g*cm/sec  
            dpzApprch_2=dpApprch_2Crrnt[2,:]                       # g*cm/sec  
         else:  
# Total distance from origin of the coordinate system to electron along the trajectory:
 	    b_2=np.concatenate((b_2,bCrrnt),axis=0)                # cm
# Total log10 of two important ratios; dimensionless :  
	    larmR_b_2=np.concatenate((larmR_b_2,larmR_bCrrnt_2),axis=0)                  
	    uPot_enrgKin_2=np.concatenate((uPot_enrgKin_2,uPot_enrgKinCrrnt_2),axis=0)        
# Total values deltaPapprch_2:
            dpxApprch_2=np.concatenate((dpxApprch_2,dpApprch_2Crrnt[0,:]),axis=0)                      # g*cm/sec  
            dpyApprch_2=np.concatenate((dpyApprch_2,dpApprch_2Crrnt[1,:]),axis=0)                      # g*cm/sec
            dpzApprch_2=np.concatenate((dpzApprch_2,dpApprch_2Crrnt[2,:]),axis=0)                      # g*cm/sec 
#          print 'trackNumb_2:%d: shapes: b=%d, larmR_b_2=%d, uPot=%d, dpx=%d, dpy=%d, dpz=%d' % \
# 	         (trackNumb_2,b.shape[0],larmR_b_2.shape[0],uPot_enrgKin_2.shape[0], \
# 		 dpxApprch_2.shape[0],dpyApprch_2.shape[0],dpzApprch_2.shape[0])  
#
         lastTrackNumber_2=trackNumb_2+1                           # quantity of tracks = trackNumb_2 + 1!     
         sumPoints_2 += pointTrack_2[trackNumb_2]
         timeEnd=os.times()
	 cpuTime_2[trackNumb_2]=1.e+6*(float(timeEnd[0])-float(timeStart[0]))        # CPU time , mks
         cpuTimeTotal += cpuTime_2[trackNumb_2] 
#
#------- End of calculations for approach_2 --------------
#

print 'Approach_1: for %d tracks number of points is %d' % (lastTrackNumber,sumPoints)
print 'Approach_2: for %d tracks number of points is %d' % (lastTrackNumber_2,sumPoints_2)

# for i in range(lastTrackNumber):
#    print 'Track %d: larmor turns=%d, cpuTime(mks)=%e, time per turn(mks)=%6.1f' % \
#          (i,larmorNumber[i],cpuTime[i],cpuTime[i]/larmorNumber[i])

print 'cpuTimeTotal(mksec) = %e' % cpuTimeTotal

###################################################
#
#       Data processing of the first approach:
#: layout of arrays xA=uPot_enrgKin, yB=larmR_b to nBins channels
# and arrays zApprch1dpx, zApprch1dpy, zApprch1dpz to (nBins x nBins) channels
#
###################################################   

nBins=80

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

timeStart=os.times()
for nPoint in range(int(sumPoints)):
   for iA in range(nBins):
      searchAflag=0
      if (xAedges[iA] <= uPot_enrgKin[nPoint] < xAedges[iA+1]):
         if xAnumb[iA] == 0:
	    xA[iA]=uPot_enrgKin[nPoint]                            # log10(Upot/Ekin)
	 else:
	    xA[iA]=(xA[iA]*xAnumb[iA]+uPot_enrgKin[nPoint])/(xAnumb[iA]+1)           # averaging inside bin iA
         xAnumb[iA] += 1
	 searchAflag=1
	 break
   if searchAflag == 0:
      xA[nBins-1]=(xA[nBins-1]*xAnumb[nBins-1]+uPot_enrgKin[nPoint])/(xAnumb[nBins-1]+1)          # averaging inside bin iA 
      xAnumb[nBins-1] += 1
   for iB in range(nBins):
      searchBflag=0
      if (yBedges[iB] <= larmR_b[nPoint] < yBedges[iB+1]):
         if yBnumb[iB] == 0:
	    yB[iB]=larmR_b[nPoint]                                 # log10(Rlarm/b)
	 else:
	    yB[iB]=(yB[iB]*yBnumb[iB]+larmR_b[nPoint])/(yBnumb[iB]+1)                # averaging inside bin iB  
         yBnumb[iB] += 1
	 searchBflag=1
	 break
   if searchBflag == 0:
      yB[nBins-1]=(yB[nBins-1]*yBnumb[nBins-1]+larmR_b[nPoint])/(yBnumb[nBins-1]+1)  # averaging inside bin iB 
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
#      print 'iA=%d, iB=%d: xA=%f (%d), xB=%f (%d), zApprch1dpx=%e' % \
#            (iA,iB,xA[iA],xAnumb[iA],yB[iB],yBnumb[iB],zApprch1dpx[iA,iB])
   
print 'sumA=%d, sumB=%d, sumC=%d; sumPoints=%d' % (sumAnumb,sumBnumb,sumCnumb,int(sumPoints)) 

timeEnd=os.times()
runTime=1.e+6*(float(timeEnd[0])-float(timeStart[0]))              # CPU time , mks
print 'runTime(mksec) = %e' % runTime

#
# Checking of the first trajectories:
#
pointsTot=len(prtclCoor[0,:])
turns=10                                                           # Number of larmorturns for drawing 
points=turns*stepsNumberOnGyro                                     # Number of points for drawing

lengthArrowElc=4
fig10=plt.figure(10)
ax10=fig10.gca(projection='3d')
ax10.plot(1.e+4*prtclCoor[0,0:points],1.e+4*prtclCoor[2,0:points],1.e+4*prtclCoor[4,0:points],'-r',linewidth=2)
ax10.plot(1.e+4*prtclCoor[0,0:lengthArrowElc],1.e+4*prtclCoor[2,0:lengthArrowElc],1.e+4*prtclCoor[4,0:lengthArrowElc],'-b',linewidth=2)
ax10.plot(1.e+4*prtclCoor[0,points-lengthArrowElc:points],1.e+4*prtclCoor[2,points-lengthArrowElc:points], \
          1.e+4*prtclCoor[4,points-lengthArrowElc:points],'-b',linewidth=2)
plt.xlabel('x, $\mu m$',color='m',fontsize=16)
plt.ylabel('y, $\mu m$',color='m',fontsize=16)
ax10.set_zlabel('z, $\mu m$',color='m',fontsize=16)
plt.title(('First Electron Trajectory (Start; $N_L=$%d):\nImpact Parameter=%5.2f $\mu$m, $R_L$=%5.2f $\mu$m' \
           % (larmorNumber[0],1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)),color='m',fontsize=16)

# Data to traw the arrow in the direction of motion:
#
arrwBegxIon=1.e+7*prtclCoor[6,points/2+50]
arrwBegyIon=1.e+7*prtclCoor[8,points/2+50]
arrwBegzIon=1.e+7*prtclCoor[10,points/2+50]
arrwEndxIon=1.e+7*prtclCoor[6,points/2]
arrwEndyIon=1.e+7*prtclCoor[8,points/2]
arrwEndzIon=1.e+7*prtclCoor[10,points/2]
fig15=plt.figure(15)
ax15=fig15.gca(projection='3d')
ax15.plot(1.e+7*prtclCoor[6,0:points],1.e+7*prtclCoor[8,0:points],1.e+7*prtclCoor[10,0:points],'-b',linewidth=2)
ax15.plot([arrwBegxIon,arrwEndxIon],[arrwBegyIon,arrwEndyIon],[arrwBegzIon,arrwEndzIon],color='r',alpha=0.8,lw=2)
plt.xlabel('x, $nm$',color='m',fontsize=16)
plt.ylabel('y, $nm$',color='m',fontsize=16)
ax15.set_zlabel('z, $nm$',color='m',fontsize=16)
plt.title('First Ion Trajectory (Start)',color='m',fontsize=16)

fig20=plt.figure(20)
ax20=fig20.gca(projection='3d')
ax20.plot(1.e+4*prtclCoor[0,pointsTot-points:pointsTot],1.e+4*prtclCoor[2,pointsTot-points:pointsTot], \
          1.e+4*prtclCoor[4,pointsTot-points:pointsTot],'-r',linewidth=2)
ax20.plot(1.e+4*prtclCoor[0,pointsTot-lengthArrowElc:pointsTot],1.e+4*prtclCoor[2,pointsTot-lengthArrowElc:pointsTot], \
          1.e+4*prtclCoor[4,pointsTot-lengthArrowElc:pointsTot],'-b',linewidth=2)
ax20.plot(1.e+4*prtclCoor[0,pointsTot-points:pointsTot-points+lengthArrowElc], \
          1.e+4*prtclCoor[2,pointsTot-points:pointsTot-points+lengthArrowElc], \
          1.e+4*prtclCoor[4,pointsTot-points:pointsTot-points+lengthArrowElc],'-b',linewidth=2)
plt.xlabel('x, $\mu m$',color='m',fontsize=16)
plt.ylabel('y, $\mu m$',color='m',fontsize=16)
ax20.set_zlabel('z, $\mu m$',color='m',fontsize=16)
plt.title(('First Electron Trajectory (End; $N_L=$%d):\nImpact Parameter=%5.2f $\mu$m, $R_L$=%5.2f $\mu$m' \
           % (larmorNumber[0],1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)),color='m',fontsize=16)

# Data to traw the arrow in the direction of motion:
#
arrwBegxIon=1.e+7*prtclCoor[6,pointsTot-points/2]
arrwBegyIon=1.e+7*prtclCoor[8,pointsTot-points/2]
arrwBegzIon=1.e+7*prtclCoor[10,pointsTot-points/2]
arrwEndxIon=1.e+7*prtclCoor[6,pointsTot-points/2-50]
arrwEndyIon=1.e+7*prtclCoor[8,pointsTot-points/2-50]
arrwEndzIon=1.e+7*prtclCoor[10,pointsTot-points/2-50]
fig25=plt.figure(25)
ax25=fig25.gca(projection='3d')
ax25.plot(1.e+7*prtclCoor[6,pointsTot-points:pointsTot],1.e+7*prtclCoor[8,pointsTot-points:pointsTot], \
          1.e+7*prtclCoor[10,pointsTot-points:pointsTot],'-b',linewidth=2)
ax25.plot([arrwBegxIon,arrwEndxIon],[arrwBegyIon,arrwEndyIon],[arrwBegzIon,arrwEndzIon],color='r',alpha=0.8,lw=2)
plt.xlabel('x, $nm$',color='m',fontsize=16)
plt.ylabel('y, $nm$',color='m',fontsize=16)
ax25.set_zlabel('z, $nm$',color='m',fontsize=16)
plt.title('First Ion Trajectory (End)',color='m',fontsize=16)

plt.figure(30)
plt.plot(range(pointsTot),1.e4*prtclCoor[4,0:pointsTot],'-r',linewidth=2)
plt.xlabel('Points',color='m',fontsize=16)
plt.ylabel('z, $\mu$m',color='m',fontsize=16)
plt.title(('First Trajectory ($N_L=$%d): Impact Parameter=%5.2f $\mu$m, $R_L$=%5.2f $\mu$m' % \
           (larmorNumber[0],1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)),color='m',fontsize=16)
plt.grid(True)

plt.figure(40)
plt.plot(1.e4*prtclCoor[4,0:pointsTot],1.e4*prtclCoor[12,0:pointsTot],'-r',linewidth=2)
plt.xlabel('z, $\mu$m',color='m',fontsize=16)
plt.ylabel('Distance from Motionless Ion, $\mu$m',color='m',fontsize=16)
plt.title(('First Trajectory ($N_L=$%d): Impact Parameter=%5.2f $\mu$m, $R_L$=%5.2f $\mu$m' % \
           (larmorNumber[0],1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)),color='m',fontsize=16)
plt.grid(True)

plt.figure(55)
plt.plot(uPot_enrgKin,larmR_b,'.r')
plt.xlim([minUpot_enrgKin_1-.1,maxUpot_enrgKin_1+.1])
plt.ylim([minLarmR_b_1-.1,maxLarmR_b_1+.1])
plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
plt.ylabel('$B=log_{10}(R_L/b)$',color='m',fontsize=16)
plt.title('Map for Transfered Momenta $dP_x,dP_y,dP_z$ Calculations', color='m',fontsize=20)
plt.grid(True)

X,Y=np.meshgrid(xA,yB) 
fig70=plt.figure(70)
ax70=fig70.gca(projection='3d')
# surf=ax70.plot_surface(X,Y,zApprch1dpx,cmap=cm.coolwarm,linewidth=0,antialiased=False)
surf=ax70.plot_surface(X,Y,zApprch1dpx,cmap=cm.jet,linewidth=0,antialiased=False)
plt.title('Transfered Momentum $dP_x$:\n$dP_x=q_e^2/b \cdot C_x$', color='m',fontsize=20)
plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
plt.ylabel('$B=log_{10}(R_L/b)$',color='m',fontsize=16)
ax70.set_zlabel('$C_x$, $cm^{-2}$',color='m',fontsize=16)
cb = fig70.colorbar(surf)
# cbar=fig70.colorbar(surf,ticks=[0,1000,2000,3000,4000,5000,6000,7000])  # Next 2 commands not work
# cbar.ax.set_yticklabels(['0','1000','2000','3000','4000','5000','6000','7000'])
# labels=np.arange(0,8000,1000)               # Next 4 commands not work
# location=labels
# cb.set_ticks(location)
# cb.set_ticklabels(labels)
# tick_locator = ticker.MaxNLocator(nbins=10) # Next 3 commands not work
# cb.locator = tick_locator
# cb.update_ticks()
plt.grid(True)


fig75=plt.figure(75)
ax=fig75.add_subplot(111)         # for contours poltting
X,Y=np.meshgrid(xA,yB) 
mapDpx=ax.contourf(X,Y,zApprch1dpx)   
# mapDpx=ax.contourf(X,Y,dpxApprch_1,levels)   
# Contourrange=[int(NlarmCutofDown+1)]
# mapTurnCutoff=ax.contour(X,Y,Nlarm,Contourrange,format='%d',colors=('w'),linewidths=(2)) 
# plt.clabel(mapTurnCutoff,inline=1,fontsize=14,manual=[(-3,-1.5)])  
plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
plt.ylabel('$B=log_{10}(R_L/b)$',color='m',fontsize=16)
plt.title('$C_x (cm^{-2})$ for Transf. Momntm. $dP_x$: $dP_x=q_e^2/b\cdot C_x$', color='m',fontsize=20)
fig75.colorbar(mapDpx)

'''
fig80=plt.figure(80)
ax80=fig80.gca(projection='3d')
surf=ax80.plot_surface(X,Y,zApprch1dpy,cmap=cm.coolwarm,linewidth=0,antialiased=False)
plt.title('Transfered Momentum $dP_y$:\n$dP_y=q_e^2/b \cdot C_y$', color='m',fontsize=20)
plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
plt.ylabel('$B=log_{10}(R_L/b)$',color='m',fontsize=16)
ax80.set_zlabel('$C_y$, $cm^{-2}$',color='m',fontsize=16)
fig80.colorbar(surf)
plt.grid(True)


fig90=plt.figure(90)
ax90=fig90.gca(projection='3d')
surf=ax90.plot_surface(X,Y,zApprch1dpz,cmap=cm.coolwarm,linewidth=0,antialiased=False)
plt.title('Transfered Momentum $dP_z$:\n$dP_z=q_e^2/b \cdot C_z$', color='m',fontsize=20)
plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
plt.ylabel('$B=log_{10}(R_L/b)$',color='m',fontsize=16)
ax90.set_zlabel('$C_z$, $cm^{-2}$',color='m',fontsize=16)
fig90.colorbar(surf)
plt.grid(True)
'''

plt.show()   
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

outfile.write ('\n                              First trajectory of electron (total points %d)' % pointsTot)
outfile.write ('\n      x, cm            y, cm            z, cm         px, g*cm/sec     py, g*cm/sec     pz, g*cm/sec\n\n')

for it in range(pointsTot):
   strLine=        '  {: 12.6e}'.format(prtclCoor[0,it])+',   '+'{: 12.6e}'.format(prtclCoor[2,it])+',   '+'{: 12.6e}'.format(prtclCoor[4,it])+',   '
   strLine=strLine+'{: 12.6e}'.format(prtclCoor[1,it])+',   '+'{: 12.6e}'.format(prtclCoor[3,it])+',   '+'{: 12.6e}'.format(prtclCoor[5,it])
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

outfile.write ('\n                                First trajectory of ion (total points %d)' % pointsTot)
outfile.write ('\n      x, cm            y, cm            z, cm         px, g*cm/sec     py, g*cm/sec     pz, g*cm/sec\n\n')

for it in range(pointsTot):
   strLine=        '  {: 12.6e}'.format(prtclCoor[6,it])+',   '+'{: 12.6e}'.format(prtclCoor[8,it])+',   '+'{: 12.6e}'.format(prtclCoor[10,it])+',   '
   strLine=strLine+'{: 12.6e}'.format(prtclCoor[7,it])+',   '+'{: 12.6e}'.format(prtclCoor[9,it])+',   '+'{: 12.6e}'.format(prtclCoor[11,it])
   outfile.write ('%s\n' % strLine)

outfile.close()
print 'Close the written output file "%s"' % outputFile

sys.exit()   

###################################################
#
# Writing the results to output file for later processing and vizualization: 
#
###################################################   
#
#
# Opening the output file: 
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
         if valPrev != 0. and valCrrnt == 0:
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
if k == 0:
   if nZero !=0:   
      zRes_line=zRes_line+'('+'{:d}'.format(nZero)+')'
      countsRes += nZero
      outfile.write ('%s\n' % zRes_line)
if k != 0:
   if nZero !=0:   
      zRes_line=zRes_line+'('+'{:d}'.format(nZero)+')'
      countsRes += nZero
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
         if valPrev != 0. and valCrrnt == 0:
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
if k == 0:
   if nZero !=0:   
      zRes_line=zRes_line+'('+'{:d}'.format(nZero)+')'
      countsRes += nZero
      outfile.write ('%s\n' % zRes_line)
if k != 0:
   if nZero !=0:   
      zRes_line=zRes_line+'('+'{:d}'.format(nZero)+')'
      countsRes += nZero
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
         if valPrev != 0. and valCrrnt == 0:
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
if k == 0:
   if nZero !=0:   
      zRes_line=zRes_line+'('+'{:d}'.format(nZero)+')'
      countsRes += nZero
      outfile.write ('%s\n' % zRes_line)
if k != 0:
   if nZero !=0:   
      zRes_line=zRes_line+'('+'{:d}'.format(nZero)+')'
      countsRes += nZero
      outfile.write ('%s\n' % zRes_line)

outfile.write ('\n')
            

print 'Total counts for zApprch1dpz data =%d' % countsRes
if countsRes != nBins*nBins:
   print 'Something wrong in the writing of zApprch1dpz data to output file!'

outfile.close()
print 'Close the written output file "%s"' % outputFile

# sys.exit()

###################################################
#
# Reading the data from output file for later processing and vizualization: 
#
###################################################   
#
#
# Opening the input file: 
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
# Header for xA-Data:
   if lines == 2: 
      words=lineData.split()
      indx1=words.index('Entries:')
      xAentries=int(words[indx1+1])
      indx2=words.index('with')
      perLine=int(words[indx2+1])
#       print 'xAdata-Header from %d: words =%s, index1=%d, entries=%d, index2=%d, perLine=%d' % (lines,words,indx1,xAentries,indx2,perLine)
      xAdata=np.zeros(xAentries)
      linesFull=xAentries//perLine
      entriesRem=xAentries-perLine*linesFull
      yBheaderLineNumber=linesFull+5
      if entriesRem > 0:
         yBheaderLineNumber += 1
#      print 'yBheaderLineNumber=%d' % yBheaderLineNumber
   if xAdataFlag == 0:
# xA-Data:
      if lines > 3 and lines <= yBheaderLineNumber-2:
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
	 
sys.exit()

