# from __future__ import division

#-------------------------------------
#
#        Started at 07/25/2017 (YuE)
# 
#-------------------------------------

import os, sys
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import matplotlib as mpl

from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D

import scipy.integrate as integrate
from scipy.integrate import quad, nquad, dblquad

from scipy.constants import pi
from scipy.constants import speed_of_light as clight
from scipy.constants import epsilon_0 as eps0
from scipy.constants import mu_0 as mu0
from scipy.constants import elementary_charge as qe
from scipy.constants import electron_mass as me
from scipy.constants import proton_mass as mp
from scipy.constants import Boltzmann as kB


eVtoErg=1.602e-12                # energy from eV to erg (from CI to CGS)
# Indices:
(Ix, Ipx, Iy, Ipy, Iz, Ipz) = range(6)

#
# Initial parameters:
#
Z_ion = qe*2.997e+9                                                # charge of ion (proton), CGSE units of the charge
M_ion = mp*1.e+3                                                   # mass of ion (proton), g
q_elec = qe*2.997e+9                                     # charge of electron, CGSE units of the charge (without sign!)
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
# For each electron z_e=(x_e,px_e,y_e,py_e,z_e,pz_e) --> zgc_e=(phi,p_phi,y_gc,py_gc,z_e,pz_e);
# z_c and zgc_e are 2D arrays with dimension (6,n_elec) 
#
def toGuidingCenter(z_e):
    mOmega=m_elec*omega_L                                          # g/sec
    zgc_e=z_e.copy()                                               # 2D array with dimension (6,n_elec)
    zgc_e[Ix,:] = np.arctan2(z_e[Ipx,:]+mOmega*z_e[Iy,:],z_e[Ipy,:])             # radians
    zgc_e[Ipx,:]= (((z_e[Ipx,:]+mOmega*z_e[Iy,:])**2+z_e[Ipy,:]**2)/(2.*mOmega)) # g*cm**2/sec
    zgc_e[Iy,:] =-z_e[Ipx,:]/mOmega                                # cm
    zgc_e[Ipy,:]= z_e[Ipy,:]+mOmega*z_e[Ix,:]                      # g/sec
    return zgc_e

#
# Convertion from guiding-center coordinates to electron's "coordinates":
# For each electron zgc_e=(phi,p_phi,y_gc,py_gc,z_e,pz_e) --> z_e=(x_e,px_e,y_e,py_e,z_e,pz_e);
# zgc_c and z_e are 2D arrays with dimension (6,n_elec) 
#
def fromGuidingCenter(zgc_e):
    mOmega=m_elec*omega_L                                          # g/sec
    rho_larm=np.sqrt(2.*zgc_e[Ipx,:]/mOmega)                       # cm
    z_e = zgc.copy()                                               # 2D array with dimension (6,n_elec)
    z_e[Ix,:] = zgc_e[Ipy,:]/mOmega-rho_larm*np.cos(zgc_e[Ix,:])   # cm
    z_e[Ipx,:]=-mOmega*zgc_e[Iy,:]                                 # g*cm/sec
    z_e[Iy,:] = zgc_e[Iy,:]+rho_larm*np.sin(zgc_e[Ix,:])           # cm
    z_e[Ipy,:]= mOmega*rho_larm*np.cos(zgc_e[Ix,:])                # g*cm/sec
    return z_e

#
# Matrix to dragg electron through the solenoid with field 'B_mag' during time interval 'deltaT':
#
def solenoid_eMatrix(B_mag,deltaT):
    import numpy as np   # Next np.identity does not work without that! Why?!
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
# Calculation of the ranges for the ratio larmR_b=R_larmor/b and for the 
# ratio eUpot_enrgKin=potential_energy/kinetic_energy of the electron
# as functions on ratio impctParMin_Rlarm=min_impact_parameter/R_larmor and 
#
minAlpha=-pi/2
maxAlpha=pi/2

stepsAlpha=100
stepAlpha=(maxAlpha-minAlpha)/stepsAlpha

minImpctParMin_Rlarm=1.1                                           # dimensionless
maxImpctParMin_Rlarm=5                                             # dimensionless
stepsImpctParMin_Rlarm=2
stepImpctParMin_Rlarm=(maxImpctParMin_Rlarm-minImpctParMin_Rlarm)/(stepsImpctParMin_Rlarm-1)    # dimensionless

minVlong=.1*eVrmsLong                                              # cm/sec
maxVlong=5.*eVrmsLong                                              # cm/sec
stepsVlong=100
stepVlong=(maxVlong-minVlong)/(stepsVlong-1)                       # cm/sec

minVtran=.1*eVrmsTran                                              # cm/sec
maxVtran=5.*eVrmsTran                                              # cm/sec
stepsVtran=2
stepVtran=(maxVtran-minVtran)/(stepsVtran-1)                       # cm/sec

z_elecCrrnt=np.zeros(6)                                            # 6-vector for electron
begTracks=np.zeros(stepsVtran*stepsImpctParMin_Rlarm+1)
endTracks=np.zeros(stepsVtran*stepsImpctParMin_Rlarm)
minLarmR_b=np.zeros(stepsImpctParMin_Rlarm)                        # dimensionless     
maxLarmR_b=np.zeros(stepsImpctParMin_Rlarm)                        # dimensionless
minUpot_enrgKin=np.zeros(stepsImpctParMin_Rlarm)                   # dimensionless      
maxUpot_enrgKin=np.zeros(stepsImpctParMin_Rlarm)                   # dimensionless      
'''
#========================== Obsoleted (begin) ===========================
#
# Case: the same transversal (eVrmsTran) and different longidudinal velocities: 
#
for i in range(1):    # stepsVlong):
   eVlong=maxVlong-stepVlong*i                                     # cm/sec
   kinEnergy=m_elec*(eVrmsTran**2+eVlong**2)/2.                    # erg
   for j in range(stepsImpctParMin_Rlarm):
      impctParMin_Rlarm=maxImpctParMin_Rlarm-stepImpctParMin_Rlarm*j         # dimensionless
      minImpctPar=impctParMin_Rlarm*ro_larmRMS                     # cm
      halfIntrcnLength=minImpctPar*math.tan(pi/2-stepAlpha)        # cm
      timePath=halfIntrcnLength/eVlong
      numbLarmor=int(timePath/T_larm)
      timePoints=numbLarmor*stepsNumberOnGyro
      timestep=timePath/timePoints                                 # sec
      matr_elec=solenoid_eMatrix(B_mag,timeStep)       # matrix for electron for half of the time step (magnetic field)
      print 'minImpctPar(mkm)=%e, halfIntrcnLength=%e, timePath(sec)=%e,  numbLarmor=%d, timePoints=%d, timestep=%e' % \
            (1.e+4*minImpctPar,1.e+4*halfIntrcnLength,timePath,numbLarmor,timePoints,timestep)
      matr_elec=solenoid_eMatrix(B_mag,timeStep)                   # matrix for electron for time step (magnetic field)
      if j == 0:
         elecCoor=np.zeros((3,timePoints,2))                       # points of trajectory; cm
# Distance from origin of the coordinate system to electron along the trajectory; cm:
         b=np.zeros((stepsImpctParMin_Rlarm,timePoints))           # cm
# log10 of two important ratios; dimensionless:
         larmR_b=np.zeros((stepsImpctParMin_Rlarm,timePoints))     # ratio R_larmor/b
         upot_enrgKin=np.zeros((stepsImpctParMin_Rlarm,timePoints))          # ratio potential_energy/kinetic_energy
      for m in range(6):
          z_elecCrrnt[m]=0.                                        # Initial zero-vector for electron
      z_elecCrrnt[Ix]=minImpctPar                                  # x, cm
      z_elecCrrnt[Iz]=-halfIntrcnLength                            # z, cm
      z_elecCrrnt[Ipy]=m_elec*eVrmsTran                            # py, g*cm/sec
      z_elecCrrnt[Ipz]=m_elec*eVlong                               # pz, g*cm/sec
      for k in range(timePoints):
 	  z_elecCrrnt=matr_elec.dot(z_elecCrrnt)
          for ic in (Ix,Iy,Iz):
             elecCoor[ic//2,pointTrack[j],j]=z_elecCrrnt[ic]            # cm
# Distance from origin of the coordinate system to electron along the trajectory; cm:
 	  b[j,k]=np.sqrt(elecCoor[0,pointTrack[j],j]**2+elecCoor[1,pointTrack[j],j]**2+elecCoor[2,pointTrack[j],j]**2)
# log10 of two important ratios; dimensionless:
	  larmR_b[j,k]=math.log10(ro_larmRMS/b[j,k])               # dimensionless 
	  upot_enrgKin[j,k]=math.log10((q_elec**2/b[j,k])/kinEnergy)         # dimensionless   
          pointTrack[j] += 1
      minLarmR_b[j]=min(larmR_b[j,:])      
      maxLarmR_b[j]=max(larmR_b[j,:])
      print 'For j=%d larmR_b: min=%e, max=%e' % (j,minLarmR_b[j],maxLarmR_b[j])      
      minUpot_enrgKin[j]=min(upot_enrgKin[j,:])
      maxUpot_enrgKin[j]=max(upot_enrgKin[j,:])
      print 'For j=%d upot_enrgKin: min=%e, max=%e' % (j,minUpot_enrgKin[j],maxUpot_enrgKin[j])      
#========================== Obsoleted (end) ===========================
'''
#===============================================================================
#
# Case: the same longidudinal (eVrmsLong) and different transversal velocities: 
#

'''
#
# To define the range for b (distance to the origin of coordinates):
#
roMin=1.e+8
roMax=0.
roLarmMin=1.e+8
roLarmMax=0.
halfIntrcnLengthMax=0.
timeStart=os.times()
for i in range(stepsVtran):
#    pointTrack=np.zeros(stepsImpctParMin_Rlarm)             
   eVtran=maxVtran-stepVtran*i                                     # cm/sec
   ro_larm=eVtran/omega_L                                          # cm
   if roLarmMin > ro_larm:
      roLarmMin=ro_larm
   if roLarmMax < ro_larm:
      roLarmMax=ro_larm	 
   kinEnergy=m_elec*(eVtran**2+eVrmsLong**2)/2.                    # erg
   for j in range(stepsImpctParMin_Rlarm):
      impctParMin_Rlarm=maxImpctParMin_Rlarm-stepImpctParMin_Rlarm*j         # dimensionless
      impctPar=impctParMin_Rlarm*ro_larm                        # cm
      if roMin > impctPar:
         roMin=impctPar
      if roMax < impctPar:
         roMax=impctPar	 
      halfIntrcnLength=impctPar*math.tan(pi/2-.25*stepAlpha)    # cm
      timePath=halfIntrcnLength/eVrmsLong
      numbLarmor=int(timePath/T_larm)
#       print 'impctPar(mkm)=%e, rhoLarmour (mkm)=%e, halfIntrcnLength=%e, numbLarmor=%d' % \
#             (1.e+4*impctPar,1.e+4*ro_larm,1.e+4*halfIntrcnLength,numbLarmor)
      if halfIntrcnLengthMax < halfIntrcnLength:
         halfIntrcnLengthMax=halfIntrcnLength
print 'roMin=%e, roMax=%e, roLarmMin=%e, roLarmMax=%e, halfIntrcnLengthMax=%e (all sizes in mkm)' % \
      (1.e+4*roMin,1.e+4*roMax,1.e+4*roLarmMin,1.e+4*roLarmMax,1.e+4*halfIntrcnLengthMax)
timeEnd=os.times()
cpuTime=1.e+6*float(timeEnd[0])-1.e+6*float(timeStart[0])   # CPU time , mks
print 'cpuTime = %e' % cpuTime

sys.exit()
'''
#
# Results of the commented above block:
#
roMin=2.623566e-4                                                  # cm
roMax=5.962650e-2                                                  # cm
roLarmMin=2.385060e-4                                              # cm
roLarmMax=1.192530e-2                                              # cm
halfIntrcnLengthMax=7.591726e+0                                    # cm

print 'roMin=%e, roMax=%e, roLarmMin=%e, roLarmMax=%e, halfIntrcnLengthMax=%e (all sizes in mkm)' % \
      (1.e+4*roMin,1.e+4*roMax,1.e+4*roLarmMin,1.e+4*roLarmMax,1.e+4*halfIntrcnLengthMax)

bMin=roMin-roLarmMin                                               # cm
bMax=.5*np.sqrt(roMax**2+halfIntrcnLengthMax**2)+roLarmMax            # cm

nTotal=2000
bEdges=np.zeros(nTotal+1)
bStep=(bMax-bMin)/nTotal                                           # cm

print 'bMin(mkm)=%e, bMax(mkm)=%e, bStep(mkm)=%e' % (1.e+4*bMin,1.e+4*bMax,1.e+4*bStep)

for i in range(nTotal+1):
   bEdges[i]=bMin+bStep*i                                          # cm
   
print 'bEdges[0]=%e, bEdges[nTotal]=%e (all sizes in mkm)' % (1.e+4*bEdges[0],1.e+4*bEdges[nTotal])
   
#
# All data will be placed in nTotal bins (for range of impact parameters):
#
bTot=np.zeros(nTotal)                                              # cm
larmR_bTot=np.zeros(nTotal)                                        # ratio R_larmor/b
uPot_enrgKinTot=np.zeros(nTotal)                                   # ratio eVca/Ekin
population=np.zeros(nTotal)                                        # number of data in each bin

minLarmorTurns=5000
mPrev=0
bPrev=0.
for i in range(stepsVtran):
   timeStart=os.times()
#    pointTrack=np.zeros(stepsImpctParMin_Rlarm)             
   eVtran=maxVtran-stepVtran*i                                     # cm/sec
   ro_larm=eVtran/omega_L                                          # cm
   kinEnergy=m_elec*(eVtran**2+eVrmsLong**2)/2.                    # erg
   for j in range(stepsImpctParMin_Rlarm):
      impctParMin_Rlarm=maxImpctParMin_Rlarm-stepImpctParMin_Rlarm*j         # cm
      minImpctPar=impctParMin_Rlarm*ro_larm                        # cm
      alphaStep=stepAlpha*np.random.uniform(low=0.25,high=1.,size=1)
      halfIntrcnLength=minImpctPar*math.tan(pi/2-alphaStep)        # cm
#      halfIntrcnLength=50.e-4                                      # cm
      timePath=halfIntrcnLength/eVrmsLong
      numbLarmor=int(timePath/T_larm)
      numbLarmorCalc=numbLarmor
      if numbLarmorCalc < minLarmorTurns:
         numbLarmorCalc=minLarmorTurns
      timePoints=numbLarmor*stepsNumberOnGyro
      timePoints=int(timePath/timeStep)
      pointTrack=np.zeros(timePoints)             
#      timeStep=timePath/timePoints                                 # sec
      if numbLarmorCalc == minLarmorTurns:
         halfIntrcnLength=timePoints*timeStep*eVrmsLong            # cm
      print 'minImpctPar(mkm)=%e, rhoLarmour (mkm)=%e, halfIntrcnLength (mkm)=%e, numbLarmor=%d (calculated %d)' % \
            (1.e+4*minImpctPar,1.e+4*ro_larm,1.e+4*halfIntrcnLength,numbLarmor,numbLarmorCalc)
      matr_elec=solenoid_eMatrix(B_mag,timeStep)                   # matrix for electron for timeStep in magnetic field
# To draw only first trajectory (for checking only):
      if i == 0 and j == 0:
         elecCoor=np.zeros((4,timePoints))                         # points of trajectory; cm
# Current distance from origin of the coordinate system to electron along the trajectory; cm:
      bCrrnt=np.zeros(timePoints)                                  # cm
# Current log10 of two important ratios; dimensionless:
      larmR_bCrrnt=np.zeros(timePoints)                            # ratio R_larmor/b
      uPot_enrgKinCrrnt=np.zeros(timePoints)                       # ratio potential_energy/kinetic_energy
      for m in range(6):
          z_elecCrrnt[m]=0.                                        # Initial zero-vector for electron
      z_elecCrrnt[Ix]=minImpctPar                                  # x, cm
      z_elecCrrnt[Iz]=-halfIntrcnLength                            # z, cm
      z_elecCrrnt[Ipy]=-m_elec*eVtran                              # py, g*cm/sec
      z_elecCrrnt[Ipz]=m_elec*eVrmsLong                            # pz, g*cm/sec
#-----------------------------------------------
# Main action - dragging of the current trajectory (for given i and j):
#
      for k in range(timePoints):
 	 z_elecCrrnt=matr_elec.dot(z_elecCrrnt)
# Current distance from origin of the coordinate system to electron along the trajectory; cm:
 	 bCrrnt[k]=np.sqrt(z_elecCrrnt[0]**2+z_elecCrrnt[2]**2+z_elecCrrnt[4]**2)
# Current log10 of two important ratios:  
	 larmR_bCrrnt[k]=math.log10(ro_larm/bCrrnt[k])             # dimensionless 
	 uPot_enrgKinCrrnt[k]=math.log10((q_elec**2/bCrrnt[k])/kinEnergy)      # dimensionless   
# To draw only first trajectory (for checking only):
         if i==0 and j==0:
            for ic in (Ix,Iy,Iz):
               elecCoor[ic//2,pointTrack[j]]=z_elecCrrnt[ic]       # cm
	       elecCoor[3,pointTrack[j]]=bCrrnt[k]                 # cm
#++++++++++++++++++++++++++++++++++
# "Depopulating":
#
         if bCrrnt[k] > bPrev:
 	    mBeg=mPrev-1
	    if mBeg < 0:
	       mBeg=0
  	    mEnd=nTotal
 	    mIncr=1
	 else:
 	    mBeg=mPrev+1
	    if mBeg > nTotal:
	       mBeg=nTotal
 	    mEnd=-1	
 	    mIncr=-1
         for m in range(mBeg,mEnd,mIncr):
	    if bEdges[m] < bCrrnt[k] <= bEdges[m+1]:
	       if population[m] == 0:
	          bTot[m]=bCrrnt[k]
	          larmR_bTot[m]=larmR_bCrrnt[k]
		  uPot_enrgKinTot[m]=uPot_enrgKinCrrnt[k]
	          mPrev=m
	          bPrev=bCrrnt[k]
	       population[m] += 1
	       break
#++++++++++++++++++++++++++++++++++
#  	 print 'i=%d, j=%d, k=%d: bPrev=%e (m=%d), bCrrnt=%e, x=%e, y=%e, z=%e' % \
#  	       (i,j,k,1.e+4*bPrev,mPrev,1.e+4*bCrrnt[k],1.e+4*z_elecCrrnt[0],1.e+4*z_elecCrrnt[2],1.e+4*z_elecCrrnt[4])       
         pointTrack[j] += 1
# End of gragging of the current trajectory	  
#-----------------------------------------------
      minLarmR_b[j]=min(larmR_bCrrnt)      
      maxLarmR_b[j]=max(larmR_bCrrnt)
#       print 'For j=%d larmR_bCrrnt: min=%e, max=%e' % (j,minLarmR_b[j],maxLarmR_b[j])      
      minUpot_enrgKin[j]=min(uPot_enrgKinCrrnt)
      maxUpot_enrgKin[j]=max(uPot_enrgKinCrrnt)
#       print 'For j=%d upot_enrgKinCrrnt: min=%e, max=%e' % (j,minUpot_enrgKin[j],maxUpot_enrgKin[j])      
#========================== Obsoleted (begin) ===========================
      if i == 0 and j == 0: 
# First definition of the total distance from origin of the coordinate system to electron along the trajectory; cm:
 	 b=bCrrnt
# First definition of the total log10 of two important ratios; dimensionless:  
	 larmR_b=larmR_bCrrnt                                     # dimensionless 
	 uPot_enrgKin=uPot_enrgKinCrrnt                           # dimensionless 
      else:  
# Total distance from origin of the coordinate system to electron along the trajectory:
 	 b=np.concatenate((b,bCrrnt),axis=0)                      # cm
# Total log10 of two important ratios; dimensionless :  
	 larmR_b=np.concatenate((larmR_b,larmR_bCrrnt),axis=0)    # dimensionless               
	 uPot_enrgKin=np.concatenate((uPot_enrgKin,uPot_enrgKinCrrnt),axis=0)          # cm
#========================== Obsoleted (end) ===========================
      trackNumb=stepsImpctParMin_Rlarm*i+j 
      endTracks[trackNumb]=begTracks[trackNumb]+pointTrack[j]     # index of the start of the current track     
      print 'i=%d, j=%d: beg=%d, end=%d' % (i,j,begTracks[trackNumb],endTracks[trackNumb])# index of the end of current track
      begTracks[trackNumb+1]=endTracks[trackNumb]+1      
bTotalPoints=b.shape
bDataSize=bTotalPoints[0]
print 'totalPoints: bDataSize=%d' % bDataSize

nonZeroRslts=0
for m in range(nTotal):
   if population[m] != 0:
      nonZeroRslts += 1
#       print '%d: left=%e, value(mkm)=%e (%d), right=%e' % (m,1.e+4*bEdges[m],1.e+4*bTot[m],population[m],1.e+4*bEdges[m+1])
print 'Summa of population=%d; nonZeroRslts=%d' % (sum(population),nonZeroRslts)

timeEnd=os.times()
cpuTime=1.e+6*float(timeEnd[0])-1.e+6*float(timeStart[0])   # CPU time , mks
print 'cpuTime(mksec) = %e' % cpuTime

#
# Checking of the first trajectory:
#
pointsTot=len(elecCoor[0,:])
turns=10                                                           # Number of larmorturns for drawing 
points=turns*stepsNumberOnGyro                                     # Number of points for drawing

fig10=plt.figure(10)
ax10=fig10.gca(projection='3d')
ax10.plot(1.e+4*elecCoor[0,0:points],1.e+4*elecCoor[1,0:points],1.e+4*elecCoor[2,0:points],'-r',linewidth=2)
plt.xlabel('x, $\mu m$',color='m',fontsize=16)
plt.ylabel('y, $\mu m$',color='m',fontsize=16)
plt.title(('First %d Larmor Turns of the First Trajectory' % turns),color='m',fontsize=16)

fig20=plt.figure(20)
ax20=fig20.gca(projection='3d')
ax20.plot(1.e+4*elecCoor[0,pointsTot-points:pointsTot],1.e+4*elecCoor[1,pointsTot-points:pointsTot], \
          1.e+4*elecCoor[2,pointsTot-points:pointsTot],'-r',linewidth=2)
plt.xlabel('x, $\mu m$',color='m',fontsize=16)
plt.ylabel('y, $\mu m$',color='m',fontsize=16)
plt.title(('Last %d Larmor Turns of the First Trajectory' % turns),color='m',fontsize=16)

# plt.show()   


plt.figure(30)
plt.plot(range(pointsTot),10*elecCoor[2,0:pointsTot],'-r',linewidth=2)
plt.xlabel('Points',color='m',fontsize=16)
plt.ylabel('z, mm',color='m',fontsize=16)
plt.title('First Trajectory',color='m',fontsize=16)
plt.grid(True)

plt.figure(40)
plt.plot(10*elecCoor[2,0:pointsTot],10*elecCoor[3,0:pointsTot],'-r',linewidth=2)
plt.xlabel('z, mm',color='m',fontsize=16)
plt.ylabel('Distance from Motionless Ion, mm',color='m',fontsize=16)
plt.title('First Trajectory',color='m',fontsize=16)
plt.grid(True)

plt.show()   

# sys.exit()

#
# Opening the output file: 
#
apprch1_file='apprch1_res.dat'
print 'Open output file "%s"...' % apprch1_file
apprch1_flag=0
try:
   outfile = open(apprch1_file,'w')
   apprch1_flag=1
except:
   print 'Problem to open output file "%s"' % apprch1_file

#
# Writing the results to output file: 
#
nInLine=10
totEntries=nonZeroRslts

outfile.write ('\n     Distance from Motionless Ion, mm ( Entries %d )\n\n' % totEntries)
k=0
for m in range(totEntries):   
   valCrrnt=1.e+3*bTot[m]
   strVal='{:f}'.format(valCrrnt)
   if k == 0:
      bLine=strVal
   else:
      bLine=bLine+', '+strVal
   k += 1
   if k == nInLine or m == totEntries-1:
      outfile.write ('%s\n' % bLine)
      k=0


outfile.write ('\n             Population ( Entries %d )\n\n' % totEntries)
k=0
for m in range(totEntries):   
   valCrrnt=int(population[m])
   strVal='{:d}'.format(valCrrnt)
   if k == 0:
      bLine=strVal
   else:
      bLine=bLine+', '+strVal
   k += 1
   if k == nInLine or m == totEntries-1:
      outfile.write ('%s\n' % bLine)
      k=0

outfile.write ('\n      Log10 of Ratio eVca/Ekin ( Entries %d )\n\n' % totEntries)
k=0
for m in range(totEntries):   
   valCrrnt=uPot_enrgKinTot[m]
   strVal='{:f}'.format(valCrrnt)
   if k == 0:
      bLine=strVal
   else:
      bLine=bLine+', '+strVal
   k += 1
   if k == nInLine or m == totEntries-1:
      outfile.write ('%s\n' % bLine)
      k=0

outfile.write ('\n       Log10 of Ratio R_larmor/b ( Entries %d )\n\n' % totEntries)
k=0
for m in range(totEntries):  
   valCrrnt=larmR_bTot[m]
   strVal='{:f}'.format(valCrrnt)
   if k == 0:
      bLine=strVal
   else:
      bLine=bLine+', '+strVal
   k += 1
   if k == nInLine or m == totEntries-1:
      outfile.write ('%s\n' % bLine)
      k=0

outfile.write ('\n        Approach1 Results:  deltaP ( Entries %d )\n\n' % totEntries)
k=0
for m in range(totEntries):   
   valCrrnt=larmR_b[m]
   strVal='{:f}'.format(valCrrnt)
   if k == 0:
      bLine=strVal
   else:
      bLine=bLine+', '+strVal
   k += 1
   if k == nInLine or m == totEntries-1:
      outfile.write ('%s\n' % bLine)
      k=0

outfile.write ('\n   Edges of the Bins, mm ( Entries %d )\n\n' % totEntries)
k=0
for m in range(totEntries+1):   
   valCrrnt=1.e+3*bEdges[m]
   strVal='{:f}'.format(valCrrnt)
   if k == 0:
      bLine=strVal
   else:
      bLine=bLine+', '+strVal
   k += 1
   if k == nInLine or m == totEntries-1:
      outfile.write ('%s\n' % bLine)
      k=0

outfile.close()
print 'Close the writtenn output file "%s"...' % apprch1_file

# sys.exit()

#
# Opening the input file: 
#
apprch1_file='apprch1_res.dat'
print 'Open input file "%s"...' % apprch1_file
apprch1_flag=0
try:
   inpfile = open(apprch1_file,'r')
   apprch1_flag=1
except:
   print 'Problem to open input file "%s"' % apprch1_file

#
# Reading the results from input file: 
#
lines=0                                                            # Number of current line from input file   
linesFull=0                                                        # Number of fully filled rows with each type of data
poplHeaderLineNumber=0                                             # Serial number of line with header for population-Data
xHeaderLineNumber=0                                                # Serial number of line with header for x-Data
yHeaderLineNumber=0                                                # Serial number of line with header for y-Data
dataHeaderLineNumber=0                                             # Serial number of line with header result-Data
edgesHeaderLineNumber=0                                            # Serial number of line with header edges-Data
lastLineNumber=0                                                   # Number of the last line 
distDataFlag=0                                                     # =1 when distance-Data already read
poplDataFlag=0                                                     # =1 when population-Data already read
xDataFlag=0                                                        # =1 when x-Data already read
yDataFlag=0                                                        # =1 when y-Data already read 
dataFlag=0                                                         # =1 when result-Data already read
edgesFlag=0                                                        # =1 when edges-Data already read
dataNumber=0                                                       # Number of current value of each type of Data
distData=[]                                                        # Array of distance-Data
poplData=[]                                                        # Array of population-Data
xData=[]                                                           # Array of x-Data
yData=[]                                                           # Array of y-Data
data=[]                                                            # Array of result-Data
edgesData=[]                                                       # Array of edges-Data

while True:
   lineData=inpfile.readline()
   if not lineData:
      break
   lines += 1
#    print '%d: %s' % (lines,lineData)
# Header for distance-Data:
   if lines == 2: 
      words=lineData.split()
      indx=words.index('Entries')
      entries=int(words[indx+1])
#       print 'distance-Header from %d: words =%s, index=%d, entries=%d' % (lines,words,indx,entries)
      distData=np.zeros(entries)
      poplData=np.zeros(entries)
      xData=np.zeros(entries)
      yData=np.zeros(entries)
      data=np.zeros(entries)
      edgesData=np.zeros(entries)
      linesFull=entries//10
      entriesRem=entries-10*linesFull
      poplHeaderLineNumber=linesFull+5
      if entriesRem > 0:
         poplHeaderLineNumber += 1
      print 'poplHeaderLineNumber=%d' % poplHeaderLineNumber
   if distDataFlag == 0:
# distance-Data:
      if lines > 3 and lines <= poplHeaderLineNumber-2:
         words=lineData.split()
         nWords=len(words)
#          print 'distance-Data from %d: words=%s, number of entries = %d' % (lines,words,nWords)
         for m in range(nWords):
            wordCrrnt=words[m].split(",")
            distData[dataNumber]=float(wordCrrnt[0])
	    dataNumber += 1
      if lines == poplHeaderLineNumber-2:
         distDataFlag=1   
         print 'distance-Data already read'  
# Header for population-Data:
   if lines == poplHeaderLineNumber:
      xHeaderLineNumber=poplHeaderLineNumber+linesFull+3
      if entriesRem > 0:
         xHeaderLineNumber += 1
      print 'xHeaderLineNumber=%d' % xHeaderLineNumber
      dataNumber=0
   if distDataFlag == 1 and poplDataFlag == 0:
# population-Data:
      if lines >  poplHeaderLineNumber+1 and lines <= xHeaderLineNumber-2:
         words=lineData.split()
         nWords=len(words)
#          print 'population-Data from %d: words=%s, number of entries = %d' % (lines,words,nWords)
         for m in range(nWords):
            wordCrrnt=words[m].split(",")
            poplData[dataNumber]=int(wordCrrnt[0])
#            print 'poplData(%d)=%d' % (dataNumber,poplData[dataNumber])
	    dataNumber += 1
      if lines == xHeaderLineNumber-2:
         poplDataFlag=1   
         print 'population-Data already read'  
# Header for x-Data:
   if lines == xHeaderLineNumber:
      yHeaderLineNumber=xHeaderLineNumber+linesFull+3
      if entriesRem > 0:
         yHeaderLineNumber += 1
      print 'yHeaderLineNumber=%d' % yHeaderLineNumber
      dataNumber=0
   if poplDataFlag == 1 and xDataFlag == 0:
# x-Data:
      if lines >  xHeaderLineNumber+1 and lines <= yHeaderLineNumber-2:
         words=lineData.split()
         nWords=len(words)
#          print 'x-Data from %d: words=%s, number of entries = %d' % (lines,words,nWords)
         for m in range(nWords):
            wordCrrnt=words[m].split(",")
            xData[dataNumber]=float(wordCrrnt[0])
	    dataNumber += 1
      if lines == yHeaderLineNumber-2:
         xDataFlag=1   
         print 'x-Data already read'  
# Header for y-Data:
   if lines == yHeaderLineNumber:
      dataHeaderLineNumber=yHeaderLineNumber+linesFull+3
      if entriesRem > 0:
         dataHeaderLineNumber += 1
      print 'dataHeaderLineNumber=%d' % dataHeaderLineNumber
      dataNumber=0
   if xDataFlag == 1 and yDataFlag == 0:
# y-Data:
      if lines >  yHeaderLineNumber+1 and lines <= dataHeaderLineNumber-2:
         words=lineData.split()
         nWords=len(words)
#          print 'y-Data from %d: words=%s, number of entries = %d' % (lines,words,nWords)
         for m in range(nWords):
            wordCrrnt=words[m].split(",")
            yData[dataNumber]=float(wordCrrnt[0])
	    dataNumber += 1
      if lines == dataHeaderLineNumber-2:
         yDataFlag=1   
         print 'y-Data already read'  
# Header for result-Data:
   if lines == dataHeaderLineNumber:
      edgesHeaderLineNumber=dataHeaderLineNumber+linesFull+3
      if entriesRem > 0:
         edgesHeaderLineNumber += 1
      print 'edgesHeaderLineNumber=%d' % edgesHeaderLineNumber
      dataNumber=0
   if yDataFlag == 1 and dataFlag == 0:    
# result-Data:
      if lines >  dataHeaderLineNumber+1 and lines <= edgesHeaderLineNumber-2:
         words=lineData.split()
         nWords=len(words)
#          print 'Data from %d: words=%s, number of entries = %d' % (lines,words,nWords)
         for m in range(nWords):
            wordCrrnt=words[m].split(",")
            data[dataNumber]=float(wordCrrnt[0])
	    dataNumber += 1
      if lines == edgesHeaderLineNumber-2:
         dataFlag=1   
         print 'Data already read'  
# Header for edges-Data:
   if lines == edgesHeaderLineNumber:
      lastLineNumber=edgesHeaderLineNumber+1+linesFull
      if entriesRem > 0:
         lastLineNumber += 1
      print 'lastLineNumber=%d' % lastLineNumber
      dataNumber=0
   if dataFlag == 1 and edgesFlag == 0:    
# edges-Data:
      if lines >  edgesHeaderLineNumber+1 and lines <= lastLineNumber:
         words=lineData.split()
         nWords=len(words)
#          print 'Data from %d: words=%s, number of entries = %d' % (lines,words,nWords)
         for m in range(nWords):
            wordCrrnt=words[m].split(",")
            data[dataNumber]=float(wordCrrnt[0])
	    dataNumber += 1
      if lines == lastLineNumber:
         edgesFlag=1   
         print 'Edges already read'  
inpfile.close()
print 'Close input file "%s"' % apprch1_file
print 'Length of any type of Data =%d' % entries     
	 
timeEnd=os.times()
cpuTime=1.e+6*float(timeEnd[0])-1.e+6*float(timeStart[0])   # CPU time , mks
print 'cpuTime = %e' % cpuTime

# 
#
#
'''
for i in range(stepsVtran):
   for j in range(stepsImpctParMin_Rlarm):
      trackNumb=stepsImpctParMin_Rlarm*i+j 
      if trackNumb < 10:
         plt.figure(trackNumb)
         plt.plot(uPot_enrgKin[begTracks[trackNumb]:endTracks[trackNumb]], \
                       larmR_b[begTracks[trackNumb]:endTracks[trackNumb]],'.r')
lastTrack=stepsImpctParMin_Rlarm*stepsVtran-1
plt.figure(100)
plt.plot(uPot_enrgKin[0:endTracks[lastTrack]],larmR_b[0:endTracks[lastTrack]],'.r')

plt.show()   
'''

#
# Checking of writing/reading data from the output file:
#
flagError=0
for m in range(entries):
   if population[m] != poplData[m]:
      print '%d: out=%d, inp=%d ==> delta=%d' % (m,population[m],poplData[m],population[m]-poplData[m])
      flagError=1  
if flagError==0:
   print 'No writing/reading errors'
sys.exit()
