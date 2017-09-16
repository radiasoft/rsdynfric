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

z_elecCrrnt=np.zeros(6)                                            # 6-vector for electron

#===============================================================================
#
# Case: the same longidudinal (eVrmsLong).
# The selection of main parameters is sescribed in the
# document "notesToChoiceCodeParameters.docx
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
max_bCrrnt=0

matr_elec=solenoid_eMatrix(B_mag,timeStep)        # matrix for electron for timeStep in magnetic field

#
# Electrons are magnetized for impact parameter >> rhoCrit:
#
rhoCrit=math.pow(q_elec**2/(m_elec*omega_L**2),1./3)               # cm
alpha=2*q_elec**2*omega_L/(m_elec*eVrmsLong**3)                    # dimensionless
print 'rhoCrit(mkm)=%f, alpha=%f' % (1.e+4*rhoCrit, alpha)

nA=35
crrntA=np.zeros(nA)
minA=-5.
maxA=0.
stepA=(maxA-minA)/(nA-1)

nB=35
crrntB=np.zeros(nB)
minB=-3.
maxB=-.5
stepB=(maxB-minB)/(nB-1)

pointTrack=np.zeros(nA*nB)  
larmorNumber=np.zeros(nA*nB)
cpuTime=np.zeros(nA*nB)

# minLarmR_b=np.zeros(stepsRho)                                      # dimensionless     
# maxLarmR_b=np.zeros(stepsRho)                                      # dimensionless
# minUpot_enrgKin=np.zeros(stepsRho)                                 # dimensionless      
# maxUpot_enrgKin=np.zeros(stepsRho)                                 # dimensionless      

## rhoLarm=np.zeros((nA,nB))
## dist=np.zeros((nA,nB))
## rho=np.zeros((nA,nB))
## Lint=np.zeros((nA,nB))
## Nlarm=np.zeros((nA,nB))

minLarmR_b=1.e8                                                   # dimensionless     
maxLarmR_b=-1.e8                                                    # dimensionless
minUpot_enrgKin=1.e8                                              # dimensionless      
maxUpot_enrgKin=-1.e8                                               # dimensionless      

cpuTimeTotal=0
trackNumb=-1                                        # Tracks will be enumerated from 0!
for iA in range(nA):
   crrntA[iA]=minA+stepA*iA
   for iB in range(nB):
      crrntB[iB]=minB+stepB*iB
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
         larmorNumber[trackNumb]=numbLarmor
	 print 'iA=%d, iB=%d, trackNumber=%d, numbLarmor=%d' % (iA,iB,trackNumb,numbLarmor)   
         timePoints=int(numbLarmor*stepsNumberOnGyro)              # dimensionless
# To draw only first trajectory (for checking only):
         if trackNumb == 0:
	    rhoFirstTurn=rhoCrrnt
	    rhoLarmorFirstTurn=rho_larm
# Points of the first trajectory (x,y,z for indices=0,1,2 and b for index=3; all values in cm):
            elecCoor=np.zeros((4,timePoints))                    
# Current distance from origin of the coordinate system to electron along the trajectory; cm
         bCrrnt=np.zeros(timePoints)                               # cm
# Current log10 of two important ratios; dimensionless:
         larmR_bCrrnt=np.zeros(timePoints)                         # ratio R_larmor/b; dimensionless
         uPot_enrgKinCrrnt=np.zeros(timePoints)                    # ratio potential_energy/kinetic_energy; dimensionless
         dpApprch_1Crrnt=np.zeros((3,timePoints))                  # 1/cm^2; deltaPapprch_1=q_e^2*timeStep*|dpApprch_1Crrnt|
#---------------------------------------------
         for m in range(6): 
            z_elecCrrnt[m]=0.                                      # Initial zero-vector for electron
         z_elecCrrnt[Ix]=rhoCrrnt                                  # x, cm
         z_elecCrrnt[Iz]=-halfLintr                                # z, cm
         z_elecCrrnt[Ipy]=m_elec*eVtran                            # py, g*cm/sec
         z_elecCrrnt[Ipz]=m_elec*eVrmsLong                         # pz, g*cm/sec
#-----------------------------------------------
# Main action - dragging of the current trajectory (for given i and j):
#
         for k in range(timePoints):
 	    z_elecCrrnt=matr_elec.dot(z_elecCrrnt)
# Current distance from origin of the coordinate system to electron along the trajectory; cm:
 	    bCrrnt[k]=np.sqrt(z_elecCrrnt[0]**2+z_elecCrrnt[2]**2+z_elecCrrnt[4]**2)
	    if max_bCrrnt < bCrrnt[k]:
	       max_bCrrnt=bCrrnt[k]
# Current log10 of two important ratios:  
	    larmR_bCrrnt[k]=math.log10(rho_larm/bCrrnt[k])         # dimensionless 
	    if maxLarmR_b < larmR_bCrrnt[k]:
	       maxLarmR_b=larmR_bCrrnt[k]
	    if minLarmR_b > larmR_bCrrnt[k]:
	       minLarmR_b=larmR_bCrrnt[k]
	    uPot_enrgKinCrrnt[k]=math.log10((q_elec**2/bCrrnt[k])/kinEnergy)         # dimensionless 
	    if maxUpot_enrgKin < uPot_enrgKinCrrnt[k]:
	       maxUpot_enrgKin=uPot_enrgKinCrrnt[k]
	    if minUpot_enrgKin > uPot_enrgKinCrrnt[k]:
	       minUpot_enrgKin=uPot_enrgKinCrrnt[k]
# Current values to calculate deltaPapprch_1=q_e^2*timeStep*|dpApprch_1Crrnt|:  
            for ic in range(3):
	       dpApprch_1Crrnt[ic,k]= abs(z_elecCrrnt[2*ic])/ bCrrnt[k]**2           # 1/cm^2;  
# To draw only first trajectory (for checking only):
            if trackNumb == 0:
               for ic in (Ix,Iy,Iz):
                  elecCoor[ic//2,pointTrack[trackNumb]]=z_elecCrrnt[ic]              # cm
	          elecCoor[3,pointTrack[trackNumb]]=bCrrnt[k]      # cm
###################################################
#
# Here was be placed block A): not in used now
#
###################################################   
#  	 print 'i=%d, j=%d, k=%d: bPrev=%e (m=%d), bCrrnt=%e, x=%e, y=%e, z=%e' % \
#  	       (i,j,k,1.e+4*bPrev,mPrev,1.e+4*bCrrnt[k],1.e+4*z_elecCrrnt[0],1.e+4*z_elecCrrnt[2],1.e+4*z_elecCrrnt[4])       
            pointTrack[trackNumb] += 1
# End of gragging of the current trajectory	  
# 
         if trackNumb == 0: 
# First definition of the total distance from origin of the coordinate system to electron along the trajectory; cm:
 	    b=bCrrnt
# First definition of the total log10 of two important ratios; dimensionless:  
	    larmR_b=larmR_bCrrnt                                # dimensionless 
	    uPot_enrgKin=uPot_enrgKinCrrnt                      # dimensionless 
# First definition of the values to calculate deltaPapprch_1=q_e^2*timeStep*|dpApprch_1Crrnt|:
            dpxApprch_1=dpApprch_1Crrnt[0,:]                    # 1/cm^2;  
            dpyApprch_1=dpApprch_1Crrnt[1,:]                    # 1/cm^2;  
            dpzApprch_1=dpApprch_1Crrnt[2,:]                    # 1/cm^2;  
         else:  
# Total distance from origin of the coordinate system to electron along the trajectory:
 	    b=np.concatenate((b,bCrrnt),axis=0)                 # cm
# Total log10 of two important ratios; dimensionless :  
	    larmR_b=np.concatenate((larmR_b,larmR_bCrrnt),axis=0)                  
	    uPot_enrgKin=np.concatenate((uPot_enrgKin,uPot_enrgKinCrrnt),axis=0)        
# Total values to calculate deltaPapprch_1=q_e^2*timeStep*|dpApprch_1Crrnt|:
            dpxApprch_1=np.concatenate((dpxApprch_1,dpApprch_1Crrnt[0,:]),axis=0)      # 1/cm^2;  
            dpyApprch_1=np.concatenate((dpyApprch_1,dpApprch_1Crrnt[1,:]),axis=0)      # 1/cm^2;  
            dpzApprch_1=np.concatenate((dpzApprch_1,dpApprch_1Crrnt[2,:]),axis=0)      # 1/cm^2;  
#          print 'trackNumb:%d: shapes: b=%d, larmR_b=%d, uPot=%d, dpx=%d, dpy=%d, dpz=%d' % \
# 	          (trackNumb,b.shape[0],larmR_b.shape[0],uPot_enrgKin.shape[0], \
# 		   dpxApprch_1.shape[0],dpyApprch_1.shape[0],dpzApprch_1.shape[0])  
#
#--------------------------------------------------------------------
         timeEnd=os.times()
	 cpuTime[trackNumb]=1.e+6*(float(timeEnd[0])-float(timeStart[0]))       # CPU time , mks
	 cpuTimeTotal += cpuTime[trackNumb] 
         lastTrackNumber=trackNumb+1                                  # quantity of tracks = trackNumber + 1!     
         sumPoints += pointTrack[trackNumb]

# for i in range(lastTrackNumber):
#    print 'Track %d: larmor turns=%d, cpuTime(mks)=%e, time per turn(mks)=%6.1f' % \
#          (i,larmorNumber[i],cpuTime[i],cpuTime[i]/larmorNumber[i])
print 'For %d tracks number of points is %d' % (lastTrackNumber,sumPoints)

###################################################
#
# Here was be placed block B): not in used now
#
###################################################   
print 'cpuTimeTotal(mksec) = %e' % cpuTimeTotal

nBins=80

# xA=np.array(1.*np.zeros(nBins))
xA=np.zeros(nBins)
xAedges=np.zeros(nBins+1)
xAnumb=np.zeros(nBins)
xAstep=(maxUpot_enrgKin-minUpot_enrgKin)/nBins
for i in range(nBins+1):
   xAedges[i]=minUpot_enrgKin+xAstep*i

# yB=np.array(1.*np.zeros(nBins))
yB=np.zeros(nBins)
yBedges=np.zeros(nBins+1)
yBnumb=np.zeros(nBins)
yBstep=(maxLarmR_b-minLarmR_b)/nBins
for i in range(nBins+1):
   yBedges[i]=minLarmR_b+yBstep*i

# zApprch1dpx=np.array(np.zeros((nBins,nBins)))
zApprch1dpx=np.zeros((nBins,nBins))
zApprch1dpxNumb=np.zeros((nBins,nBins))

timeStart=os.times()

for nPoint in range(int(sumPoints)):
   for iA in range(nBins):
      searchAflag=0
      if (xAedges[iA] <= uPot_enrgKin[nPoint] < xAedges[iA+1]):
         if xAnumb[iA] == 0:
	    xA[iA]=uPot_enrgKin[nPoint]
	 else:
	    xA[iA]=(xA[iA]*xAnumb[iA]+uPot_enrgKin[nPoint])/(xAnumb[iA]+1)          # averaging
         xAnumb[iA] += 1
	 searchAflag=1
	 break
   if searchAflag == 0:
      xA[nBins-1]=(xA[nBins-1]*xAnumb[nBins-1]+uPot_enrgKin[nPoint])/(xAnumb[nBins-1]+1)          # averaging 
      xAnumb[nBins-1] += 1
   for iB in range(nBins):
      searchBflag=0
      if (yBedges[iB] <= larmR_b[nPoint] < yBedges[iB+1]):
         if yBnumb[iB] == 0:
	    yB[iB]=larmR_b[nPoint]
	 else:
	    yB[iB]=(yB[iB]*yBnumb[iB]+larmR_b[nPoint])/(yBnumb[iB]+1)               # averaging  
         yBnumb[iB] += 1
	 searchBflag=1
	 break
   if searchBflag == 0:
      yB[nBins-1]=(yB[nBins-1]*yBnumb[nBins-1]+larmR_b[nPoint])/(yBnumb[nBins-1]+1) # averaging 
      yBnumb[nBins-1] += 1
   if zApprch1dpxNumb[iA,iB] == 0:
      zApprch1dpx[iA,iB]=dpxApprch_1[nPoint]
   else:
      zApprch1dpx[iA,iB]= \
      (zApprch1dpx[iA,iB]*zApprch1dpxNumb[iA,iB]+dpxApprch_1[nPoint])/(zApprch1dpxNumb[iA,iB]+1)  # averaging
   zApprch1dpxNumb[iA,iB] += 1

# Some checkings:
sumAnumb=0
sumBnumb=0
sumCnumb=0
for iA in range(nBins):
   sumAnumb += xAnumb[iA]
   sumBnumb += yBnumb[iA]
   for iB in range(nBins):
      sumCnumb += zApprch1dpxNumb[iA,iB]
#      print 'iA=%d, iB=%d: xA=%f (%d), xB=%f (%d), zApprch1dpx=%e' % \
#            (iA,iB,xA[iA],xAnumb[iA],yB[iB],yBnumb[iB],zApprch1dpx[iA,iB])
   
print 'sumA=%d, sumB=%d, sumC=%d; sumPoints=%d' % (sumAnumb,sumBnumb,sumCnumb,int(sumPoints)) 

timeEnd=os.times()
runTime=1.e+6*(float(timeEnd[0])-float(timeStart[0]))       # CPU time , mks
print 'runTime(mksec) = %e' % runTime

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
ax10.set_zlabel('z, $\mu m$',color='m',fontsize=16)
plt.title(('First Trajectory (Start; $N_L=$%d):\nImpact Parameter=%5.2f $\mu$m, $R_L$=%5.2f $\mu$m' \
           % (larmorNumber[0],1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)),color='m',fontsize=16)

fig20=plt.figure(20)
ax20=fig20.gca(projection='3d')
ax20.plot(1.e+4*elecCoor[0,pointsTot-points:pointsTot],1.e+4*elecCoor[1,pointsTot-points:pointsTot], \
          1.e+4*elecCoor[2,pointsTot-points:pointsTot],'-r',linewidth=2)
plt.xlabel('x, $\mu m$',color='m',fontsize=16)
plt.ylabel('y, $\mu m$',color='m',fontsize=16)
ax20.set_zlabel('z, $\mu m$',color='m',fontsize=16)
plt.title(('First Trajectory (End; $N_L=$%d):\nImpact Parameter=%5.2f $\mu$m, $R_L$=%5.2f $\mu$m' \
           % (larmorNumber[0],1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)),color='m',fontsize=16)

plt.figure(30)
plt.plot(range(pointsTot),1.e4*elecCoor[2,0:pointsTot],'-r',linewidth=2)
plt.xlabel('Points',color='m',fontsize=16)
plt.ylabel('z, $\mu$m',color='m',fontsize=16)
plt.title(('First Trajectory ($N_L=$%d): Impact Parameter=%5.2f $\mu$m, $R_L$=%5.2f $\mu$m' % \
           (larmorNumber[0],1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)),color='m',fontsize=16)
plt.grid(True)

plt.figure(40)
plt.plot(1.e4*elecCoor[2,0:pointsTot],1.e4*elecCoor[3,0:pointsTot],'-r',linewidth=2)
plt.xlabel('z, $\mu$m',color='m',fontsize=16)
plt.ylabel('Distance from Motionless Ion, $\mu$m',color='m',fontsize=16)
plt.title(('First Trajectory ($N_L=$%d): Impact Parameter=%5.2f $\mu$m, $R_L$=%5.2f $\mu$m' % \
           (larmorNumber[0],1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)),color='m',fontsize=16)
plt.grid(True)

# plt.show()   

# sys.exit()

plt.figure(55)
plt.plot(uPot_enrgKin,larmR_b,'.r')
plt.xlim([minUpot_enrgKin-.1,maxUpot_enrgKin+.1])
plt.ylim([minLarmR_b-.1,maxLarmR_b+.1])
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
plt.title('$C_x (cm^{-2})$ for Transf. Momntm $dP_x$: $dP_x=q_e^2/b\cdot C_x$', color='m',fontsize=20)
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

sys.exit()

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

outfile.write ('\n     Distance from Motionless Ion, mkm ( Entries %d )\n\n' % totEntries)

k=0
for m in range(nTotal): 
   if population[m] != 0:  
      valCrrnt=1.e+4*bTot[m]
      strVal='{:f}'.format(valCrrnt)
      if k == 0:
         bLine=strVal
      else:
         bLine=bLine+', '+strVal
      k += 1
      if k == nInLine:
         outfile.write ('%s\n' % bLine)
         k=0
if k != 0:
   outfile.write ('%s\n' % bLine)

outfile.write ('\n             Population ( Entries %d )\n\n' % totEntries)

k=0
for m in range(nTotal):   
   if population[m] != 0:  
      valCrrnt=int(population[m])
      strVal='{:d}'.format(valCrrnt)
      if k == 0:
         bLine=strVal
      else:
         bLine=bLine+', '+strVal
      k += 1
      if k == nInLine:
         outfile.write ('%s\n' % bLine)
         k=0
if k != 0:
   outfile.write ('%s\n' % bLine)

outfile.write ('\n      Log10 of Ratio eVca/Ekin ( Entries %d )\n\n' % totEntries)

k=0
for m in range(nTotal):   
   if population[m] != 0:  
      valCrrnt=uPot_enrgKinTot[m]
      strVal='{:f}'.format(valCrrnt)
      if k == 0:
         bLine=strVal
      else:
         bLine=bLine+', '+strVal
      k += 1
      if k == nInLine:
         outfile.write ('%s\n' % bLine)
         k=0
if k != 0:
   outfile.write ('%s\n' % bLine)

outfile.write ('\n       Log10 of Ratio R_larmor/b ( Entries %d )\n\n' % totEntries)

k=0
for m in range(nTotal):  
   if population[m] != 0:  
      valCrrnt=larmR_bTot[m]
      strVal='{:f}'.format(valCrrnt)
      if k == 0:
         bLine=strVal
      else:
         bLine=bLine+', '+strVal
      k += 1
      if k == nInLine:
         outfile.write ('%s\n' % bLine)
         k=0
if k != 0:
   outfile.write ('%s\n' % bLine)

outfile.write ('\n        Approach1 Results:  deltaPx ( Entries %d )\n\n' % totEntries)

k=0
for m in range(nTotal):   
   if population[m] != 0:  
      valCrrnt=dpApprch_1Tot[0,m]
      strVal='{:e}'.format(valCrrnt)
      if k == 0:
         bLine=strVal
      else:
         bLine=bLine+', '+strVal
      k += 1
      if k == nInLine:
         outfile.write ('%s\n' % bLine)
         k=0
if k != 0:
   outfile.write ('%s\n' % bLine)

outfile.write ('\n        Approach1 Results:  deltaPy ( Entries %d )\n\n' % totEntries)

k=0
for m in range(nTotal):   
   if population[m] != 0:  
      valCrrnt=dpApprch_1Tot[1,m]
      strVal='{:e}'.format(valCrrnt)
      if k == 0:
         bLine=strVal
      else:
         bLine=bLine+', '+strVal
      k += 1
      if k == nInLine:
         outfile.write ('%s\n' % bLine)
         k=0
if k != 0:
   outfile.write ('%s\n' % bLine)

outfile.write ('\n        Approach1 Results:  deltaPz ( Entries %d )\n\n' % totEntries)

k=0
for m in range(nTotal):   
   if population[m] != 0:  
      valCrrnt=dpApprch_1Tot[2,m]
      strVal='{:e}'.format(valCrrnt)
      if k == 0:
         bLine=strVal
      else:
         bLine=bLine+', '+strVal
      k += 1
      if k == nInLine:
         outfile.write ('%s\n' % bLine)
         k=0
if k != 0:
   outfile.write ('%s\n' % bLine)

outfile.write ('\n   Edges of the Bins, mkm ( Entries %d )\n\n' % totEntries)

k=0
for m in range(nTotal):   
   if population[m] != 0:  
      valCrrnt=1.e+4*bEdges[m]
      strVal='{:f}'.format(valCrrnt)
      if k == 0:
         bLine=strVal
      else:
         bLine=bLine+', '+strVal
      k += 1
      if k == nInLine:
         outfile.write ('%s\n' % bLine)
         k=0
if k != 0:
   outfile.write ('%s\n' % bLine)

outfile.close()
print 'Close the written output file "%s"...' % apprch1_file

# sys.exit()

#####################################################
#
# Reading the output file for checking:
#
#####################################################
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
poplHeaderLineNumber=0                                             # Serial number of line with header for population-Data
xHeaderLineNumber=0                                                # Serial number of line with header for x-Data
yHeaderLineNumber=0                                                # Serial number of line with header for y-Data
dataDpxHeaderLineNumber=0                                          # Serial number of line with header for dataDpx
dataDpyHeaderLineNumber=0                                          # Serial number of line with header for dataDpy
dataDpzHeaderLineNumber=0                                          # Serial number of line with header for dataDpz
edgesHeaderLineNumber=0                                            # Serial number of line with header edges-Data
lastLineNumber=0                                                   # Number of the last line 
distDataFlag=0                                                     # =1 when distance-Data already read
poplDataFlag=0                                                     # =1 when population-Data already read
xDataFlag=0                                                        # =1 when x-Data already read
yDataFlag=0                                                        # =1 when y-Data already read 
dataDpxFlag=0                                                      # =1 when dataDpx already read
dataDpyFlag=0                                                      # =1 when dataDpy already read
dataDpzFlag=0                                                      # =1 when dataDpz already read
edgesFlag=0                                                        # =1 when edges-Data already read
distData=[]                                                        # Array of distance-Data
poplData=[]                                                        # Array of population-Data
xData=[]                                                           # Array of x-Data
yData=[]                                                           # Array of y-Data
dataDpx=[]                                                         # Array of dataDpx
dataDpy=[]                                                         # Array of dataDpy
dataDpz=[]                                                         # Array of dataDpz
edgesData=[]                                                       # Array of edges-Data

lines=0                                                            # Number of current line from input file   
linesFull=0                                                        # Number of fully filled rows with each type of data
dataNumber=0                                                       # Number of current value of any types of Data
while True:
   lineData=inpfile.readline()
   if not lineData:
      break
   lines += 1
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
      dataDpx=np.zeros(entries)
      dataDpy=np.zeros(entries)
      dataDpz=np.zeros(entries)
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
      dataDpxHeaderLineNumber=yHeaderLineNumber+linesFull+3
      if entriesRem > 0:
         dataDpxHeaderLineNumber += 1
      print 'dataDpxHeaderLineNumber=%d' % dataDpxHeaderLineNumber
      dataNumber=0
   if xDataFlag == 1 and yDataFlag == 0:
# y-Data:
      if lines >  yHeaderLineNumber+1 and lines <= dataDpxHeaderLineNumber-2:
         words=lineData.split()
         nWords=len(words)
#          print 'y-Data from %d: words=%s, number of entries = %d' % (lines,words,nWords)
         for m in range(nWords):
            wordCrrnt=words[m].split(",")
            yData[dataNumber]=float(wordCrrnt[0])
	    dataNumber += 1
      if lines == dataDpxHeaderLineNumber-2:
         yDataFlag=1   
         print 'y-Data already read'  
# Header for dataDpx:
   if lines == dataDpxHeaderLineNumber:
      dataDpyHeaderLineNumber=dataDpxHeaderLineNumber+linesFull+3
      if entriesRem > 0:
         dataDpyHeaderLineNumber += 1
      print 'dataDpysHeaderLineNumber=%d' % dataDpyHeaderLineNumber
      dataNumber=0
   if yDataFlag == 1 and dataDpxFlag == 0:    
# dataDpx:
      if lines >  dataDpxHeaderLineNumber+1 and lines <= dataDpyHeaderLineNumber-2:
         words=lineData.split()
         nWords=len(words)
#          print 'Data from %d: words=%s, number of entries = %d' % (lines,words,nWords)
         for m in range(nWords):
            dataDpx[dataNumber]=float(wordCrrnt[0])
	    dataNumber += 1
      if lines == dataDpyHeaderLineNumber-2:
         dataDpxFlag=1   
         print 'dataDpx already read'  
# Header for dataDpy:
   if lines == dataDpyHeaderLineNumber:
      dataDpzHeaderLineNumber=dataDpyHeaderLineNumber+linesFull+3
      if entriesRem > 0:
         dataDpzHeaderLineNumber += 1
      print 'dataDpzHeaderLineNumber=%d' % dataDpzHeaderLineNumber
      dataNumber=0
   if dataDpxFlag == 1 and dataDpyFlag == 0:    
# dataDpy:
      if lines >  dataDpyHeaderLineNumber+1 and lines <= dataDpzHeaderLineNumber-2:
         words=lineData.split()
         nWords=len(words)
#          print 'Data from %d: words=%s, number of entries = %d' % (lines,words,nWords)
         for m in range(nWords):
            wordCrrnt=words[m].split(",")
            dataDpy[dataNumber]=float(wordCrrnt[0])
	    dataNumber += 1
      if lines == dataDpzHeaderLineNumber-2:
         dataDpyFlag=1   
         print 'dataDpy already read'  
# Header for dataDpz:
   if lines == dataDpzHeaderLineNumber:
      edgesHeaderLineNumber=dataDpzHeaderLineNumber+linesFull+3
      if entriesRem > 0:
         edgesHeaderLineNumber += 1
      print 'edgesHeaderLineNumber=%d' % edgesHeaderLineNumber
      dataNumber=0
   if dataDpyFlag == 1 and dataDpzFlag == 0:    
# dataDpz:
      if lines >  dataDpzHeaderLineNumber+1 and lines <= edgesHeaderLineNumber-2:
         words=lineData.split()
         nWords=len(words)
#          print 'Data from %d: words=%s, number of entries = %d' % (lines,words,nWords)
         for m in range(nWords):
            wordCrrnt=words[m].split(",")
            dataDpz[dataNumber]=float(wordCrrnt[0])
	    dataNumber += 1
      if lines == edgesHeaderLineNumber-2:
         dataDpzFlag=1   
         print 'dataDpz already read'  
# Header for edges-Data:
   if lines == edgesHeaderLineNumber:
      lastLineNumber=edgesHeaderLineNumber+1+linesFull
      if entriesRem > 0:
         lastLineNumber += 1
      print 'lastLineNumber=%d' % lastLineNumber
      dataNumber=0
   if dataDpzFlag == 1 and edgesFlag == 0:    
# edges-Data:
      if lines >  edgesHeaderLineNumber+1 and lines <= lastLineNumber:
         words=lineData.split()
         nWords=len(words)
#          print 'Data from %d: words=%s, number of entries = %d' % (lines,words,nWords)
         for m in range(nWords):
            wordCrrnt=words[m].split(",")
            edgesData[dataNumber]=float(wordCrrnt[0])
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
# Checking of writing/reading data from the output file:
#
flagError=0
k=0
for m in range(nTotal):
   if population[m] != 0:
      if population[m] != poplData[k]:
#         print 'm=%d, k=%d: out=%d, inp=%d ==> delta=%d' % (m,k,population[m],poplData[k],population[m]-poplData[k])
         flagError=1  
      k +=1
if flagError==0:
   print 'No writing/reading errors'

sys.exit()

#++++++++++++++++++++++++++++++++++
###################################################
#
# Block A): not in used now
#
###################################################   
#++++++++++++++++++++++++++++++++++
# "Depopulating":
#
#            if bCrrnt[k] > bPrev:
# 	       mBeg=mPrev-1
#	       if mBeg < 0:
#	          mBeg=0
#  	       mEnd=nTotal
# 	       mIncr=1
#	    else:
# 	       mBeg=mPrev+1
#	       if mBeg > nTotal:
#	          mBeg=nTotal
# 	       mEnd=-1	
# 	       mIncr=-1
#            for m in range(mBeg,mEnd,mIncr):
#	       if bEdges[m] < bCrrnt[k] <= bEdges[m+1]:
#	          if population[m] == 0:
#	             bTot[m]=bCrrnt[k]                             # cm
#	             larmR_bTot[m]=larmR_bCrrnt[k]                 # dimensionless
#		     uPot_enrgKinTot[m]=uPot_enrgKinCrrnt[k]       # dimensionless
#                     for ic in range(3):
#	                dpApprch_1Tot[ic,m]=dpApprch_1Crrnt[ic,k]  # 1/cm^2;  
#	             mPrev=m
#	             bPrev=bCrrnt[k]
#		     depopFlag=1
#	          population[m] += 1
#	          break
#++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++
###################################################
#
# Nlock B): not in used now
#
###################################################   
#nonZeroRslts=0
#for m in range(nTotal):
#   if population[m] != 0:
#      nonZeroRslts += 1
#
#print 'Summa of population=%d; nonZeroRslts=%d' % (sum(population),nonZeroRslts)
#++++++++++++++++++++++++++++++++++



###################################################
#
# These plots can be made, when Blocks A) and B) in in used!
#
###################################################   

# plt.figure(50)
# plt.plot(uPot_enrgKinTot,larmR_bTot,'.r')
# plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
# plt.ylabel('$B=log_{10}(R_L/b)$',color='m',fontsize=16)
# plt.title('Map for Transfered Momenta $dP_x,dP_y,dP_z$ Calculations', color='m',fontsize=20)
# plt.grid(True)

# X,Y=np.meshgrid(uPot_enrgKinTot,larmR_bTot)      
# fig60=plt.figure(60)
# ax60=fig60.gca(projection='3d')
# surf=ax60.plot_surface(X,Y,dpApprch_1Tot[0,:],cmap=cm.coolwarm,linewidth=0,antialiased=False)
# plt.title('Transfered Momentum $dP_x$:\n$dP_x=q_e^2/b \cdot C_x$', color='m',fontsize=20)
# plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
# plt.ylabel('$B=log_{10}(R_L/b)$',color='m',fontsize=16)
# ax60.set_zlabel('$C_x$, $cm^{-2}$',color='m',fontsize=16)
# fig60.colorbar(surf)
# plt.grid(True)

