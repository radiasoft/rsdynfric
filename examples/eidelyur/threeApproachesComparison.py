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
# Case: the same longidudinal (eVrmsLong) and different transversal velocities: 
#

rhoMin=1.e-4                                                       # R_crit=1 mkm; cm
rShield=100.e-4                                                    # max value of the R_shield=50 mkm; cm
rhoMax=rShield                                                     # cm
stepsRho=5                                                         # 33
stepRho=(rhoMax-rhoMin)/(stepsRho-1)                               # cm

stepsVtran=5                                                       # 33

pointTrack=np.zeros((stepsVtran-1)*(stepsRho-1))  
larmorNumber=np.zeros((stepsVtran-1)*(stepsRho-1))
# minLarmR_b=np.zeros(stepsRho)                                      # dimensionless     
# maxLarmR_b=np.zeros(stepsRho)                                      # dimensionless
# minUpot_enrgKin=np.zeros(stepsRho)                                 # dimensionless      
# maxUpot_enrgKin=np.zeros(stepsRho)                                 # dimensionless      

bMin=rhoMin                                                        # cm
bMax=2*rhoMax                                                      # cm (In case with L_int=R_shield sqrt(2) instead 2 is  needed)

nTotal=500
bEdges=np.zeros(nTotal+1)
bStep=(bMax-bMin)/nTotal                                           # cm

print 'bMin(mkm)=%e, bMax(mkm)=%e, bStep(mkm)=%e' % (1.e+4*bMin,1.e+4*bMax,1.e+4*bStep)

for i in range(nTotal+1):
   bEdges[i]=bMin+bStep*i                                          # cm
   
print 'bEdges[0]=%e, bEdges[nTotal]=%e (all sizes in mkm)' % (1.e+4*bEdges[0],1.e+4*bEdges[nTotal])

'''
minRhoLarm=1.e+8
maxRhoLarm=0.

minHalfLintr=1.e+8
maxHalfLintr=0.

minNumbLarm=100000000
maxNumbLarm=0
'''
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
dpApprch_1Tot=np.zeros((3,nTotal))                                   # 1/cm^2
dpApprch_2Tot=np.zeros((3,nTotal))                                   # 1/cm^2
dpApprch_3Tot=np.zeros((3,nTotal))                                   # 1/cm^2

###################################################
#
# Here was be placed block A): not in used now
#
###################################################   

mPrev=0
bPrev=0.
sumPoints=0
timeStart=os.times()
max_bCrrnt=0

matr_elec=solenoid_eMatrix(B_mag,timeStep)                         # matrix for electron for timeStep in magnetic field
for j in range(1,stepsRho):
   rhoCrrnt=rhoMin+stepRho*j                                       # cm
   eVtranMax=(rhoCrrnt-rhoCrit)*omega_L                            # cm/sec
   stepVtran=eVtranMax/(stepsVtran-1)                              # cm/sec
   for i in range(1,stepsVtran):                                   # 
      trackNumb=(stepsVtran-1)*(j-1)+(i-1)                         # Tracks are enumerated from 0!     
      eVtran=stepVtran*i                                           # cm/sec
      rho_larm=eVtran/omega_L                                      # cm
      kinEnergy=m_elec*(eVtran**2+eVrmsLong**2)/2.                 # erg
      halfLintr=np.sqrt(rShield**2+rhoCrrnt**2)                    # cm
      timePath=halfLintr/eVrmsLong                                 # sec
      numbLarmor=int(timePath/T_larm)
      larmorNumber[trackNumb]=numbLarmor
      timePoints=int(numbLarmor*stepsNumberOnGyro)
# To draw only first trajectory (for checking only):
      if i == 1 and j == 1:
	 rhoFirstTurn=rhoCrrnt
	 rhoLarmorFirstTurn=rho_larm
#---------------------------------------------
# Definition of the current arrays:
#	 
# Points of the first trajectory (x,y,z for first index=0,1,2 and b for first index=3; cm):
         elecCoor=np.zeros((4,timePoints))                    
# Current distance from origin of the coordinate system to electron along the trajectory; cm
      bCrrnt=np.zeros(timePoints)                                  # cm
# Current log10 of two important ratios; dimensionless:
      larmR_bCrrnt=np.zeros(timePoints)                            # ratio R_larmor/b; dimensionless
      uPot_enrgKinCrrnt=np.zeros(timePoints)                       # ratio potential_energy/kinetic_energy; dimensionless
      dpApprch_1Crrnt=np.zeros((3,timePoints))                     # 1/cm^2; deltaPapprch_1=q_e^2*timeStep*|coeffApprch_1|
#---------------------------------------------
      for m in range(6): 
         z_elecCrrnt[m]=0.                                         # Initial zero-vector for electron
      z_elecCrrnt[Ix]=rhoCrrnt                                     # x, cm
      z_elecCrrnt[Iz]=-halfLintr                                   # z, cm
      z_elecCrrnt[Ipy]=-m_elec*eVtran                              # py, g*cm/sec
      z_elecCrrnt[Ipz]=m_elec*eVrmsLong                            # pz, g*cm/sec
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
	 larmR_bCrrnt[k]=math.log10(rho_larm/bCrrnt[k])            # dimensionless 
	 uPot_enrgKinCrrnt[k]=math.log10((q_elec**2/bCrrnt[k])/kinEnergy)            # dimensionless 
# Current values to calculate deltaPapprch_1=q_e^2*timeStep*|coeffApprch_1|:  
         for ic in range(3):
	    dpApprch_1Crrnt[ic,k]= abs(z_elecCrrnt[2*ic])/ bCrrnt[k]**2              # 1/cm^2;  
# To draw only first trajectory (for checking only):
         if i==1 and j==1:
            for ic in (Ix,Iy,Iz):
               elecCoor[ic//2,pointTrack[trackNumb]]=z_elecCrrnt[ic]                 # cm
	       elecCoor[3,pointTrack[trackNumb]]=bCrrnt[k]         # cm
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
	          bTot[m]=bCrrnt[k]                                # cm
	          larmR_bTot[m]=larmR_bCrrnt[k]                    # dimensionless
		  uPot_enrgKinTot[m]=uPot_enrgKinCrrnt[k]          # dimensionless
                  for ic in range(3):
	             dpApprch_1Tot[ic,m]=dpApprch_1Crrnt[ic,k]     # 1/cm^2;  
	          mPrev=m
	          bPrev=bCrrnt[k]
		  depopFlag=1
	       population[m] += 1
	       break
#++++++++++++++++++++++++++++++++++
#  	 print 'i=%d, j=%d, k=%d: bPrev=%e (m=%d), bCrrnt=%e, x=%e, y=%e, z=%e' % \
#  	       (i,j,k,1.e+4*bPrev,mPrev,1.e+4*bCrrnt[k],1.e+4*z_elecCrrnt[0],1.e+4*z_elecCrrnt[2],1.e+4*z_elecCrrnt[4])       
         pointTrack[trackNumb] += 1
# End of gragging of the current trajectory	  
#-----------------------------------------------
###################################################
#
# Here was be placed block B): not in used now
#
###################################################   
      lastTrackNumber=trackNumb+1                                  # quantity of tracks = trackNumber + 1!     
      sumPoints += pointTrack[trackNumb]
print 'For %d tracks number of points is %d' % (lastTrackNumber,sumPoints)
# print 'larmorNumber: ', larmorNumber
# print 'pointTrack: ', pointTrack

nonZeroRslts=0
for m in range(nTotal):
   if population[m] != 0:
      nonZeroRslts += 1

timeEnd=os.times()
cpuTime=1.e+6*float(timeEnd[0])-1.e+6*float(timeStart[0])          # CPU time , mks
print 'Summa of population=%d; nonZeroRslts=%d' % (sum(population),nonZeroRslts)
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
ax10.set_zlabel('z, $\mu m$',color='m',fontsize=16)
plt.title(('First %d Larmor Turns of the First Trajectory:\nImpact Parameter=%6.3f $\mu$m, $R_L$=%6.3f $\mu$m' \
           % (turns,1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)),color='m',fontsize=16)

fig20=plt.figure(20)
ax20=fig20.gca(projection='3d')
ax20.plot(1.e+4*elecCoor[0,pointsTot-points:pointsTot],1.e+4*elecCoor[1,pointsTot-points:pointsTot], \
          1.e+4*elecCoor[2,pointsTot-points:pointsTot],'-r',linewidth=2)
plt.xlabel('x, $\mu m$',color='m',fontsize=16)
plt.ylabel('y, $\mu m$',color='m',fontsize=16)
ax20.set_zlabel('z, $\mu m$',color='m',fontsize=16)
plt.title(('Last %d Larmor Turns of the First Trajectory:\nImpact Parameter=%6.3f $\mu$m, $R_L$=%6.3f $\mu$m' \
           % (turns,1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)),color='m',fontsize=16)

# plt.show()   

plt.figure(30)
plt.plot(range(pointsTot),1.e4*elecCoor[2,0:pointsTot],'-r',linewidth=2)
plt.xlabel('Points',color='m',fontsize=16)
plt.ylabel('z, $\mu$m',color='m',fontsize=16)
plt.title(('First Trajectory: Impact Parameter=%6.3f $\mu$m, $R_L$=%6.3f $\mu$m' % \
           (1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)),color='m',fontsize=16)
plt.grid(True)

plt.figure(40)
plt.plot(1.e4*elecCoor[2,0:pointsTot],1.e4*elecCoor[3,0:pointsTot],'-r',linewidth=2)
plt.xlabel('z, $\mu$m',color='m',fontsize=16)
plt.ylabel('Distance from Motionless Ion, $\mu$m',color='m',fontsize=16)
plt.title(('First Trajectory: Impact Parameter=%6.3f $\mu$m, $R_L$=%6.3f $\mu$m' % \
           (1.e+4*rhoFirstTurn,1.e+4*rhoLarmorFirstTurn)),color='m',fontsize=16)
plt.grid(True)

plt.figure(50)
plt.plot(larmR_bTot,uPot_enrgKinTot,'.r')
plt.xlabel('$log_{10}(R_L/b)$',color='m',fontsize=16)
plt.ylabel('$log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
plt.title('Map for Transfered Momenta $dP_x,dP_y,dP_z$ Calculations', color='m',fontsize=20)
plt.grid(True)

X,Y=np.meshgrid(larmR_bTot,uPot_enrgKinTot)      
fig60=plt.figure(60)
ax60=fig60.gca(projection='3d')
surf=ax60.plot_surface(X,Y,dpApprch_1Tot[0,:],cmap=cm.coolwarm,linewidth=0,antialiased=False)
plt.title('Transfered Momentum $dP_x$:\n$dP_x=q_e^2/b \cdot C_x$', color='m',fontsize=20)
plt.xlabel('$log_{10}(R_L/b)$',color='m',fontsize=16)
plt.ylabel('$log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
ax60.set_zlabel('$C_x$, $cm^{-2}$',color='m',fontsize=16)
fig60.colorbar(surf, shrink=0.5, aspect=5)
plt.grid(True)

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

#####################################################################
#
# Block A)
#
#--------------------------------------------------------------------
# 
#       if i == 1 and j == 1: 
# # First definition of the total distance from origin of the coordinate system to electron along the trajectory; cm:
#  	 b=bCrrnt
# # First definition of the total log10 of two important ratios; dimensionless:  
# 	 larmR_b=larmR_bCrrnt                                      # dimensionless 
# 	 uPot_enrgKin=uPot_enrgKinCrrnt                            # dimensionless 
# # First definition of the values to calculate deltaPapprch_1=q_e^2*timeStep*|coeffApprch_1|:
# #         for ic in range(3):
# #            for k in range(timePoints):
# #	       apprch_1[ic,k]=coeffApprch_1[ic,k]                     # 1/cm^2;  
#       else:  
# # Total distance from origin of the coordinate system to electron along the trajectory:
#  	 b=np.concatenate((b,bCrrnt),axis=0)                       # cm
# # Total log10 of two important ratios; dimensionless :  
# 	 larmR_b=np.concatenate((larmR_b,larmR_bCrrnt),axis=0)                  
# 	 uPot_enrgKin=np.concatenate((uPot_enrgKin,uPot_enrgKinCrrnt),axis=0)        
# # Total values to calculate deltaPapprch_1=q_e^2*timeStep*|coeffApprch_1|:
# #         for ic in range(3):
# #	    apprch_1[ic,:]=np.concatenate((apprch_1[ic,:],coeffApprch_1[ic,:]),axis=0)      # 1/cm^2;  
#
##################### End of block A) ###############################
