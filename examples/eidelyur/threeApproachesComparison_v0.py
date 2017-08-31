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
stepsVtran=100
stepVtran=(maxVtran-minVtran)/(stepsVtran-1)                       # cm/sec

z_elecCrrnt=np.zeros(6)                                            # 6-vector for electron
pointTrack=np.zeros(2)                                             # to draw only 2 trajectories
minLarmR_b=np.zeros(stepsImpctParMin_Rlarm)                        # dimensionless     
maxLarmR_b=np.zeros(stepsImpctParMin_Rlarm)                        # dimensionless
minUpot_enrgKin=np.zeros(stepsImpctParMin_Rlarm)                   # dimensionless      
maxUpot_enrgKin=np.zeros(stepsImpctParMin_Rlarm)                   # dimensionless      
turns=10                                                           # Number of larmorturns for drawing 
points=turns*stepsNumberOnGyro                                     # Number of points for drawing
'''
#
# Case: the same transversal (eVrmsTran) and different longidudinal velocity: 
#
for i in range(1):    # stepsVlong):
   eVlong=maxVlong-stepVlong*i                                     # cm/sec
   kinEnergy=m_elec*(eVrmsTran**2+eVlong**2)/2.                    # erg
   for j in range(stepsImpctParMin_Rlarm):
      impctParMin_Rlarm=maxImpctParMin_Rlarm-stepImpctParMin_Rlarm*j         # cm
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
'''
#
# Case: the same longidudinal (eVrmsLong) and different transversal velocity: 
#
for i in range(1):    # stepsVtran):
   eVtran=maxVtran-stepVtran*i                                     # cm/sec
   ro_larm = eVtran/omega_L                                        # cm
   kinEnergy=m_elec*(eVtran**2+eVrmsLong**2)/2.                    # erg
   for j in range(stepsImpctParMin_Rlarm):
      impctParMin_Rlarm=maxImpctParMin_Rlarm-stepImpctParMin_Rlarm*j         # cm
      minImpctPar=impctParMin_Rlarm*ro_larm                        # cm
      halfIntrcnLength=minImpctPar*math.tan(pi/2-stepAlpha)        # cm
      timePath=halfIntrcnLength/eVrmsLong
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
      z_elecCrrnt[Ipy]=m_elec*eVtran                               # py, g*cm/sec
      z_elecCrrnt[Ipz]=m_elec*eVrmsLong                            # pz, g*cm/sec
      for k in range(timePoints):
 	  z_elecCrrnt=matr_elec.dot(z_elecCrrnt)
          for ic in (Ix,Iy,Iz):
             elecCoor[ic//2,pointTrack[j],j]=z_elecCrrnt[ic]            # cm
# Distance from origin of the coordinate system to electron along the trajectory; cm:
 	  b[j,k]=np.sqrt(elecCoor[0,pointTrack[j],j]**2+elecCoor[1,pointTrack[j],j]**2+elecCoor[2,pointTrack[j],j]**2)
# log10 of two important ratios; dimensionless:  
	  larmR_b[j,k]=math.log10(ro_larm/b[j,k])                  # dimensionless 
	  upot_enrgKin[j,k]=math.log10((q_elec**2/b[j,k])/kinEnergy)         # dimensionless   
          pointTrack[j] += 1
      minLarmR_b[j]=min(larmR_b[j,:])      
      maxLarmR_b[j]=max(larmR_b[j,:])
      print 'For j=%d larmR_b: min=%e, max=%e' % (j,minLarmR_b[j],maxLarmR_b[j])      
      minUpot_enrgKin[j]=min(upot_enrgKin[j,:])
      maxUpot_enrgKin[j]=max(upot_enrgKin[j,:])
      print 'For j=%d upot_enrgKin: min=%e, max=%e' % (j,minUpot_enrgKin[j],maxUpot_enrgKin[j])      

#
# Checking of trajectories:
#
turns=10
points=turns*stepsNumberOnGyro

fig10=plt.figure(10)
ax10=fig10.gca(projection='3d')
ax10.plot(1.e+4*elecCoor[0,0:points,0],1.e+4*elecCoor[1,0:points,0],1.e+4*elecCoor[2,0:points,0],'-r',linewidth=2)
plt.xlabel('x, $\mu m$',color='m',fontsize=16)
plt.xlabel('x, $\mu m$',color='m',fontsize=16)
plt.title(('First %d Larmor Turns' % turns),color='m',fontsize=16)

fig20=plt.figure(20)
ax20=fig20.gca(projection='3d')
ax20.plot(1.e+4*elecCoor[0,0:points,1],1.e+4*elecCoor[1,0:points,1],1.e+4*elecCoor[2,0:points,1],'-b',linewidth=2)
plt.xlabel('x, $\mu m$',color='m',fontsize=16)
plt.xlabel('x, $\mu m$',color='m',fontsize=16)
plt.title(('First %d Larmor Turns' % turns),color='m',fontsize=16)

# 
#
#
fig30=plt.figure(30)
plt.plot(upot_enrgKin[0,0:pointTrack[0]],larmR_b[0,0:pointTrack[0]],'.r', \
         upot_enrgKin[1,0:pointTrack[1]],larmR_b[1,0:pointTrack[1]],'.b',linewidth=2)

fig40=plt.figure(40)
plt.plot(range(int(pointTrack[0])),upot_enrgKin[0,0:pointTrack[0]],'.r')


fig50=plt.figure(50)
plt.plot(range(int(pointTrack[1])),upot_enrgKin[0,0:pointTrack[1]],'.b')





plt.show()   




sys.exit()
