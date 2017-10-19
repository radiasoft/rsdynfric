# from __future__ import division

#-------------------------------------
#
#        Started at 10/18/2017 (YuE)
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
beta=2*q_elec**2/(m_elec*eVrmsLong**2)                             # cm
print 'Rshield(mkm)=%e, rhoCrit(mkm)=%f, alpha=%f' % (1.e+4*Rshield[0],rhoCrit, alpha)

#
# Array A=log10(Upot/Ekin)=log10([q_e^2/rho]/[m_e((V_transverse^2+V_longitudinal^2)/2]): 
#
nA=75
crrntA=np.zeros(nA)
minA=-5.
maxA=0.
stepA=(maxA-minA)/(nA-1)

#
# Array B=log10(R_larm/rho): 
#
nB=75
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
rrhoCrrntMax_B=0
trackNumb=-1                                                     # Tracks are enumerated from 0!  
track_1=-1                                                       # Tracks are enumerated from 0!
track_2=-1                                                       # Tracks are enumerated from 0!
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
         print 'iA=%d, iB=%d: track=%d, numbLarmor[iA,iB]=%d: Lintr=%e, rhoLarm=%e, rhoCrrnt=%e' % \
	 (iA,iB,trackNumb,numbLarmor[iA,iB],1.e4*halfLintr[iA,iB],1.e4*rho_larm[iA,iB],1.e4*rhoCrrnt[iA,iB])  
      else: 
         print 'iA=%d, iB=%d, rhoLarm=%e, rhoCrrnt=%e' % (iA,iB,1.e4*rho_larm[iA,iB],1.e4*rhoCrrnt[iA,iB]) 
#-----
# To plot area of the impact parameter b:
      if Rshield[0] <= rhoCrrnt[iA,iB] < Rshield[1]:
         track_1 += 1  
	 mapAimpctParmtr[1,track_1]=crrntA[iA]                     # Tracks are enumerated from 0!  
	 mapBimpctParmtr[1,track_1]=crrntB[iB]                     # Tracks are enumerated from 0!  
      if Rshield[1] <= rhoCrrnt[iA,iB] < Rshield[2]:
         track_2 += 1  
	 mapAimpctParmtr[2,track_2]=crrntA[iA]                     # Tracks are enumerated from 0!  
	 mapBimpctParmtr[2,track_2]=crrntB[iB]                     # Tracks are enumerated from 0!  
#-----
	  
print 'rhoCrrntMax=%e for iA=%e, iB=%e' % (1.e4*rhoCrrntMax,rhoCrrntMax_A,rhoCrrntMax_B)

plt.figure (10)
plt.plot(mapAimpctParmtr[0,0:trackNumb],mapBimpctParmtr[0,0:trackNumb],'.r',markersize=10)
plt.plot(mapAimpctParmtr[1,0:track_1],mapBimpctParmtr[1,0:track_1],'.b',markersize=10)
plt.plot(mapAimpctParmtr[2,0:track_2],mapBimpctParmtr[2,0:track_2],'.m',markersize=10)
plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
Rshield_mkm=1.e4*Rshield[0]
Rshield_150_mkm=1.e4*Rshield[1]
Rshield_200_mkm=1.e4*Rshield[2]
plt.title('Area of Impact Parameter $b$', color='m',fontsize=16)
plt.xlim([minA,maxA])
plt.ylim([minB,maxB])
plt.legend([('$b$ < %6.1f $\mu$m' % Rshield_mkm),('%6.1f $\mu$m < $b$ < %6.1f $\mu$m' % (Rshield_mkm,Rshield_150_mkm)),('%6.1f $\mu$m < $b$ < %6.1f $\mu$m' % (Rshield_150_mkm,Rshield_200_mkm))],loc='lower left',fontsize=16)
plt.grid(True)

X,Y=np.meshgrid(crrntA,crrntB)      

fig20=plt.figure(20)
ax20=fig20.gca(projection='3d')
surf=ax20.plot_surface(X,Y,1.e+4*halfLintr,cmap=cm.coolwarm,linewidth=0,antialiased=False)
plt.title('Half Length of Interaction $L_{interaction}$, $\mu$m', color='m',fontsize=16)
plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
ax20.set_zlabel('$1/2 \cdot L_{interaction}$',color='m',fontsize=16)
fig20.colorbar(surf, shrink=1., aspect=10)
plt.grid(True)


fig30=plt.figure(30)
ax30=fig30.gca(projection='3d')
surf=ax30.plot_surface(X,Y,numbLarmor,cmap=cm.coolwarm,linewidth=0,antialiased=False)
plt.title('Number of Larmor Turns $N_{Larmor}$', color='m',fontsize=16)
plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
plt.ylabel('$B=log_{10}(R_{Larm}/b)$',color='m',fontsize=16)
ax30.set_zlabel('$N_{Larmor}$',color='m',fontsize=16)
fig30.colorbar(surf, shrink=1., aspect=10)
plt.grid(True)

plt.show()

sys.exit()

