# from __future__ import division

#-------------------------------------
#
#        Started at 09/08/2017 (YuE)
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
q_elec = qe*2.997e+9                                               # charge of electron, CGSE unit (without sign, i.e. q_elec > 0!)
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


alpha=2*q_elec**2*omega_L/(m_elec*eVrmsLong**3)                    # dimensionless
print 'alpha = %f' % alpha

# A=log_10(U_pot/E_kin):
#
nA=150
crrntA=np.zeros(nA)
minA=-5.
maxA=0.
stepA=(maxA-minA)/(nA-1)

# B=log_10(Rlarm/b):
#
nB=150
crrntB=np.zeros(nB)
minB=-3.
maxB=-.45
stepB=(maxB-minB)/(nB-1)

C=np.zeros((nA,nB))
logC=np.zeros((nA,nB))
root=np.zeros((nA,nB))
vTransv=np.zeros((nA,nB))
rhoLarm=np.zeros((nA,nB))
dist=np.zeros((nA,nB))
rho=np.zeros((nA,nB))
Lint=np.zeros((nA,nB))
Nlarm=np.zeros((nA,nB))
smallNlarm=np.zeros((nA*nB,5))
largeNlarm=np.zeros((nA*nB,3))

mapA=np.zeros(nA*nB)
mapB=np.zeros(nA*nB)
mapNlarm=np.zeros((nA*nB,nA*nB))

minC=1.e8
maxC=0.
minS=1.e8
maxS=0.

minR=1.e8
maxR=0.
minR=1.e8
maxR=0.

minLeftPartFirst=1.e8
maxLeftPartFirst=-1.e8

minLeftPartLast=1.e8
maxLeftPartLast=-1.e8

vTransvMin=3.e10
vTransvMax=0.

minRhoLarm=1.e8
maxRhoLarm=0.

minDist=1.e8
maxDist=0.

minRatio=1.e8
maxRatio=0.

minRho=1.e8
maxRho=0.

minLint=1.e8
maxLint=0.

minNlarm=10000000
maxNlarm=0

lenMap=0
NlarmCutofDown=39
NlarmCutofUp=2000
indxSmallNlarm=0
indxLargeNlarm=0

for iA in range(nA):
   crrntA[iA]=minA+stepA*iA
   for iB in range(nB):
      crrntB[iB]=minB+stepB*iB
      C[iA,iB]=alpha*math.pow(10.,crrntB[iB]-crrntA[iA])
      if minC > C[iA,iB]:
         minC=C[iA,iB]
      if maxC < C[iA,iB]:
         maxC=C[iA,iB]
      logC[iA,iB]=math.log10(C[iA,iB])
#      print 'A=%f, B=%f: C=%e, logC=%e' % (crrntA[iA],crrntB[iB],C[iA,iB],logC[iA,iB])
      q=-C[iA,iB]
#
# Equation: y^3+y+q=0:
#
      D=1./27+(q/2)**2
#       print 'A=%f, B=%f: D=%e' % (crrntA[iA],crrntB[iB],D)
      root[iA,iB]=math.pow(np.sqrt(D)-q/2,1./3)-math.pow(np.sqrt(D)+q/2,1./3)	 
      if minR > root[iA,iB]:
         minR=root[iA,iB]
      if maxR < root[iA,iB]:
         maxR=root[iA,iB]
# Checking, that root^3+root+q=0:
      leftPart=root[iA,iB]**3+root[iA,iB]+q  
      if minLeftPartFirst > leftPart:
         minLeftPartFirst=leftPart
      if maxLeftPartFirst < leftPart:
         maxLeftPartFirst=leftPart
# Corrections of the root:
      mm=0
      while (abs(leftPart) > 1.e-8) and (mm  < 5):
         deltaRoot=-leftPart/(root[iA,iB]*(3.*root[iA,iB]+1))
         root[iA,iB] += deltaRoot
         leftPart=root[iA,iB]**3+root[iA,iB]+q  
	 mm += 1
      if minLeftPartLast > leftPart:
         minLeftPartLast=leftPart
      if maxLeftPartLast < leftPart:
         maxLeftPartLast=leftPart

      vTransv[iA,iB]=eVrmsLong*root[iA,iB] 
      if vTransvMin > vTransv[iA,iB]:
         vTransvMin=vTransv[iA,iB]
      if vTransvMax < vTransv[iA,iB]:
         vTransvMax=vTransv[iA,iB]

      rhoLarm[iA,iB]=vTransv[iA,iB]/omega_L	 
      if minRhoLarm > rhoLarm[iA,iB]:
         minRhoLarm=rhoLarm[iA,iB]
      if maxRhoLarm < rhoLarm[iA,iB]:
         maxRhoLarm=rhoLarm[iA,iB]

      dist[iA,iB]=rhoLarm[iA,iB]/math.pow(10.,crrntB[iB])
      if minDist > dist[iA,iB]:
         minDist=dist[iA,iB]
      if maxDist < dist[iA,iB]:
         maxDist=dist[iA,iB]

      ratio=dist[iA,iB]/(rhoLarm[iA,iB]+rhoCrit)
      if minRatio > ratio:
         minRatio=ratio
      if maxRatio < ratio:
         maxRatio=ratio
      if ratio < 1.:
         print 'A=%f, B=%f: vTransv=%e, rhoLarm=%e, dist=%e, dist/rhoLarm=%e' \
               (crrntA[iA],crrntB[iB],vTransv[iA,iB],1e+4*rhoLarm[iA,iB],1e+4*dist[iA,iB],ratio)  
# print 'minC=%e, maxC=%e' % (minC,maxC)

      rho[iA,iB]=rhoCrit+rhoLarm[iA,iB]
      if minRho > rho[iA,iB]:
         minRho=rho[iA,iB]
      if maxRho < rho[iA,iB]:
         maxRho=rho[iA,iB]

      Lint[iA,iB]=np.sqrt(dist[iA,iB]**2-rho[iA,iB]**2)
      if minLint > Lint[iA,iB]:
         minLint=Lint[iA,iB]
      if maxLint < Lint[iA,iB]:
         maxLint=Lint[iA,iB]

      Nlarm[iA,iB]=int(Lint[iA,iB]/eVrmsLong/T_larm)      
      if minNlarm > Nlarm[iA,iB]:
         minNlarm=Nlarm[iA,iB]
      if maxNlarm < Nlarm[iA,iB]:
         maxNlarm=Nlarm[iA,iB]

#       print 'A=%f, B=%f: vTransv=%e, rhoLarm=%e, dist=%e, dist/rhoLarm=%e, rho=%e, Lint=%e, Nlarm=%d' % \
#             (crrntA[iA],crrntB[iB],vTransv[iA,iB],1e+4*rhoLarm[iA,iB],1e+4*dist[iA,iB],ratio, \
# 	     1e+4*rho[iA,iB],1e+4*Lint[iA,iB],Nlarm[iA,iB]) 
      
      if (Nlarm[iA,iB] >= NlarmCutofDown) and (Nlarm[iA,iB] <= NlarmCutofUp):
         mapA[lenMap]=crrntA[iA]
         mapB[lenMap]=crrntB[iB]
         mapNlarm[lenMap,lenMap]=Nlarm[iA,iB]
         lenMap += 1
      if Nlarm[iA,iB] <= NlarmCutofDown:
         smallNlarm[indxSmallNlarm,0]=crrntA[iA]	 
         smallNlarm[indxSmallNlarm,1]=crrntB[iB]	 
         smallNlarm[indxSmallNlarm,2]=Nlarm[iA,iB]
         smallNlarm[indxSmallNlarm,3]=dist[iA,iB]
         smallNlarm[indxSmallNlarm,4]=ratio
	 indxSmallNlarm += 1	 
      if Nlarm[iA,iB] >= NlarmCutofUp:
         largeNlarm[indxLargeNlarm,0]=crrntA[iA]	 
         largeNlarm[indxLargeNlarm,1]=crrntB[iB]	 
         largeNlarm[indxLargeNlarm,2]=Nlarm[iA,iB]	 
	 indxLargeNlarm += 1	 
print 'Root: min=%e, max=%e, LeftPart: minFirst=%e, maxFirst=%e, minLast=%e, maxLast=%e' % \
      (minR,maxR,minLeftPartFirst,maxLeftPartFirst,minLeftPartLast,maxLeftPartLast)
print 'Transverse velocity (cm/sec): min=%e, max=%e' % (vTransvMin,vTransvMax)
print 'Larmor radius: (mkm): min=%e, max=%e' % (1e+4*minRhoLarm,1e+4*maxRhoLarm)
print 'Distance (mkm): min=%e, max=%e; dist/rholarm: min=%e, max=%e' % (1e+4*minDist,1e+4*maxDist,minRatio,maxRatio)
print 'minRho=%e, maxRho=%e, minLint=%e, maxLint=%e' % (1e+4*minRho,1e+4*maxRho,1e+4*minLint,1e+4*maxLint)
print ' minNlarm=%d, maxNlarm=%d' % (minNlarm,maxNlarm)
print 'lenMap=%d' % lenMap

print '%d Cases with Small Larmor Radius (from Total %d)' % (indxSmallNlarm, nA*nB)
# for i in range(indxSmallNlarm):
#    print '%d): A=%f, B=%f, Nlarm=%d, dist (mkm)=%f, dist/Rlarm=%f' % \
#          (i,smallNlarm[i,0],smallNlarm[i,1],smallNlarm[i,2],1e+4*smallNlarm[i,3],smallNlarm[i,4])

print '%d Cases with Large Larmor Radius (from Total %d)' % (indxLargeNlarm, nA*nB)
# for i in range(indxLargeNlarm):
#    print '%d): A=%f, B=%f, Nlarm=%d' % (i,largeNlarm[i,0],largeNlarm[i,1],largeNlarm[i,2])

# A=log_10(U_pot/E_kin), B=log_10(Rlarm/b):
#
plt.figure(10)
plt.plot(mapA[0:lenMap-1],mapB[0:lenMap-1],'.r')
plt.xlabel('$B=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
plt.ylabel('$A=log_{10}(R_L/b)$',color='m',fontsize=16)
plt.title('Map for Transfered Momenta $dP_x,dP_y,dP_z$ Calculations', color='m',fontsize=20)
plt.text(-4.0,-.65,'Magnetized Electrons:', color='m',fontsize=20)
plt.text(-4.3,-.8,('Number of Larmor Turns > %d' % NlarmCutofDown), color='m',fontsize=20)
# plt.grid(True)

# A=log_10(U_pot/E_kin), B=log_10(Rlarm/b):
#
stepLevel=int(maxNlarm/20)
levels=np.arange(0,maxNlarm,stepLevel)
fig15=plt.figure(15)
ax=fig15.add_subplot(111)
X,Y=np.meshgrid(crrntA,crrntB) 
mapTurns=ax.contourf(X,Y,Nlarm,levels)   
# X,Y=np.meshgrid(mapA[0:lenMap-1],mapB[0:lenMap-1]) 
# aa=ax.contourf(X,Y,mapNlarm[0:lenMap-1,0:lenMap-1])   
plt.xlabel('$B=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
plt.ylabel('$A=log_{10}(R_L/b)$',color='m',fontsize=16)
plt.title('Magnetized Electrons: Number of Larmor Turns', color='m',fontsize=20)
# plt.text(-4.0,-.65,'Magnetized Electrons:', color='m',fontsize=20)
# plt.text(-4.3,-.8,('Number of Larmor Turns > %d' % NlarmCutof), color='m',fontsize=20)
fig15.colorbar(mapTurns)

'''
X=np.zeros(nA)
for iA in range(nA):
   X=crrntA[iA]*np.ones(nA)
   plt.figure(15)
   plt.hold(True)
   plt.plot(X,crrntB,'.r')
plt.xlabel('$log_{10}(R_L/b)$',color='m',fontsize=16)
plt.ylabel('$log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
plt.title('Map for Transfered Momenta $dP_x,dP_y,dP_z$ Calculations', color='m',fontsize=20)
plt.grid(True)
'''

'''      
fig10=plt.figure(10)
plt.plot(crrntA,logC[0,:],'-b',linewidth=1.5)
plt.hold(True)
for iB in range(1,nB):
   plt.plot(crrntA,logC[iB,:],'-b',linewidth=1.5)
plt.title('$C=log_{10}[2q_e^2 \cdot \omega_L/(m_e \cdot V_{e||}^3)] + B - A$ $(B=log_{10}(q_e^2/b/E_{kin}))$', color='m',fontsize=16)
plt.xlabel('$A=log_{10}(R_L/b)$',color='m',fontsize=16)
plt.ylabel('$C$',color='m',fontsize=16)
plt.xlim([minA-.1,maxA+.1])
for iB in range(nB):
#    plt.text(-3.5,-.9+.555*iB,('B%f' % crrntB[iB]),color='m')
   plt.text(-3.5,4.095-.555*iB,('B%f' % crrntB[iB]),color='m')
plt.grid(True)
'''

# A=log_10(U_pot/E_kin), B=log_10(Rlarm/b):
#
X,Y=np.meshgrid(crrntA,crrntB)      
fig20=plt.figure(20)
ax20=fig20.gca(projection='3d')
surf=ax20.plot_surface(X,Y,logC,cmap=cm.coolwarm,linewidth=0,antialiased=False)
plt.title('Map C:\n$C=log_{10}[2q_e^2 \cdot \omega_L/(m_e \cdot V_{e||}^3)] + B - A$', color='m',fontsize=16)
plt.xlabel('$B=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
plt.ylabel('$A=log_{10}(R_L/b)$',color='m',fontsize=16)
ax20.set_zlabel('$C$',color='m',fontsize=16)
# fig20.colorbar(surf, shrink=0.5, aspect=5)
fig20.colorbar(surf)
plt.grid(True)

# A=log_10(U_pot/E_kin), B=log_10(Rlarm/b):
#
X,Y=np.meshgrid(crrntA,crrntB)      
fig30=plt.figure(30)
ax30=fig30.gca(projection='3d')
surf=ax30.plot_surface(X,Y,1.e+4*rhoLarm,cmap=cm.coolwarm,linewidth=0,antialiased=False)
plt.title('Larmor Radius $R_L$', color='m',fontsize=16)
plt.xlabel('$B=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
plt.ylabel('$A=log_{10}(R_L/b)$',color='m',fontsize=16)
ax30.set_zlabel('$R_L$, $\mu m$',color='m',fontsize=16)
# fig20.colorbar(surf, shrink=0.5, aspect=5)
fig30.colorbar(surf)
plt.grid(True)

# A=log_10(U_pot/E_kin), B=log_10(Rlarm/b):
#
stepLevel=int(1.e+4*maxRhoLarm/20)
levels=np.arange(0,1.e+4*maxRhoLarm,stepLevel)
fig35=plt.figure(35)
ax=fig35.add_subplot(111)
X,Y=np.meshgrid(crrntA,crrntB) 
mapTurns=ax.contourf(X,Y,1.e+4*rhoLarm,levels)   
# X,Y=np.meshgrid(mapA[0:lenMap-1],mapB[0:lenMap-1]) 
# aa=ax.contourf(X,Y,mapNlarm[0:lenMap-1,0:lenMap-1])   
plt.xlabel('$B=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
plt.ylabel('$A=log_{10}(R_L/b)$',color='m',fontsize=16)
plt.title('Larmor Radius $R_L$, $\mu$m', color='m',fontsize=16)
# plt.text(-4.0,-.65,'Magnetized Electrons:', color='m',fontsize=20)
# plt.text(-4.3,-.8,('Number of Larmor Turns > %d' % NlarmCutof), color='m',fontsize=20)
fig35.colorbar(mapTurns)

# A=log_10(U_pot/E_kin), B=log_10(Rlarm/b):
#
X,Y=np.meshgrid(crrntA,crrntB)      
fig40=plt.figure(40)
ax40=fig40.gca(projection='3d')
surf=ax40.plot_surface(X,Y,1.e+4*dist,cmap=cm.coolwarm,linewidth=0,antialiased=False)
plt.title('Distance $b$', color='m',fontsize=16)
plt.xlabel('$B=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
plt.ylabel('$A=log_{10}(R_L/b)$',color='m',fontsize=16)
ax40.set_zlabel('$b$, $\mu$m',color='m',fontsize=16)
# fig20.colorbar(surf, shrink=0.5, aspect=5)
fig40.colorbar(surf)
plt.grid(True)

# A=log_10(U_pot/E_kin), B=log_10(Rlarm/b):
#
stepLevel=int(1.e+4*maxDist/20)
levels=np.arange(0,1.e+4*maxDist,stepLevel)
fig45=plt.figure(45)
ax=fig45.add_subplot(111)
X,Y=np.meshgrid(crrntA,crrntB) 
mapTurns=ax.contourf(X,Y,1.e+4*dist,levels)   
# X,Y=np.meshgrid(mapA[0:lenMap-1],mapB[0:lenMap-1]) 
# aa=ax.contourf(X,Y,mapNlarm[0:lenMap-1,0:lenMap-1])   
plt.xlabel('$B=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
plt.ylabel('$A=log_{10}(R_L/b)$',color='m',fontsize=16)
plt.title('Distance $b$, $\mu$m', color='m',fontsize=16)
# plt.text(-4.0,-.65,'Magnetized Electrons:', color='m',fontsize=20)
# plt.text(-4.3,-.8,('Number of Larmor Turns > %d' % NlarmCutof), color='m',fontsize=20)
fig45.colorbar(mapTurns)

plt.show()   

sys.exit()

