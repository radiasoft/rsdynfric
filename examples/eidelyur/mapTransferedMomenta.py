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


alpha=2*q_elec**2*omega_L/(m_elec*eVrmsLong**3)
print 'alpha = %f' % alpha

nA=150
crrntA=np.zeros(nA)
minA=-5.
maxA=0.
stepA=(maxA-minA)/(nA-1)

nB=150
crrntB=np.zeros(nB)
minB=-3.
maxB=-.5
stepB=(maxB-minB)/(nB-1)

C=np.zeros((nA,nB))
logC=np.zeros((nA,nB))
root=np.zeros((nA,nB))
vTransv=np.zeros((nA,nB))
rhoLarm=np.zeros((nA,nB))
dist=np.zeros((nA,nB))

minC=1.e8
maxC=0.
minS=1.e8
maxS=0.

minR=1.e8
maxR=0.
minR=1.e8
maxR=0.

minLeftPart=1.e8
maxLeftPart=-1.e8

vTransvMin=3.e10
vTransvMax=0.

minRhoLarm=1.e8
maxRhoLarm=0.

minDist=1.e8
maxDist=0.

minRatio=1.e8
maxRatio=0.

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
# Equation: x^3+x+q=0:
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
      if minLeftPart > leftPart:
         minLeftPart=leftPart
      if maxLeftPart < leftPart:
         maxLeftPart=leftPart
#       if abs(leftPart) > 1.e-8:
#          print 'A=%f, B=%f: C=%e, root=%e, leftPart=%e' % (crrntA[iA],crrntB[iB],C[iA,iB],root[iA,iB],leftPart)  
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
      ratio=dist[iA,iB]/rhoLarm[iA,iB]
      if minRatio > ratio:
         minRatio=ratio
      if maxRatio < ratio:
         maxRatio=ratio
      if ratio < 1.:
         print 'A=%f, B=%f: vTransv=%e, rhoLarm=%e, dist=%e, dist/rhoLarm=%e' \
               (crrntA[iA],crrntB[iB],vTransv[iA,iB],1e+4*rhoLarm[iA,iB],1e+4*dist[iA,iB],ratio)  
#       print 'A=%f, B=%f: vTransv=%e, rhoLarm=%e, dist=%e, dist/rhoLarm=%e' % \
#             (crrntA[iA],crrntB[iB],vTransv[iA,iB],1e+4*rhoLarm[iA,iB],1e+4*dist[iA,iB],ratio)  
# print 'minC=%e, maxC=%e' % (minC,maxC)
print 'minRoot=%e, maxroot=%e, minLeftPart=%e, maxLeftPart=%e' % (minR,maxR,minLeftPart,maxLeftPart)
print 'Transverse velocity (cm/sec): min=%e, max=%e' % (vTransvMin,vTransvMax)
print 'Larmor radius: (mkm): min=%e, max=%e' % (1e+4*minRhoLarm,1e+4*maxRhoLarm)
print 'Distance (mkm): min=%e, max=%e; dist/rholarm: min=%e, max=%e' % (1e+4*minDist,1e+4*maxDist,minRatio,maxRatio)


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

X,Y=np.meshgrid(crrntA,crrntB)      
fig20=plt.figure(20)
ax20=fig20.gca(projection='3d')
surf=ax20.plot_surface(X,Y,logC,cmap=cm.coolwarm,linewidth=0,antialiased=False)
plt.title('Map C:\n$C=log_{10}[2q_e^2 \cdot \omega_L/(m_e \cdot V_{e||}^3)] + B - A$', color='m',fontsize=16)
plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
plt.ylabel('$B=log_{10}(R_L/b)$',color='m',fontsize=16)
ax20.set_zlabel('$C$',color='m',fontsize=16)
# fig20.colorbar(surf, shrink=0.5, aspect=5)
fig20.colorbar(surf)
plt.grid(True)

X,Y=np.meshgrid(crrntA,crrntB)      
fig30=plt.figure(30)
ax30=fig30.gca(projection='3d')
surf=ax30.plot_surface(X,Y,1.e+4*rhoLarm,cmap=cm.coolwarm,linewidth=0,antialiased=False)
plt.title('Larmor Radius $R_L$', color='m',fontsize=16)
plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
plt.ylabel('$B=log_{10}(R_L/b)$',color='m',fontsize=16)
ax30.set_zlabel('$R_L$, $\mu m$',color='m',fontsize=16)
# fig20.colorbar(surf, shrink=0.5, aspect=5)
fig30.colorbar(surf)
plt.grid(True)

X,Y=np.meshgrid(crrntA,crrntB)      
fig40=plt.figure(40)
ax40=fig40.gca(projection='3d')
surf=ax40.plot_surface(X,Y,1.e+4*dist,cmap=cm.coolwarm,linewidth=0,antialiased=False)
plt.title('Distance $b$', color='m',fontsize=16)
plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
plt.ylabel('$B=log_{10}(R_L/b)$',color='m',fontsize=16)
ax40.set_zlabel('$b$, $\mu m$',color='m',fontsize=16)
# fig20.colorbar(surf, shrink=0.5, aspect=5)
fig40.colorbar(surf)
plt.grid(True)



plt.show()   

      
      





sys.exit()





