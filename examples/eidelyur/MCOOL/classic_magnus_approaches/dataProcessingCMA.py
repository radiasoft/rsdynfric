# from __future__ import division

#-------------------------------------
#
#        Started at 07/04/2018 (YuE)
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

#
# Opening the input file: 
#
inputFile='dependenceDeltaEonImpactParameter.data'
print 'Open input file "%s"...' % inputFile
inpfileFlag=0
try:
   inpfile = open(inputFile,'r')
   inpfileFlag=1
except:
   print 'Problem to open input file "%s"' % inputFile
if inpfileFlag == 1:
   print 'No problem to open input file "%s"' % inputFile


lines = 0                                # Number of current line from input file   
lineData = inpfile.readline()
lineData = inpfile.readline()
words = lineData.split()
nWords = len(words)
# print 'Data from %d: words=%s, number of entries = %d' % (lines+1,words,nWords)
nImpctPrmtr=int(words[1])
print 'nImpctPrmtr = %d' % nImpctPrmtr
lineData = inpfile.readline()
words = lineData.split()
nWords = len(words)
# print 'Data from %d: words=%s, number of entries = %d' % (lines+1,words,nWords)
nVion = int(nWords-1)
print 'nVion = %d' % nVion

vIon = np.zeros(nVion)
for i in range(nVion):
   vIon[i] = float(words[i+1])
   print 'vIon[%d] = %10.4e' % (i,vIon[i])


rhoInit=np.zeros((nImpctPrmtr,nVion))
deltaEnrgIon_c=np.zeros((nImpctPrmtr,nVion))

lineData = inpfile.readline()
lineData = inpfile.readline()
lineData = inpfile.readline()
lines += 7
while True:
   lineData=inpfile.readline()
#   print 'line=%d: %s' % (lines,lineData)
   if not lineData:
      break
   lines += 1
   words = lineData.split()
   nWords = len(words)
   for i in range(nVion):
      n = int(words[0])
      rhoInit[n,i] = float(words[1+2*i])
      deltaEnrgIon_c[n,i] = float(words[2+2*i])
#      print 'rho[%2d,%2d]=%10.4e, deltaE[%2d,%2d] = %10.4e' % (n,i,rhoInit[n,i],n,i,deltaEnrgIon_c[n,i])

inpfile.close()
print 'Close input file "%s"' % inputFile

'''
for n in range(nImpctPrmtr):
   print '%2d' % n,
   for i in range(nVion):
      print '  %10.4e  %10.4e' % (rhoInit[n,i],deltaEnrgIon_c[n,i]),
   print ''
'''

#
# Fitting:
#
fitA = np.zeros(nVion)
fitB = np.zeros(nVion)
maxIndx = np.zeros(nVion)
log10rhoInit = np.zeros((nImpctPrmtr,nVion))
log10deltaEnrgIon_c = np.zeros((nImpctPrmtr,nVion))

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
   for n in range(int(maxIndx[i])+1):
      sumRho[i] += log10rhoInit[n,i]
      sumRho2[i] += log10rhoInit[n,i]**2
      sumEnrg[i] += log10deltaEnrgIon_c[n,i]
      sumRhoEnrg[i] += log10rhoInit[n,i]*log10deltaEnrgIon_c[n,i]

   delta = maxIndx[i]*sumRho2[i]-sumRho[i]**2
   fitA[i] = (sumRho2[i]*sumEnrg[i]-sumRho[i]*sumRhoEnrg[i])/delta
   fitB[i] = (maxIndx[i]*sumRhoEnrg[i]-sumRho[i]*sumEnrg[i])/delta
#   print 'fitA(%d) = %e, fitB(%d) = %e' % (i,fitA[i],i,fitB[i])

rhoInitFit = np.zeros((maxIndx[0]+1,nVion))
deltaEnrgIon_c_Fit = np.zeros((maxIndx[0]+1,nVion))
funcHi2 = np.zeros(nVion)
for i in range(nVion):
   factorA = math.pow(10.,fitA[i])
   for n in range(int(maxIndx[i])+1):
      rhoInitFit[n,i] = math.pow(10.,log10rhoInit[n,i])
      deltaEnrgIon_c_Fit[n,i] = factorA*math.pow(rhoInitFit[n,i],fitB[i])
#      funcHi2[i] += (deltaEnrgIon_c[n,i]-deltaEnrgIon_c_Fit[n,i])**2  
      funcHi2[i] += (1.-deltaEnrgIon_c_Fit[n,i]/deltaEnrgIon_c[n,i])**2  
   print 'i=%2d: fitA = %e, fitB = %e, hi2 = %e' % \
         (i,fitA[i],fitB[i],funcHi2[i])

'''
#
# Optimization of the fitting:
#
fitAopt = np.zeros(nVion)
fitBopt = np.zeros(nVion)
deltaEnrgIon_c_FitOpt = np.zeros((maxIndx[0],nVion))
funcHi2opt = np.zeros(nVion)

indxCrrnt = 5
fitAopt[indxCrrnt] = fitA[indxCrrnt]
fitBopt[indxCrrnt] = fitB[indxCrrnt]

val1 = fitA[indxCrrnt]
val2 = fitB[indxCrrnt]
endInputFlag = 0
while (endInputFlag == 0):
   print ('Old values: fitAopt[%d] = %10.4e, fitBopt[%d] = %10.4e' % \
          (indxCrrnt,val1,indxCrrnt,val2))
   val1,val2 = raw_input('New values: ').split(',') 
   val1 = float(val1) 
   val2 = float(val2) 
   print ('New values: fitAopt[%d] = %10.4e, fitBopt[%d] = %10.4e' % \
          (indxCrrnt,val1,indxCrrnt,val2))
   if ((val1 == 0.) or (val2 == 0.)):
      print 'Script is terminated' 
      sys.exit()
   else:
      fitAopt[indxCrrnt] = val1
      fitBopt[indxCrrnt] = val2

   rho1 = rhoInit[0,indxCrrnt]
   rho2 = rhoInit[nImpctPrmtr-1,indxCrrnt]
   dE1 = deltaEnrgIon_c[0,indxCrrnt]
   dE2 = deltaEnrgIon_c[nImpctPrmtr-1,indxCrrnt]

   print 'rho1 = %10.4e, dE1 = %10.4e, rho2 = %10.4e, dE2 = %10.4e' % \
         (rho1,dE1,rho2,dE2),
   fitBopt[indxCrrnt] = np.log10(dE2/dE1)/np.log10(rho2/rho1)
   fitAopt[indxCrrnt] = np.log10(dE2) - fitBopt[indxCrrnt]*np.log10(rho2) 
   print 'fiAopt = %10.4e, fitBopt = %10.4e' % (fitAopt[indxCrrnt],fitBopt[indxCrrnt])


   factorA = math.pow(10.,fitAopt[indxCrrnt])
   funcHi2opt[indxCrrnt] = 0.
   for n in range(int(maxIndx[indxCrrnt])):
      rhoInitFit[n,indxCrrnt] = math.pow(10.,log10rhoInit[n,i])
      deltaEnrgIon_c_FitOpt[n,indxCrrnt] = factorA*math.pow(rhoInitFit[n,indxCrrnt],fitBopt[indxCrrnt])
#         funcHi2[i] += (deltaEnrgIon_c[n,i]-deltaEnrgIon_c_Fit[n,i])**2  
      funcHi2opt[indxCrrnt] += (1.-deltaEnrgIon_c_FitOpt[n,indxCrrnt]/deltaEnrgIon_c[n,indxCrrnt])**2  
   print 'index=%2d: fitA = %e, fitB = %e, hi2 = %e' % \
         (indxCrrnt,fitAopt[indxCrrnt],fitBopt[indxCrrnt],funcHi2opt[indxCrrnt])

   yLimit=[1.5e-10,2.e-8]

   indxPlot=5
   VionCrrnt = vIon[indxPlot]
   powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
   mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
   fig1000=plt.figure (1000)
   plt.loglog(rhoInit[0:nImpctPrmtr,indxPlot], \
              deltaEnrgIon_c[0:nImpctPrmtr,indxPlot],'-xr', \
              rhoInitFit[0:maxIndx[indxPlot]+1,indxPlot], \
   	   deltaEnrgIon_c_Fit[0:maxIndx[indxPlot]+1,indxPlot],'ob', \
              rhoInitFit[0:maxIndx[indxPlot]+1,indxPlot], \
   	   deltaEnrgIon_c_FitOpt[0:maxIndx[indxPlot]+1,indxPlot],'om',linewidth=2)
#    plt.loglog(rhoInit[0:nImpctPrmtr-1,indxPlot], \
#               deltaEnrgIon_c[0:nImpctPrmtr-1,indxPlot],'-xr', \
#               rhoInitFit[0:maxIndx[indxPlot],indxPlot], \
#    	   deltaEnrgIon_c_Fit[0:maxIndx[indxPlot],indxPlot],'ob', \
#    	   linewidth=2)
   plt.xlabel('Track Initial Impact Parameter $rho_{Init}$', \
              color='m',fontsize=16)
   plt.ylabel('$\Delta E_{ion}$, $eV$', color='m',fontsize=16)
   titleHeader = 'Transferred Energy $\Delta E_{ion}$ to Single Ion for '
   titleHeader += ' $V_{ion}=%5.3f\cdot10^{%2d}$ $cm/s$'
   plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=16)
#    plt.xlim([.95*rhoInit[0,0],1.05*rhoInit[nImpctPrmtr-1,0]])
#    plt.ylim(yLimit)  # for "indxPlot=49"
   plt.legend(['Calculated Data','Fitted Data'],loc='lower center',fontsize=16)
#    plt.text(2.5e-3,1.e-8, \
#             ('Fitted $\Delta E_{ion}$ are proportional to $1/rho_{Init}^{%5.3f}$' % \
#    	  abs(fitBopt[indxPlot])),color='m',fontsize=16)
#    plt.text(2.5e-3,1.e-8, \
#             ('Fitted $\Delta E_{ion}$ are proportional to $1/rho_{Init}^{%5.3f}$' % \
#    	  abs(fitB[indxPlot])),color='m',fontsize=16)
   plt.grid(True)

   plt.show()

'''
indxPlot=0
VionCrrnt = vIon[indxPlot]
powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
fig500=plt.figure (500)
# plt.loglog(rhoInit[0:nImpctPrmtr,indxPlot], \
#            deltaEnrgIon_c[0:nImpctPrmtr,indxPlot],'-xr', \
#            rhoInitFit[0:maxIndx[indxPlot]+1,indxPlot], \
# 	   deltaEnrgIon_c_Fit[0:maxIndx[indxPlot]+1,indxPlot],'ob',linewidth=2)
# plt.loglog(rhoInit[0:nImpctPrmtr,indxPlot], \
#            deltaEnrgIon_c[0:nImpctPrmtr,indxPlot],'-xr', \
#            rhoInitFit[0:maxIndx[indxPlot]+1,indxPlot], \
# 	   deltaEnrgIon_c_Fit[0:maxIndx[indxPlot]+1,indxPlot],'ob', \
#            rhoInitFit[0:maxIndx[indxPlot]+1,indxPlot], \
# 	   deltaEnrgIon_c_FitOpt[0:maxIndx[indxPlot]+1,indxPlot],'om',linewidth=2)
plt.loglog(rhoInit[0:nImpctPrmtr,indxPlot], \
           deltaEnrgIon_c[0:nImpctPrmtr,indxPlot],'-xr',linewidth=2)
plt.xlabel('Track Initial Impact Parameter $rho_{Init}$', \
           color='m',fontsize=16)
plt.ylabel('$\Delta E_{ion}$, $eV$', color='m',fontsize=16)
titleHeader = 'Transferred Energy $\Delta E_{ion}$ to Single Ion for '
titleHeader += ' $V_{ion}=%5.3f\cdot10^{%2d}$ $cm/s$'
plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=16)
plt.xlim([.95*rhoInit[0,0],1.05*rhoInit[nImpctPrmtr-1,0]])
# plt.ylim([3.e-11,2.e-8])
plt.legend(['Calculated Data','Fitted Data'],loc='lower center',fontsize=16)
# plt.text(2.5e-3,1.e-8, \
#          ('Fitted $\Delta E_{ion}$ are proportional to $1/rho_{Init}^{%5.3f}$' % \
# 	  abs(fitBopt[indxPlot])),color='m',fontsize=16)
# plt.text(2.5e-3,1.e-8, \
#          ('Fitted $\Delta E_{ion}$ are proportional to $1/rho_{Init}^{%5.3f}$' % \
# 	  abs(fitB[indxPlot])),color='m',fontsize=16)
plt.grid(True)

indxPlot=1
VionCrrnt = vIon[indxPlot]
powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
fig600=plt.figure (600)
# plt.loglog(rhoInit[0:nImpctPrmtr,indxPlot], \
#            deltaEnrgIon_c[0:nImpctPrmtr,indxPlot],'-xr', \
#            rhoInitFit[0:maxIndx[indxPlot]+1,indxPlot], \
# 	   deltaEnrgIon_c_Fit[0:maxIndx[indxPlot]+1,indxPlot],'ob',linewidth=2)
# plt.loglog(rhoInit[0:nImpctPrmtr,indxPlot], \
#            deltaEnrgIon_c[0:nImpctPrmtr,indxPlot],'-xr', \
#            rhoInitFit[0:maxIndx[indxPlot]+1,indxPlot], \
# 	   deltaEnrgIon_c_Fit[0:maxIndx[indxPlot]+1,indxPlot],'ob', \
#            rhoInitFit[0:maxIndx[indxPlot]+1,indxPlot], \
# 	   deltaEnrgIon_c_FitOpt[0:maxIndx[indxPlot]+1,indxPlot],'om',linewidth=2)
plt.loglog(rhoInit[0:nImpctPrmtr,indxPlot], \
           deltaEnrgIon_c[0:nImpctPrmtr,indxPlot],'-xr',linewidth=2)
plt.xlabel('Track Initial Impact Parameter $rho_{Init}$', \
           color='m',fontsize=16)
plt.ylabel('$\Delta E_{ion}$, $eV$', color='m',fontsize=16)
titleHeader = 'Transferred Energy $\Delta E_{ion}$ to Single Ion for '
titleHeader += ' $V_{ion}=%5.3f\cdot10^{%2d}$ $cm/s$'
plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=16)
# plt.xlim([.95*rhoInit[0,0],1.05*rhoInit[nImpctPrmtr,0]])
# plt.ylim([8.e-11,2.e-8])
plt.legend(['Calculated Data','Fitted Data'],loc='lower center',fontsize=16)
# plt.text(2.5e-3,1.e-8, \
#          ('Fitted $\Delta E_{ion}$ are proportional to $1/rho_{Init}^{%5.3f}$' % \
# 	  abs(fitBopt[indxPlot])),color='m',fontsize=16)
# plt.text(2.5e-3,1.e-8, \
#          ('Fitted $\Delta E_{ion}$ are proportional to $1/rho_{Init}^{%5.3f}$' % \
# 	  abs(fitB[indxPlot])),color='m',fontsize=16)
plt.grid(True)

yLimit=[1.5e-10,2.e-8]

indxPlot=1
VionCrrnt = vIon[indxPlot]
powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
fig700=plt.figure (700)
# plt.loglog(rhoInit[0:nImpctPrmtr,indxPlot], \
#            deltaEnrgIon_c[0:nImpctPrmtr,indxPlot],'-xr', \
#            rhoInitFit[0:maxIndx[indxPlot]+1,indxPlot], \
# 	   deltaEnrgIon_c_Fit[0:maxIndx[indxPlot]+1,indxPlot],'ob',linewidth=2)
# plt.loglog(rhoInit[0:nImpctPrmtr,indxPlot], \
#            deltaEnrgIon_c[0:nImpctPrmtr,indxPlot],'-xr', \
#            rhoInitFit[0:maxIndx[indxPlot]+1,indxPlot], \
# 	   deltaEnrgIon_c_Fit[0:maxIndx[indxPlot]+1,indxPlot],'ob', \
#            rhoInitFit[0:maxIndx[indxPlot]+1,indxPlot], \
# 	   deltaEnrgIon_c_FitOpt[0:maxIndx[indxPlot]+1,indxPlot],'om',linewidth=2)
plt.loglog(rhoInit[0:nImpctPrmtr,indxPlot], \
           deltaEnrgIon_c[0:nImpctPrmtr,indxPlot],'-xr',linewidth=2)
plt.xlabel('Track Initial Impact Parameter $rho_{Init}$', \
           color='m',fontsize=16)
plt.ylabel('$\Delta E_{ion}$, $eV$', color='m',fontsize=16)
titleHeader = 'Transferred Energy $\Delta E_{ion}$ to Single Ion for '
titleHeader += ' $V_{ion}=%5.3f\cdot10^{%2d}$ $cm/s$'
plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=16)
# plt.xlim([.95*rhoInit[0,0],1.05*rhoInit[nImpctPrmtr,0]])
# plt.ylim(yLimit)  # for "indxPlot=19"
plt.legend(['Calculated Data','Fitted Data'],loc='lower center',fontsize=16)
# plt.text(2.5e-3,1.e-8, \
#          ('Fitted $\Delta E_{ion}$ are proportional to $1/rho_{Init}^{%5.3f}$' % \
# 	  abs(fitBopt[indxPlot])),color='m',fontsize=16)
# plt.text(2.5e-3,1.e-8, \
#          ('Fitted $\Delta E_{ion}$ are proportional to $1/rho_{Init}^{%5.3f}$' % \
# 	  abs(fitB[indxPlot])),color='m',fontsize=16)
plt.grid(True)

indxPlot=3
VionCrrnt = vIon[indxPlot]
powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
fig800=plt.figure (800)
# plt.loglog(rhoInit[0:nImpctPrmtr,indxPlot], \
#            deltaEnrgIon_c[0:nImpctPrmtr,indxPlot],'-xr', \
#            rhoInitFit[0:maxIndx[indxPlot]+1,indxPlot], \
# 	   deltaEnrgIon_c_Fit[0:maxIndx[indxPlot]+1,indxPlot],'ob',linewidth=2)
# plt.loglog(rhoInit[0:nImpctPrmtr,indxPlot], \
#            deltaEnrgIon_c[0:nImpctPrmtr,indxPlot],'-xr', \
#            rhoInitFit[0:maxIndx[indxPlot]+1,indxPlot], \
# 	   deltaEnrgIon_c_Fit[0:maxIndx[indxPlot]+1,indxPlot],'ob', \
#            rhoInitFit[0:maxIndx[indxPlot]+1,indxPlot], \
# 	   deltaEnrgIon_c_FitOpt[0:maxIndx[indxPlot]+1,indxPlot],'om',linewidth=2)
plt.loglog(rhoInit[0:nImpctPrmtr,indxPlot], \
           deltaEnrgIon_c[0:nImpctPrmtr,indxPlot],'-xr',linewidth=2)
plt.xlabel('Track Initial Impact Parameter $rho_{Init}$', \
           color='m',fontsize=16)
plt.ylabel('$\Delta E_{ion}$, $eV$', color='m',fontsize=16)
titleHeader = 'Transferred Energy $\Delta E_{ion}$ to Single Ion for '
titleHeader += ' $V_{ion}=%5.3f\cdot10^{%2d}$ $cm/s$'
plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=16)
# plt.xlim([.95*rhoInit[0,0],1.05*rhoInit[nImpctPrmtr,0]])
# plt.ylim(yLimit)  # for "indxPlot=29"
plt.legend(['Calculated Data','Fitted Data'],loc='lower center',fontsize=16)
# plt.text(2.5e-3,1.e-8, \
#          ('Fitted $\Delta E_{ion}$ are proportional to $1/rho_{Init}^{%5.3f}$' % \
# 	  abs(fitBopt[indxPlot])),color='m',fontsize=16)
# plt.text(2.5e-3,1.e-8, \
#          ('Fitted $\Delta E_{ion}$ are proportional to $1/rho_{Init}^{%5.3f}$' % \
# 	  abs(fitB[indxPlot])),color='m',fontsize=16)
plt.grid(True)

indxPlot=4
VionCrrnt = vIon[indxPlot]
powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
fig900=plt.figure (900)
# plt.loglog(rhoInit[0:nImpctPrmtr,indxPlot], \
#            deltaEnrgIon_c[0:nImpctPrmtr,indxPlot],'-xr', \
#            rhoInitFit[0:maxIndx[indxPlot]+1,indxPlot], \
# 	   deltaEnrgIon_c_Fit[0:maxIndx[indxPlot]+1,indxPlot],'ob',linewidth=2)
# plt.loglog(rhoInit[0:nImpctPrmtr,indxPlot], \
#            deltaEnrgIon_c[0:nImpctPrmtr,indxPlot],'-xr', \
#            rhoInitFit[0:maxIndx[indxPlot]+1,indxPlot], \
# 	   deltaEnrgIon_c_Fit[0:maxIndx[indxPlot]+1,indxPlot],'ob', \
#            rhoInitFit[0:maxIndx[indxPlot]+1,indxPlot], \
# 	   deltaEnrgIon_c_FitOpt[0:maxIndx[indxPlot]+1,indxPlot],'om',linewidth=2)
plt.loglog(rhoInit[0:nImpctPrmtr,indxPlot], \
           deltaEnrgIon_c[0:nImpctPrmtr,indxPlot],'-xr',linewidth=2)
plt.xlabel('Track Initial Impact Parameter $rho_{Init}$', \
           color='m',fontsize=16)
plt.ylabel('$\Delta E_{ion}$, $eV$', color='m',fontsize=16)
titleHeader = 'Transferred Energy $\Delta E_{ion}$ to Single Ion for '
titleHeader += ' $V_{ion}=%5.3f\cdot10^{%2d}$ $cm/s$'
plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=16)
# plt.xlim([.95*rhoInit[0,0],1.05*rhoInit[nImpctPrmtr,0]])
# plt.ylim(yLimit)  # for "indxPlot=39"
plt.legend(['Calculated Data','Fitted Data'],loc='lower center',fontsize=16)
# plt.text(2.5e-3,1.e-8, \
#          ('Fitted $\Delta E_{ion}$ are proportional to $1/rho_{Init}^{%5.3f}$' % \
# 	  abs(fitBopt[indxPlot])),color='m',fontsize=16)
# plt.text(2.5e-3,1.e-8, \
#          ('Fitted $\Delta E_{ion}$ are proportional to $1/rho_{Init}^{%5.3f}$' % \
# 	  abs(fitB[indxPlot])),color='m',fontsize=16)
plt.grid(True)

yLimit=[1.5e-10,2.e-8]

indxPlot=5
VionCrrnt = vIon[indxPlot]
powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
fig1000=plt.figure (1000)
# plt.loglog(rhoInit[0:nImpctPrmtr,indxPlot], \
#            deltaEnrgIon_c[0:nImpctPrmtr,indxPlot],'-xr', \
#            rhoInitFit[0:maxIndx[indxPlot]+1,indxPlot], \
# 	   deltaEnrgIon_c_Fit[0:maxIndx[indxPlot]+1,indxPlot],'ob',linewidth=2)
# plt.loglog(rhoInit[0:nImpctPrmtr,indxPlot], \
#            deltaEnrgIon_c[0:nImpctPrmtr,indxPlot],'-xr', \
#            rhoInitFit[0:maxIndx[indxPlot]+1,indxPlot], \
# 	   deltaEnrgIon_c_Fit[0:maxIndx[indxPlot]+1,indxPlot],'ob', \
#            rhoInitFit[0:maxIndx[indxPlot]+1,indxPlot], \
# 	   deltaEnrgIon_c_FitOpt[0:maxIndx[indxPlot]+1,indxPlot],'om',linewidth=2)
# plt.loglog(rhoInit[0:nImpctPrmtr-1,indxPlot], \
#            deltaEnrgIon_c[0:nImpctPrmtr-1,indxPlot],'-xr', \
#            rhoInitFit[0:maxIndx[indxPlot],indxPlot], \
# 	   deltaEnrgIon_c_Fit[0:maxIndx[indxPlot],indxPlot],'ob', \
# 	   linewidth=2)
plt.loglog(rhoInit[0:nImpctPrmtr,indxPlot], \
           deltaEnrgIon_c[0:nImpctPrmtr,indxPlot],'-xr',linewidth=2)
plt.xlabel('Track Initial Impact Parameter $rho_{Init}$', \
           color='m',fontsize=16)
plt.ylabel('$\Delta E_{ion}$, $eV$', color='m',fontsize=16)
titleHeader = 'Transferred Energy $\Delta E_{ion}$ to Single Ion for '
titleHeader += ' $V_{ion}=%5.3f\cdot10^{%2d}$ $cm/s$'
plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=16)
# plt.xlim([.95*rhoInit[0,0],1.05*rhoInit[nImpctPrmtr,0]])
# plt.ylim(yLimit)  # for "indxPlot=49"
plt.legend(['Calculated Data','Fitted Data'],loc='lower center',fontsize=16)
# plt.text(2.5e-3,1.e-8, \
#          ('Fitted $\Delta E_{ion}$ are proportional to $1/rho_{Init}^{%5.3f}$' % \
# 	  abs(fitBopt[indxPlot])),color='m',fontsize=16)
# plt.text(2.5e-3,1.e-8, \
#          ('Fitted $\Delta E_{ion}$ are proportional to $1/rho_{Init}^{%5.3f}$' % \
# 	  abs(fitB[indxPlot])),color='m',fontsize=16)
plt.grid(True)

plt.show()



sys.exit()
