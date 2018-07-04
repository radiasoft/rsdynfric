# from __future__ import division

#-------------------------------------
#
#        Started at 07/04/2018 (YuE)
#
#import os, sys
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
indxPlot=0
VionCrrnt = vIon[indxPlot]
powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
fig500=plt.figure (500)
# plt.loglog(rhoInit[0:nImpctPrmtr-1,indxPlot], \
#            deltaEnrgIon_c[0:nImpctPrmtr-1,indxPlot],'-xr', \
#            rhoInitFit[0:maxIndx[indxPlot],indxPlot], \
# 	   deltaEnrgIon_c_Fit[0:maxIndx[indxPlot],indxPlot],'ob', \
#            rhoInitFit[0:maxIndx[indxPlot],indxPlot], \
# 	   deltaEnrgIon_c_FitOpt[0:maxIndx[indxPlot],indxPlot],'om',linewidth=2)
plt.loglog(rhoInit[0:nImpctPrmtr-1,indxPlot], \
           deltaEnrgIon_c[0:nImpctPrmtr-1,indxPlot],'-xr',linewidth=2)
plt.xlabel('Track Initial Impact Parameter $rho_{Init}$', \
           color='m',fontsize=16)
plt.ylabel('$\Delta E_{ion}$, $eV$', color='m',fontsize=16)
titleHeader = 'Transferred Energy $\Delta E_{ion}$ to Single Ion for '
titleHeader += ' $V_{ion}=%5.3f\cdot10^{%2d}$ $cm/s$'
plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=16)
plt.xlim([.95*rhoInit[0,0],1.05*rhoInit[nImpctPrmtr-1,0]])
plt.ylim([3.e-11,2.e-8])
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
# plt.loglog(rhoInit[0:nImpctPrmtr-1,indxPlot], \
#            deltaEnrgIon_c[0:nImpctPrmtr-1,indxPlot],'-xr', \
#            rhoInitFit[0:maxIndx[indxPlot],indxPlot], \
# 	   deltaEnrgIon_c_Fit[0:maxIndx[indxPlot],indxPlot],'ob', \
#            rhoInitFit[0:maxIndx[indxPlot],indxPlot], \
# 	   deltaEnrgIon_c_FitOpt[0:maxIndx[indxPlot],indxPlot],'om',linewidth=2)
plt.loglog(rhoInit[0:nImpctPrmtr-1,indxPlot], \
           deltaEnrgIon_c[0:nImpctPrmtr-1,indxPlot],'-xr',linewidth=2)
plt.xlabel('Track Initial Impact Parameter $rho_{Init}$', \
           color='m',fontsize=16)
plt.ylabel('$\Delta E_{ion}$, $eV$', color='m',fontsize=16)
titleHeader = 'Transferred Energy $\Delta E_{ion}$ to Single Ion for '
titleHeader += ' $V_{ion}=%5.3f\cdot10^{%2d}$ $cm/s$'
plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=16)
plt.xlim([.95*rhoInit[0,0],1.05*rhoInit[nImpctPrmtr-1,0]])
plt.ylim([8.e-11,2.e-8])
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
# plt.loglog(rhoInit[0:nImpctPrmtr-1,indxPlot], \
#            deltaEnrgIon_c[0:nImpctPrmtr-1,indxPlot],'-xr', \
#            rhoInitFit[0:maxIndx[indxPlot],indxPlot], \
# 	   deltaEnrgIon_c_Fit[0:maxIndx[indxPlot],indxPlot],'ob', \
#            rhoInitFit[0:maxIndx[indxPlot],indxPlot], \
# 	   deltaEnrgIon_c_FitOpt[0:maxIndx[indxPlot],indxPlot],'om',linewidth=2)
plt.loglog(rhoInit[0:nImpctPrmtr-1,indxPlot], \
           deltaEnrgIon_c[0:nImpctPrmtr-1,indxPlot],'-xr',linewidth=2)
plt.xlabel('Track Initial Impact Parameter $rho_{Init}$', \
           color='m',fontsize=16)
plt.ylabel('$\Delta E_{ion}$, $eV$', color='m',fontsize=16)
titleHeader = 'Transferred Energy $\Delta E_{ion}$ to Single Ion for '
titleHeader += ' $V_{ion}=%5.3f\cdot10^{%2d}$ $cm/s$'
plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=16)
plt.xlim([.95*rhoInit[0,0],1.05*rhoInit[nImpctPrmtr-1,0]])
plt.ylim(yLimit)  # for "indxPlot=19"
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
# plt.loglog(rhoInit[0:nImpctPrmtr-1,indxPlot], \
#            deltaEnrgIon_c[0:nImpctPrmtr-1,indxPlot],'-xr', \
#            rhoInitFit[0:maxIndx[indxPlot],indxPlot], \
# 	   deltaEnrgIon_c_Fit[0:maxIndx[indxPlot],indxPlot],'ob', \
#            rhoInitFit[0:maxIndx[indxPlot],indxPlot], \
# 	   deltaEnrgIon_c_FitOpt[0:maxIndx[indxPlot],indxPlot],'om',linewidth=2)
plt.loglog(rhoInit[0:nImpctPrmtr-1,indxPlot], \
           deltaEnrgIon_c[0:nImpctPrmtr-1,indxPlot],'-xr',linewidth=2)
plt.xlabel('Track Initial Impact Parameter $rho_{Init}$', \
           color='m',fontsize=16)
plt.ylabel('$\Delta E_{ion}$, $eV$', color='m',fontsize=16)
titleHeader = 'Transferred Energy $\Delta E_{ion}$ to Single Ion for '
titleHeader += ' $V_{ion}=%5.3f\cdot10^{%2d}$ $cm/s$'
plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=16)
plt.xlim([.95*rhoInit[0,0],1.05*rhoInit[nImpctPrmtr-1,0]])
plt.ylim(yLimit)  # for "indxPlot=29"
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
# plt.loglog(rhoInit[0:nImpctPrmtr-1,indxPlot], \
#            deltaEnrgIon_c[0:nImpctPrmtr-1,indxPlot],'-xr', \
#            rhoInitFit[0:maxIndx[indxPlot],indxPlot], \
# 	   deltaEnrgIon_c_Fit[0:maxIndx[indxPlot],indxPlot],'ob', \
#            rhoInitFit[0:maxIndx[indxPlot],indxPlot], \
# 	   deltaEnrgIon_c_FitOpt[0:maxIndx[indxPlot],indxPlot],'om',linewidth=2)
plt.loglog(rhoInit[0:nImpctPrmtr-1,indxPlot], \
           deltaEnrgIon_c[0:nImpctPrmtr-1,indxPlot],'-xr',linewidth=2)
plt.xlabel('Track Initial Impact Parameter $rho_{Init}$', \
           color='m',fontsize=16)
plt.ylabel('$\Delta E_{ion}$, $eV$', color='m',fontsize=16)
titleHeader = 'Transferred Energy $\Delta E_{ion}$ to Single Ion for '
titleHeader += ' $V_{ion}=%5.3f\cdot10^{%2d}$ $cm/s$'
plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=16)
plt.xlim([.95*rhoInit[0,0],1.05*rhoInit[nImpctPrmtr-1,0]])
plt.ylim(yLimit)  # for "indxPlot=39"
plt.legend(['Calculated Data','Fitted Data'],loc='lower center',fontsize=16)
# plt.text(2.5e-3,1.e-8, \
#          ('Fitted $\Delta E_{ion}$ are proportional to $1/rho_{Init}^{%5.3f}$' % \
# 	  abs(fitBopt[indxPlot])),color='m',fontsize=16)
# plt.text(2.5e-3,1.e-8, \
#          ('Fitted $\Delta E_{ion}$ are proportional to $1/rho_{Init}^{%5.3f}$' % \
# 	  abs(fitB[indxPlot])),color='m',fontsize=16)
plt.grid(True)

indxPlot=5
VionCrrnt = vIon[indxPlot]
powVionCrrnt = math.floor(np.log10(VionCrrnt)) 
mantVionCrrnt = VionCrrnt/(10**powVionCrrnt) 
fig1000=plt.figure (1000)
# plt.loglog(rhoInit[0:nImpctPrmtr-1,indxPlot], \
#            deltaEnrgIon_c[0:nImpctPrmtr-1,indxPlot],'-xr', \
#            rhoInitFit[0:maxIndx[indxPlot],indxPlot], \
# 	   deltaEnrgIon_c_Fit[0:maxIndx[indxPlot],indxPlot],'ob', \
#            rhoInitFit[0:maxIndx[indxPlot],indxPlot], \
# 	   deltaEnrgIon_c_FitOpt[0:maxIndx[indxPlot],indxPlot],'om',linewidth=2)
plt.loglog(rhoInit[0:nImpctPrmtr-1,indxPlot], \
           deltaEnrgIon_c[0:nImpctPrmtr-1,indxPlot],'-xr',linewidth=2)
plt.xlabel('Track Initial Impact Parameter $rho_{Init}$', \
           color='m',fontsize=16)
plt.ylabel('$\Delta E_{ion}$, $eV$', color='m',fontsize=16)
titleHeader = 'Transferred Energy $\Delta E_{ion}$ to Single Ion for '
titleHeader += ' $V_{ion}=%5.3f\cdot10^{%2d}$ $cm/s$'
plt.title(titleHeader % (mantVionCrrnt,powVionCrrnt),color='m',fontsize=16)
plt.xlim([.95*rhoInit[0,0],1.05*rhoInit[nImpctPrmtr-1,0]])
plt.ylim(yLimit)  # for "indxPlot=49"
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
