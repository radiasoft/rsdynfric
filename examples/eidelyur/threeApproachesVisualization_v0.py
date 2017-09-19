# from __future__ import division

#-------------------------------------
#
#        Started at 09/19/2017 (YuE)
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
dpxDataHeaderLineNumber=0                                          # Serial number of line with header for dpxData
dpyDataHeaderLineNumber=0                                          # Serial number of line with header for dpyData
dpzDataHeaderLineNumber=0                                          # Serial number of line with header for dpzData
lastLineNumber=0                                                   # Number of the last line 
xAdataFlag=0                                                       # =1 when xA-Data already read
yBdataFlag=0                                                       # =1 when yB-Data already read 
dpxDataFlag=0                                                      # =1 when dpxData already read
dpyDataFlag=0                                                      # =1 when dpyData already read
dpzDataFlag=0                                                      # =1 when dpzData already read
xAdata=[]                                                          # Array of xA-Data
yBdata=[]                                                          # Array of yB-Data
dpxData=[]                                                         # Array of dpxData
dpyData=[]                                                         # Array of dpyData
dpzData=[]                                                         # Array of dpzData

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
#   
# xAdata=log10(Upot/Ekin): 
#
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
      dpxDataHeaderLineNumber=yBheaderLineNumber+linesFull+3
      if entriesRem > 0:
         dpxDataHeaderLineNumber += 1
#      print 'dpxDataHeaderLineNumber=%d' % dpxDataHeaderLineNumber
      dataNumber=0
   if xAdataFlag == 1 and yBdataFlag == 0:
#   
# yBdata=log10(R_larm/b): 
#
      if lines >  yBheaderLineNumber+1 and lines <= dpxDataHeaderLineNumber-2:
         words=lineData.split()
         nWords=len(words)
#          print 'yB-Data from %d: words=%s, number of entries = %d' % (lines,words,nWords)
         for m in range(nWords):
            wordCrrnt=words[m].split(",")
            yBdata[dataNumber]=float(wordCrrnt[0])
	    dataNumber += 1
      if lines == dpxDataHeaderLineNumber-2:
         yBdataFlag=1   
         print 'yB-Data(%d entries) already read' % yBentries  
# Header for dpxData:
   if lines == dpxDataHeaderLineNumber:
      words=lineData.split()
      indx1=words.index('Entries:')
      entriesA=int(words[indx1+1])
      indx2=words.index('x')
      entriesB=int(words[indx2+1])
      indx3=words.index('with')
      perLine=int(words[indx3+1])
#       print 'dpxData-Header from %d: words =%s, index1=%d, entriesA=%d, index2=%d, entriesB=%d, index3=%d, perLine=%d' % \
#             (lines,words,indx1,entriesA,indx2,entriesB,indx3,perLine)
      dpxData=np.zeros((entriesA,entriesB))
# If data are written without skiping of the repeated zero values:      
#       linesFull=entriesA*entriesB//perLine
#       entriesRem=entriesA*entriesB-perLine*linesFull
      mA=0
      mB=0
      total_mAmB=0
   if yBdataFlag == 1 and dpxDataFlag == 0:    
#   
# dpxData=zApprch1dpx 
# Format in input file: at the begining, the index mA for "xA-direction" increases and only after that index mB for "yB-direction")
#   
      if lines >  dpxDataHeaderLineNumber+2:
#            print 'line %d: "%s' % (lines,lineData)
         words=lineData.split()
         nWords=len(words)
#          print 'Data from %d: words=%s, number of entries = %d' % (lines,words,nWords)
	 if nWords == 0:
	    dpxDataFlag=1
	    dpyDataHeaderLineNumber=lines+1   
            print 'dpyDataHeaderLineNumber=%d' % dpyDataHeaderLineNumber
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
                  dpxData[mA,mB]=float(wordCrrnt[0])
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
		     dpxData[mA,mB]=0.
                     total_mAmB += 1
                     mA += 1
	             if mA  == entriesA:
	                mA=0
		        mB += 1
         if dpxDataFlag == 1:
            print 'dpxData(%d x %d entries) already read (total %d)' % (entriesA,entriesB,total_mAmB) 
# Header for dpyData:
   if lines == dpyDataHeaderLineNumber:
      words=lineData.split()
      indx1=words.index('Entries:')
      entriesA=int(words[indx1+1])
      indx2=words.index('x')
      entriesB=int(words[indx2+1])
      indx3=words.index('with')
      perLine=int(words[indx3+1])
#       print 'dpyData-Header from %d: words =%s, index1=%d, entriesA=%d, index2=%d, entriesB=%d, index3=%d, perLine=%d' % \
#             (lines,words,indx1,entriesA,indx2,entriesB,indx3,perLine)
      dpyData=np.zeros((entriesA,entriesB))
# If data are written without skiping of the repeated zero values:      
#       linesFull=entriesA*entriesB//perLine
#       entriesRem=entriesA*entriesB-perLine*linesFull
      mA=0
      mB=0
      total_mAmB=0
   if dpxDataFlag == 1 and dpyDataFlag == 0:    
#   
# dpyData=zApprch1dpy 
# Format in input file: at the begining, the index mA for "xA-direction" increases and only after that index mB for "yB-direction")
#   
      if lines >  dpyDataHeaderLineNumber+2:
#            print 'line %d: "%s' % (lines,lineData)
         words=lineData.split()
         nWords=len(words)
#          print 'Data from %d: words=%s, number of entries = %d' % (lines,words,nWords)
	 if nWords == 0:
	    dpyDataFlag=1
	    dpzDataHeaderLineNumber=lines+1   
            print 'dpzDataHeaderLineNumber=%d' % dpzDataHeaderLineNumber
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
                  dpyData[mA,mB]=float(wordCrrnt[0])
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
		     dpyData[mA,mB]=0.
                     total_mAmB += 1
                     mA += 1
	             if mA  == entriesA:
	                mA=0
		        mB += 1
         if dpyDataFlag == 1:
            print 'dpyData(%d x %d entries) already read (total %d)' % (entriesA,entriesB,total_mAmB) 
# Header for dpzData:
   if lines == dpzDataHeaderLineNumber:
      words=lineData.split()
      indx1=words.index('Entries:')
      entriesA=int(words[indx1+1])
      indx2=words.index('x')
      entriesB=int(words[indx2+1])
      indx3=words.index('with')
      perLine=int(words[indx3+1])
#       print 'dpzData-Header from %d: words =%s, index1=%d, entriesA=%d, index2=%d, entriesB=%d, index3=%d, perLine=%d' % \
#             (lines,words,indx1,entriesA,indx2,entriesB,indx3,perLine)
      dpzData=np.zeros((entriesA,entriesB))
# If data are written without skiping of the repeated zero values:      
#       linesFull=entriesA*entriesB//perLine
#       entriesRem=entriesA*entriesB-perLine*linesFull
      mA=0
      mB=0
      total_mAmB=0
   if dpyDataFlag == 1 and dpzDataFlag == 0:    
#   
# dpzData=zApprch1dpz 
# Format in input file: at the begining, the index mA for "xA-direction" increases and only after that index mB for "yB-direction")
#   
      if lines >  dpzDataHeaderLineNumber+2:
#            print 'line %d: "%s' % (lines,lineData)
         words=lineData.split()
         nWords=len(words)
#          print 'Data from %d: words=%s, number of entries = %d' % (lines,words,nWords)
	 if nWords == 0:
	    dpzDataFlag=1
# Not necessary:	    
# 	    dpzDataHeaderLineNumber=lines+1   
#             print 'dpzDataHeaderLineNumber=%d' % dpzDataHeaderLineNumber
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
                  dpzData[mA,mB]=float(wordCrrnt[0])
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
		     dpzData[mA,mB]=0.
                     total_mAmB += 1
                     mA += 1
	             if mA  == entriesA:
	                mA=0
		        mB += 1
         if dpzDataFlag == 1:
            print 'dpyData(%d x %d entries) already read (total %d)' % (entriesA,entriesB,total_mAmB)
	    break 
   
inpfile.close()
print 'Close input file "%s"' % inputFile

#
# Visualization of results
#

xAmin=min(xAdata)
xAmax=np.max(xAdata)

yBmin=min(yBdata)
yBmax=max(yBdata)


plt.figure(55)
plt.plot(xAdata,yBdata,'.r')
plt.xlim([xAmin-.1,xAmax+.1])
plt.ylim([yBmin-.1,yBmax+.1])
plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
plt.ylabel('$B=log_{10}(R_L/b)$',color='m',fontsize=16)
plt.title('Map for Transfered Momenta $dP_x,dP_y,dP_z$ Calculations', color='m',fontsize=20)
plt.grid(True)

X,Y=np.meshgrid(xAdata,yBdata) 
fig70=plt.figure(70)
ax70=fig70.gca(projection='3d')
# surf=ax70.plot_surface(X,Y,zApprch1dpx,cmap=cm.coolwarm,linewidth=0,antialiased=False)
surf=ax70.plot_surface(X,Y,dpxData,cmap=cm.jet,linewidth=0,antialiased=False)
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
X,Y=np.meshgrid(xAdata,yBdata) 
mapDpx=ax.contourf(X,Y,dpxData)   
# mapDpx=ax.contourf(X,Y,dpxApprch_1,levels)   
# Contourrange=[int(NlarmCutofDown+1)]
# mapTurnCutoff=ax.contour(X,Y,Nlarm,Contourrange,format='%d',colors=('w'),linewidths=(2)) 
# plt.clabel(mapTurnCutoff,inline=1,fontsize=14,manual=[(-3,-1.5)])  
plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
plt.ylabel('$B=log_{10}(R_L/b)$',color='m',fontsize=16)
plt.title('$C_x (cm^{-2})$ for Transf. Momntm. $dP_x$: $dP_x=q_e^2/b\cdot C_x$', color='m',fontsize=20)
fig75.colorbar(mapDpx)

X,Y=np.meshgrid(xAdata,yBdata) 
fig80=plt.figure(80)
ax80=fig80.gca(projection='3d')
# surf=ax80.plot_surface(X,Y,zApprch1dpy,cmap=cm.coolwarm,linewidth=0,antialiased=False)
surf=ax80.plot_surface(X,Y,dpyData,cmap=cm.jet,linewidth=0,antialiased=False)
plt.title('Transfered Momentum $dP_y$:\n$dP_y=q_e^2/b \cdot C_y$', color='m',fontsize=20)
plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
plt.ylabel('$B=log_{10}(R_L/b)$',color='m',fontsize=16)
ax80.set_zlabel('$C_y$, $cm^{-2}$',color='m',fontsize=16)
cb = fig80.colorbar(surf)
# cbar=fig80.colorbar(surf,ticks=[0,1000,2000,3000,4000,5000,6000,7000])  # Next 2 commands not work
# cbar.ax.set_yticklabels(['0','1000','2000','3000','4000','5000','6000','7000'])
# labels=np.arange(0,8000,1000)               # Next 4 commands not work
# location=labels
# cb.set_ticks(location)
# cb.set_ticklabels(labels)
# tick_locator = ticker.MaxNLocator(nbins=10) # Next 3 commands not work
# cb.locator = tick_locator
# cb.update_ticks()
plt.grid(True)

fig85=plt.figure(85)
ax=fig85.add_subplot(111)         # for contours poltting
X,Y=np.meshgrid(xAdata,yBdata) 
mapDpy=ax.contourf(X,Y,dpyData)   
# mapDpy=ax.contourf(X,Y,dpyApprch_1,levels)   
# Contourrange=[int(NlarmCutofDown+1)]
# mapTurnCutoff=ax.contour(X,Y,Nlarm,Contourrange,format='%d',colors=('w'),linewidths=(2)) 
# plt.clabel(mapTurnCutoff,inline=1,fontsize=14,manual=[(-3,-1.5)])  
plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
plt.ylabel('$B=log_{10}(R_L/b)$',color='m',fontsize=16)
plt.title('$C_y (cm^{-2})$ for Transf. Momntm. $dP_y$: $dP_y=q_e^2/b\cdot C_y$', color='m',fontsize=20)
fig85.colorbar(mapDpy)

X,Y=np.meshgrid(xAdata,yBdata) 
fig90=plt.figure(90)
ax90=fig90.gca(projection='3d')
# surf=ax90.plot_surface(X,Y,zApprch1dpy,cmap=cm.coolwarm,linewidth=0,antialiased=False)
surf=ax90.plot_surface(X,Y,dpzData,cmap=cm.jet,linewidth=0,antialiased=False)
plt.title('Transfered Momentum $dP_z$:\n$dP_z=q_e^2/b \cdot C_z$', color='m',fontsize=20)
plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
plt.ylabel('$B=log_{10}(R_L/b)$',color='m',fontsize=16)
ax80.set_zlabel('$C_z$, $cm^{-2}$',color='m',fontsize=16)
cb = fig90.colorbar(surf)
# cbar=fig90.colorbar(surf,ticks=[0,1000,2000,3000,4000,5000,6000,7000])  # Next 2 commands not work
# cbar.ax.set_yticklabels(['0','1000','2000','3000','4000','5000','6000','7000'])
# labels=np.arange(0,8000,1000)               # Next 4 commands not work
# location=labels
# cb.set_ticks(location)
# cb.set_ticklabels(labels)
# tick_locator = ticker.MaxNLocator(nbins=10) # Next 3 commands not work
# cb.locator = tick_locator
# cb.update_ticks()
plt.grid(True)

fig95=plt.figure(95)
ax=fig95.add_subplot(111)         # for contours poltting
X,Y=np.meshgrid(xAdata,yBdata) 
mapDpz=ax.contourf(X,Y,dpzData)   
# mapDpz=ax.contourf(X,Y,dpzApprch_1,levels)   
# Contourrange=[int(NlarmCutofDown+1)]
# mapTurnCutoff=ax.contour(X,Y,Nlarm,Contourrange,format='%d',colors=('w'),linewidths=(2)) 
# plt.clabel(mapTurnCutoff,inline=1,fontsize=14,manual=[(-3,-1.5)])  
plt.xlabel('$A=log_{10}(q_e^2/b/E_{kin})$',color='m',fontsize=16)
plt.ylabel('$B=log_{10}(R_L/b)$',color='m',fontsize=16)
plt.title('$C_z (cm^{-2})$ for Transf. Momntm. $dP_z$: $dP_z=q_e^2/b\cdot C_z$', color='m',fontsize=20)
fig95.colorbar(mapDpz)

plt.show()   

sys.exit()

