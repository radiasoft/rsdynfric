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
#
# Case: the same transversal (eVrmsTran) and different longidudinal velocities: 
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
# Case: the same longidudinal (eVrmsLong) and different transversal velocities: 
#
minLarmorTurns=5000

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
      timePath=halfIntrcnLength/eVrmsLong
      numbLarmor=int(timePath/T_larm)
      if numbLarmor < minLarmorTurns:
         numbLarmor=minLarmorTurns
      timePoints=numbLarmor*stepsNumberOnGyro
      pointTrack=np.zeros(timePoints)             
#      timeStep=timePath/timePoints                                 # sec
      if numbLarmor == minLarmorTurns:
         halfIntrcnLength=timePoints*timeStep*eVrmsLong            # cm
      print 'minImpctPar(mkm)=%e, rhoLarmour (mkm)=%e, halfIntrcnLength=%e, numbLarmor=%d' % \
            (1.e+4*minImpctPar,1.e+4*ro_larm,1.e+4*halfIntrcnLength,numbLarmor)
      matr_elec=solenoid_eMatrix(B_mag,timeStep)                   # matrix for electron for timeStep in magnetic field
# To draw only first trajectory (for checking only):
      if i == 0 and j == 0:
         elecCoor=np.zeros((3,timePoints))                         # points of trajectory; cm
# Current distance from origin of the coordinate system to electron along the trajectory; cm:
      bCrrnt=np.zeros(timePoints)                                  # cm
# Current log10 of two important ratios; dimensionless:
      larmR_bCrrnt=np.zeros(timePoints)                            # ratio R_larmor/b
      uPot_enrgKinCrrnt=np.zeros(timePoints)                       # ratio potential_energy/kinetic_energy
      for m in range(6):
          z_elecCrrnt[m]=0.                                        # Initial zero-vector for electron
      z_elecCrrnt[Ix]=minImpctPar                                  # x, cm
      z_elecCrrnt[Iz]=-halfIntrcnLength                            # z, cm
      z_elecCrrnt[Ipy]=m_elec*eVtran                               # py, g*cm/sec
      z_elecCrrnt[Ipz]=m_elec*eVrmsLong                            # pz, g*cm/sec
#-----------------------------------------------
# Main action - dragging of the current trajectory (for given i and j):
#
      for k in range(timePoints):
 	  z_elecCrrnt=matr_elec.dot(z_elecCrrnt)
# To draw only first trajectory (for checking only):
          if i==0 and j==0:
             for ic in (Ix,Iy,Iz):
                elecCoor[ic//2,pointTrack[j]]=z_elecCrrnt[ic] # cm
# Current distance from origin of the coordinate system to electron along the trajectory; cm:
 	  bCrrnt[k]=np.sqrt(z_elecCrrnt[0]**2+z_elecCrrnt[2]**2+z_elecCrrnt[4]**2)
# Current log10 of two important ratios:  
	  larmR_bCrrnt[k]=math.log10(ro_larm/bCrrnt[k])            # dimensionless 
	  uPot_enrgKinCrrnt[k]=math.log10((q_elec**2/bCrrnt[k])/kinEnergy)      # dimensionless   
          pointTrack[j] += 1
# End of gragging of the current trajectory	  
#-----------------------------------------------
      minLarmR_b[j]=min(larmR_bCrrnt)      
      maxLarmR_b[j]=max(larmR_bCrrnt)
#       print 'For j=%d larmR_bCrrnt: min=%e, max=%e' % (j,minLarmR_b[j],maxLarmR_b[j])      
      minUpot_enrgKin[j]=min(uPot_enrgKinCrrnt)
      maxUpot_enrgKin[j]=max(uPot_enrgKinCrrnt)
#       print 'For j=%d upot_enrgKinCrrnt: min=%e, max=%e' % (j,minUpot_enrgKin[j],maxUpot_enrgKin[j])      
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
      trackNumb=stepsImpctParMin_Rlarm*i+j 
      endTracks[trackNumb]=begTracks[trackNumb]+pointTrack[j]     # index of the start of the current track     
      print 'i=%d, j=%d: beg=%d, end=%d' % (i,j,begTracks[trackNumb],endTracks[trackNumb])# index of the end of current track
      begTracks[trackNumb+1]=endTracks[trackNumb]+1      
totalPoints=b.shape
ySize=totalPoints[0]
uPot_enrgKinSize=uPot_enrgKin.shape 
xSize=uPot_enrgKinSize[0]     
print 'totalPoints=ySize=%d, xSize=%d' % (ySize,xSize)

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
xEntries=100   # xSize
outfile.write ('\n        Ratio: eVca/Ekin ( x-Entries %d )\n\n' % xEntries)
k=0
for m in range(xEntries):   
   valCrrnt=uPot_enrgKin[m]
   strVal='{:f}'.format(valCrrnt)
   if k == 0:
      bLine=strVal
   else:
      bLine=bLine+', '+strVal
   k += 1
   if k == nInLine or m == xEntries-1:
      outfile.write ('%s\n' % bLine)
      k=0

yEntries=100   # ySize
outfile.write ('\n        Ratio: R_larmor/b ( y-Entries %d )\n\n' % yEntries)
k=0
for m in range(yEntries):  
   valCrrnt=larmR_b[m]
   strVal='{:f}'.format(valCrrnt)
   if k == 0:
      bLine=strVal
   else:
      bLine=bLine+', '+strVal
   k += 1
   if k == nInLine or m == yEntries-1:
      outfile.write ('%s\n' % bLine)
      k=0

entries=xEntries*yEntries
outfile.write ('\n        Approach1 Results:  deltaP ( Entries %d x %d )\n\n' % (xEntries,yEntries))
k=0
for m in range(entries):   
   valCrrnt=larmR_b[m]
   strVal='{:f}'.format(valCrrnt)
   if k == 0:
      bLine=strVal
   else:
      bLine=bLine+', '+strVal
   k += 1
   if k == nInLine or m == entries-1:
      outfile.write ('%s\n' % bLine)
      k=0
outfile.close()

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
xLinesFull=0                                                       # Number of fully filled rows with x-Data
yLinesFull=0                                                       # Number of fully filled rows with y-Data
dataLinesFull=0                                                    # Number of fully filled rows with data
yHeaderLineNumber=0                                                # Number of line with header for y-Data
dataHeaderLineNumber=0                                             # Number of line with header for data
xDataFlag=0                                                        # =1 when x-Data already read
yDataFlag=0                                                        # =1 when y-Data already read 
dataFlag=0                                                         # =1 when data already read
xNumber=0                                                          # number of current value of x-Data
yNumber=0                                                          # number of current value of y-Data
dataNumber=0                                                       # number of current value of data
lastDataLineNumber=0                                               # number of last line with data
xData=[]                                                           # array of x-Data
yData=[]                                                           # array of y-Data
data=[]                                                            # array of data
while True:
   lineData=inpfile.readline()
   if not lineData:
      break
   lines += 1
#    print '%d: %s' % (lines,lineData)
# Header for x-Data:
   if lines == 2: 
      words=lineData.split()
      indx=words.index('x-Entries')
      xEntries=int(words[indx+1])
#       print 'x-Header from %d: words =%s, index=%d, xEntries=%d' % (lines,words,indx,xEntries)
      xData=np.zeros(xEntries)
      xLinesFull=xEntries//10
      xEntriesRem=xEntries-10*xLinesFull
      yHeaderLineNumber=xLinesFull+5
      if xEntriesRem > 0:
         yHeaderLineNumber += 1
   if xDataFlag == 0:
# x-Data:
      if lines > 3 and lines <= yHeaderLineNumber-2:
         words=lineData.split()
         nWords=len(words)
#          print 'x-Data from %d: words=%s, number of entries = %d' % (lines,words,nWords)
         for m in range(nWords):
            wordCrrnt=words[m].split(",")
            xData[xNumber]=float(wordCrrnt[0])
	    xNumber += 1
      if lines == yHeaderLineNumber-2:
         xDataFlag=1   
         print 'x-Data already read'  
# Header for y-Data:
   if lines == yHeaderLineNumber:
      words=lineData.split()
      indx=words.index('y-Entries')
      yEntries=int(words[indx+1])
#      print 'y-Header from  %d: words =%s, index=%d, yEntries=%d' % (lines,words,indx,yEntries)
      yData=np.zeros(yEntries)
      yLinesFull=yEntries//10
      yEntriesRem=yEntries-10*yLinesFull
      dataHeaderLineNumber=yHeaderLineNumber+yLinesFull+3
      if yEntriesRem > 0:
         dataHeaderLineNumber += 1
   if xDataFlag == 1 and yDataFlag == 0:
# y-Data:
      if lines >  yHeaderLineNumber+1 and lines <= dataHeaderLineNumber-2:
         words=lineData.split()
         nWords=len(words)
#          print 'y-Data from %d: words=%s, number of entries = %d' % (lines,words,nWords)
         for m in range(nWords):
            wordCrrnt=words[m].split(",")
            yData[yNumber]=float(wordCrrnt[0])
	    yNumber += 1
      if lines == dataHeaderLineNumber-2:
         yDataFlag=1   
         print 'y-Data already read'  
# Header for data:
   if lines == dataHeaderLineNumber:
      words=lineData.split()
      indx=words.index('Entries')
      xEntries=int(words[indx+1])
      yEntries=int(words[indx+3])
      entries=xEntries*yEntries
#       print 'data-Header from %d: words =%s, index=%d, xEntries=%d, yEntries=%d' % (lines,words,indx,xEntries,yEntries)
      data=np.zeros(entries)
      dataLinesFull=entries//10
      dataEntriesRem=entries-10*dataLinesFull
      lastDataLineNumber=dataHeaderLineNumber+1+dataLinesFull
      if dataEntriesRem > 0:
         lastDataLineNumber += 1
# Data:
   if yDataFlag == 1 and dataFlag == 0:    
      if lines >  dataHeaderLineNumber+1 and lines <= lastDataLineNumber:
         words=lineData.split()
         nWords=len(words)
#          print 'Data from %d: words=%s, number of entries = %d' % (lines,words,nWords)
         for m in range(nWords):
            wordCrrnt=words[m].split(",")
            data[dataNumber]=float(wordCrrnt[0])
	    dataNumber += 1
      if lines == lastDataLineNumber:
         dataFlag=1   
         print 'Data already read'  
inpfile.close()
# print 'x-Data, total=', (xData,xEntries)     
# print 'y-Data, total=', (yData,yEntries)      
# print 'data, total=', (data,entries)      
print 'Length of x-Data =%d' % xEntries     
print 'Length of y-Data =%d' % yEntries      
	 
         



timeEnd=os.times()
cpuTime=1.e+6*float(timeEnd[0])-1.e+6*float(timeStart[0])   # CPU time , mks
print 'cpuTime = %e' % cpuTime
#
# Checking of the first trajectory:
#
turns=10                                                           # Number of larmorturns for drawing 
points=turns*stepsNumberOnGyro                                     # Number of points for drawing

fig10=plt.figure(10)
ax10=fig10.gca(projection='3d')
ax10.plot(1.e+4*elecCoor[0,0:points],1.e+4*elecCoor[1,0:points],1.e+4*elecCoor[2,0:points],'-r',linewidth=2)
plt.xlabel('x, $\mu m$',color='m',fontsize=16)
plt.xlabel('x, $\mu m$',color='m',fontsize=16)
plt.title(('First %d Larmor Turns' % turns),color='m',fontsize=16)

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
'''
lastTrack=stepsImpctParMin_Rlarm*stepsVtran-1
plt.figure(100)
plt.plot(uPot_enrgKin[0:endTracks[lastTrack]],larmR_b[0:endTracks[lastTrack]],'.r')

#
# Checking of data from the output/input file:
#
plt.figure(110)
plt.plot(xData,yData,'.r')

plt.show()   

sys.exit()
