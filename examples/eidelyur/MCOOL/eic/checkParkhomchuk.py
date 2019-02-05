def main(): 

   velEff = 1.
   nTeta = 100
   nVel = 100
   teta = np.zeros(nTeta)
   velRel = np.zeros(nVel)
   Plong2 = np.zeros(nTeta)
   Ptrnsv3 = np.zeros(nTeta)
   Plong4 = np.zeros((nTeta,nVel))
   Ptrnsv4 = np.zeros((nTeta,nVel))

   for i in range(nTeta):
      teta[i] = .5*np.pi*i/(nTeta-1)

   for i in range(nVel):
      velRelLog =-1.+2.*i/(nVel-1)
      velRel[i] = np.power(10.,velRelLog) 

   for i in range(nTeta):
      Plong2[i] = 3.*math.sin(teta[i])**2*math.cos(teta[i])
      Ptrnsv3[i] = (math.sin(teta[i])**2-2.*math.cos(teta[i])**2)*math.sin(teta[i])
      for k in range(nVel):
         Plong4[i,k]  = velRel[k]*np.cos(teta[i])/np.power(1.+velRel[k]**2,1.5)
         Ptrnsv4[i,k] = velRel[k]*np.sin(teta[i])/np.power(1.+velRel[k]**2,1.5)

   plt.figure(10)
   plt.plot(teta,Plong2,'-r',teta,Ptrnsv3,'-b')
   titleHeader = \
             'Components $FF_{\parallel}$ and $FF_{\perp}$ of Friction Force "Parjhomchuk"'
   plt.title(titleHeader,color='m',fontsize=14)
   plt.xlabel('Angle Relative Beam Axis, rad',color='m',fontsize=14)
   plt.ylabel('$FF_{\parallel}$ and $FF_{\perp}$',color='m',fontsize=14)
   plt.legend(['$FF_{\parallel}$ (from (2))','$FF_{\perp}$ (from (3))'], \
              loc='lower right',fontsize=14)
   plt.grid(True)

   plt.figure(20)
   plt.plot(np.log10(velRel),Plong4[0,:],'-r',np.log10(velRel),Ptrnsv4[0,:],'-b')
   titleHeader = \
             'Components $FF_{\parallel}$ and $FF_{\perp}$ of Friction Force "Parjhomchuk"'
   plt.title(titleHeader,color='m',fontsize=14)
   plt.xlabel('$log_{10}(V/V_{eff})$',color='m',fontsize=14)
   plt.ylabel('$FF_{\parallel}$ and $FF_{\perp}$',color='m',fontsize=14)
   plt.legend(['$FF_{\parallel}$ (from (4))','$FF_{\perp}$ (from (4))'], \
              loc='upper right',fontsize=14)
   plt.text(-.8,.05,'Angle Rlative Beam Axis = 0',color='m',fontsize=14)
   plt.grid(True)

   plt.figure(30)
   plt.plot(np.log10(velRel),Plong4[nTeta-1,:],'-r',np.log10(velRel),Ptrnsv4[nTeta-1,:],'-b')
   titleHeader = \
             'Components $FF_{\parallel}$ and $FF_{\perp}$ of Friction Force "Parjhomchuk"'
   plt.title(titleHeader,color='m',fontsize=14)
   plt.xlabel('$log_{10}(V/V_{eff})$',color='m',fontsize=14)
   plt.ylabel('$FF_{\parallel}$ and $FF_{\perp}$',color='m',fontsize=14)
   plt.legend(['$FF_{\parallel}$ (from (4))','$FF_{\perp}$ (from (4))'], \
              loc='upper right',fontsize=14)
   plt.text(-.9,.05,'Angle Rlative Beam Axis = $\pi/2$',color='m',fontsize=14)
   plt.grid(True)

   X_teta = np.zeros((nTeta,nVel))
   Y_vel = np.zeros((nTeta,nVel))
   for i in range(nTeta):
      for k in range(nVel):
         X_teta[i,k] = teta[i]
         Y_vel[i,k] = np.log10(velRel[k])

   fig40=plt.figure (40)
   ax = fig40.add_subplot(111)                    # for contours plotting
   mapPlong4 = ax.contourf(X_teta,Y_vel,Plong4,cmap='jet') 
   plt.ylabel('Relative Ion Velocity, $log_{10}(V/V_{eff})$',color='m',fontsize=14)
   plt.xlabel('Angle Relative Beam Axis, rad',color='m',fontsize=14)
   titleHeader = 'Longitudinal Friction Force "Parjhomchuk"'
   plt.title(titleHeader,color='m',fontsize=14)
   plt.text(.25,.8,'Angle Rlative Beam Axis = 0',color='w',fontsize=14)
   fig40.colorbar(mapPlong4)
   plt.grid(True)
#   fig40.savefig('mapDeltaEnrgFit_fig40.png')    
#   print ('File "mapDeltaEnrgFit_fig40.png" is written')

   fig50=plt.figure (50)
   ax = fig50.add_subplot(111)                    # for contours plotting
   mapPlong4 = ax.contourf(X_teta,Y_vel,Ptrnsv4,cmap='jet') 
   plt.ylabel('Relative Ion Velocity, $log_{10}(V/V_{eff})$',color='m',fontsize=14)
   plt.xlabel('Angle Relative Beam Axis, rad',color='m',fontsize=14)
   titleHeader = 'Transversal Friction Force "Parjhomchuk"'
   plt.title(titleHeader,color='m',fontsize=14)
   plt.text(.25,.8,'Angle Rlative Beam Axis = $\pi/2$',color='w',fontsize=14)
   fig50.colorbar(mapPlong4)
   plt.grid(True)
#   fig50.savefig('mapDeltaEnrgFit_fig50.png')    
#   print ('File "mapDeltaEnrgFit_fig50.png" is written')

   X,Y=np.meshgrid(teta,np.log10(velRel))      

   fig60=plt.figure(60)
   ax60=fig60.gca(projection='3d')
   surf=ax60.plot_surface(X,Y,Plong4,cmap=cm.jet,linewidth=0,antialiased=False)
   plt.ylabel('$log_{10}(V/V_{eff})$',color='m',fontsize=14)
   plt.xlabel('Angle, rad',color='m',fontsize=14)
   titleHeader = 'Longitudinal Friction Force "Parjhomchuk"'
   plt.title(titleHeader,color='m',fontsize=14)
   ax60.set_zlabel('$F_{\parallel}$',color='m',fontsize=14)
   fig60.colorbar(surf, shrink=1., aspect=10)
   plt.grid(True)
#   fig60.savefig('mapDeltaEnrgFit_fig60.png')    
#   print ('File "mapDeltaEnrgFit_fig60.png" is written')

   fig70=plt.figure(70)
   ax70=fig70.gca(projection='3d')
   surf=ax70.plot_surface(X,Y,Ptrnsv4,cmap=cm.jet,linewidth=0,antialiased=False)
   plt.ylabel('$log_{10}(V/V_{eff})$',color='m',fontsize=14)
   plt.xlabel('Angle, rad',color='m',fontsize=14)
   titleHeader = 'Transversal Friction Force "Parjhomchuk"'
   plt.title(titleHeader,color='m',fontsize=14)
   ax70.set_zlabel('$F_{\perp}$',color='m',fontsize=14)
   fig70.colorbar(surf, shrink=1., aspect=10)
   plt.grid(True)
#   fig70.savefig('mapDeltaEnrgFit_fig70.png')    
#   print ('File "mapDeltaEnrgFit_fig70.png" is written')

   plt.show()

   sys.exit()
  
if __name__=="__main__":
  import os, sys
  import numpy as np
  import math
  import matplotlib as mpl
  import matplotlib.pyplot as plt 
  from matplotlib.legend_handler import HandlerLine2D
  from matplotlib.colors import LogNorm
  from matplotlib import ticker
  from matplotlib import markers
  from mpl_toolkits.mplot3d import Axes3D
  from matplotlib import cm

  main() 
  
