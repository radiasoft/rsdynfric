def main(): 

# Yury's results (for sourceFlag = 1), othersise Ilya's results (for sourceFlag = 2): 
   sourceFlag = 1   
   
#
# Opening the input file: 
#
   inputFile='frictionForceLong_Nsplit_5_ve_3.dat'
   inputFile='HESR_coldBeam.dat'
   print ('Open input file "%s"...' % inputFile)
   inpfileFlag=0
   try:
      inpfile = open(inputFile,'r')
      inpfileFlag=1
   except:
      print ('Problem to open input file "%s"' % inputFile)
   if inpfileFlag == 1:
      print ('No problem to open input file "%s"' % inputFile)

   lines=0                           # Number of current line from input file   
   dataNumber=0                      # Number of current value of any types of Data

   Delta_e_par_arr=np.zeros(3)
   Vion=np.zeros(100)
   frctnFrcLong=np.zeros((100,3))
   frctnFrcLong_P=np.zeros(100)

   while True:
      lineData=inpfile.readline()
#      print ('line=%d: %s' % (lines,lineData))
      if not lineData:
         break
      lines += 1
##      if lines == 11:
##         words=lineData.split()
##         nWords=len(words)
###         print ('Data from %d: words=%s, number of entries = %d' % (lines,words,nWords))
##         Delta_e_par_arr[0] = float(words[1])
##         Delta_e_par_arr[1] = float(words[2])
##         Delta_e_par_arr[2] = float(words[3])
##         print ('Delta_e_par_arr: %e, %e, %e' % \
##                (Delta_e_par_arr[0],Delta_e_par_arr[1],Delta_e_par_arr[2]))
      if lines > 13:
         words=lineData.split()
         nWords=len(words)
#         print ('Data from %d: words=%s, number of entries = %d' % (lines,words,nWords))
         Vion[dataNumber]=float(words[0])
         frctnFrcLong[dataNumber,0]=float(words[1])
         frctnFrcLong_P[dataNumber]=float(words[2])
##         frctnFrcLong[dataNumber,1]=float(words[2])
##         frctnFrcLong[dataNumber,2]=float(words[3])
         dataNumber += 1

   inpfile.close()
   print ('Close input file "%s"' % inputFile)

#   sys.exit()

   v_i_numb = dataNumber
   v_e_numb = 3
#    print ('v_i_numb=%d' % v_i_numb)

   Nsplit = 5
#    
#   Delta_e_par = 1.e5
   Delta_e_par = .6e5

#   powDelta_e_par = np.zeros(v_e_numb)
#   mantDelta_e_par = np.zeros(v_e_numb)
#   for k in range(v_e_numb):
#      powDelta_e_par[k] = round(np.log10(Delta_e_par_arr[k])) 
#      mantDelta_e_par[k] = Delta_e_par_arr[k]/(10**powDelta_e_par[k]) 

#   powV0 = powDelta_e_par[1]
#   mantV0 = mantDelta_e_par[1] 
   powV0 =  round(np.log10(Delta_e_par))
   mantV0 =  Delta_e_par/(10**powV0)

   m_e=9.10938356e-31          # electron mass, kg
   eVtoJ=1.6021766208e-19    # 1 eV = 1.6...e-12 erg

   v0 = mantV0*np.power(10.,powV0)
   eVrms = [.1*v0,1.*v0,4.5*v0]
   eT_eV = .5*m_e*np.power(eVrms,2)/eVtoJ
   for k in range(3):
      print ('eVrms = %e m/s, T = %e eV' % (eVrms[k],eT_eV[k])) 
   
#
# Picture for checking of Yury's input data:
#
   fig40 = plt.figure(40)
   plt.plot(Vion[0:dataNumber], -frctnFrcLong[0:dataNumber,0],'-xr') 
   plt.plot(Vion[0:dataNumber], -frctnFrcLong_P[0:dataNumber],'-xb') 
   plt.xlabel('$V_{i\parallel}/v_{e\parallel}}$',color='m',fontsize=14) 
   plt.ylabel('$-F_{\parallel} (eV / m)$',color='m',fontsize=14) 
   titleHeader = ('HESR Longitudinal Friction Force: $V_{e||}=%3.1f\cdot10^{%2d}$ m/s' % \
#                   (mantDelta_e_par[1],powDelta_e_par[1]))
                  (mantV0,powV0))
   plt.title(titleHeader,color='m',fontsize=14)
   plt.legend(['Simulation','"Parkhomchuk"'],loc='upper right',fontsize=14)
   plt.grid(True)

#   plt.show()
   
#   sys.exit()

# =========================     Ponomarev's data:    ==================
#
#    V_i = np.array([0., 5.e+3, 15.e+3, 25.e+3, 35.e+3, 45.e+3, 55.e+3, 65.e+3, 75.e+3, 85.e+3, 95.e+3, 
#                    105.e+3, 155.e+3, 205.e+3, 255.e+3, 305.e+3])  # m/s 
#    Fpar = np.array([0.0, 822.823619243 +17.9896920134, 2312.80389435 +53.3197380211, 3399.66269332 \
#                     +85.5985841022, 3964.28354122 +115.916869842, 4040.93182278 +142.790579008, \
#    		 3779.38657397 +157.326074971, 3328.16981233 +168.937215016, 2881.84622455 \
#    		 +183.098201561, 2477.14550915 +193.697487947, 2181.01256832 +118.453878358, \
#    		 1854.1191578 +123.625249815, 908.83863436 +129.464534718, 502.035082987 \
#    		 +120.367640755,  303.212901962 +105.948844318, 197.61917079 \
#    		 +91.3366767275])                                       # eV / m 


#
# Picture for comparison of input Yury's and Ilya's tabled data:
#
#    fig45 = plt.figure(45)
#    plt.plot(Vion[0:dataNumber], -frctnFrcLong[0:dataNumber,1],'-b',1.e-5*V_i,Fpar,'xr') 
#    plt.xlabel('$V_i/V_0$',color='m',fontsize=14) 
#    plt.ylabel('$-F_{\parallel} (eV / m)$',color='m',fontsize=14) 
#    titleHeader = ('Longitudinal Friction Force: $V_0=%3.1f\cdot10^{%2d}$ m/s' % \
#                   (mantDelta_e_par[1],powDelta_e_par[1]))
#    plt.title(titleHeader,color='m',fontsize=14)
#    plt.legend(['Yury','Ilya'],loc='upper right',fontsize=14)
#    plt.grid(True)

   ffLong_intrpln_rslt = np.vectorize(ffLong_intrpln) 
   V_smln  = np.linspace(-9.e+5, 9.e+5, 3601) # note that the step in v is the same as in el. distribution 
   F_smln = ffLong_intrpln_rslt(V_smln) 

#
# Picture for checking of Ilya's input data:
#
#   plt.figure(47)
#   plt.plot(V_smln/v0, F_smln, '-b') 
#   plt.xlabel('Ion Velocity $V_{ion\parallel}$ (m/s)', fontsize=14) 
#   plt.ylabel('$F_{\parallel}$ (eV/m) ', fontsize=14) 
#   plt.title('Longitudinal Friction Force $-F_{\parallel}$ for $T_e = 0$',color='m',fontsize=14)
#   plt.grid(True)

#=======================================================================

#
# Auxiliary picture:
#
#   fig50 = plt.figure(50)
#   plt.plot(Vion[0:v_i_numb], -frctnFrcLong[0:v_i_numb,0],'-b', \
#            Vion[0:v_i_numb], -frctnFrcLong[0:v_i_numb,1],'-m', \
#            Vion[0:v_i_numb], -frctnFrcLong[0:v_i_numb,2],'-g', \
#            Vion[0:v_i_numb], -frctnFrcLong[0:v_i_numb,0],'xr', \
#            Vion[0:v_i_numb], -frctnFrcLong[0:v_i_numb,1],'xr', \
#            Vion[0:v_i_numb], -frctnFrcLong[0:v_i_numb,2],'xr') 
#   plt.xlabel('$V_{i\parallel}/v_{e\parallel}}$',color='m',fontsize=14) 
#   plt.ylabel('$-F_{\parallel} (eV / m)$',color='m',fontsize=14) 
#   titleHeader = ('Longitudinal Friction Force: $N_{split}$ = %d' % (int(Nsplit)))
#   plt.title(titleHeader,color='m',fontsize=14)
#   plt.legend([(('$v_{e\parallel}=%3.1f \cdot 10^{%2d}$ m/s' % (mantDelta_e_par[0],powDelta_e_par[0]))), \
#               (('$v_{e\parallel}=%3.1f \cdot 10^{%2d}$ m/s' % (mantDelta_e_par[1],powDelta_e_par[1]))), \
#               (('$v_{e\parallel}=%3.1f \cdot 10^{%2d}$ m/s' % (mantDelta_e_par[2],powDelta_e_par[2])))], \
#       	     loc='upper right',fontsize=14)
#   plt.grid(True)
#   fig50.savefig('frictionForceLong_fig50.png')    
#   print ('File "frictionForceLong_fig50.png" is written')   

#
# Auxiliary picture:
#
#   fig55 = plt.figure(55)
#   plt.semilogx(Delta_e_par_arr[0]*Vion[0:v_i_numb], -frctnFrcLong[0:v_i_numb,0],'-b', \
#                Delta_e_par_arr[1]*Vion[0:v_i_numb], -frctnFrcLong[0:v_i_numb,1],'-m', \
#                Delta_e_par_arr[2]*Vion[0:v_i_numb], -frctnFrcLong[0:v_i_numb,2],'-g', \
#                Delta_e_par_arr[0]*Vion[0:v_i_numb], -frctnFrcLong[0:v_i_numb,0],'xr', \
#                Delta_e_par_arr[1]*Vion[0:v_i_numb], -frctnFrcLong[0:v_i_numb,1],'xr', \
#                Delta_e_par_arr[2]*Vion[0:v_i_numb], -frctnFrcLong[0:v_i_numb,2],'xr') 
#   plt.xlabel('$V_{i\parallel}$, m/s',color='m',fontsize=14) 
#   plt.ylabel('$-F_{\parallel} (eV / m)$',color='m',fontsize=14) 
#   titleHeader = ('Longitudinal Friction Force: $N_{split}$ = %d' % (int(Nsplit)))
#   plt.title(titleHeader,color='m',fontsize=14)
#   plt.legend([(('$v_{e\parallel}=%3.1f \cdot 10^{%2d}$ m/s' % (mantDelta_e_par[0],powDelta_e_par[0]))), \
#               (('$v_{e\parallel}=%3.1f \cdot 10^{%2d}$ m/s' % (mantDelta_e_par[1],powDelta_e_par[1]))), \
#               (('$v_{e\parallel}=%3.1f \cdot 10^{%2d}$ m/s' % (mantDelta_e_par[2],powDelta_e_par[2])))], \
#       	     loc='upper right',fontsize=14)
#   plt.grid(True)
#   fig55.savefig('frictionForceLong_fig55.png')    
#   print ('File "frictionForceLong_fig55.png" is written')   

#
# Auxiliary picture:
#
#   fig60 = plt.figure(60)
#   plt.semilogx(Delta_e_par_arr[0]*Vion[0:v_i_numb], -frctnFrcLong[0:v_i_numb,1],'-b', \
#                Delta_e_par_arr[0]*Vion[0:v_i_numb], -frctnFrcLong[0:v_i_numb,1],'xr') 
#   plt.xlabel("Electron's velocity $V_{e\parallel}$ in the Ion Beam Frame, m/s", \
#              color='m',fontsize=14) 
#   plt.ylabel('$-F_{\parallel} (eV / m)$',color='m',fontsize=14) 
#   titleHeader = ('Longitudinal Friction Force: $N_{split}$ = %d, $T_e = 0$' % (int(Nsplit)))
#   plt.title(titleHeader,color='m',fontsize=14)
#   plt.grid(True)
#    # fig60.savefig('frictionForceLong_fig60.png')    
#    # print ('File "frictionForceLong_fig60.png" is written')   
   
   if (sourceFlag == 1):
      n = 2*v_i_numb+1
      ffLong_fft = np.zeros(n)
      vi_arg = np.zeros(n)
      for i in range(v_i_numb):
         vi_arg[v_i_numb+1+i] = Vion[i] 
         vi_arg[i] = -Vion[v_i_numb-i-1] 
         ffLong_fft[v_i_numb+1+i] = -frctnFrcLong[i,0]
         ffLong_fft[i] = frctnFrcLong[v_i_numb-i-1,0]

#       for k in range(n):
#          print ('%d) v_i = %e,  ffLong = %e' % (k,vi_arg[k],ffLong_fft[k]))

   if (sourceFlag == 2):
      n = np.size(V_smln)
      vi_arg = V_smln/v0
      ffLong_fft = F_smln


   Fk_ffLong = fft.fft(ffLong_fft)/n         # Fourier coefficients (divided by n)
#    Fk = fft.fftshift(Fk)                     # Shift zero freq to center
#    for i in range(n):
#       print i,') v_i = ',vi_arg[i],'  Fk_ffLong =',  Fk_ffLong[i]

   ffLong_check = n*fft.ifft(Fk_ffLong)

#
# Picture for checking of FFT and inverse FFT:
#
#   fig70 = plt.figure(70)
#   plt.plot(vi_arg, ffLong_fft,'-b',vi_arg, np.real(ffLong_check),'xr') 
#   plt.xlabel("Electron's velocity $V_{e||}/V_0$ in the Ion Beam Frame", \
#              color='m',fontsize=14) 
#   plt.ylabel('$-F_{\parallel} (eV / m)$',color='m',fontsize=14) 
#   titleHeader = ('Longitudinal Friction Force: $T_e = 0, V_0=%3.1f \cdot 10^{%2d}$ m/s' % \
#                  (mantV0,powV0))
#   plt.title(titleHeader,color='m',fontsize=14)
#   plt.legend(['Before FFT','After FFT'],loc='lower right',fontsize=14)
#   plt.grid(True)

   results= np.zeros((n,3))

   for k in range(3):
      gauss_fft = np.zeros(n)
      for i in range(n):
         gauss_fft[i] = np.exp(-.5*m_e*(vi_arg[i]*v0)**2/(eVtoJ*eT_eV[k]))
      gauss_fft = gauss_fft/sum(gauss_fft)   
#      print ('%d) v_i = %e,  gauss = %e' % (i,vi_arg[i],gauss_fft[i]))

      Fk_gauss = fft.fft(gauss_fft)/n           # Fourier coefficients (divided by n)
#       Fk = fft.fftshift(Fk)                     # Shift zero freq to center
#       for i in range(n):
#          print i,') v_i = ',vi_arg[i],'  Fk_gauss =',  Fk_gauss[i]

      gauss_check = n*fft.ifft(Fk_gauss)
#
# Picture for checking of FFT and inverse FFT:
#
#      fig90 = plt.figure(90+k)
#      plt.plot(vi_arg,gauss_fft,'-b',vi_arg,np.real(gauss_check),'xr') 
#      plt.xlabel("Electron's velocity $V_{e\parallel}}/V_0$ in the Ion Beam Frame", \
#                 color='m',fontsize=14) 
#      plt.ylabel('Gauss Distribution',color='m',fontsize=14) 
#      titleHeader = ('Gauss Distribution: $T_e = %5.3f eV, V_0=%3.1f \cdot 10^{%2d}$ m/s' % \
#                     (eT_eV[k],mantV0,powV0))
#      plt.title(titleHeader,color='m',fontsize=14)
#      plt.legend(['Before FFT','After FFT'],loc='upper right',fontsize=14)
#      plt.grid(True)

# --------------- My convolution: --------------------
#
#      Fk_ffLong_vs_vi = Fk_ffLong*Fk_gauss 
#       for i in range(n):
#          print i,') v_i = ',vi_arg[i],'  Fk_ffLong_vs_vi =',  Fk_ffLong_vs_vi[i]

#      ffLong_vs_vi = n**2*fft.ifft(Fk_ffLong_vs_vi) 
#      result = np.real(ffLong_vs_vi)
#      for i in range((n-1)/2):
#         results[i,k] = result[(n-1)/2-i]
#         results[(n-1)/2+1+i,k] = result[n-1-i]
#
# --------------- End of convolution: ----------------

# --------------- Python's convolution: --------------------
#
      ffLong_vs_vi = np.convolve(gauss_fft, ffLong_fft, mode='same') 

      results[:,k] = np.real(ffLong_vs_vi)

      if (sourceFlag == 1):
         xLimit = [-.2,3.]
      
      if (sourceFlag == 2):
         xLimit = [-.2,6.]

   fig100 = plt.figure(100+sourceFlag)
   plt.plot(vi_arg[(n-1)/2:n],ffLong_fft[(n-1)/2:n],'-b', \
            vi_arg[(n-1)/2:n],   results[(n-1)/2:n,0],'-r', \
   	    vi_arg[(n-1)/2:n],   results[(n-1)/2:n,1],'-g', \
   	    vi_arg[(n-1)/2:n],   results[(n-1)/2:n,2],'-m', \
	    Vion[0:dataNumber],-frctnFrcLong_P[0:dataNumber],'-xr')
   plt.xlim(xLimit) 
   plt.xlabel('Ion velocity $V_i/V_0$', color='m',fontsize=14) 
   plt.ylabel('$-F_{\parallel} (eV / m)$',color='m',fontsize=14) 
   titleHeader = ('HESR Longitudinal Friction Force $F_{\parallel}$: $V_0=%3.1f \cdot 10^{%2d}$ m/s' % \
                  (mantV0,powV0))
   plt.title(titleHeader,color='m',fontsize=14)
   plt.legend(['$eV_{rms}=0$ $(T_e = 0)$', \
               ('$eV_{rms}=%3.1f \cdot V_0$ $(T_e =%6.4f$ eV)' % (eVrms[0]/v0,eT_eV[0])), \
               ('$eV_{rms}=%3.1f \cdot V_0$ $(T_e =%5.3f$ eV)' % (eVrms[1]/v0,eT_eV[1])), \
   	       ('$eV_{rms}=%3.1f \cdot V_0$ $(T_e =%5.3f$ eV)' % (eVrms[2]/v0,eT_eV[2])), \
	       ('"Parkhomchuk"')], \
   	      loc='upper right',fontsize=14)
   plt.grid(True)
   figureFile = 'HESRsimulation_fig10'+str(sourceFlag)+'.png'
   print 'File  with figure is "',figureFile,'"'   
   fig100.savefig(figureFile)    
   print 'File "', figureFile,'" is written'  
   
#   someChecking() 

   plt.show()

   sys.exit()

def ffLong_intrpln(v): 
  
  s = np.sign(v) 
  v = np.abs(v) 
  
  V_i = np.array([0., 5.e+3, 15.e+3, 25.e+3, 35.e+3, 45.e+3, 55.e+3, 65.e+3, 75.e+3, \
                  85.e+3, 95.e+3, 105.e+3, 155.e+3, 205.e+3, 255.e+3, 305.e+3])    # m/s
  Fpar = np.array([0.0, 822.823619243+17.9896920134, 2312.80389435+53.3197380211, \
                       3399.66269332+85.5985841022,  3964.28354122+115.916869842, \
		       4040.93182278+142.790579008,  3779.38657397+157.326074971, \
		       3328.16981233+168.937215016,  2881.84622455+183.098201561, \
		       2477.14550915+193.697487947,  2181.01256832+118.453878358, \
		       1854.1191578+123.625249815,    908.83863436+129.464534718, \
		        502.035082987+120.367640755,  303.212901962+105.948844318, \
		        197.61917079+91.3366767275])                               # eV/m 
  
  F_interp = interp.interp1d(V_i, Fpar, kind='cubic') 

#  plt.figure(10)
#  plt.plot(V_i, Fpar, 'b-', label='cold electrons') 
#  plt.plot(V_i, F_interp, 'xr', label='After interpolation') 


  
  p105 = np.log(Fpar[12]/Fpar[11]) / np.log(V_i[12]/V_i[11]) 
  p155 = np.log(Fpar[13]/Fpar[12]) / np.log(V_i[13]/V_i[12]) 
  p205 = np.log(Fpar[14]/Fpar[13]) / np.log(V_i[14]/V_i[13]) 
  p255 = np.log(Fpar[15]/Fpar[14]) / np.log(V_i[15]/V_i[14]) 
  p305 = -2.0 
  
  if (0. <= v < 5.0e+3): 
    return s *v *Fpar[1] /V_i[1] 
  
  if (5.0e+3 <= v < 105.e+3): 
    return s *F_interp(v) 
  
  if (105.e+3 <= v < 155.e+3): 
    return s *Fpar[11] *np.power(v /V_i[11], p105) 
  
  
  if (155.e+3 <= v < 205.e+3): 
    return s *Fpar[12] *np.power(v /V_i[12], p155) 
  
  if (205.e+3 <= v < 255.e+3): 
    return s *Fpar[13] *np.power(v /V_i[13], p205) 
  
  if (255.e+3 <= v < 305.e+3): 
    return s *Fpar[14] *np.power(v /V_i[14], p255) 
  
  if (305.e+3 <= v ): 
    return s *Fpar[15] *np.power(v /V_i[15], p305)

def someChecking():

   I2 = 1./3.
   I4 = 1./5.
   I6 = 1./7.
   x10 = np.sqrt( (np.sqrt(5.)+2.)/(3.*np.sqrt(5.))) # -.7946544
   x20 = np.sqrt( (np.sqrt(5.)-2.)/(3.*np.sqrt(5.))) # -.1875924
   print 'x10 =',x10,', x20 = ',x20
   
   delta = 24.*x10*x20*(x10**2-x20**2)*(I2*(x10**2+x20**2)-(x10*x20)**2-I2**2)
   
   d1 = 2.*x20*(x20**2-I2)*(3.*(x20**2+I2)*(I4-I2**2)-2.*(I6-I2**3))/delta
   d2 = 2.*x10*(x10**2-I2)*(2.*(I6-I2**3)-3.*(x10**2+I2)*(I4-I2**2))/delta
   x3_2 = I2-2.*x10**2*d1-2.*x20**2*d2
   
   print 'd1 = ',d1,', d2 = ',d2,', x3^2 = ',x3_2
   
   x10 = x10+d1
   x20 = x20+x3_2  
   
   delta = 24.*x10*x20*(x10**2-x20**2)*(I2*(x10**2+x20**2)-(x10*x20)**2-I2**2)
   
   d1 = 2.*x20*(x20**2-I2)*(3.*(x20**2+I2)*(I4-I2**2)-2.*(I6-I2**3))/delta
   d2 = 2.*x10*(x10**2-I2)*(2.*(I6-I2**3)-3.*(x10**2+I2)*(I4-I2**2))/delta
   x3_2 = I2-2.*x10**2*d1-2.*x20**2*d2
   
   print 'd1 = ',d1,', d2 = ',d2,', x3^2 = ',x3_2
   
  
if __name__=="__main__": 
   import os, sys
   from numpy import fft
   import numpy as np
   import math
   import matplotlib.pyplot as plt 
   from matplotlib.legend_handler import HandlerLine2D
   from scipy import interpolate as interp 
   ffLong_intrpln_rslt = np.vectorize(ffLong_intrpln) 
   main() 
