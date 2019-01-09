import os, sys
from numpy import fft
import numpy as np
import math
import matplotlib.pyplot as plt 
from matplotlib.legend_handler import HandlerLine2D
#
# Opening the input file: 
#
inputFile='frictionForceLong_Nsplit_5_ve_3.dat'
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

while True:
   lineData=inpfile.readline()
#   print ('line=%d: %s' % (lines,lineData))
   if not lineData:
      break
   lines += 1
   if lines == 11:
      words=lineData.split()
      nWords=len(words)
#      print ('Data from %d: words=%s, number of entries = %d' % (lines,words,nWords))
      Delta_e_par_arr[0] = float(words[1])
      Delta_e_par_arr[1] = float(words[2])
      Delta_e_par_arr[2] = float(words[3])
      print ('Delta_e_par_arr: %e, %e, %e' % \
             (Delta_e_par_arr[0],Delta_e_par_arr[1],Delta_e_par_arr[2]))
   if lines > 14:
      words=lineData.split()
      nWords=len(words)
#      print ('Data from %d: words=%s, number of entries = %d' % (lines,words,nWords))
      Vion[dataNumber]=float(words[0])
      frctnFrcLong[dataNumber,0]=float(words[1])
      frctnFrcLong[dataNumber,1]=float(words[2])
      frctnFrcLong[dataNumber,2]=float(words[3])
      dataNumber += 1

inpfile.close()
print ('Close input file "%s"' % inputFile)

v_i_numb = dataNumber
v_e_numb = 3
# print ('v_i_numb=%d' % v_i_numb)

Nsplit = 5
# Delta_e_par = 1.e5

powDelta_e_par = np.zeros(v_e_numb)
mantDelta_e_par = np.zeros(v_e_numb)
for k in range(v_e_numb):
   powDelta_e_par[k] = round(np.log10(Delta_e_par_arr[k])) 
   mantDelta_e_par[k] = Delta_e_par_arr[k]/(10**powDelta_e_par[k]) 

# fig40 = plt.figure(40)
# plt.plot(Vion[0:dataNumber], -frctnFrcLong[0:dataNumber],'-b', \
#          Vion[0:dataNumber], -frctnFrcLong[0:dataNumber],'xr') 
# plt.xlabel('$V_{i||}/v_{e||}}$',color='m',fontsize=14) 
# plt.ylabel('$-F_{\parallel} (eV / m)$',color='m',fontsize=14) 
# titleHeader = ('Longitudinal Friction Force: $V_{e||}=%3.1f\cdot10^{%2d}$ m/s' % \
#                (mantDelta_e_par,powDelta_e_par))
# plt.title(titleHeader,color='m',fontsize=14)
# plt.legend(['$N_{split}$ = %d' % Nsplit],loc='upper right',fontsize=14)
# plt.grid(True)

# fig50 = plt.figure(50)
# plt.plot(Vion[0:v_i_numb], -frctnFrcLong[0:v_i_numb,0],'-b', \
#          Vion[0:v_i_numb], -frctnFrcLong[0:v_i_numb,1],'-m', \
#          Vion[0:v_i_numb], -frctnFrcLong[0:v_i_numb,2],'-g', \
#          Vion[0:v_i_numb], -frctnFrcLong[0:v_i_numb,0],'xr', \
#          Vion[0:v_i_numb], -frctnFrcLong[0:v_i_numb,1],'xr', \
#          Vion[0:v_i_numb], -frctnFrcLong[0:v_i_numb,2],'xr') 
# plt.xlabel('$V_{i||}/v_{e||}}$',color='m',fontsize=14) 
# plt.ylabel('$-F_{\parallel} (eV / m)$',color='m',fontsize=14) 
# titleHeader = ('Longitudinal Friction Force: $N_{split}$ = %d' % (int(Nsplit)))
# plt.title(titleHeader,color='m',fontsize=14)
# plt.legend([(('$v_{e||}=%3.1f \cdot 10^{%2d}$ m/s' % (mantDelta_e_par[0],powDelta_e_par[0]))), \
#             (('$v_{e||}=%3.1f \cdot 10^{%2d}$ m/s' % (mantDelta_e_par[1],powDelta_e_par[1]))), \
#             (('$v_{e||}=%3.1f \cdot 10^{%2d}$ m/s' % (mantDelta_e_par[2],powDelta_e_par[2])))], \
#     	     loc='upper right',fontsize=14)
# plt.grid(True)
# # fig50.savefig('frictionForceLong_fig50.png')    
# # print ('File "frictionForceLong_fig50.png" is written')   

# fig55 = plt.figure(55)
# plt.semilogx(Delta_e_par_arr[0]*Vion[0:v_i_numb], -frctnFrcLong[0:v_i_numb,0],'-b', \
#              Delta_e_par_arr[1]*Vion[0:v_i_numb], -frctnFrcLong[0:v_i_numb,1],'-m', \
#              Delta_e_par_arr[2]*Vion[0:v_i_numb], -frctnFrcLong[0:v_i_numb,2],'-g', \
#              Delta_e_par_arr[0]*Vion[0:v_i_numb], -frctnFrcLong[0:v_i_numb,0],'xr', \
#              Delta_e_par_arr[1]*Vion[0:v_i_numb], -frctnFrcLong[0:v_i_numb,1],'xr', \
#              Delta_e_par_arr[2]*Vion[0:v_i_numb], -frctnFrcLong[0:v_i_numb,2],'xr') 
# plt.xlabel('$V_{i||}$, m/s',color='m',fontsize=14) 
# plt.ylabel('$-F_{\parallel} (eV / m)$',color='m',fontsize=14) 
# titleHeader = ('Longitudinal Friction Force: $N_{split}$ = %d' % (int(Nsplit)))
# plt.title(titleHeader,color='m',fontsize=14)
# plt.legend([(('$v_{e||}=%3.1f \cdot 10^{%2d}$ m/s' % (mantDelta_e_par[0],powDelta_e_par[0]))), \
#             (('$v_{e||}=%3.1f \cdot 10^{%2d}$ m/s' % (mantDelta_e_par[1],powDelta_e_par[1]))), \
#             (('$v_{e||}=%3.1f \cdot 10^{%2d}$ m/s' % (mantDelta_e_par[2],powDelta_e_par[2])))], \
#     	     loc='upper right',fontsize=14)
# plt.grid(True)
# # fig55.savefig('frictionForceLong_fig55.png')    
# # print ('File "frictionForceLong_fig55.png" is written')   

# fig60 = plt.figure(60)
# plt.semilogx(Delta_e_par_arr[0]*Vion[0:v_i_numb], -frctnFrcLong[0:v_i_numb,1],'-b', \
#              Delta_e_par_arr[0]*Vion[0:v_i_numb], -frctnFrcLong[0:v_i_numb,1],'xr') 
# plt.xlabel("Electron's velocity $V_{e||}$ in the Ion Beam Frame, m/s", \
#            color='m',fontsize=14) 
# plt.ylabel('$-F_{\parallel} (eV / m)$',color='m',fontsize=14) 
# titleHeader = ('Longitudinal Friction Force: $N_{split}$ = %d, $T_e = 0$' % (int(Nsplit)))
# plt.title(titleHeader,color='m',fontsize=14)
# plt.grid(True)
# # fig60.savefig('frictionForceLong_fig60.png')    
# # print ('File "frictionForceLong_fig60.png" is written')   

n = 2*v_i_numb+1
ffLong_fft = np.zeros(n)
vi_arg = np.zeros(n)
for i in range(v_i_numb):
   vi_arg[v_i_numb+1+i] = Vion[i] 
   vi_arg[i] = -Vion[v_i_numb-i-1] 
   ffLong_fft[v_i_numb+1+i] = -frctnFrcLong[i,1]
   ffLong_fft[i] = frctnFrcLong[v_i_numb-i-1,1]

# for k in range(n):
#    print ('%d) v_i = %e,  ffLong = %e' % (k,vi_arg[k],ffLong_fft[k]))


Fk_ffLong = fft.fft(ffLong_fft)/n         # Fourier coefficients (divided by n)
# Fk = fft.fftshift(Fk)                     # Shift zero freq to center
# for i in range(n):
#    print i,') v_i = ',vi_arg[i],'  Fk_ffLong =',  Fk_ffLong[i]

ffLong_check = n*fft.ifft(Fk_ffLong)

powV0 = powDelta_e_par[1]
mantV0 = mantDelta_e_par[1] 

fig70 = plt.figure(70)
plt.plot(vi_arg, ffLong_fft,'-b',vi_arg, np.real(ffLong_check),'xr') 
plt.xlabel("Electron's velocity $V_{e||}/V_0$ in the Ion Beam Frame", \
           color='m',fontsize=14) 
plt.ylabel('$-F_{\parallel} (eV / m)$',color='m',fontsize=14) 
titleHeader = ('Longitudinal Friction Force: $T_e = 0, V_0=%3.1f \cdot 10^{%2d}$ m/s' % \
               (mantV0,powV0))
plt.title(titleHeader,color='m',fontsize=14)
plt.legend(['Before FFT','After FFT'],loc='lower right',fontsize=14)
plt.grid(True)
# fig70.savefig('frictionForceLong_fig70.png')    
# print ('File "frictionForceLong_fig70.png" is written')   


m_e=9.10938356e-31          # electron mass, kg
eVtoJ=1.6021766208e-19    # 1 eV = 1.6...e-12 erg

v0 = mantV0*np.power(10.,powV0)
eVrms = v0
eT_eV = .5*m_e*eVrms**2/eVtoJ
# print ('eVrms = %e m/s, T = %e eV' % (eVrms,eT_eV)) 

gauss_fft = np.zeros(n)
for i in range(n):
   gauss_fft[i] = np.exp(-.5*m_e*(vi_arg[i]*v0)**2/(eVtoJ*eT_eV))
gauss_fft = gauss_fft/sum(gauss_fft)   
#    print ('%d) v_i = %e,  gauss = %e' % (i,vi_arg[i],gauss_fft[i]))

# fig80 = plt.figure(80)
# plt.plot(vi_arg,gauss_fft,'-b') 
# plt.xlabel("Electron's velocity $V_{e||}/V_0$ in the Ion Beam Frame", \
#            color='m',fontsize=14) 
# plt.ylabel('Gauss Distribution',color='m',fontsize=14) 
# titleHeader = ('Gauss Distribution: $T_e = %5.3f eV, V_0=%3.1f \cdot 10^{%2d}$ m/s' % \
#                (eT_eV,mantV0,powV0))
# plt.title(titleHeader,color='m',fontsize=14)
# plt.grid(True)
# # fig80.savefig('frictionForceLong_fig80.png')    
# # print ('File "frictionForceLong_fig80.png" is written')   

Fk_gauss = fft.fft(gauss_fft)/n           # Fourier coefficients (divided by n)
# Fk = fft.fftshift(Fk)                     # Shift zero freq to center
# for i in range(n):
#    print i,') v_i = ',vi_arg[i],'  Fk_gauss =',  Fk_gauss[i]

gauss_check = n*fft.ifft(Fk_gauss)

fig90 = plt.figure(90)
plt.plot(vi_arg,gauss_fft,'-b',vi_arg,np.real(gauss_check),'xr') 
plt.xlabel("Electron's velocity $V_{e||}/V_0$ in the Ion Beam Frame", \
           color='m',fontsize=14) 
plt.ylabel('Gauss Distribution',color='m',fontsize=14) 
titleHeader = ('Gauss Distribution: $T_e = %5.3f eV, V_0=%3.1f \cdot 10^{%2d}$ m/s' % \
               (eT_eV,mantV0,powV0))
plt.title(titleHeader,color='m',fontsize=14)
plt.legend(['Before FFT','After FFT'],loc='upper right',fontsize=14)
plt.grid(True)
# fig90.savefig('frictionForceLong_fig90.png')    
# print ('File "frictionForceLong_fig90.png" is written')   


Fk_ffLong_vs_vi = Fk_ffLong*Fk_gauss 
# for i in range(n):
#    print i,') v_i = ',vi_arg[i],'  Fk_ffLong_vs_vi =',  Fk_ffLong_vs_vi[i]

ffLong_vs_vi = n**2*fft.ifft(Fk_ffLong_vs_vi) 

result = -np.real(ffLong_vs_vi)

fig100 = plt.figure(100)
plt.plot(Vion[0:v_i_numb], -frctnFrcLong[0:v_i_numb,1],'-b', \
         vi_arg[v_i_numb+1:n],result[v_i_numb+1:n],'xr') 
plt.xlabel("Electron's velocity $V_{e||}/V_0$ in the Ion Beam Frame", \
           color='m',fontsize=14) 
plt.ylabel('$-F_{\parallel} (eV / m)$',color='m',fontsize=14) 
titleHeader = ('Longitudinal Friction Force: $V_0=%3.1f \cdot 10^{%2d}$ m/s' % \
               (mantV0,powV0))
plt.title(titleHeader,color='m',fontsize=14)
plt.legend(['$T_e = 0$',('$T_e = %5.3f$ eV' % (eT_eV))],loc='upper right',fontsize=14)
plt.grid(True)
fig100.savefig('frictionForceLong_fig100.png')    
print ('File "frictionForceLong_fig100.png" is written')   

plt.show()

sys.exit()
