# from __future__ import division

#####################################################
#
# Comparison of different methods of integrations
#
#####################################################

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
from scipy.integrate import quad, nquad, dblquad, fixed_quad, quadrature, romberg
                            

from scipy.constants import pi

#----------------------------------------------------
#
# Will be calculated integral from function x*exp(-x^2) over 
# interval (a,b); results will be compare with analytical one.
#
#----------------------------------------------------

a=0.1
b=2.15

#
#  Analytical result:
#
rsltAnltcl = .5*(math.exp(-a**2)-math.exp(-b**2))
print 'rsltAnltcl = %20.14e' % rsltAnltcl 

#
#  Python method "QUAD":
#
rsltQUAD = integrate.quad(lambda x: x*np.exp(-x**2), a, b)
print '\nrsltQUAD = %20.14e, error = %20.14e' % (rsltQUAD[0],rsltQUAD[1]) 

#
#  Python method "FIXED_QUAD":
#
rsltFIXED_QUAD,dummy = integrate.fixed_quad(lambda x: x*np.exp(-x**2),a,b,n=5)
print '\nrsltFIXED_QUAD = %20.14e' % rsltFIXED_QUAD

#
#  Python method "QUADRATURE":
#
rsltQUADRATURE = integrate.quadrature(lambda x: x*np.exp(-x**2),a,b, \
                                      tol=1.e-8,maxiter=20)
print '\nrsltQUADRATURE = %20.14e, error = %20.14e' % \
      (rsltQUADRATURE[0],rsltQUADRATURE[1])

#
#  Python method "ROMBERG":
#
rsltROMBERG = integrate.romberg(lambda x: x*np.exp(-x**2),a,b,tol=1.e-8)
print '\nrsltROMBERG = %20.14e' % rsltROMBERG

#----------------------------------------------------
#
#  My "Simpson" (trapeze) method (SM):
#
#----------------------------------------------------
pointsSM = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, \
            1000, 2000, 10000, 50000]

intgrlSM = np.zeros(18)
diffSM = np.zeros(18)

print 's'
for i in range(18):
   stepSM = (b-a)/(pointsSM[i]-1)
   xLeft = a
   fLeft = xLeft*math.exp(-xLeft**2)
   for n in range(pointsSM[i]):
      xRight = xLeft+stepSM
      fRight = xRight*math.exp(-xRight**2)
      intgrlSM[i] += .5*(fLeft+fRight)*stepSM
      xLeft = xRight
      fLeft = fRight
   diffSM[i]=100.*(intgrlSM[i]/rsltAnltcl-1.)
   print 'intgrlSM(%5d) = %20.14e' % (pointsSM[i],intgrlSM[i])
   
#----------------------------------------------------
#
#  My "Gauss-Kronrod" method (GK)
#
#----------------------------------------------------
#
# Data for GK:
#
#----------------------------------------------------
#
# Points (psi_i) and weigths (w_i) to integrate for interval from -1 to 1;
# These data are from William H. Beyer. "Handbook of Mathematical Science".
# 5th Edition, CRC Press, Inc, 1978.
#
# To integrate for interval from 0 to 1 it is necessary to change points
# psi_i with points ksi_i=(1+psi_i)/2;
#
# For method with order N for function f(x):
#      int_(-1)^1 = sum_1^N [w_i* f(psi_i)]
#
#----------------------------------------------------

psi5 =np.array([-0.9061798, -0.5384093,  0.0,       0.5384093, 0.9061798])
w5   =np.array([ 0.2369269,  0.4786287,  0.5688889, 0.4786287, 0.2369269])

psi6 =np.array([-0.9324695, -0.6612094, -0.2386192, 0.2386192, 0.6612094, \
                0.9324695])
w6   =np.array([ 0.1713245,  0.3607616,  0.4679139, 0.4679139, 0.3607616, \
                0.1713245])

psi7 =np.array([-0.9491079, -0.7415312, -0.4058452, 0.0,       0.4058452, \
                0.7415312,  0.9491079])
w7   =np.array([ 0.1294850,  0.2797054,  0.3818301, 0.4179592, 0.3818301, \
                0.2797054,  0.1294850])

psi8 =np.array([-0.9602899, -0.7966665, -0.5255324, -0.1834346, 0.1834346, \
                 0.5255324,  0.7966665,  0.9602899])
w8   =np.array([ 0.1012285,  0.2223810,  0.3137066,  0.3626838, 0.3626838, \
                 0.3137066,  0.2223810,  0.1012285])

psi9 =np.array([-0.9681602, -0.8360311, -0.6133714, -0.3242534, 0.0,       \
                 0.3242534,  0.6133714,   0.8360311, 0.9681602])
w9   =np.array([ 0.0812744,  0.1806482,  0.2606107,  0.3123471, 0.3302394, \
                 0.3123471,  0.2606107,  0.1806482,  0.0812744])

psi10=np.array([-0.9739065, -0.8650634, -0.6794096, -0.4333954, -0.1488743, \
                 0.1488743,  0.4333954,  0.6794096,  0.8650634,  0.9739065])
w10  =np.array([ 0.0666713,  0.1494513,  0.2190864,  0.2692602,  0.2955242, \
                 0.2955242,  0.2692602,  0.2190864,  0.1494513,  0.0666713])

psi11=np.array([-0.9782287, -0.8870626, -0.7301520, -0.5190961, -0.2695432, \
                 0.0,        0.2695432,  0.5190961,  0.7301520,  0.8870626, \
		 0.9782287])
w11  =np.array([ 0.0556686,  0.1255804,  0.1862902,  0.2331938,  0.2628045, \
                 0.2729251,  0.2628045,  0.2331938,  0.1862902,  0.1255804, \
		 0.0556686])

psi12=np.array([-0.9815606, -0.9041173, -0.7699027, -0.5873180, -0.3678315, \
                -0.1253334,  0.1253334,  0.3678315,  0.5873180,  0.7699027, \
		 0.9041173,  0.9815606])
w12  =np.array([ 0.0471753,  0.1069393,  0.1600783,  0.2031674,  0.2334925, \
                 0.2491470,  0.2491470,  0.2334925,  0.2031674,  0.1600783, \
		 0.1069393,  0.0471753])

psi16=np.array([-0.9894009, -0.9445750, -0.8656312, -0.7554044, -0.6178762, \
                -0.4580168, -0.2816036, -0.0950125,  0.0950125,  0.2816036, \
		 0.4580168,  0.6178762,  0.7554044,  0.8656312,  0.9445750, \
		 0.9894009])
w16  =np.array([ 0.0271525,  0.0622535,  0.0951585,  0.1246290,  0.1495960, \
                 0.1691565,  0.1826034,  0.1894506,  0.1894506,  0.1826034, \
		 0.1691565,  0.1495960,  0.1246290,  0.0951585,  0.0622535, \
		 0.0271525])

#----------------------------------------------------
#
#                I =int_a^b [x*F(x)]dx  
#
# Let y_i(x_i) = 2*x_i/(b-a) - (b+a)/(b-a), then
#    x_i(y_i) = (b-a)/2*y_i+(a+b)/2 and y_i = psi_i
# or x_i(y_i) = (b-a)/2*psi_i+(a+b)/2, so that 
# I = (b-a)/2 * int_(-1)^1 {x(y)*F[x(y)]}dy =
#   = (b-a)/2 * sum over i_1^16 w_i*x_i*F(x_i)  
#
#----------------------------------------------------
pointsGK = [5, 6, 7, 8, 9, 10, 11, 12, 16]

psi = np.zeros((16,9))
w   = np.zeros((16,9))

for n in range(pointsGK[0]):
   psi[n,0] = psi5[n]
   w[n,0]   = w5[n]

for n in range(pointsGK[1]):
   psi[n,1] = psi6[n]
   w[n,1]   = w6[n]

for n in range(pointsGK[2]):
   psi[n,2] = psi7[n]
   w[n,2]   = w7[n]

for n in range(pointsGK[3]):
   psi[n,3] = psi8[n]
   w[n,3]   = w8[n]

for n in range(pointsGK[4]):
   psi[n,4] = psi9[n]
   w[n,4]   = w9[n]

for n in range(pointsGK[5]):
   psi[n,5] = psi10[n]
   w[n,5]   = w10[n]

for n in range(pointsGK[6]):
   psi[n,6] = psi11[n]
   w[n,6]   = w11[n]

for n in range(pointsGK[7]):
   psi[n,7] = psi12[n]
   w[n,7]   = w12[n]

for n in range(pointsGK[8]):
   psi[n,8] = psi16[n]
   w[n,8]   = w16[n]


intgrlGK = np.zeros(9)
diffGK = np.zeros(9)

print ''
for i in range(9):
   for n in range(pointsGK[i]):
      y = psi[n,i]*(b-a)/2 + (b+a)/2
      intgrlGK[i] += .5*(b-a)*w[n,i]*y*math.exp(-y**2)
   diffGK[i]=1.e4*(intgrlGK[i]/rsltAnltcl-1.)
   print 'intgrlGK(%2d) = %20.14e; GK-anltcl=%20.14e' % \
         (pointsGK[i],intgrlGK[i],diffGK[i])
print ''

'''
fig10=plt.figure(10)
# plt.semilogx([pointsSM[0],pointsSM[9]],[rsltAnltcl,rsltAnltcl],'-r', \
#          pointsSM,intgrlSM,'-xb',linewidth=2)
plt.semilogx(pointsSM,diffSM,'-xr',linewidth=2)
plt.grid(True)
xLimit=[.9*pointsSM[0],1.1*pointsSM[17]]
plt.xlim(xLimit)
# yLimit=[0.,1.05*max(intgrlSM)]
if (min(diffSM) > 0):
   yLimit=[-.1,1.1*max(diffSM)]
else: 
   yLimit=[1.1*min(diffSM),1.1*max(diffSM)]
plt.ylim(yLimit)
plt.xlabel('Intervals',color='m',fontsize=16)
plt.ylabel('"Simpson"/Analytical-1, %',color='m',fontsize=16)
plt.title( \
    ('Integration: $\int_{a=%4.2f}^{b=%4.2f} x \cdot exp(-x^2) dx$' % (a,b)), \
    color='m',fontsize=16)

fig20=plt.figure(20)
plt.plot(pointsGK,diffGK,'-xr',linewidth=2)
plt.grid(True)
xLimit=[.9*pointsGK[0],1.1*pointsGK[8]]
plt.xlim(xLimit)
# yLimit=[-.1,1.1*max(diffSM)]
# plt.ylim(yLimit)
plt.xlabel('Intervals',color='m',fontsize=16)
plt.ylabel('"GK"/Analytical-1, %%',color='m',fontsize=16)
plt.title( \
   ('Integration: $\int_{a=%4.2f}^{b=%4.2f} x \cdot exp(-x^2) dx$' % (a,b)), \
   color='m',fontsize=16)

plt.show()

fig10.savefig('myTest_simpson_fig10.jpg')    
fig20.savefig('myTest_gauss-kronrod_fig20.jpg')    
'''

sys.exit()

'''
#----------------------------------------------------
#
# Results of integration on interval (a=0.1, b=2.15):
#
#----------------------------------------------------
rsltAnltcl = 4.90110819456894e-01

rsltQUAD = 4.90110819456894e-01, error = 6.52910891241611e-15

rsltFIXED_QUAD = 4.90087869004433e-01

rsltQUADRATURE = 4.90110819482838e-01, error = 2.70023226001115e-10

rsltROMBERG = 4.90110819457543e-01

intgrlSM(   10) = 4.88901849396376e-01
intgrlSM(   20) = 4.90974001226727e-01
intgrlSM(   30) = 4.90986106657393e-01
intgrlSM(   40) = 4.90876650177943e-01
intgrlSM(   50) = 4.90775245578819e-01
intgrlSM(   60) = 4.90693018838555e-01
intgrlSM(   70) = 4.90627202524878e-01
intgrlSM(   80) = 4.90573996490510e-01
intgrlSM(   90) = 4.90530348989691e-01
intgrlSM(  100) = 4.90494010141370e-01
intgrlSM(  200) = 4.90314971519685e-01
intgrlSM(  300) = 4.90249693323094e-01
intgrlSM(  400) = 4.90216012063713e-01
intgrlSM(  500) = 4.90195470771605e-01
intgrlSM( 1000) = 4.90153641637933e-01
intgrlSM( 2000) = 4.90132354502638e-01
intgrlSM(10000) = 4.90115146282058e-01
intgrlSM(50000) = 4.90111685614262e-01

intgrlGK( 5) = 4.90103742297861e-01; GK-anltcl=-1.44399159374142e-01
intgrlGK( 6) = 4.90111567055930e-01; GK-anltcl=1.52536733777353e-02
intgrlGK( 7) = 4.90110842416122e-01; GK-anltcl=4.68449719015496e-04
intgrlGK( 8) = 4.90110782013873e-01; GK-anltcl=-7.63970512407397e-04
intgrlGK( 9) = 4.90110875772042e-01; GK-anltcl=1.14902886982904e-03
intgrlGK(10) = 4.90107036669213e-01; GK-anltcl=-7.71822928846699e-02
intgrlGK(11) = 4.90110822415593e-01; GK-anltcl=6.03679595201356e-05
intgrlGK(12) = 4.90108630924974e-01; GK-anltcl=-4.46538177378830e-02
intgrlGK(16) = 4.90110811033566e-01; GK-anltcl=-1.71865782716552e-04
'''
