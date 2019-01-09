import os, sys
from numpy import fft
import numpy as np
import matplotlib.pyplot as plt

n = 1000                                           # Number of data points
dx = 5.0                                           # Sampling period (in meters)
x = dx*np.arange(0,n)                              # x coordinates
w1 = 100.0                                         # wavelength (meters)
w2 = 20.0                                          # wavelength (meters)

fx = np.sin(2*np.pi*x/w1) + 2*np.cos(2*np.pi*x/w2) # signal

Fk = fft.fft(fx)/n                                 # Fourier coefficients (divided by n)
nu = fft.fftfreq(n,dx)                             # Natural frequencies
Fk = fft.fftshift(Fk)                              # Shift zero freq to center
nu = fft.fftshift(nu)                              # Shift zero freq to center

# f, ax = plt.subplots(3,1,sharex=True)
# ax[0].plot(nu,np.real(Fk))                         # Plot Cosine terms
# ax[0].set_ylabel(r'$Re[F_k]$', size = 'x-large')
# ax[1].plot(nu,np.imag(Fk))                         # Plot Sine terms
# ax[1].set_ylabel(r'$Im[F_k]$', size = 'x-large')
# ax[2].plot(nu,np.absolute(Fk)**2)                  # Plot spectral power
# ax[2].set_ylabel(r'$| F_k | ^2$', size = 'x-large')
# ax[2].set_xlabel(r'$\widetilde{\nu}$', size = 'x-large')


fx_check = n*fft.ifft(Fk)

plt.figure(1)
plt.plot(x[0:100],fx[0:100],'-r',x[0:100],np.real(fx_check[0:100]),'xb')

plt.figure(2)
plt.plot(x[0:100],np.real(fx_check[0:100]),'-xb')

plt.figure(3)
plt.plot(x[0:100],np.imag(fx_check[0:100]),'-xb')




plt.show()

sys.exit()
