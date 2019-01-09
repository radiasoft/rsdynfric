# This differs from the Fpar_1D_01(a).py in that the average of the force on the ion is computed over several 
# nearest-neighbour particles, matching for each physical particle either the lowest 3 or lowest 5 moments of 
# the (spherical) distribution's projection onto the z axis.  Note that this code is not optimized for performance: 
# A run on a laptop with 27k nearest neighbors takes ~ 20 mins, with 64k nearest neighbors ~ 50 mins (Nstep = 1000). 
# For small v, results converge for k_max ~ 1k.  

def main(): 
  
  Tsim = 4.0e-10  # s 
  Nstep = 1000  # number of timesteps 
  
  Z = 79 
  Bz = 5.0  # T (working in MKS units)
  n_e = 2.0e+15  # m^-3, in the beam frame 
  Delta_e_par = 1.0e+5  # m/s 
  Delta_e_tr = 4.2e+5  # m/s 
  
  frac_delta_e_par = 0.75  # ratio of the ion velocity to the thermal electron velocity (abs value) 
  # The fact that Delta_e_par is the rms thermal velocity is not used here; it really is only setting the vel. scale 
  V_i_par = frac_delta_e_par *Delta_e_par  # m/s 
  
  #k = 1  # k-th nearest neighbor, if contribs from only one particle are computed (Not done here; better to introduce k_min?) 
  k_max = 1000  # average the force on the ion from up to (and including) the k_max-th nearest neighbour 
  Nsplit = 5  # split each nearest neighbor into Nsplit "sub-particles", matching Nsplit+1 lowest z-moments of the isotropic distribution if Nsplit is even, Nsplit moments if Nsplit is odd; Nsplit must be <= 5, for now. 
   
  q_el = -1.6021766208e-19  # Coulombs 
  m_el = 9.10938356e-31  # kg 
  c0 = 299792458.  # m/s 
  r_cl_el = 2.8179403227e-15  # m, classical radius of electron 
		
  Zrcc = Z *r_cl_el *c0 *c0  #  Z * 253.2638458397924 m^3 / s^2 == Z*e*e /([4*pi*eps0]*m_e) 
		
  w_L = np.abs(q_el*Bz) /m_el  # magnitude of gyration frequency in MKS 
  T_L = 2.*np.pi /w_L 
  r_L = Delta_e_tr /w_L  # a "typical" gyro-radius 
  
  w_p = 56.4146*np.sqrt(n_e) # s^-1, electrom plasma frequency; n_e in m^-3 
  T_p = 2.*np.pi /w_p 
  
  
  t_arr, F_ave_arr  = aveOverNearestNeibs01(n_e, -1. *V_i_par, Zrcc, Tsim, Nstep, k_max, Nsplit) 
  
  F_ave_arr *= m_el /1.60218e-19  # eV / m = 1.60218e-19 N 
  
  F_time_aved = (Tsim /np.float64(Nstep)) *np.sum( F_ave_arr[1:-1] ) # time-integrated acceleration of electron(s) 
  F_time_aved /= Tsim  # average force on the ion during the interaction time (redundant *Tsim/Tsim for clarity) 
  
  print "Nsplit = ", Nsplit, ", n-r of nearest neighbors = ", k_max
  print "Vz_ion = ", V_i_par, " m/s, time-averaged longitudinal force = ", F_time_aved, ' eV / m' 
  
  plt.plot(t_arr, F_ave_arr, 'r-') 
  plt.plot(t_arr, 0.0*t_arr, 'g--') 
  plt.xlabel('t (s)') 
  plt.ylabel('$F_{\parallel} (eV / m)$') 
  plt.xlim(0., Tsim) 
  plt.show()
  


# For a sphere, dS(z)/dz = 2*pi*r = const; but instead of uniformly sampling the [-r,r], we are choosing the z 
# coordinates of the initial conditions so as to match the lowest-order z-moments of the spherical shell distribution:
# for Nsplit odd, z-moments up to order Nsplit are matched, for Nsplit even, z-moments up to order Nsplit+1 are matched 

def aveOverNearestNeibs01(n_e, v0, C0, T, Nstep, k_max, Nsplit=4): 
  
  dt = T /np.float64(Nstep) 
  t_arr = np.zeros(Nstep+2, dtype=np.float64) 
  t_arr[-1] = T 
  t_arr[1:Nstep+1] = 0.5*dt +np.linspace(0., T -dt, Nstep) 
  F_cumul_arr = 0.0 *t_arr 
  
  z_split = np.zeros(Nsplit, dtype=np.float64) 
  
  if Nsplit == 5: 
    z_split[0] = -np.sqrt( (5. +np.sqrt(11.)) /12. ) # -0.8 if uniform in z & not matching moments, -0.8324974870 as is 
    z_split[1] = -np.sqrt( (5. -np.sqrt(11.)) /12. ) # -.4  --"--, -0.37454141 as is 
    z_split[3] =  np.sqrt( (5. -np.sqrt(11.)) /12. ) #  .4  --"--,  0.37454141 as is 
    z_split[4] =  np.sqrt( (5. +np.sqrt(11.)) /12. ) #  .8  --"--,  0.8324974870 as is 
    
  if Nsplit == 4: 
    z_split[0] = -np.sqrt( (np.sqrt(5.)+2.)/(3.*np.sqrt(5.)) ) 
    z_split[1] = -np.sqrt( (np.sqrt(5.)-2.)/(3.*np.sqrt(5.)) ) 
    z_split[2] =  np.sqrt( (np.sqrt(5.)-2.)/(3.*np.sqrt(5.)) ) 
    z_split[3] =  np.sqrt( (np.sqrt(5.)+2.)/(3.*np.sqrt(5.)) ) 
  
  if Nsplit == 3: 
    z_split[0] = -np.sqrt(2.) /2. 
    z_split[2] =  np.sqrt(2.) /2. 
  
  if Nsplit == 2: 
    z_split[0] = -1. /np.sqrt(3.) 
    z_split[1] =  1. /np.sqrt(3.) 
  
  if Nsplit == 1: 
    print "Nsplit = 1? Really?"
  
  for ik in np.arange(k_max): 
    k = ik +1 
    rk = np.power(3.*k/(4.*np.pi*n_e), 1./3.)  # expectation value of distance to the k-th nearest neighbor, assuming isotropic const density 
    z_k = rk *z_split  
    
    for i_sp in np.arange(Nsplit): 
      D0 = np.sqrt(rk*rk -z_k[i_sp]*z_k[i_sp]) 
      a1, a2, a3, a4, Fz_arr, a5, a6 = magDyn1D_leapfr02(z_k[i_sp], v0, C0, D0, T, Nstep)
      F_free_arr = unpert1D(z_k[i_sp], v0, C0, D0, T, Nstep)
      F_cumul_arr += (Fz_arr -F_free_arr) /np.float64(Nsplit)
    
  
  #F_cumul_arr /= np.float64(k_max) 
  
  return t_arr, F_cumul_arr 




def unpert1D(z_ini, v0, C0, D0, T, Nstep): 
  
  dt = T /np.float64(Nstep) 
  t_arr = np.zeros(Nstep+2, dtype=np.float64) 
  t_arr[-1] = T 
  t_arr[1:Nstep+1] = 0.5*dt +np.linspace(0., T -dt, Nstep) 
  z_free_arr = z_ini +v0 *t_arr 
  r = np.sqrt(D0*D0 +z_free_arr*z_free_arr) 
  F_free_arr = C0 *z_free_arr /(r*r*r)  # from the electron on the ion, hence the "+" sign 
  
  return F_free_arr 




def magDyn1D_leapfr02(z_ini, v0, C0, D0, T, Nstep): 
  dt = T /np.float64(Nstep) 
  t_arr = np.zeros(Nstep+2, dtype=np.float64) 
  z_arr = 1.0*t_arr  
  v_arr = np.zeros(Nstep+1, dtype=np.float64) 
  t_v_arr = np.linspace(0.0, T, Nstep+1) 
  Fz_arr = 1.0*t_arr  # from electron on the ion 
  
  z = z_ini 
  vz = v0 
  #print v0
  
  z_arr[0] = z 
  v_arr[0] = vz 
  r = np.sqrt(D0*D0 +z*z) 
  r3 = r*r*r
  Fz_arr[0] = C0 *z / r3  #  "+", because the force from electron on the ion 
  
  z += vz * dt/2. 
  
  t_arr[1] = dt/2. 
  z_arr[1] = z 
  r = np.sqrt(D0*D0 +z*z) 
  r3 = r*r*r
  Fz_arr[1] = C0 *z / r3
  
  for istep in np.arange(Nstep-1):
    vz += -dt *C0 *z / r3  # C0 > 0 
    z += vz *dt 
    r = np.sqrt(D0*D0 +z*z) 
    r3 = r*r*r 
    
    t_arr[2+istep] = (1.5 +istep) *dt 
    z_arr[2+istep] = z 
    Fz_arr[2+istep] = C0 *z / r3 
    v_arr[1+istep] = vz 
  
  vz += -dt *C0 *z / r3  # C0 > 0 
  z += vz *dt/2. 
  
  t_arr[Nstep+1] = T  
  z_arr[Nstep+1] = z 
  r = np.sqrt(D0*D0 +z*z) 
  r3 = r*r*r
  Fz_arr[Nstep+1] = C0 *z / r3 
  v_arr[Nstep] = vz 
  
  return z, vz, t_arr, z_arr, Fz_arr, t_v_arr, v_arr   



if __name__=="__main__":
  import numpy as np 
  import matplotlib.pyplot as plt 
  main() 
  
