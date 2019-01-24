def main(): 
  
  # Parkhomchuk: 
  
  Z = 92 
  Bz = 0.1  # T (working in MKS units)
  n_e = 8.323e+12  # m^-3, in the beam frame 
  Delta_e_par = 0.0  #4.2e+4  # m/s 
  Delta_e_tr = 1.88e+5  # m/s 
  T_inter = 2.26e-9  # s, in the beam frame 
  
  q_el = -1.6021766208e-19  # Coulombs 
  m_el = 9.10938356e-31  # kg 
  c0 = 299792458.  # m/s 
  r_cl_el = 2.8179403227e-15  # m, classical radius of electron 
#  Z * 253.2638458397924 m^3 / s^2 == Z*e*e /([4*pi*eps0]*m_e): 
  Zrcc = Z *r_cl_el *c0 *c0  
		
  w_L = q_el *Bz /m_el  # gyration frequency, in MKS; can be positive or negative  
  T_L = 2.*np.pi /np.abs(w_L) 
  r_L = Delta_e_tr /np.abs(w_L)  # a "typical" gyro-radius 
  w_p = 56.4146e+0*np.sqrt(n_e) # s^-1, electrom plasma frequency; n_e in m^-3 
  T_p = 2.*np.pi /w_p 
  
# expectation value of distance to the nearest neighbor, 
# assuming isotropic const density: 
  r1 = np.power(3./(4.*np.pi*n_e), 1./3.)  
  
  print ' ' 
  print "omega_L = ", w_L, "s^-1, T_L = ", T_L, "s, r_L = ", r_L, "m, \
        r_nearest (in 3D) = ", r1, 'm; w_p = ', w_p,'1/s, T_p =',T_p,' s\n'
  
  Nv =1000 
  v_i_pos = np.linspace(0., 6.e+5, Nv+1)[1:] 
  rho_min = Zrcc / (v_i_pos*v_i_pos) 
  rho_max = v_i_pos / max(w_p, 1./T_inter) 
  L_P = np.log( (rho_max +rho_min +r_L) /(rho_min +r_L) ) # velocity-dependent Coulomb log 
  
  plt.plot(1.e-5*v_i_pos, L_P, 'b-'); 
  plt.xlabel('$V_{ion,\parallel} (10^5 m/s)$', fontsize=18); 
  plt.ylabel('$L_P$', fontsize=18);  
  plt.show() 
  
  V_e_eff = Delta_e_par 
  #V_e_eff = np.sqrt(Delta_e_par*Delta_e_par +Delta_e_tr*Delta_e_tr) 
  
  F_Parkhom_cold = 4. *m_el *Zrcc *Zrcc *n_e *L_P *v_i_pos / \
                   np.power(v_i_pos*v_i_pos +V_e_eff*V_e_eff, 1.5) 
  F_Parkhom_cold /= 1.60218e-19  # eV / m = 1.60218e-19 N 
  
  V_e_eff = 1.88e+5  # m/s 
  F_Parkhom_warm = 4. *m_el *Zrcc *Zrcc *n_e *L_P *v_i_pos / \
                   np.power(v_i_pos*v_i_pos +V_e_eff*V_e_eff, 1.5) 
  F_Parkhom_warm /= 1.60218e-19  # eV / m = 1.60218e-19 N 
  
  V_e_eff = Delta_e_par # reset 
  #V_e_eff = np.sqrt(Delta_e_par*Delta_e_par +Delta_e_tr*Delta_e_tr) 
		
  plt.plot(1.e-5*v_i_pos, F_Parkhom_cold, 'r-', \
           label='Parkhomchuk ($\Delta_{e \parallel} = 0$)') 
  plt.plot(1.e-5*v_i_pos, F_Parkhom_warm, 'r--', \
           label='Parkhomchuk ($\Delta_{e \parallel} = 1.88 \cdot 10^5 m/s$)') 
		
  #plt.ticklabel_format(axis='both', style='sci', scilimits=(-3,3)) 
  plt.ticklabel_format(axis='x', style='sci', scilimits=(-3,3)) 
  plt.legend(loc='upper right') 
  plt.xlabel('$V_{ion,\parallel} (10^5 m/s)$', fontsize=18) 
  plt.ylabel('$- F_{\parallel} (eV / m)$ ', fontsize=18) 
  plt.xlim(0., 6.0)  #6.05e+5)
  #plt.ylim(0., 1.1 *np.max(F_Parkh_cold)) 
  
  #plt.savefig('compare_3models.pdf') 
  plt.show() 
  

if __name__=="__main__": 
  import numpy as np 
  import matplotlib.pyplot as plt 
  #from scipy import interpolate as interp 
  #vec_piece_fit_Z79_01 = np.vectorize(piece_fit_Z79_01) 
  main() 





