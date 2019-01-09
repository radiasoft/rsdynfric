def main(): 
  
  V_sm  = np.linspace(-9.e+5, 9.e+5, 3601) # note that the step in v is the same as in el. distribution 
  F_sm = vec_piece_fit_Z79_01(V_sm) 
  plt.plot(V_sm, F_sm, 'b-') 
  plt.plot(V_sm, 0.0 *F_sm, 'k--') 
  plt.xlim(-6.e+5, 6.e+5) 
  plt.xlabel('$V_{ion,\parallel} (m/s)$', fontsize=18) 
  plt.ylabel('$- F_{\parallel} (eV / m)$ ', fontsize=18) 
  #plt.plot(V_i, Fpar, 'ro') 
  #plt.show()  
  
  V_i = np.array([0., 5.e+3, 15.e+3, 25.e+3, 35.e+3, 45.e+3, 55.e+3, 65.e+3, 75.e+3, 85.e+3, 95.e+3, 
                    105.e+3, 155.e+3, 205.e+3, 255.e+3, 305.e+3])  # m/s 
  Fpar = np.array([0.0, 822.823619243 +17.9896920134, 2312.80389435 +53.3197380211, 3399.66269332 +85.5985841022, 3964.28354122 +115.916869842, 4040.93182278 +142.790579008, 3779.38657397 +157.326074971, 3328.16981233 +168.937215016, 2881.84622455 +183.098201561, 2477.14550915 +193.697487947, 2181.01256832 +118.453878358, 1854.1191578 +123.625249815, 908.83863436 +129.464534718, 502.035082987 +120.367640755,  303.212901962 +105.948844318, 197.61917079 +91.3366767275  ]) # eV / m 
  
  plt.plot(V_i, Fpar, 'ro') 
  plt.ticklabel_format(axis='x', style='sci', scilimits=(-2,2)) 
  plt.show() 
  
  #p105 = np.log(Fpar[12]/Fpar[11]) / np.log(V_i[12]/V_i[11]) 
  #p155 = np.log(Fpar[13]/Fpar[12]) / np.log(V_i[13]/V_i[12]) 
  #p205 = np.log(Fpar[14]/Fpar[13]) / np.log(V_i[14]/V_i[13]) 
  #p255 = np.log(Fpar[15]/Fpar[14]) / np.log(V_i[15]/V_i[14]) 
  #p305 = -2.0
  
  #print p105, p155, p205, p255, p305 
  
  v_el = np.linspace(-3.e+5, 3.e+5, 1201) # note that the step in v is the same as in initial F(v) above (F_sm)
  sig_el = 1.0e+5  # m/s 
  f_el = np.exp( -v_el*v_el /(2. *sig_el *sig_el) ) 
  f_el /= np.sum(f_el) 
  #plt.plot(v_el, f_el, 'y-') 
  #plt.show() 
  
  Fwarm = np.convolve(f_el, F_sm, mode='same') 
  plt.plot(V_sm, F_sm, 'b-', label='cold electrons') 
  plt.plot(V_sm, Fwarm, 'g-', label='$\Delta_{e \parallel} = 10^5 m/s$')        #'warm electrons') 
  plt.plot(V_sm, 0.0 *Fwarm, 'k--') 
  plt.xlim(-3.e+5, 3.e+5) 
  plt.ticklabel_format(axis='x', style='sci', scilimits=(-2,2)) 
  plt.xlabel('$V_{ion,\parallel} (m/s)$', fontsize=18) 
  plt.ylabel('$- F_{\parallel} (eV / m)$ ', fontsize=18) 
  plt.legend() 
  plt.show() 



def piece_fit_Z79_01(v): 
  
  s = np.sign(v) 
  v = np.abs(v) 
  
  V_i = np.array([0., 5.e+3, 15.e+3, 25.e+3, 35.e+3, 45.e+3, 55.e+3, 65.e+3, 75.e+3, 85.e+3, 95.e+3, 
                    105.e+3, 155.e+3, 205.e+3, 255.e+3, 305.e+3])  # m/s 
  Fpar = np.array([0.0, 822.823619243 +17.9896920134, 2312.80389435 +53.3197380211, 3399.66269332 +85.5985841022, 3964.28354122 +115.916869842, 4040.93182278 +142.790579008, 3779.38657397 +157.326074971, 3328.16981233 +168.937215016, 2881.84622455 +183.098201561, 2477.14550915 +193.697487947, 2181.01256832 +118.453878358, 1854.1191578 +123.625249815, 908.83863436 +129.464534718, 502.035082987 +120.367640755,  303.212901962 +105.948844318, 197.61917079 +91.3366767275  ]) # eV / m 
  
  F_interp = interp.interp1d(V_i, Fpar, kind='cubic') 
  
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
  

if __name__=="__main__": 
  import numpy as np 
  import matplotlib.pyplot as plt 
  from scipy import interpolate as interp 
  vec_piece_fit_Z79_01 = np.vectorize(piece_fit_Z79_01) 
  main() 

