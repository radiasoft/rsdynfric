#============================================================================
#
# Firstly (12/27/19), this is my copy of Ponomarev's script Fpar_1D_ave_-2.py
#
#============================================================================

#
# This differs from the Fpar_1D_01(a).py in that the average of the force on the ion is computed over 
# several nearest-neighbour particles, matching for each physical particle either the lowest 3 or
# lowest 5 moments of  the (spherical) distribution's projection onto the z axis.  Note that this 
# code is not optimized for performance: A run on a laptop with 27k nearest neighbors takes ~ 20 mins,
# with 64k nearest neighbors ~ 50 mins (Nstep = 1000). For small v, results converge for k_max ~ 1k.  
#

def main(): 
  
  eVtoJ=1.6021766208e-19      # 1 eV = 1.6...e-12 erg

  q_el = -1.6021766208e-19    # Coulombs 
  m_el = 9.10938356e-31       # kg 
  m_p = 1.672622e-27          # kg
  c0 = 299792458.             # m/s 
  r_cl_el = 2.8179403227e-15  # m, classical radius of electron 
  Z_p = 1.                    # Z of proton
  Z_Au = 79.                  # Z of gold
  Z_U = 92.                   # Z of uranium
  eps0 = 8.854188e-12         # permitivity of the space, F/m
		
  Nstep = 1000  # number of timesteps 

#-------------------------------------------------------
#
#                Input data for EIC cooler:
#
  Tsim = 4.0e-10  # s 
  
  Z = Z_Au 
		
  Bz = 5.0                    # T (working in MKS units)
  n_e = 2.0e+15               # m^-3, in the beam frame 
  Delta_e_par = 1.0e+5        # m/s 
  Delta_e_tr = 4.2e+5         # m/s 
  
#-------------------------------------------------------
#
#                Input data for HESR cooler:
#
#  Z * 253.2638458397924 m^3 / s^2 == Z*e*e /([4*pi*eps0]*m_e): 
  Zrcc = Z *r_cl_el *c0 *c0  
  L_cooling = 2.7                    # m

  Ekin_p = [200.,353.,580.,1670.]    # MeV
  m_p_MeV = m_p*c0**2/(1.e6*eVtoJ)
  Ekin_el = [109.,192.,316.,908.]    # keV
  m_el_keV = m_el*c0**2/(1.e3*eVtoJ)
  for i in range(np.size(Ekin_p)):
     v_p = c0*np.sqrt(Ekin_p[i]/m_p_MeV/(Ekin_p[i]/m_p_MeV+1.))
     v_el = c0*np.sqrt(Ekin_el[i]/m_el_keV/(Ekin_el[i]/m_el_keV+1.))
     T_cool = L_cooling/v_p
     print 'v_p = ', v_p, 'm/s, v_el = ', v_el, 'm/s; T_cool = ',T_cool,'s'


  Tsim = T_cool               # s;   for Ekin_p  = 1670 MeV 
  T_inter = 2.26e-9  # s, in the beam frame 
  T_inter = T_cool  # s, in the beam frame 
  v0_beam = v_el              # m/s; for Ekin_el =  908 keV 
#  Z = Z_p 
  Z = Z_U 


  Bz = .1                     # T (working in MKS units)
  I_beam = 0.5                # A
  a_beam = .01                # m
  S_beam = np.pi*a_beam**2    # m^2
  n_e = I_beam/S_beam/(abs(q_el)*v0_beam)   # m^-3, in the beam frame
  print 'v0_beam = ',v0_beam,'m/s, n_e = ', n_e,' 1/m^3' 
  Delta_e_par = 6.0e+4        # m/s 
  Delta_e_tr = 2.7e+5         # m/s 
  
#
#                     End of input data
#
#-------------------------------------------------------
  
  eT_par_eV = .5*m_el*np.power(Delta_e_par,2)/eVtoJ
  eT_tr_eV  = .5*m_el*np.power(Delta_e_tr,2)/eVtoJ
  print 'Delta_e_par = ',Delta_e_par,' m/s  ==> eT_par_eV=',eT_par_eV,' eV'
  print 'Delta_e_tr = ',Delta_e_tr,' m/s  ==> eT_tr_eV=',eT_tr_eV,' eV'

# k-th nearest neighbor, if contribs from only one particle are computed (Not done here;
# better to introduce k_min?):
  #k = 1   

# average the force on the ion from up to (and including) the k_max-th nearest neighbour 
  k_max = 1000

# split each nearest neighbor into Nsplit "sub-particles", matching Nsplit+1 lowest 
# z-moments of the isotropic distribution if Nsplit is even, Nsplit moments if Nsplit
# is odd; Nsplit must be <= 5,
# for now: 
  Nsplit = 5  
   
  w_L = np.abs(q_el*Bz) /m_el # Larmor frequency in MKS, 1/s 
  T_L = 2.*np.pi /w_L         # Larmor period, s
  r_L = Delta_e_tr /w_L       # Larmor radius, m 
  
  w_p = 56.4146*np.sqrt(n_e)  # electrom plasma frequency; n_e in m^-3, 1/s 
  T_p = 2.*np.pi /w_p         # period of electrom plasma oscillation, s
 
  print ('w_L (1/s) = %e, T_L (s) =%e, r_L (m) =%e; w_p (1/s) = %e, T_p (s) = %e' % \
         (w_L,T_L,r_L,w_p,T_p))
 
#  sys.exit()

#  v_i_rel = [.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9, \
#             2.]
  v_i_rel = [.01,.02,.03,.04,.05,.06,.07,.08,.09,.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6, \
             .65,.7,.75,.8,.85,.9,.95,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2., 2.2, \
	     2.4,2.6,2.8,3.,3.5,4.,4.5,5.,5.5,6.]
  v_i_numb = np.size(v_i_rel)
  print ('v_i_numb = %d' % v_i_numb)
#  v_i_numb = 5                                                   # For debugging
#  Delta_e_par_arr = [.5e+5,1.0e+5,2.0e5]      # m/s 
#  Delta_e_par_arr = [1.0e+5]      # m/s 

  Delta_e_par_arr = [.6e+5]      # m/s 

  v_e_numb = np.size(Delta_e_par_arr)
  print ('v_e_numb = %d' % v_e_numb)

  v_i_abs = np.zeros((v_i_numb,v_e_numb)) 
# For longitudinal friction forse from Parchomchuk formulae:
  v_i_abs = np.zeros((v_i_numb,v_e_numb)) 
  log_par = np.zeros((v_i_numb,v_e_numb))         # Coulomb log  
  F_P_par_w = np.zeros((v_i_numb,v_e_numb))       # without Coulomb log       
  F_P_par_t = np.zeros((v_i_numb,v_e_numb))       # wit Coulomb log            
  for k in range(v_e_numb):
     v_eff = Delta_e_par_arr[k]
     eT_eff_eV  = .5*m_el*np.power(v_eff,2)/eVtoJ
     print ('v_eff = %e m/s, T_eff = %e' % (v_eff,eT_eff_eV))
     for i in range(v_i_numb):
        v_i_abs[i,k] = v_i_rel[i] *Delta_e_par_arr[k]  # m/s
        rho_min = Zrcc/v_i_abs[i,k]**2
        rho_max = v_i_abs[i,k]/max(w_p,1./T_inter) 
        log_par[i,k] = np.log((rho_max+rho_min+r_L)/(rho_min+r_L) ) 
        F_P_par_w[i,k] = -4.*m_el*Zrcc*Zrcc*n_e*v_i_abs[i,k]/ \
                       np.power(v_i_abs[i,k]*v_i_abs[i,k] +v_eff*v_eff, 1.5)/eVtoJ 
        F_P_par_t[i,k] = F_P_par_w[i,k]* log_par[i,k]
	
  powDelta_e_par_p = np.zeros(v_e_numb)
  mantDelta_e_par_p = np.zeros(v_e_numb)
  for k in range(v_e_numb):
     powDelta_e_par_p[k] = round(np.log10(Delta_e_par_arr[k])) 
     mantDelta_e_par_p[k] = Delta_e_par_arr[k]/np.power(10,powDelta_e_par_p[k]) 

  fig1 = plt.figure(1)
  plt.plot(v_i_rel[0:v_i_numb], log_par[0:v_i_numb,0],'-xr')
  plt.xlabel('$V_{i||}/v_{e||}}$',color='m',fontsize=14) 
  plt.ylabel('$\Lambda$',color='m',fontsize=14) 
  titleHeader = ('Coulomb Logarithm $\Lambda$: "Parkhomchuk"')
  plt.title(titleHeader,color='m',fontsize=14)
  plt.grid(True)
  fig1.savefig('HESR_CoulombLog_Parkhomchuk_fig1.png')    
  print ('File "HESR_CoulombLog_Parkhomchuk_fig1.png" is written')   

  fig5 = plt.figure(5)
  plt.plot(v_i_rel[0:v_i_numb], -F_P_par_t[0:v_i_numb,0],'-xr', \
           v_i_rel[0:v_i_numb], -F_P_par_w[0:v_i_numb,0],'-xb')
  plt.xlabel('$V_{i||}/v_{e||}}$',color='m',fontsize=14) 
  plt.ylabel('$-F_{\parallel} (eV / m)$',color='m',fontsize=14) 
  titleHeader = ('HESR Longitudinal Friction Force: "Parkhomchuk"')
  plt.title(titleHeader,color='m',fontsize=14)
  plt.text(2.,75.,('"Cold" Beam: $V_{e||}=%3.1f\cdot10^{%2d}$ m/s' % \
                       (mantDelta_e_par_p[0],powDelta_e_par_p[0])), \
           color='m',fontsize=14)
  plt.legend([('With Coulomb Logaruthm'),('Without Coulomb Logaruthm')], \
             loc='lower center',fontsize=14)
  plt.grid(True)
  fig5.savefig('HESR_coldBeam_Parkhomchuk_fig5.png')    
  print ('File "HESR_coldBeam_Parkhomchuk_fig5.png" is written')   

#  plt.show()
  
#  sys.exit()




#  Nsplit_arr = [9,7,6,5]
#  Nsplit_arr = [5,4,3,2]
  Nsplit_arr = [5]
  numbNsplit = np.size(Nsplit_arr)
  ffLong = np.zeros((v_i_numb,numbNsplit,v_e_numb)) 
  ffLong_p = np.zeros((v_i_numb,numbNsplit,v_e_numb)) 
  z_traject = np.zeros((Nstep+3,2,v_i_numb)) 
  vz_traject = np.zeros((Nstep+3,2,v_i_numb)) 

  indx_i = [0,8,15,17,19,21,22,23,24,25,30,35,40]
  for j in range(numbNsplit):
     Nsplit = int(Nsplit_arr[j])
     z_split = spliting(Nsplit)
     for k in range(v_e_numb):
        for i in range(v_i_numb):
#        for ii in range(np.size(indx_i)):
#	   i = indx_i[ii]
           timeStart=os.times()
# The fact that Delta_e_par is the rms thermal velocity is not used here; it really
# is only 

#           t_arr, F_ave_arr, F_free_arr, zFirst, zLast, vzFirst, vzLast  = \
#	   aveOverNearestNeibs01(n_e, -1.*v_i_abs_p[i,k],Zrcc,Tsim,Nstep,k_max,Nsplit,z_split) 
           t_arr, F_ave_arr, F_free_arr, zFirst, zLast, vzFirst, vzLast  = \
	   aveOverNearestNeibs01(n_e,-1.*v_i_abs[i,k],Zrcc,Tsim,Nstep,k_max,Nsplit,z_split) 
  
           F_ave_arr *= m_el /1.60218e-19  # eV / m = 1.60218e-19 N 
# time-integrated acceleration of electron(s):
           F_time_aved = (Tsim /np.float64(Nstep)) *np.sum( F_ave_arr[1:-1] ) 

# average force on the ion during the interaction time (redundant *Tsim/Tsim for clarity):
           F_time_aved /= Tsim   
#           ffLong[i,j,k] = F_time_aved
           ffLong_p[i,j,k] = F_time_aved
# Trajectories data:
	   z_length = np.size(zFirst)
	   z_traject[0:z_length,0,i] = zFirst
	   z_traject[0:z_length,1,i] = zLast
	   vz_length = np.size(vzFirst)
	   vz_traject[0:vz_length,0,i] = vzFirst
	   vz_traject[0:vz_length,1,i] = vzLast
     
           timeEnd=os.times()
           cpuTime=float(timeEnd[0])-float(timeStart[0])  # CPU time , s
  
           print "   cpuTime (s) = ", cpuTime, " for v_i_rel = ", v_i_rel[i], \
	         ", Nsplit = ",Nsplit," and Delta_e_par_arr = ",Delta_e_par_arr[k], \
		 "ffLong_p = ",ffLong_p[i,j,k]," eV/m"
  
     #  plt.plot(t_arr, F_ave_arr, 'r-') 
     #  plt.plot(t_arr, 0.0*t_arr, 'g--') 
     #  plt.xlabel('t (s)') 
     #  plt.ylabel('$F_{\parallel} (eV / m)$') 
     #  plt.xlim(0., Tsim) 
     #  plt.show()

  powDelta_e_par = np.zeros(v_e_numb)
  mantDelta_e_par = np.zeros(v_e_numb)
    
  for k in range(v_e_numb):
     powDelta_e_par[k] = round(np.log10(Delta_e_par_arr[k])) 
     mantDelta_e_par[k] = Delta_e_par_arr[k]/(10**powDelta_e_par[k]) 

#
# Opening the output file (for 4 values of Nsplit): 
#
###  results_file='frictionForceLong_Nsplit_5-6-7-9.dat'
###  results_file='frictionForceLong_Nsplit_2-3-4-5.dat'
##  results_file='frictionForceLong_Nsplit_5-6-7-9_noBckgrnd.dat'
###  results_file='frictionForceLong_Nsplit_2-3-4-5_noBckgrnd.dat'
##  print ('Open output file "%s"...' % results_file)
##  results_file_flag=0
##  try:
##     outfile = open(results_file,'w')
##     results_file_flag=1
##  except:
##     print ('Problem to open output file "%s"' % results_file)
#
# Writing the results to output file: 
#
##  outfile.write ('\n                       Initial data:\n\n')
##  outfile.write ('Velocities of Electrons: rmsTran = %5.3e m/s, rmsLong = %5.3e m/s' % \
##                 (Delta_e_tr,Delta_e_par))
##  outfile.write ('\nTime of Cooling = %5.3e s During %d Steps' % (Tsim,Nstep))
##  outfile.write ('\nBfield=%4.2f T, Zion = %d, Plasma Density = %5.3e 1/m^3' % (Bz,Z,n_e))
##  outfile.write ('\nLarmor Frequency = %e 1/s, Plasma Frequency = %e 1/s' % (w_L,w_p))
##  outfile.write ('\nLarmor Period = %e s, Plasma Period = %e s;' % (T_L,T_p))
##  outfile.write ('\nTypical Larmor Radius = %e m' % (r_L))  
##  outfile.write ('\n\n    Nsplit        %d              %d               %d              %d' % \
##                 (int(Nsplit_arr[0]),int(Nsplit_arr[1]),int(Nsplit_arr[2]),int(Nsplit_arr[3])))
##  outfile.write ('\nVion/VeLong               F r i c t i o n   F o r c e   (eV/m)\n')

##  for i in range(v_i_numb):
##     outfile.write ('\n %5.3f       %5.3e     %5.3e     %5.3e     %5.3e'  % \
##                    (v_i_rel[i],ffLong[i,0],ffLong[i,1],ffLong[i,2],ffLong[i,3]))
##  outfile.close()
##  print ('File "%s" is written' % results_file)
#----------------------------------------------------
#
# Opening the output file (for one valuee each of Nsplit and v_e_numb): 
#
#  results_file='frictionForceLong_Nsplit_5_ve_3.dat'
#  print ('Open output file "%s"...' % results_file)
#  results_file_flag=0
#  try:
#     outfile = open(results_file,'w')
#     results_file_flag=1
#  except:
#     print ('Problem to open output file "%s"' % results_file)
#
# Writing the results to output file: 
#
#  outfile.write ('\n                       Initial data:\n\n')
#  outfile.write ('Velocity of Electrons: rmsTran = %5.3e m/s' % Delta_e_tr)
#  outfile.write ('\nTime of Cooling = %5.3e s During %d Steps' % (Tsim,Nstep))
#  outfile.write ('\nBfield=%4.2f T, Zion = %d, Plasma Density = %5.3e 1/m^3' % (Bz,Z,n_e))
#  outfile.write ('\nLarmor Frequency = %e 1/s, Plasma Frequency = %e 1/s' % (w_L,w_p))
#  outfile.write ('\nLarmor Period = %e s, Plasma Period = %e s;' % (T_L,T_p))
#  outfile.write ('\nTypical Larmor Radius = %e m; Nsplit = %d' % (r_L,int(Nsplit_arr[0])))  
#  outfile.write ('\n                  %e      %e       %e' % \
#                 (Delta_e_par_arr[0],Delta_e_par_arr[1],Delta_e_par_arr[2]))  
#  
#  outfile.write ('\n\nVion/VeLong            F r i c t i o n    F o r c e     (eV/m)\n')
#
#  for i in range(v_i_numb):
#     outfile.write ('\n    %5.3f         %5.3e         %5.3e         %5.3e'  % \
#                    (v_i_rel[i],ffLong[i,0,0],ffLong[i,0,1],ffLong[i,0,2]))
#  outfile.close()
#  print ('File "%s" is written' % results_file)
#
#----------------------------------------------------

#----------------------------------------------------
#
# Opening the output file (for one values of Nsplit and v_e_numb value of 
# electron velocity): 
#

  results_file='HESR_coldBeam.dat'
  print ('Open output file "%s"...' % results_file)
  results_file_flag=0
  try:
     outfile = open(results_file,'w')
     results_file_flag=1
  except:
     print ('Problem to open output file "%s"' % results_file)
#
# Writing the results to output file: 
#
  outfile.write ('\n                       Initial data:\n\n')
  outfile.write ('Velocity of Electrons: rmsTran = %5.3e m/s (T = %5.3f eV)' % \
                 (Delta_e_tr,eT_tr_eV))
  outfile.write ('\nVelocity of Electrons: rmsLong = %5.3e m/s (T = %5.3f eV)' % \
                 (Delta_e_par,eT_par_eV))
  outfile.write ('\nTime of Cooling = %5.3e s During %d Steps' % (Tsim,Nstep))
  outfile.write ('\nBfield=%4.2f T, Zion = %d, Plasma Density = %5.3e 1/m^3' % (Bz,Z,n_e))
  outfile.write ('\nLarmor Frequency = %e 1/s, Plasma Frequency = %e 1/s' % (w_L,w_p))
  outfile.write ('\nLarmor Period = %e s, Plasma Period = %e s;' % (T_L,T_p))
  outfile.write ('\nTypical Larmor Radius = %e m; Nsplit = %d' % (r_L,int(Nsplit_arr[0])))  
  
  outfile.write ('\n\nVion/VeLong   Friction Force, eV/m       "Parkhomchuk", eV/m\n')

  for i in range(v_i_numb):
     outfile.write ('\n    %5.3f         %5.3e                 %5.3e'  % \
                   (v_i_rel[i],ffLong_p[i,0,0],F_P_par_t[i,0]))
  outfile.close()
  print ('File "%s" is written' % results_file)
  
#
#----------------------------------------------------


#
# Plotting of results: 
#
#  fig20 = plt.figure(20)
#  plt.plot(v_i_rel[0:v_i_numb], -ffLong[0:v_i_numb,0],'-b',v_i_rel[0:v_i_numb], -ffLong[0:v_i_numb,1],'-m', \
#           v_i_rel[0:v_i_numb], -ffLong[0:v_i_numb,2],'-g',v_i_rel[0:v_i_numb], -ffLong[0:v_i_numb,3],'-k', \
#           v_i_rel[0:v_i_numb], -ffLong[0:v_i_numb,0],'xr',v_i_rel[0:v_i_numb], -ffLong[0:v_i_numb,1],'xr', \
# 	   v_i_rel[0:v_i_numb], -ffLong[0:v_i_numb,1],'xr',v_i_rel[0:v_i_numb], -ffLong[0:v_i_numb,3],'xr') 
#  plt.xlabel('$V_{i||}/v_{e||}}$',color='m',fontsize=14) 
#  plt.ylabel('$-F_{\parallel} (eV / m)$',color='m',fontsize=14) 
#  titleHeader = ('Longitudinal Friction Force: $V_{e||}=%3.1f\cdot10^{%2d}$ m/s' % \
#                 (mantDelta_e_par,powDelta_e_par))
#  plt.title(titleHeader,color='m',fontsize=14)
#  plt.legend([('$N_{split}$ = %d' % Nsplit_arr[0]),('$N_{split}$ = %d' % Nsplit_arr[1]), \
#              ('$N_{split}$ = %d' % Nsplit_arr[2]),('$N_{split}$ = %d' % Nsplit_arr[3])], \
#    	      loc='upper right',fontsize=14)
#  plt.grid(True)
#  fig20.savefig('frictionForceLong_fig20.png')    
#  print ('File "frictionForceLong_fig20.png" is written')   


#  fig30 = plt.figure(30)
#  plt.plot(v_i_rel, -ffLong[:,0],'-b',v_i_rel, -ffLong[:,1],'-m', \
#           v_i_rel, -ffLong[:,2],'-g',v_i_rel, -ffLong[:,3],'-k', \
#           v_i_rel, -ffLong[:,0],'xr',v_i_rel, -ffLong[:,1],'xr', \
#	   v_i_rel, -ffLong[:,1],'xr',v_i_rel, -ffLong[:,3],'xr') 
#  plt.xlabel('$V_{i||}/v_{e||}}$',color='m',fontsize=14) 
#  plt.ylabel('$-F_{\parallel} (eV / m)$',color='m',fontsize=14) 
#  titleHeader = ('Longitudinal Friction Force: $V_{e||}=%3.1f\cdot10^{%2d}$ m/s' % \
#                 (mantDelta_e_par,powDelta_e_par))
#  plt.title(titleHeader,color='m',fontsize=14)
#  plt.legend([('$N_{split}$ = %d' % Nsplit_arr[0]),('$N_{split}$ = %d' % Nsplit_arr[1]), \
#              ('$N_{split}$ = %d' % Nsplit_arr[2]),('$N_{split}$ = %d' % Nsplit_arr[3])], \
#   	      loc='lower right',fontsize=14)
#  plt.grid(True)
#  fig30.savefig('frictionForceLong_fig30_noBckgrnd.png')    
#  print ('File "frictionForceLong_fig30_noBckgrnd.png" is written')   

#  fig40 = plt.figure(40)
#  plt.plot(v_i_rel[0:v_i_numb], -ffLong[0:v_i_numb,0],'-b', \
#           v_i_rel[0:v_i_numb], -ffLong[0:v_i_numb,0],'xr') 
#  plt.xlabel('$V_{i||}/v_{e||}}$',color='m',fontsize=14) 
#  plt.ylabel('$-F_{\parallel} (eV / m)$',color='m',fontsize=14) 
#  titleHeader = ('Longitudinal Friction Force: $V_{e||}=%3.1f\cdot10^{%2d}$ m/s' % \
#                 (mantDelta_e_par,powDelta_e_par))
#  plt.title(titleHeader,color='m',fontsize=14)
#  plt.legend([('$N_{split}$ = %d' % Nsplit_arr[0])],loc='upper right',fontsize=14)
#  plt.grid(True)
#  fig40.savefig('frictionForceLong_fig40.png')    
#  print ('File "frictionForceLong_fig40.png" is written')   

  fig45 = plt.figure(45)
  plt.plot(v_i_rel[0:v_i_numb], -ffLong_p[0:v_i_numb,0,0],'-xr') 
  plt.xlabel('$V_{i||}/v_{e||}}$',color='m',fontsize=14) 
  plt.ylabel('$-F_{\parallel} (eV / m)$',color='m',fontsize=14) 
  titleHeader = ('HESR Longitudinal Friction Force: $V_{e||}=%3.1f\cdot10^{%2d}$ m/s' % \
                 (mantDelta_e_par,powDelta_e_par))
  plt.title(titleHeader,color='m',fontsize=14)
#  plt.legend([('$N_{split}$ = %d' % Nsplit_arr[0])],loc='upper right',fontsize=14)
  plt.grid(True)
  fig45.savefig('HESR_coldBeam_Simulation_fig45.png')    
  print ('File "HESR_coldBeam_Simulation_fig45.png" is written')   

  fig46 = plt.figure(46)
  plt.plot(v_i_rel[0:v_i_numb], -ffLong_p[0:v_i_numb,0,0],'-xr') 
  plt.xlabel('$V_{i||}/v_{e||}}$',color='m',fontsize=14) 
  plt.ylabel('$-F_{\parallel} (eV / m)$',color='m',fontsize=14) 
  titleHeader = ('HESR Longitudinal Friction Force: $V_{e||}=%3.1f\cdot10^{%2d}$ m/s' % \
                 (mantDelta_e_par,powDelta_e_par))
  plt.title(titleHeader,color='m',fontsize=14)
  plt.xlim([0.,2.])
#  plt.legend([('$N_{split}$ = %d' % Nsplit_arr[0])],loc='upper right',fontsize=14)
  plt.grid(True)
  fig46.savefig('HESR_coldBeam_Simulation_fig46.png')    
  print ('File "HESR_coldBeam_Simulation_fig46.png" is written')   

#  fig50 = plt.figure(50)
#  plt.plot(v_i_rel[0:v_i_numb], -ffLong[0:v_i_numb,0,0],'-b', \
#           v_i_rel[0:v_i_numb], -ffLong[0:v_i_numb,0,1],'-m', \
#           v_i_rel[0:v_i_numb], -ffLong[0:v_i_numb,0,2],'-g', \
#           v_i_rel[0:v_i_numb], -ffLong[0:v_i_numb,0,0],'xr', \
#           v_i_rel[0:v_i_numb], -ffLong[0:v_i_numb,0,1],'xr', \
#           v_i_rel[0:v_i_numb], -ffLong[0:v_i_numb,0,2],'xr') 
#  plt.xlabel('$V_{i||}/v_{e||}}$',color='m',fontsize=14) 
#  plt.ylabel('$-F_{\parallel} (eV / m)$',color='m',fontsize=14) 
#  titleHeader = ('Longitudinal Friction Force: $N_{split}$ = %d' % (int(Nsplit_arr[0])))
#  plt.title(titleHeader,color='m',fontsize=14)
#  plt.legend([(('$v_{e||}=%3.1f \cdot 10^{%2d}$ m/s' % (mantDelta_e_par[0],powDelta_e_par[0]))), \
#              (('$v_{e||}=%3.1f \cdot 10^{%2d}$ m/s' % (mantDelta_e_par[1],powDelta_e_par[1]))), \
#              (('$v_{e||}=%3.1f \cdot 10^{%2d}$ m/s' % (mantDelta_e_par[2],powDelta_e_par[2])))], \
#    	     loc='upper right',fontsize=14)
#  plt.grid(True)
#  fig50.savefig('frictionForceLong_fig50.png')    
#  print ('File "frictionForceLong_fig50.png" is written') 

#------------------- Potting of thrajectories ------------
#
###  arg_x = np.linspace(start=0,stop=z_length,num=z_length,dtype=int)
  
###  typeLines = ['-r','-b','-m','-k','-g','.r','.b','.m','.k','.g', \
###               '-r','-b','-m','-k','-g','.r','.b','.m','.k','.g']
####  for i in range(0,v_i_numb,5):
###  fig60 = plt.figure(60)
###  plt.xlabel('Time Step',color='m',fontsize=14)
###  plt.ylabel('z, $\mu$m',color='m',fontsize=14)
###  titleHeader=('First Trajectories: $V_i/V_0$=%5.3f' % v_i_rel[i])
###  plt.title(titleHeader,color='m',fontsize=14)
###  plt.grid(True)
###  for ii in range(np.size(indx_i)):
###     i = indx_i[ii]
###     plt.plot(arg_x,1.e6*z_traject[0:z_length,0,i],typeLines[ii])
     

####  for i in range(0,v_i_numb,5):
###  fig160 = plt.figure(160)
###  plt.xlabel('Time Step',color='m',fontsize=14)
###  plt.ylabel('z, $\mu$m',color='m',fontsize=14)
###  titleHeader=('Last Trajectories: $V_i/V_0$=%5.3f' % v_i_rel[i])
###  plt.title(titleHeader,color='m',fontsize=14)
###  plt.grid(True)
###  for ii in range(np.size(indx_i)):
###     i = indx_i[ii]
###     plt.plot(arg_x,1.e6*z_traject[0:z_length,1,i],typeLines[ii])
#
#-------------------------------------------------------
  
  plt.show()

def spliting(Nsplit):
#  print ('Start of spliting with Nsplit = %d' % (Nsplit))
  z_split = np.zeros(Nsplit, dtype=np.float64) 
   
  if Nsplit == 9: 
    z_split[0] = -np.sqrt(0.831224489796)  # Ilya; = -.9117151
    z_split[1] = -np.sqrt(0.355918367347)  # Ilya; = -.5965889
    z_split[2] = -np.sqrt(0.285510204082)  # Ilya; = -.5343315 
    z_split[3] = -np.sqrt(0.027346938776)  # Ilya; = -.1653691
    z_split[5] =  np.sqrt(0.027346938776)  # Ilya; =  .1653691
    z_split[6] =  np.sqrt(0.285510204082)  # Ilya; = . 5343315
    z_split[7] =  np.sqrt(0.355918367347)  # Ilya; =  .5965889  
    z_split[8] =  np.sqrt(0.831224489796)  # Ilya; =  .9117151 
    
  if ((Nsplit == 7) or (Nsplit == 6)):  

     if (Nsplit == 6):
        a1 = 1.
        a2 = .6
        a3 = 3./7.

     if (Nsplit == 7):
        a1 = 7./6.
        a2 = .7
        a3 = .5

# Cubic equations for  z^2 with Nsplit = 6,7:     
     a = 3.
     b = -3.*a1
     c = 1.5*(a1**2-a2)
     d = -(.5*a1*(a1*a1-3.*a2)+a3)
#     print ('For Nsplit = %d: a = %e, b = %e, c = %e, d =%e' % (Nsplit,a,b,c,d))
     r = b/a
     s = c/a
     t = d/a
     p = s-r*r/3.
     q = 2.*r*r*r/27.-r*s/3.+t
     D = (p/3.)**3 +(q/2.)**2
  
#     print ('For Nsplit = %d: r = %e,s = %e,t = %e;\np = %e, q = %e, D = %e' % (Nsplit,r,s,t,p,q,D))

     roots=np.zeros(3)
     if (D > 0.):
        u = np.power(-.5*q+np.sqrt(D),1./3.)
        v = np.power(-.5*q-np.sqrt(D),1./3.)
        roots[2] = u+v
        print('For D > 0: u = %e, v = %e, roots[2] = %e' % (u,v,roots[2]))
   
        roots[1] = (a1-roots[2]+np.sqrt(-3.*roots[2]**2+2.*a1*roots[2]+a1**2+2.*a2))/2.
        roots[0] = a1-roots[1] -roots[2]
        print('For Nsplit = %d: roots[2] = %e, roots[1] = %e, roots[0] = %e' % \
	      (Nsplit,roots[2],roots[1],roots[0]))

# For testing: should be root2_2 = roots[0] and root1_2 = roots[1]
        root2_2 = (a1-roots[2]-np.sqrt(-3.*roots[2]**2+2.*a1*roots[2]-a1**2+2.*s2))/2.
        root1_2 = a1-root2_2 -roots[2]
#        print('Testing for Nsplit = %d: : roots[2] = %e, root2_2 = %e, root1_2 = %e' % \
#	      (Nsplit,roots[2],root2_2,root1_2))

        roots.sort()
#        print('For Nsplit = %d (final): root1 = %e, root2 = %e, root3 = %e' % \
#	      (Nsplit,roots[0],roots[1],roots[2]))
	
     if (D < 0.):
        rho = np.sqrt(-p**3/27.)
        phi = math.acos(-q/(2.*rho))

#        print('For D < 0: rho = %e, phi = %e' % (rho,phi))
  
        roots[0] = 2.*np.power(rho,1./3.)*np.cos(phi/3.)-r/3.
        roots[1] = 2.*np.power(rho,1./3.)*np.cos((phi+2.*np.pi)/3.)-r/3.
        roots[2] = 2.*np.power(rho,1./3.)*np.cos((phi+4.*np.pi)/3.)-r/3.

        roots.sort()
#        print('For Nsplit = %d: roots[0] = %e, roots[1] = %e, roots[2] = %e' % \
#	      (Nsplit,roots[0],roots[1],roots[2]))

# For testing: should be root2_2 = roots[0] and root1_2 = roots[1]
        root2_2 = (a1-roots[2]-np.sqrt(-3.*roots[2]**2+2.*a1*roots[2]-a1**2+2.*a2))/2.
        root1_2 = a1-root2_2 -roots[2]
#        print('Testing for Nsplit = %d: : roots[2] = %e, root2_2 = %e, root1_2 = %e' % \
#	      (Nsplit,roots[2],root2_2,root1_2))

  if Nsplit == 7: 
#    z_split[0] = - .8838 (Ilya)
#    z_split[1] = - .5301 (Ilya)
#    z_split[2] = - .3232 (Ilya) 
#    z_split[4] = .3232   (Ilya)
#    z_split[5] = .5301   (Ilya) 
#    z_split[6] = .8838   (Ilya)
    z_split[0] = -np.sqrt(roots[2]) # -.8858617
    z_split[1] = -np.sqrt(roots[1]) # -.5296568
    z_split[2] = -np.sqrt(roots[0]) # -.3239118
    z_split[4] =  np.sqrt(roots[0]) #  .3239118
    z_split[5] =  np.sqrt(roots[1]) #  .5296568  
    z_split[6] =  np.sqrt(roots[2]) #  .8858617 
    
  if Nsplit == 6:
     z_split[0] = -np.sqrt(roots[2]) # -.8662468
     z_split[1] = -np.sqrt(roots[1]) # -.4225187
     z_split[2] = -np.sqrt(roots[0]) # -.2666354 
     z_split[3] =  np.sqrt(roots[0]) #  .2666354 
     z_split[4] =  np.sqrt(roots[1]) #  .4225187  
     z_split[5] =  np.sqrt(roots[2]) #  .8662468 
    
  if Nsplit == 5: 
    z_split[0] = -np.sqrt( (5. +np.sqrt(11.)) /12.) # -.83249745
    z_split[1] = -np.sqrt( (5. -np.sqrt(11.)) /12.) # -.37454141 
    z_split[3] =  np.sqrt( (5. -np.sqrt(11.)) /12.) #  .37454141 
    z_split[4] =  np.sqrt( (5. +np.sqrt(11.)) /12.) #  .83249745 
    
  if Nsplit == 4: 
    z_split[0] = -np.sqrt( (np.sqrt(5.)+2.)/(3.*np.sqrt(5.))) # -.7946544
    z_split[1] = -np.sqrt( (np.sqrt(5.)-2.)/(3.*np.sqrt(5.))) # -.1875924
    z_split[2] =  np.sqrt( (np.sqrt(5.)-2.)/(3.*np.sqrt(5.))) #  .1875924
    z_split[3] =  np.sqrt( (np.sqrt(5.)+2.)/(3.*np.sqrt(5.))) #  .7946544
    # -
  if Nsplit == 3: 
    z_split[0] = -np.sqrt(2.) /2.# -.7071067 
    z_split[2] =  np.sqrt(2.) /2.#  .7071067
  
  if Nsplit == 2: 
    z_split[0] = -1. /np.sqrt(3.) # -.5773502
    z_split[1] =  1. /np.sqrt(3.) #  .5773502  
  
  if Nsplit == 1: 
    print "Nsplit = 1? Really?"

#  for k in range(Nsplit):  
#     print ('z_split[%d] = %e' % (k,z_split[k]))   
  return z_split  
  
# For a sphere, dS(z)/dz = 2*pi*r = const; but instead of uniformly sampling the [-r,r], we are
# choosing the z coordinates of the initial conditions so as to match the lowest-order z-moments of
# the spherical shell distribution: for Nsplit odd, z-moments up to order Nsplit are matched, for 
# Nsplit even, z-moments up to order Nsplit+1 are matched 

def aveOverNearestNeibs01(n_e, v0, C0, T, Nstep, k_max, Nsplit,z_split): 

  dt = T /np.float64(Nstep) 
  t_arr = np.zeros(Nstep+2, dtype=np.float64) 
  t_arr[-1] = T 
  t_arr[1:Nstep+1] = 0.5*dt +np.linspace(0., T -dt, Nstep) 
  F_cumul_arr = 0.0 *t_arr 
  
  for ik in np.arange(k_max): 
    k = ik +1 
    rk = np.power(3.*k/(4.*np.pi*n_e), 1./3.)  # expectation value of distance to the k-th nearest neighbor, assuming isotropic const density 
    z_k = rk *z_split  
    
    for i_sp in np.arange(Nsplit): 
      D0 = np.sqrt(rk*rk -z_k[i_sp]*z_k[i_sp]) 
      a1, a2, a3, a4, Fz_arr, a5, a6 = magDyn1D_leapfr02(z_k[i_sp], v0, C0, D0, T, Nstep)
      F_free_arr = unpert1D(z_k[i_sp], v0, C0, D0, T, Nstep)
      F_cumul_arr += (Fz_arr - F_free_arr) /np.float64(Nsplit)
  
  #F_cumul_arr /= np.float64(k_max) 
    if (k == 1): 
#       print 'First trajectory (for k=',k,')'
       zFirst = a4
       vzFirst = a6 
    if (k == k_max): 
#       print 'Last trajectory (for k=',k,')'
       zLast = a4
       vzLast = a6 
  return t_arr, F_cumul_arr, F_free_arr, zFirst, zLast, vzFirst, vzLast 

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
  import os, sys
  import numpy as np
  import math
  import matplotlib.pyplot as plt 
  from matplotlib.legend_handler import HandlerLine2D

  main() 
  
#######################################
#
##### Nsplit = 6:
####  a1 = 1.
####  a2 = .6
####  a3 = 3./7.

##### Nsplit = 7:
####  a1 = 7./6.
####  a2 = .7
####  a3 = .5

##### Cubic equations for Nsplit = 6,7:
####  a = 3.
####  b = -3.*a1
####  c = 1.5*(a1**2-a2)
####  d = -(.5*a1*(a1*a1-3.*a2)+a3)
####  print ('a = %e, b = %e, c = %e, d =%e' % (a,b,c,d))
####  r = b/a
####  s = c/a
####  t = d/a
####  p = s-r*r/3.
####  q = 2.*r*r*r/27.-r*s/3.+t
####  D = (p/3.)**3 +(q/2.)**2
  
####  print ('r = %e,s = %e,t = %e;\np = %e, q = %e, D = %e' % (r,s,t,p,q,D))

####  roots=np.zeros(3)
####  if (D > 0.):
####     u = np.power(-.5*q+np.sqrt(D),1./3.)
####     v = np.power(-.5*q-np.sqrt(D),1./3.)
####     roots[2] = u+v
####     print('For D > 0: u = %e, v = %e, roots[2] = %e' % (u,v,roots[2]))
   
####     roots[1] = (a1-roots[2]+np.sqrt(-3.*roots[2]**2+2.*a1*roots[2]+a1**2+2.*a2))/2.
####     roots[0] = a1-roots[1] -roots[2]
####     print('roots[2] = %e, roots[1] = %e, roots[0] = %e' % (roots[2],roots[1],roots[0]))

##### For testing: should be root2_2 = roots[0] and root1_2 = roots[1]
####     root2_2 = (a1-roots[2]-np.sqrt(-3.*roots[2]**2+2.*a1*roots[2]-a1**2+2.*s2))/2.
####     root1_2 = a1-root2_2 -roots[2]
####     print('Testing: roots[2] = %e, root2_2 = %e, root1_2 = %e' % (roots[2],root2_2,root1_2))

####     roots.sort()
####     print('root1 = %e, root2 = %e, root3 = %e' % (roots[0],roots[1],roots[2]))
	
####  if (D < 0.):
####     rho = np.sqrt(-p**3/27.)
####     phi = math.acos(-q/(2.*rho))

####     print('For D < 0: rho = %e, phi = %e' % (rho,phi))
  
####     roots[0] = 2.*np.power(rho,1./3.)*np.cos(phi/3.)-r/3.
####     roots[1] = 2.*np.power(rho,1./3.)*np.cos((phi+2.*np.pi)/3.)-r/3.
####     roots[2] = 2.*np.power(rho,1./3.)*np.cos((phi+4.*np.pi)/3.)-r/3.
####     roots.sort()
  
####     print('roots[0] = %e, roots[1] = %e, roots[2] = %e' % (roots[0],roots[1],roots[2]))

##### For testing: should be root2_2 = roots[0] and root1_2 = roots[1]
####     root2_2 = (a1-roots[2]-np.sqrt(-3.*roots[2]**2+2.*a1*roots[2]-a1**2+2.*a2))/2.
####     root1_2 = a1-root2_2 -roots[2]
####     print('Testing: roots[2] = %e, root2_2 = %e, root1_2 = %e' % (roots[2],root2_2,root1_2))

####  sys.exit()
#
#######################################

