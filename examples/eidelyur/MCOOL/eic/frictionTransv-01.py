#============================================================================
#
# Firstly (01/20/19), this is a copy of the script frictionLong-01.py
# Now it in progerss
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
		
  Nstep = 1000                # number of timesteps 
#  deviceFlag = 1              # HESR
  deviceFlag = 2              # EIC
#  deviceFlag = 3              # MOSOL

  if (deviceFlag == 1):
#-------------------------------------------------------
#
#                Input data for HESR cooler:
#

     L_cooling = 2.7                    # m

     Ekin_p = [200.,353.,580.,1670.]    # MeV
     m_p_MeV = m_p*c0**2/(1.e6*eVtoJ)
     Ekin_el = [109.,192.,316.,908.]    # keV
     m_el_keV = m_el*c0**2/(1.e3*eVtoJ)
     for i in range(np.size(Ekin_p)):
        v_p = c0*np.sqrt(Ekin_p[i]/m_p_MeV*(Ekin_p[i]/m_p_MeV+2.))/(Ekin_p[i]/m_p_MeV+1.)
        v_el = c0*np.sqrt(Ekin_el[i]/m_el_keV*(Ekin_el[i]/m_el_keV+2.))/ \
	       (Ekin_el[i]/m_el_keV+1.)
        T_cool = L_cooling/v_p
        print( 'v_p = ', v_p, 'm/s, v_el = ', v_el, 'm/s; T_cool = ',T_cool,'s')


     Tsim = T_cool               # s;   for Ekin_p  = 1670 MeV 
     print( 'Tsim = ',Tsim,'s')
     T_inter = 2.26e-9  # s, in the beam frame 
     T_inter = T_cool  # s, in the beam frame 
     v0_beam = v_el              # m/s; for Ekin_el =  908 keV 
#     Z = Z_p 
     Z = Z_U 


     Bz = .1                     # T (working in MKS units)
     I_beam = 0.5                # A
     a_beam = .01                # m
     S_beam = np.pi*a_beam**2    # m^2
     n_e = I_beam/S_beam/(abs(q_el)*v0_beam)   # m^-3, in the beam frame
     print( 'v0_beam = ',v0_beam,'m/s, S_beam = ',S_beam,'m^2, n_e = ', n_e,' 1/m^3') 
     Delta_e_par = 6.0e+4        # m/s 
     Delta_e_tr = 2.7e+5         # m/s 
     Delta_e_par_arr = [Delta_e_par]      # m/s 
     nameDevice = 'HESR'
  
  if (deviceFlag == 2):
#-------------------------------------------------------
#
#                Input data for EIC cooler:
#
     Tsim = 4.0e-10  # s 
     T_inter = Tsim  # s, in the beam frame 
  
     Z = Z_Au 
		
     Bz = 5.0                    # T (working in MKS units)
     n_e = 2.0e+15               # m^-3, in the beam frame 
     print( 'n_e = ', n_e,' 1/m^3, T_iner - ',T_inter,' s') 
     Delta_e_par = 1.0e+5        # m/s 
     Delta_e_tr = 4.2e+5         # m/s 
#     Delta_e_par_arr = [.5e+5,1.0e+5,2.0e5]      # m/s 
     Delta_e_par_arr = [Delta_e_par]      # m/s 
     nameDevice = 'EIC'

#-------------------------------------------------------
#
#                Input data for MOSOL cooler:
#
  if (deviceFlag == 3):
#  n_e=omega_p**2*m_e/(4.*pi*q_e**2)       # plasma density, 3.1421e+08 cm-3

#  n_e1=8.e7                               # plasma density, cm-3
#  omega_p1=np.sqrt(4.*pi*n_e1*q_e**2/m_e) # plasma frequency, 5.0459e+08 1/s  


     L_cooling = 2.4                       # m

     Ekin_p = [.85]                        # MeV
     m_p_MeV = m_p*c0**2/(1.e6*eVtoJ)
     Ekin_el = [.47]                       # keV
     m_el_keV = m_el*c0**2/(1.e3*eVtoJ)
     for i in range(np.size(Ekin_p)):
        v_p = c0*np.sqrt(Ekin_p[i]/m_p_MeV*(Ekin_p[i]/m_p_MeV+2.))/(Ekin_p[i]/m_p_MeV+1.)
        v_el = c0*np.sqrt(Ekin_el[i]/m_el_keV*(Ekin_el[i]/m_el_keV+2.))/ \
	       (Ekin_el[i]/m_el_keV+1.)
        T_cool = L_cooling/v_p
        print( 'v_p =%e m/s, v_el = %e m/s; T_cool =%e s' % (v_p,v_el,T_cool))

     T_cool = L_cooling/v_el 
     print ('T_cool = %e s' % T_cool)
     T_inter = T_cool            # s, in the beam frame 
     Tsim = T_cool               # s;   
     v0_beam = v_el              # m/s; for Ekin_el =  470 eV
     Z = Z_p 


     Bz = .4                     # T (working in MKS units)
     I_beam = 0.01               # A
     a_beam = .001               # m
     S_beam = np.pi*a_beam**2    # m^2
     n_e = I_beam/S_beam/(abs(q_el)*v0_beam)   # m^-3, in the beam frame
     print( 'v0_beam = %e m/s, S_beam = %e m^2, n_e = %e 1/m^3' % (v0_beam,S_beam,n_e)) 
     
     longT=2.0e-4                                      # longitudinal temperature, eV 
     Delta_e_par = np.sqrt(2.*longT*eVtoJ/m_el)        # m/s
     print ('Delta_e_par = %e m/s  ==> longT = %e eV' % (Delta_e_par,longT))

     trnsvT=0.5                                        # transversal temperature, eV
     Delta_e_tr = np.sqrt(2.*trnsvT*eVtoJ/m_el)        # m/s
     print ('Delta_e_tr = %e m/s  ==> trnsvT = %e eV' % (Delta_e_tr,trnsvT))

     Delta_e_par_arr = [Delta_e_par]                   # m/s 
     nameDevice = 'MOSOL'

#                     End of input data
#
#-------------------------------------------------------


  v_e_numb = np.size(Delta_e_par_arr)
  print ('v_e_numb = %d' % v_e_numb)

  savePictureFile = nameDevice
  saveDataFile = nameDevice
  
#  Z * 253.2638458397924 m^3 / s^2 == Z*e*e /([4*pi*eps0]*m_e): 
  Zrcc = Z *r_cl_el *c0 *c0  
  eT_par_eV = .5*m_el*np.power(Delta_e_par,2)/eVtoJ
  eT_tr_eV  = .5*m_el*np.power(Delta_e_tr,2)/eVtoJ
  print( 'Delta_e_par = ',Delta_e_par,' m/s  ==> eT_par_eV=',eT_par_eV,' eV')
  print( 'Delta_e_tr = ',Delta_e_tr,' m/s  ==> eT_tr_eV=',eT_tr_eV,' eV')

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


  vIonFlag = 1    # = 1 for simulation of dependence on v_ion velocity
#  vIonFlag = 2    # = 1 for simulation of dependence on y0


# For testing:
  if (vIonFlag == 1):
#     v_i_rel = [.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9, \
#                2.]
     if (deviceFlag == 1):
# ForHESR case:
        v_i_rel = [.01,.02,.03,.04,.05,.06,.07,.08,.09,.1,.15,.2,.25,.3,.35,.4,.45,.5,\
                   .55, .6,.65,.7,.75,.8,.85,.9,.95,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8, \
		   1.9,2.,2.2,2.4,2.6,2.8,3.,3.5,4.,4.5,5.,5.5,6.]
        indx_i = [0,4,9,10,11,12,13,14,15,17,22,27,32,37]
     if (deviceFlag == 2):
# For EIC case:
        v_i_rel = [.01,.02,.03,.04,.05,.06,.07,.08,.09,.1,.125,.15,.175,.2,.225,.25, \
	           .275,.3,.325,.35,.375,.4,.425,.45,.5,.55, .6,.65,.7,.75,.8,.85, \
                   .9,.95,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.]
# Specially for transverse force:		   
        v_i_rel = [.001,.002,.003,.005,.0075,.01,.02,.03,.04,.05,.06,.07,.08,.09,.1,.125,.15,.175,.2,.225,.25, \
	           .275,.3,.325,.35,.375,.4,.425,.45,.5,.55, .6,.65,.7,.75,.8,.85, \
                   .9,.95,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.]
        indx_i = [0,4,9,11,13,15,17,21,23,24,26,28,30,32,34,39,44]
     if (deviceFlag == 3):
# For MOSOL case:
        v_i_rel = [.01,.02,.03,.04,.05,.06,.07,.08,.09,.1,.11,.12,.13,.14,.15,.16,.17, \
                   .18,.19,.20,.21,.22,.23,.24,.25,.26,.27]
        indx_i = [0,8,15,17,19,21,22,23,24,25,30,35,40]
#  indx_i = [0,1,2]	
# For dependence fricion forces on initioal transversal coordinate y0:
  if (vIonFlag == 2):
     v_i_rel = [.35]
  v_i_numb = np.size(v_i_rel)
#  v_i_numb = 5                                                   # For debugging
  print ('v_i_numb = %d' % v_i_numb)
# For dependence fricion forces on initioal transversl coordinate y0:
  y0rel = [0.,.5,.75,.9]                	     
  y0_numb = np.size(y0rel)
  print ('y0_numb = %d' % y0_numb)

#==============================================
# For debugging:
#
##  for m in range(y0_numb):  
#
# Opening the output files for fist trajectories with different initial y0: 
#
##     saveDataFile = nameDevice
##     saveDataFile += '_yTracks_y0-'+'{:d}'.format(int(100.*y0rel[m]))+'_coldBeam.dat'
##     print ('Open output file "%s"...' % saveDataFile)
##     saveDataFile_flag=0
##     try:
##        outfile = open(saveDataFile,'w')
##        saveDataFile_flag=1
##     except:
##        print ('Problem to open output file "%s"' % saveDataFile)
#
# Writing the results to output file: 
#
##     strOut = '\n         '+nameDevice+':  First y-Tracks for y0_rel = '
##     strOut +='{:f}'.format(y0rel[m])+'\n\n' 
##     outfile.write (strOut)
##     outfile.write ('Longitudinal Velocity of Electrons: rmsTran = %5.3e m/s (T = %5.3f eV)' % \
##                    (Delta_e_tr,eT_tr_eV))
##     outfile.write ('\nTransversal Velocity of Electrons: rmsLong = %5.3e m/s (T = %5.3f eV)' % \
##                    (Delta_e_par,eT_par_eV))
##     outfile.write ('\nTime of Cooling = %5.3e s During %d Steps' % (Tsim,Nstep))
##     outfile.write ('\nBfield=%4.2f T, Zion = %d, Plasma Density = %5.3e 1/m^3' % (Bz,Z,n_e))
##     outfile.write ('\nLarmor Frequency = %e 1/s, Plasma Frequency = %e 1/s' % (w_L,w_p))
##     outfile.write ('\nLarmor Period = %e s, Plasma Period = %e s;' % (T_L,T_p))
##     outfile.write ('\nTypical Larmor Radius = %e m; Nsplit = %d' % (r_L,int(Nsplit_arr[0])))  
  
##     strOut = '\n\nVion/VeTrnsv  '
##     for ii in range(np.size(indx_i)):
##        i = indx_i[ii]
##        strOut += '   '+'{:f}'.format(v_i_rel[i])  
##     outfile.write(strOut)
##     outfile.close()

#  sys.exit()
   
#
#==============================================

     
  v_i_abs = np.zeros((v_i_numb,v_e_numb)) 
# For longitudinal friction forse from Parchomchuk formulae:
  log_par = np.zeros((v_i_numb,v_e_numb))            # Coulomb log  
  F_P_trnsv_w = np.zeros((v_i_numb,v_e_numb))        # without Coulomb log       
  F_P_trnsv_t = np.zeros((v_i_numb,v_e_numb))        # with Coulomb log            
  F_P_trnsv_w1 = np.zeros((v_i_numb,v_e_numb))       # without Coulomb log       
  F_P_trnsv_t1 = np.zeros((v_i_numb,v_e_numb))       # wit Coulomb log            
  for k in range(v_e_numb):
     v_eff = Delta_e_par_arr[k]
# For checking:
#     v_eff = Delta_e_tr
     eT_eff_eV  = .5*m_el*np.power(v_eff,2)/eVtoJ
     print ('v_eff = %e m/s, T_eff = %e' % (v_eff,eT_eff_eV))
     for i in range(v_i_numb):
        v_i_abs[i,k] = v_i_rel[i] *Delta_e_par_arr[k]  # m/s
        rho_min = Zrcc/v_i_abs[i,k]**2
        rho_max = v_i_abs[i,k]/max(w_p,1./T_inter) 
        log_par[i,k] = np.log((rho_max+rho_min+r_L)/(rho_min+r_L) ) 
        F_P_trnsv_w[i,k] = -4.*m_el*Zrcc*Zrcc*n_e*v_i_abs[i,k]/ \
                       np.power(v_i_abs[i,k]**2+v_eff**2, 1.5)/eVtoJ 
        F_P_trnsv_t[i,k] = F_P_trnsv_w[i,k]* log_par[i,k]
# For checking (dependence on Zrcc**2):
#        F_P_trnsv_w1[i,k] = -4.*m_el*Zrcc*Zrcc*n_e*v_i_abs[i,k]/ \
#                          np.power(v_i_abs[i,k]**2+Delta_e_tr**2, 1.5)/eVtoJ 
#        F_P_trnsv_t1[i,k] = F_P_trnsv_w1[i,k]* log_par[i,k]
	
  powDelta_e_par_p = np.zeros(v_e_numb)
  mantDelta_e_par_p = np.zeros(v_e_numb)
  for k in range(v_e_numb):
     powDelta_e_par_p[k] = round(np.log10(Delta_e_par_arr[k])) 
     mantDelta_e_par_p[k] = Delta_e_par_arr[k]/np.power(10,powDelta_e_par_p[k]) 

  fig1 = plt.figure(1)
  plt.plot(v_i_rel[0:v_i_numb], log_par[0:v_i_numb,0],'-xr')
  plt.xlabel('$V_{i\parallel}/v_{e\parallel}}$',color='m',fontsize=14) 
  plt.ylabel('$\Lambda$',color='m',fontsize=14) 
  titleHeader = nameDevice
  titleHeader += ' Coulomb Logarithm $\Lambda$: "Parkhomchuk"'
  plt.title(titleHeader,color='m',fontsize=14)
  plt.grid(True)
  savePictureFile = nameDevice
  savePictureFile += '_CoulombLog_Parkhomchuk_fig1.png'    
#  fig1.savefig(savePictureFile)    
#  print( 'File "',savePictureFile,' is written' )  

  fig5 = plt.figure(5)
  plt.plot(v_i_rel[0:v_i_numb], -F_P_trnsv_t[0:v_i_numb,0],'-xr', \
           v_i_rel[0:v_i_numb], -F_P_trnsv_w[0:v_i_numb,0],'-xb')
  plt.xlabel('$V_{i\parallel}/v_{e\parallel}}$',color='m',fontsize=14) 
  plt.ylabel('$-F_{\parallel} (eV / m)$',color='m',fontsize=14) 
  titleHeader = nameDevice
  titleHeader += ' Transversal Friction Force: "Parkhomchuk"'
  plt.title(titleHeader,color='m',fontsize=14)
  plt.text(.3,0.,('$V_{e\parallel}=%3.1f\cdot10^{%2d}$ m/s, $T_{eff}^e=%4.2f$ eV' % \
                       (mantDelta_e_par_p[0],powDelta_e_par_p[0],eT_par_eV)), \
           color='m',fontsize=14)
  plt.legend([('With Coulomb Logaruthm'),('Without Coulomb Logaruthm')], \
             loc='upper right',fontsize=14)
  plt.grid(True)
  savePictureFile = nameDevice
  savePictureFile += '_coldBeam_Parkhomchuk_fig5.png'    
#  fig5.savefig(savePictureFile)    
#  print( 'File "',savePictureFile,' is written')   

# For checking:
#  fig6 = plt.figure(6)
#  plt.plot(v_i_rel[0:v_i_numb], -F_P_trnsv_t[0:v_i_numb,0],'-xr', \
#           v_i_rel[0:v_i_numb], -F_P_trnsv_w[0:v_i_numb,0],'-xb', \
#	   v_i_rel[0:v_i_numb], -F_P_trnsv_w1[0:v_i_numb,0],'-xg')
#  plt.xlabel('$V_{i\parallel}/v_{e\parallel}}$',color='m',fontsize=14) 
#  plt.ylabel('$-F_{\parallel} (eV / m)$',color='m',fontsize=14) 
#  titleHeader = nameDevice
#  titleHeader += ' Transversal Friction Force: "Parkhomchuk"'
#  plt.title(titleHeader,color='m',fontsize=14)
#  plt.text(2.,75.,('$V_{e\parallel}=%3.1f\cdot10^{%2d}$ m/s' % \
#                       (mantDelta_e_par_p[0],powDelta_e_par_p[0])), \
#           color='m',fontsize=14)
#  plt.legend([('With Coulomb Logaruthm'),('Without Coulomb Logaruthm'),('"Pogorelov"')], \
#             loc='lower center',fontsize=14)
#  plt.grid(True)
#  savePictureFile = nameDevice
#  savePictureFile += '_CoulombLog_Parkhomchuk_fig6.png'    
#  fig6.savefig(savePictureFile)    
#  print( 'File "',savePictureFile,' is written' )  

#  plt.show()
#  sys.exit()

#  Nsplit_arr = [9,7,6,5]
#  Nsplit_arr = [5,4,3,2]
  Nsplit_arr = [5]
  numbNsplit = np.size(Nsplit_arr)
  ffLong = np.zeros((v_i_numb,numbNsplit,v_e_numb,y0_numb)) 
  ffTrnsvX = np.zeros((v_i_numb,numbNsplit,v_e_numb,y0_numb)) 
  ffTrnsvY = np.zeros((v_i_numb,numbNsplit,v_e_numb,y0_numb)) 
  z_traject = np.zeros((Nstep+3,2,v_i_numb,y0_numb)) 
  vz_traject = np.zeros((Nstep+3,2,v_i_numb,y0_numb)) 
  y_traject = np.zeros((Nstep+3,2,v_i_numb,y0_numb)) 
  vy_traject = np.zeros((Nstep+3,2,v_i_numb,y0_numb)) 

  rho = np.zeros(k_max)
  arg_rho = np.zeros(k_max)
  for ik in np.arange(k_max): 
    k = ik +1 
# expectation value of distance to the k-th nearest neighbor, assuming 
# isotropic const density: 
    rk = np.power(3.*k/(4.*np.pi*n_e), 1./3.)  
    rho[ik] = rk     
    arg_rho[ik] = k

  pow_ne = round(np.log10(n_e)) 
  man_ne = n_e/(10**pow_ne) 

  fig2 = plt.figure(2)
  plt.plot(arg_rho,1.e6*rho,'-b')
#  
# For Nsplit=5 max(z_init)=z_split[4]=.83249745: 
# so max(D0)=.83249745*rk:  
#  
  plt.plot(arg_rho,.8325e6*rho,'-r',linewidth=2)
#  
# For Nsplit=5 max(z_init)=z_split[4]=.83249745; 
# so min(D0)=sqrt(1.-.83249745**2)*rk=.554025*rk:  
#  
  plt.plot(arg_rho,.554e6*rho,'-m')
  plt.xlabel('Number $k$ of Nearest Neighbor',color='m',fontsize=14) 
  plt.ylabel('$r_k$ and $D_0$, $\mu$m',color='m',fontsize=14) 
  titleHeader = nameDevice
  titleHeader += ' Distance $r_k$ to $k$-th Nearest Neighbor \nand Impact Parameter $D_0$:'
  titleHeader += ' $n_e=%3.1f\cdot 10^{%d}$ m$^{-3}$' % ( man_ne,pow_ne)
  plt.title(titleHeader,color='m',fontsize=14)
#  plt.text(400.,150.,'Max',color='m',fontsize=14)
#  plt.text(100.,83.,'Min',color='m',fontsize=14)
  plt.text(100.,180.,'$r_k=(3k/4\pi n_e)^{1/3}$',color='m',fontsize=14)
  plt.legend([('$r_k$'),('Max $D_0$'),('Min $D_0$')], \
             loc='lower right',fontsize=14)
  plt.grid(True)
  savePictureFile = nameDevice
  savePictureFile += 'impactParameter_fig2.png'    
#  fig2.savefig(savePictureFile)    
#  print( 'File "',savePictureFile,' is written' )  

#  fig3 = plt.figure(3)
#  plt.plot(arg_rho,1.e6*rho,'-r',arg_rho,.554e6*rho,'-r')
#  plt.xlabel('Number of the Electron Closest to the Ion',color='m',fontsize=14) 
#  plt.ylabel('$D_0$, $\mu$m',color='m',fontsize=14) 
#  titleHeader = nameDevice
#  titleHeader += ' Impact Parameter $D_0$'
#  plt.title(titleHeader,color='m',fontsize=14)
#  plt.text(400.,150.,'Max',color='m',fontsize=14)
#  plt.text(400.,83.,'Min',color='m',fontsize=14)
#  plt.grid(True)
  savePictureFile = nameDevice
  savePictureFile += 'impactParameter_fig3.png'    
#  fig3.savefig(savePictureFile)    
#  print( 'File "',savePictureFile,' is written' )  

#  plt.show()
#  sys.exit()

  cpuTimeTotal = 0.
  for m in range(y0_numb):
     for j in range(numbNsplit):
        Nsplit = int(Nsplit_arr[j])
        z_split = spliting(Nsplit)
        for k in range(v_e_numb):
           for i in range(v_i_numb):
#           for i in range(3):
              timeStart=os.times()
# The fact that Delta_e_par is the rms thermal velocity is not used here; it really
# is only 

#              t_arr, F_ave_arr, F_free_arr, zFirst, zLast, vzFirst, vzLast  = \
#	      aveOverNearestNeibs01(n_e, -1.*v_i_abs_p[i,k],Zrcc,Tsim,Nstep,k_max,Nsplit,z_split) 

#
# Return from aveOverNearestNeibs01:
# t_arr, F_cumul_arr_Long,  F_free_arr_Long, \
#        F_cumul_arr_Trnsv, F_free_arr_Trnsv, \
#        zFirst, zLast, vzFirst, vzLast, yFirst, yLast, vyFirst, vyLast: 
#
              t_arr, F_ave_arr_Long, F_free_arr_Long, F_ave_arr_TrnsvY, F_free_arr_TrnsvY, \
	      zFirst, zLast, vzFirst, vzLast, yFirst, yLast, vyFirst, vyLast = \
	      aveOverNearestNeibs01(n_e,-1.*v_i_abs[i,k],Zrcc,Tsim,Nstep,k_max, \
	                            Nsplit,z_split,y0rel[m]) 
# 
# Longitudinal friction force:  
# 
              F_ave_arr_Long *= m_el /1.60218e-19  # eV / m = 1.60218e-19 N 
# time-integrated acceleration of electron(s):
              F_time_aved_Long = (Tsim /np.float64(Nstep)) *np.sum( F_ave_arr_Long[1:-1] ) 

# average force on the ion during the interaction time (redundant *Tsim/Tsim for clarity):
              F_time_aved_Long /= Tsim   
              ffLong[i,j,k,m] = F_time_aved_Long
# 
# Y-transversalal friction force:  
# 
              F_ave_arr_TrnsvY *= m_el /1.60218e-19  # eV / m = 1.60218e-19 N 
# time-integrated acceleration of electron(s):
              F_time_aved_TrnsvY = (Tsim /np.float64(Nstep)) *np.sum( F_ave_arr_TrnsvY[1:-1] ) 

# average force on the ion during the interaction time (redundant *Tsim/Tsim for clarity):
              F_time_aved_TrnsvY /= Tsim   
              ffTrnsvY[i,j,k,m] = F_time_aved_TrnsvY
# Trajectories data:
              track_length = np.size(zFirst)
              z_traject[0:track_length,0,i,m] = zFirst
              z_traject[0:track_length,1,i,m] = zLast
              y_traject[0:track_length,0,i,m] = yFirst
              y_traject[0:track_length,1,i,m] = yLast
              vtrack_length = np.size(vzFirst)
              vz_traject[0:vtrack_length,0,i,m] = vzFirst
              vz_traject[0:vtrack_length,1,i,m] = vzLast
              vy_traject[0:vtrack_length,0,i,m] = vyFirst
              vy_traject[0:vtrack_length,1,i,m] = vyLast
     
              timeEnd=os.times()
              cpuTime=float(timeEnd[0])-float(timeStart[0])  # CPU time , s
	      cpuTimeTotal += cpuTime
  
              print "   cpuTime (s) = ", cpuTime, " for v_i_rel = ", v_i_rel[i], \
	            ", Nsplit = ",Nsplit," and Delta_e_par_arr = ",Delta_e_par_arr[k], \
		    " y0rel = ",y0rel[m]," ffLong = ",ffLong[i,j,k,m]," eV/m", \
                    " ffTrnsvY = ",ffTrnsvY[i,j,k,m]," eV/m"
        #  plt.plot(t_arr, F_ave_arr, 'r-') 
        #  plt.plot(t_arr, 0.0*t_arr, 'g--') 
        #  plt.xlabel('t (s)') 
        #  plt.ylabel('$F_{\parallel} (eV / m)$') 
        #  plt.xlim(0., Tsim) 
        #  plt.show()
  print "cpuTimeTotal (s) = ", cpuTimeTotal
#  plt.show()
#  sys.exit()

  powDelta_e_par = np.zeros(v_e_numb)
  mantDelta_e_par = np.zeros(v_e_numb)
    
  for k in range(v_e_numb):
     powDelta_e_par[k] = round(np.log10(Delta_e_par_arr[k])) 
     mantDelta_e_par[k] = Delta_e_par_arr[k]/(10**powDelta_e_par[k]) 

#----------------------------------------------------
#
# Opening the output file (for 4 values of Nsplit): 
#
##   saveDataFile = nameDevice
###  saveDataFile += 'frictionForceTrnsv_Nsplit_5-6-7-9.dat'
###  saveDataFile += 'frictionForceTrnsv_Nsplit_2-3-4-5.dat'
##  saveDataFile += 'frictionForceTrnsv_Nsplit_5-6-7-9_noBckgrnd.dat'
###  saveDataFile += 'frictionForceTrnsv_Nsplit_2-3-4-5_noBckgrnd.dat'
##  print ('Open output file "%s"...' % saveDataFile)
##  saveDataFile_flag=0
##  try:
##     outfile = open(saveDataFile,'w')
##     saveDataFile_flag=1
##  except:
##     print ('Problem to open output file "%s"' % saveDataFile)
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
##  outfile.write ('\nVion/VeTrnsv               F r i c t i o n   F o r c e   (eV/m)\n')

##  for i in range(v_i_numb):
##     outfile.write ('\n %5.3f       %5.3e     %5.3e     %5.3e     %5.3e'  % \
##                    (v_i_rel[i],ffTrnsv[i,0],ffTrnsv[i,1],ffTrnsv[i,2],ffTrnsv[i,3]))
##  outfile.close()
##  print ('File "%s" is written' % saveDataFile)

#----------------------------------------------------
#
# Opening the output file (for one valuee each of Nsplit and v_e_numb): 
#
#  saveDataFile = nameDevice
#  saveDataFile += 'frictionForceTrnsv_Nsplit_5_ve_3.dat'
#  print ('Open output file "%s"...' % saveDataFile)
#  saveDataFile_flag=0
#  try:
#     outfile = open(saveDataFile,'w')
#     saveDataFile_flag=1
#  except:
#     print ('Problem to open output file "%s"' % saveDataFile)
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
#  outfile.write ('\n\nVion/VeTrnsv            F r i c t i o n    F o r c e     (eV/m)\n')
#
#  for i in range(v_i_numb):
#     outfile.write ('\n    %5.3f         %5.3e         %5.3e         %5.3e'  % \
#                    (v_i_rel[i],ffTrnsv[i,0,0],ffTrnsv[i,0,1],ffTrnsv[i,0,2]))
#  outfile.close()
#  print ('File "%s" is written' % saveDataFile)
#
#----------------------------------------------------

#----------------------------------------------------
#
# Opening the output file (for one values of Nsplit and v_e_numb): 
#
##  saveDataFile = nameDevice
##  saveDataFile += '_coldBeam.dat'
##  print ('Open output file "%s"...' % saveDataFile)
##  saveDataFile_flag=0
##  try:
##     outfile = open(saveDataFile,'w')
##     saveDataFile_flag=1
##  except:
##     print ('Problem to open output file "%s"' % saveDataFile)
#
# Writing the results to output file: 
#
##  outfile.write ('\n                       Initial data:\n\n')
##  outfile.write ('Velocity of Electrons: rmsTran = %5.3e m/s (T = %5.3f eV)' % \
##                 (Delta_e_tr,eT_tr_eV))
##  outfile.write ('\nVelocity of Electrons: rmsLong = %5.3e m/s (T = %5.3f eV)' % \
##                 (Delta_e_par,eT_par_eV))
##  outfile.write ('\nTime of Cooling = %5.3e s During %d Steps' % (Tsim,Nstep))
##  outfile.write ('\nBfield=%4.2f T, Zion = %d, Plasma Density = %5.3e 1/m^3' % (Bz,Z,n_e))
##  outfile.write ('\nLarmor Frequency = %e 1/s, Plasma Frequency = %e 1/s' % (w_L,w_p))
##  outfile.write ('\nLarmor Period = %e s, Plasma Period = %e s;' % (T_L,T_p))
##  outfile.write ('\nTypical Larmor Radius = %e m; Nsplit = %d' % (r_L,int(Nsplit_arr[0])))  
  
##  outfile.write ('\n\nVion/VeTrnsv   Friction Force, eV/m       "Parkhomchuk", eV/m\n')

##  for i in range(v_i_numb):
##     outfile.write ('\n    %5.3f         %5.3e                 %5.3e'  % \
##                   (v_i_rel[i],ffTrnsv[i,0,0],F_P_trnsv_t[i,0]))
##  outfile.close()
##  print ('File "%s" is written' % saveDataFile)
  
#
#----------------------------------------------------


#
# Plotting of results: 
#
#  fig20 = plt.figure(20)
#  plt.plot(v_i_rel[0:v_i_numb], -ffTrnsv[0:v_i_numb,0],'-b',v_i_rel[0:v_i_numb], -ffTrnsv[0:v_i_numb,1],'-m', \
#           v_i_rel[0:v_i_numb], -ffTrnsv[0:v_i_numb,2],'-g',v_i_rel[0:v_i_numb], -ffTrnsv[0:v_i_numb,3],'-k', \
#           v_i_rel[0:v_i_numb], -ffTrnsv[0:v_i_numb,0],'xr',v_i_rel[0:v_i_numb], -ffTrnsv[0:v_i_numb,1],'xr', \
# 	   v_i_rel[0:v_i_numb], -ffTrnsv[0:v_i_numb,1],'xr',v_i_rel[0:v_i_numb], -ffTrnsv[0:v_i_numb,3],'xr') 
#  plt.xlabel('$V_{i\parallel}/v_{e\parallel}}$',color='m',fontsize=14) 
#  plt.ylabel('$-F_{\parallel} (eV / m)$',color='m',fontsize=14) 
#  titleHeader = ('Transversal Friction Force: $V_{e\parallel}=%3.1f\cdot10^{%2d}$ m/s' % \
#                 (mantDelta_e_par,powDelta_e_par))
#  plt.title(titleHeader,color='m',fontsize=14)
#  plt.legend([('$N_{split}$ = %d' % Nsplit_arr[0]),('$N_{split}$ = %d' % Nsplit_arr[1]), \
#              ('$N_{split}$ = %d' % Nsplit_arr[2]),('$N_{split}$ = %d' % Nsplit_arr[3])], \
#    	      loc='upper right',fontsize=14)
#  plt.grid(True)
#  savePictureFile = nameDevice
#  savePictureFile += 'frictionForceTrnsv_fig20.png'     
#  fig20.savefig(savePictureFile)    
#  print 'File "',savePictureFile,' is written'   


#  fig30 = plt.figure(30)
#  plt.plot(v_i_rel, -ffTrnsv[:,0],'-b',v_i_rel, -ffTrnsv[:,1],'-m', \
#           v_i_rel, -ffTrnsv[:,2],'-g',v_i_rel, -ffTrnsv[:,3],'-k', \
#           v_i_rel, -ffTrnsv[:,0],'xr',v_i_rel, -ffTrnsv[:,1],'xr', \
#	   v_i_rel, -ffTrnsv[:,1],'xr',v_i_rel, -ffTrnsv[:,3],'xr') 
#  plt.xlabel('$V_{i\parallel}/v_{e\parallel}}$',color='m',fontsize=14) 
#  plt.ylabel('$-F_{\parallel} (eV / m)$',color='m',fontsize=14) 
#  titleHeader = ('Transversal Friction Force: $V_{e\parallel}=%3.1f\cdot10^{%2d}$ m/s' % \
#                 (mantDelta_e_par,powDelta_e_par))
#  plt.title(titleHeader,color='m',fontsize=14)
#  plt.legend([('$N_{split}$ = %d' % Nsplit_arr[0]),('$N_{split}$ = %d' % Nsplit_arr[1]), \
#              ('$N_{split}$ = %d' % Nsplit_arr[2]),('$N_{split}$ = %d' % Nsplit_arr[3])], \
#   	      loc='lower right',fontsize=14)
#  plt.grid(True)
#  savePictureFile = nameDevice
#  savePictureFile += 'frictionForceTrnsv_fig30_noBckgrnd.png'     
#  fig30.savefig(savePictureFile)    
#  print ('File "',savePictureFile,' is written')   

#  fig40 = plt.figure(40)
#  plt.plot(v_i_rel[0:v_i_numb], -ffTrnsv[0:v_i_numb,0],'-b', \
#           v_i_rel[0:v_i_numb], -ffTrnsv[0:v_i_numb,0],'xr') 
#  plt.xlabel('$V_{i\parallel}/v_{e\parallel}}$',color='m',fontsize=14) 
#  plt.ylabel('$-F_{\parallel} (eV / m)$',color='m',fontsize=14) 
#  titleHeader = ('Transversal Friction Force: $V_{e\parallel}=%3.1f\cdot10^{%2d}$ m/s' % \
#                 (mantDelta_e_par,powDelta_e_par))
#  plt.title(titleHeader,color='m',fontsize=14)
#  plt.legend([('$N_{split}$ = %d' % Nsplit_arr[0])],loc='upper right',fontsize=14)
#  plt.grid(True)
#  savePictureFile = nameDevice
#  savePictureFile += 'frictionForceTrnsv_fig40.png'     
#  fig40.savefig(savePictureFile)    
#  print( 'File "',savePictureFile,' is written' )  

  plotForcesFlag = 1    
#  plotForcesFlag = 0    # no plotting

  if (plotForcesFlag == 1):
#     fig1451 = plt.figure(145)
# The same as 145, but more points in the beginning:
     fig1451 = plt.figure(1451)
     plt.plot(v_i_rel[0:v_i_numb], -1.e-3*ffLong[0:v_i_numb,0,0,0],'-r') 
     plt.plot(v_i_rel[0:v_i_numb], -1.e-3*ffLong[0:v_i_numb,0,0,1],'-b') 
     plt.plot(v_i_rel[0:v_i_numb], -1.e-3*ffLong[0:v_i_numb,0,0,2],'-g') 
     plt.plot(v_i_rel[0:v_i_numb], -1.e-3*ffLong[0:v_i_numb,0,0,3],'-m') 
     plt.xlabel('$V_{i\parallel}/v_{e\parallel}}$',color='m',fontsize=14) 
     plt.ylabel('$-F_{\parallel} (keV / m)$',color='m',fontsize=14) 
     titleHeader = nameDevice
     titleHeader += (' Longitudinal Friction Force: \n$V_{e\parallel}=%3.1f\cdot10^{%2d}$ m/s' % \
                    (mantDelta_e_par,powDelta_e_par))
     plt.title(titleHeader,color='m',fontsize=14)
     plt.legend([('$y_{0,rel}=%4.2f$' % y0rel[0]),('$y_{0,rel}=%4.2f$' % y0rel[1]), \
                 ('$y_{0,rel}=%4.2f$' % y0rel[2]),('$y_{0,rel}=%4.2f$' % y0rel[3])], \
                loc='upper right',fontsize=14)
     plt.grid(True)
     savePictureFile = nameDevice
     savePictureFile += '_LongFF_coldBeam_Simulation_fig145.png'    
     fig1451.savefig(savePictureFile)    
     print( 'File "',savePictureFile,' is written' )  

  if (plotForcesFlag == 0):
     fig146 = plt.figure(146)
     plt.plot(v_i_rel[0:v_i_numb], -1.e-3*ffLong[0:v_i_numb,0,0,0],'-r') 
     plt.plot(v_i_rel[0:v_i_numb], -1.e-3*ffLong[0:v_i_numb,0,0,1],'-b') 
     plt.plot(v_i_rel[0:v_i_numb], -1.e-3*ffLong[0:v_i_numb,0,0,2],'-g') 
     plt.plot(v_i_rel[0:v_i_numb], -1.e-3*ffLong[0:v_i_numb,0,0,3],'-m') 
     plt.xlabel('$V_{i\parallel}/v_{e\parallel}}$',color='m',fontsize=14) 
     plt.ylabel('$-F_{\parallel} (keV / m)$',color='m',fontsize=14) 
     titleHeader = nameDevice
     titleHeader += (' Longitudinal Friction Force: \n$V_{e\parallel}=%3.1f\cdot10^{%2d}$ m/s' % \
                    (mantDelta_e_par,powDelta_e_par))
     plt.title(titleHeader,color='m',fontsize=14)
     plt.xlim([0.,2.])
     plt.legend([('$y_{0,rel}=%4.2f$' % y0rel[0]),('$y_{0,rel}=%4.2f$' % y0rel[1]), \
                 ('$y_{0,rel}=%4.2f$' % y0rel[2]),('$y_{0,rel}=%4.2f$' % y0rel[3])], \
                loc='upper right',fontsize=14)
     plt.grid(True)
     savePictureFile = nameDevice
     savePictureFile += '_LongFF_coldBeam_Simulation_fig146.png'    
##     fig146.savefig(savePictureFile)    
##     print( 'File "',savePictureFile,' is written' )  

#  fig50 = plt.figure(50)
#  plt.plot(v_i_rel[0:v_i_numb], -ffTrnsv[0:v_i_numb,0,0],'-b', \
#           v_i_rel[0:v_i_numb], -ffTrnsv[0:v_i_numb,0,1],'-m', \
#           v_i_rel[0:v_i_numb], -ffTrnsv[0:v_i_numb,0,2],'-g', \
#           v_i_rel[0:v_i_numb], -ffTrnsv[0:v_i_numb,0,0],'xr', \
#           v_i_rel[0:v_i_numb], -ffTrnsv[0:v_i_numb,0,1],'xr', \
#           v_i_rel[0:v_i_numb], -ffTrnsv[0:v_i_numb,0,2],'xr') 
#  plt.xlabel('$V_{i\parallel}/v_{e\parallel}}$',color='m',fontsize=14) 
#  plt.ylabel('$-F_{\parallel} (eV / m)$',color='m',fontsize=14) 
#  titleHeader = nameDevice
#  titleHeader += ('Transversal Friction Force: $N_{split}$ = %d' % (int(Nsplit_arr[0])))
#  plt.title(titleHeader,color='m',fontsize=14)
#  plt.legend([(('$v_{e\parallel}=%3.1f \cdot 10^{%2d}$ m/s' % (mantDelta_e_par[0],powDelta_e_par[0]))), \
#              (('$v_{e\parallel}=%3.1f \cdot 10^{%2d}$ m/s' % (mantDelta_e_par[1],powDelta_e_par[1]))), \
#              (('$v_{e\parallel}=%3.1f \cdot 10^{%2d}$ m/s' % (mantDelta_e_par[2],powDelta_e_par[2])))], \
#    	     loc='upper right',fontsize=14)
#  plt.grid(True)
#  savePictureFile = nameDevice
#  savePictureFile += 'frictionForceTrnsv_fig50.png'    
#  fig50.savefig(savePictureFile)    
#  print( 'File "',savePictureFile,' is written')   


  if (plotForcesFlag == 1):
#     fig1551 = plt.figure(155)
# The same as 155, but more points in the beginning:
     fig1551 = plt.figure(1551)
     plt.plot(v_i_rel[0:v_i_numb], -1.e-3*ffTrnsvY[0:v_i_numb,0,0,0],'-r') 
     plt.plot(v_i_rel[0:v_i_numb], -1.e-3*ffTrnsvY[0:v_i_numb,0,0,1],'-b') 
     plt.plot(v_i_rel[0:v_i_numb], -1.e-3*ffTrnsvY[0:v_i_numb,0,0,2],'-g') 
     plt.plot(v_i_rel[0:v_i_numb], -1.e-3*ffTrnsvY[0:v_i_numb,0,0,3],'-m') 
     plt.xlabel('$V_{i\parallel}/v_{e\parallel}}$',color='m',fontsize=14) 
     plt.ylabel('$-F_{\perp} (keV / m)$',color='m',fontsize=14) 
     titleHeader = nameDevice
     titleHeader += (' Transversal Friction Force: \n$V_{e\parallel}=%3.1f\cdot10^{%2d}$ m/s' % \
                    (mantDelta_e_par,powDelta_e_par))
     plt.title(titleHeader,color='m',fontsize=14)
     plt.legend([('$y_{0,rel}=%4.2f$' % y0rel[0]),('$y_{0,rel}=%4.2f$' % y0rel[1]), \
                 ('$y_{0,rel}=%4.2f$' % y0rel[2]),('$y_{0,rel}=%4.2f$' % y0rel[3])], \
                loc='upper right',fontsize=14)
     plt.grid(True)
     savePictureFile = nameDevice
     savePictureFile += '_TrnsvFF_coldBeam_Simulation_fig1551.png'    
     fig1551.savefig(savePictureFile)    
     print( 'File "',savePictureFile,' is written' )  

  if (plotForcesFlag == 0):
     fig156 = plt.figure(156)
     plt.plot(v_i_rel[0:v_i_numb], -1.e-3*ffTrnsvY[0:v_i_numb,0,0,0],'-r') 
     plt.plot(v_i_rel[0:v_i_numb], -1.e-3*ffTrnsvY[0:v_i_numb,0,0,1],'-b') 
     plt.plot(v_i_rel[0:v_i_numb], -1.e-3*ffTrnsvY[0:v_i_numb,0,0,2],'-g') 
     plt.plot(v_i_rel[0:v_i_numb], 1.e-3*-ffTrnsvY[0:v_i_numb,0,0,3],'-m') 
     plt.xlabel('$V_{i\parallel}/v_{e\parallel}}$',color='m',fontsize=14) 
     plt.ylabel('$-F_{\perp} (keV / m)$',color='m',fontsize=14) 
     titleHeader = nameDevice
     titleHeader += (' Transversal Friction Force: \n$V_{e\parallel}=%3.1f\cdot10^{%2d}$ m/s' % \
                    (mantDelta_e_par,powDelta_e_par))
     plt.title(titleHeader,color='m',fontsize=14)
     plt.xlim([0.,2.])
     plt.legend([('$y_{0,rel}=%4.2f$' % y0rel[0]),('$y_{0,rel}=%4.2f$' % y0rel[1]), \
                 ('$y_{0,rel}=%4.2f$' % y0rel[2]),('$y_{0,rel}=%4.2f$' % y0rel[3])], \
                loc='upper right',fontsize=14)
     plt.grid(True)
     savePictureFile = nameDevice
     savePictureFile += '_TrnsvFF_coldBeam_Simulation_fig156.png'    
##     fig156.savefig(savePictureFile)    
##     print( 'File "',savePictureFile,' is written' )  

#  plt.show()
#  sys.exit()

#------------------- Potting of thrajectories ------------
#
  arg_x = np.linspace(start=0,stop=track_length,num=track_length,dtype=int)
  
  typeLines = ['-r','-b','-m','-k','-g','.r','.b','.m','.k','.g', \
               '-r','-b','-m','-k','-g','.r','.b','.m','.k','.g']
  plotTracksFlag = 1    
  plotTracksFlag = 0    # no plotting
  
  
  if (plotTracksFlag == 1):
     fig160 = plt.figure(160)
     plt.xlabel('Time Step',color='m',fontsize=14)
     plt.ylabel('z, $m\mu$m',color='m',fontsize=14)
     titleHeader = nameDevice
     titleHeader += (' First Trajectories: $y_{0,rel}$=%5.3f' % y0rel[0])
     plt.title(titleHeader,color='m',fontsize=14)
     plt.grid(True)
     for ii in range(np.size(indx_i)):
        i = indx_i[ii]
        plt.plot(arg_x,1.e9*z_traject[0:track_length,0,i,0],typeLines[ii])
     savePictureFile = nameDevice
     savePictureFile += '_FirstPrtcl_z-track_indxY0-0_fig160.png'    
##     fig160.savefig(savePictureFile)    
##     print( 'File "',savePictureFile,' is written' )  
     
  if (plotTracksFlag == 1):
     fig161 = plt.figure(161)
     plt.xlabel('Time Step',color='m',fontsize=14)
     plt.ylabel('z, $m\mu$m',color='m',fontsize=14)
     titleHeader = nameDevice
     titleHeader += (' First Trajectories: $y_{0,rel}$=%5.3f' % y0rel[1])
     plt.title(titleHeader,color='m',fontsize=14)
     plt.grid(True)
     for ii in range(np.size(indx_i)):
        i = indx_i[ii]
        plt.plot(arg_x,1.e9*z_traject[0:track_length,0,i,1],typeLines[ii])
     savePictureFile = nameDevice
     savePictureFile += '_FirstPrtcl_z-track_indxY0-1_fig161.png'    
##     fig161.savefig(savePictureFile)    
##     print( 'File "',savePictureFile,' is written' )  
     
  if (plotTracksFlag == 1):
     fig162 = plt.figure(162)
     plt.xlabel('Time Step',color='m',fontsize=14)
     plt.ylabel('z, $m\mu$m',color='m',fontsize=14)
     titleHeader = nameDevice
     titleHeader += (' First Trajectories: $y_{0,rel}$=%5.3f' % y0rel[2])
     plt.title(titleHeader,color='m',fontsize=14)
     plt.grid(True)
     for ii in range(np.size(indx_i)):
        i = indx_i[ii]
        plt.plot(arg_x,1.e9*z_traject[0:track_length,0,i,3],typeLines[ii])
     savePictureFile = nameDevice
     savePictureFile += '_FirstPrtcl_z-track_indxY0-2_fig162.png'    
##     fig162.savefig(savePictureFile)    
##     print( 'File "',savePictureFile,' is written' )  
     
  if (plotTracksFlag == 1):
     fig163 = plt.figure(163)
     plt.xlabel('Time Step',color='m',fontsize=14)
     plt.ylabel('z, $m\mu$m',color='m',fontsize=14)
     titleHeader = nameDevice
     titleHeader += (' First Trajectories: $y_{0,rel}$=%5.3f' % y0rel[3])
     plt.title(titleHeader,color='m',fontsize=14)
     plt.grid(True)
     for ii in range(np.size(indx_i)):
        i = indx_i[ii]
        plt.plot(arg_x,1.e9*z_traject[0:track_length,0,i,3],typeLines[ii])
     savePictureFile = nameDevice
     savePictureFile += '_FirstPrtcl_z-track_indxY0-3_fig163.png'    
##     fig163.savefig(savePictureFile)    
##     print( 'File "',savePictureFile,' is written' )  
     

#  for i in range(0,v_i_numb,5):

  if (plotTracksFlag == 1):
     fig170 = plt.figure(170)
     plt.xlabel('Time Step',color='m',fontsize=14)
     plt.ylabel('z, $m\mu$m',color='m',fontsize=14)
     titleHeader = nameDevice
     titleHeader += (' Last Trajectories: $y_{0,rel}$=%5.3f' % y0rel[0])
     plt.title(titleHeader,color='m',fontsize=14)
     plt.grid(True)
     for ii in range(np.size(indx_i)):
        i = indx_i[ii]
        plt.plot(arg_x,1.e9*z_traject[0:track_length,1,i,0],typeLines[ii])
     savePictureFile = nameDevice
     savePictureFile += '_LastPrtcl_z-track_indxY0-0_fig170.png'    
##     fig170.savefig(savePictureFile)    
##     print( 'File "',savePictureFile,' is written' )  

  if (plotTracksFlag == 1):
     fig171 = plt.figure(171)
     plt.xlabel('Time Step',color='m',fontsize=14)
     plt.ylabel('z, $m\mu$m',color='m',fontsize=14)
     titleHeader = nameDevice
     titleHeader += (' Last Trajectories: $y_{0,rel}$=%5.3f' % y0rel[1])
     plt.title(titleHeader,color='m',fontsize=14)
     plt.grid(True)
     for ii in range(np.size(indx_i)):
        i = indx_i[ii]
        plt.plot(arg_x,1.e9*z_traject[0:track_length,1,i,1],typeLines[ii])
     savePictureFile = nameDevice
     savePictureFile += '_LastPrtcl_z-track_indxY0-1_fig171.png'    
##     fig171.savefig(savePictureFile)    
##     print( 'File "',savePictureFile,' is written' )  

  if (plotTracksFlag == 1):
     fig172 = plt.figure(172)
     plt.xlabel('Time Step',color='m',fontsize=14)
     plt.ylabel('z, $m\mu$m',color='m',fontsize=14)
     titleHeader = nameDevice
     titleHeader += (' Last Trajectories: $y_{0,rel}$=%5.3f' % y0rel[2])
     plt.title(titleHeader,color='m',fontsize=14)
     plt.grid(True)
     for ii in range(np.size(indx_i)):
        i = indx_i[ii]
        plt.plot(arg_x,1.e9*z_traject[0:track_length,1,i,2],typeLines[ii])
     savePictureFile = nameDevice
     savePictureFile += '_LastPrtcl_z-track_indxY0-2_fig172.png'    
##     fig172.savefig(savePictureFile)    
##     print( 'File "',savePictureFile,' is written' )  

  if (plotTracksFlag == 1):
     fig173 = plt.figure(173)
     plt.xlabel('Time Step',color='m',fontsize=14)
     plt.ylabel('x, $m\mu$m',color='m',fontsize=14)
     titleHeader = nameDevice
     titleHeader += (' Last Trajectories: $y_{0,rel}$=%5.3f' % y0rel[3])
     plt.title(titleHeader,color='m',fontsize=14)
     plt.grid(True)
     for ii in range(np.size(indx_i)):
        i = indx_i[ii]
        plt.plot(arg_x,1.e9*z_traject[0:track_length,1,i,3],typeLines[ii])
     savePictureFile = nameDevice
     savePictureFile += '_LastPrtcl_z-track_indxY0-3_fig173.png'    
##     fig173.savefig(savePictureFile)    
##     print( 'File "',savePictureFile,' is written' )  

  if (plotTracksFlag == 1):
     fig180 = plt.figure(180)
     plt.xlabel('Time Step',color='m',fontsize=14)
     plt.ylabel('y, $\mu$m',color='m',fontsize=14)
     titleHeader = nameDevice
     titleHeader += (' First Trajectories: $y_{0,rel}$=%5.3f' % y0rel[0])
     plt.title(titleHeader,color='m',fontsize=14)
     plt.grid(True)
     for ii in range(np.size(indx_i)):
        i = indx_i[ii]
        plt.plot(arg_x,1.e6*y_traject[0:track_length,0,i,0],typeLines[ii])
     savePictureFile = nameDevice
     savePictureFile += '_FirstPrtcl_y-track_indxY0-0_fig180.png'    
##     fig180.savefig(savePictureFile)    
##     print( 'File "',savePictureFile,' is written' )  
     
  if (plotTracksFlag == 1):
     fig181 = plt.figure(181)
     plt.xlabel('Time Step',color='m',fontsize=14)
     plt.ylabel('y, $\mu$m',color='m',fontsize=14)
     titleHeader = nameDevice
     titleHeader += (' First Trajectories: $y_{0,rel}$=%5.3f' % y0rel[1])
     plt.title(titleHeader,color='m',fontsize=14)
     plt.grid(True)
     for ii in range(np.size(indx_i)):
        i = indx_i[ii]
        plt.plot(arg_x,1.e6*y_traject[0:track_length,0,i,1],typeLines[ii])
     savePictureFile = nameDevice
     savePictureFile += '_FirstPrtcl_y-track_indxY0-1_fig181.png'    
##     fig181.savefig(savePictureFile)    
##     print( 'File "',savePictureFile,' is written' )  
     
  if (plotTracksFlag == 1):
     fig182 = plt.figure(182)
     plt.xlabel('Time Step',color='m',fontsize=14)
     plt.ylabel('y, $\mu$m',color='m',fontsize=14)
     titleHeader = nameDevice
     titleHeader += (' First Trajectories: $y_{0,rel}$=%5.3f' % y0rel[2])
     plt.title(titleHeader,color='m',fontsize=14)
     plt.grid(True)
     for ii in range(np.size(indx_i)):
        i = indx_i[ii]
        plt.plot(arg_x,1.e6*y_traject[0:track_length,0,i,3],typeLines[ii])
     savePictureFile = nameDevice
     savePictureFile += '_FirstPrtcl_y-track_indxY0-2_fig182.png'    
##     fig182.savefig(savePictureFile)    
##     print( 'File "',savePictureFile,' is written' )  
     
  if (plotTracksFlag == 1):
     fig183 = plt.figure(183)
     plt.xlabel('Time Step',color='m',fontsize=14)
     plt.ylabel('y, $\mu$m',color='m',fontsize=14)
     titleHeader = nameDevice
     titleHeader += (' First Trajectories: $y_{0,rel}$=%5.3f' % y0rel[3])
     plt.title(titleHeader,color='m',fontsize=14)
     plt.grid(True)
     for ii in range(np.size(indx_i)):
        i = indx_i[ii]
        plt.plot(arg_x,1.e6*y_traject[0:track_length,0,i,3],typeLines[ii])
     savePictureFile = nameDevice
     savePictureFile += '_FirstPrtcl_y-track_indxY0-3_fig183.png'    
##     fig183.savefig(savePictureFile)    
##     print( 'File "',savePictureFile,' is written' )  
     

  if (plotTracksFlag == 1):
     fig190 = plt.figure(190)
     plt.xlabel('Time Step',color='m',fontsize=14)
     plt.ylabel('y, $\mu$m',color='m',fontsize=14)
     titleHeader = nameDevice
     titleHeader += (' Last Trajectories: $y_{0,rel}$=%5.3f' % y0rel[0])
     plt.title(titleHeader,color='m',fontsize=14)
     plt.grid(True)
     for ii in range(np.size(indx_i)):
        i = indx_i[ii]
        plt.plot(arg_x,1.e6*y_traject[0:track_length,1,i,0],typeLines[ii])
     savePictureFile = nameDevice
     savePictureFile += '_LastPrtcl_y-track_indxY0-0_fig190.png'    
##     fig190.savefig(savePictureFile)    
##     print( 'File "',savePictureFile,' is written' )  

  if (plotTracksFlag == 1):
     fig191 = plt.figure(191)
     plt.xlabel('Time Step',color='m',fontsize=14)
     plt.ylabel('y, $\mu$m',color='m',fontsize=14)
     titleHeader = nameDevice
     titleHeader += (' Last Trajectories: $y_{0,rel}$=%5.3f' % y0rel[1])
     plt.title(titleHeader,color='m',fontsize=14)
     plt.grid(True)
     for ii in range(np.size(indx_i)):
        i = indx_i[ii]
        plt.plot(arg_x,1.e6*y_traject[0:track_length,1,i,1],typeLines[ii])
     savePictureFile = nameDevice
     savePictureFile += '_LastPrtcl_y-track_indxY0-1_fig191.png'    
##     fig191.savefig(savePictureFile)    
##     print( 'File "',savePictureFile,' is written' )  

  if (plotTracksFlag == 1):
     fig192 = plt.figure(192)
     plt.xlabel('Time Step',color='m',fontsize=14)
     plt.ylabel('y, $\mu$m',color='m',fontsize=14)
     titleHeader = nameDevice
     titleHeader += (' Last Trajectories: $y_{0,rel}$=%5.3f' % y0rel[2])
     plt.title(titleHeader,color='m',fontsize=14)
     plt.grid(True)
     for ii in range(np.size(indx_i)):
        i = indx_i[ii]
        plt.plot(arg_x,1.e6*y_traject[0:track_length,1,i,2],typeLines[ii])
     savePictureFile = nameDevice
     savePictureFile += '_LastPrtcl_y-track_indxY0-2_fig192.png'    
##     fig192.savefig(savePictureFile)    
##     print( 'File "',savePictureFile,' is written' )  

  if (plotTracksFlag == 1):
     fig193 = plt.figure(193)
     plt.xlabel('Time Step',color='m',fontsize=14)
     plt.ylabel('y, $\mu$m',color='m',fontsize=14)
     titleHeader = nameDevice
     titleHeader += (' Last Trajectories: $y_{0,rel}$=%5.3f' % y0rel[3])
     plt.title(titleHeader,color='m',fontsize=14)
     plt.grid(True)
     for ii in range(np.size(indx_i)):
        i = indx_i[ii]
        plt.plot(arg_x,1.e6*y_traject[0:track_length,1,i,3],typeLines[ii])
     savePictureFile = nameDevice
     savePictureFile += '_LastPrtcl_y-track_indxY0-3_fig193.png'    
##     fig193.savefig(savePictureFile)    
##     print( 'File "',savePictureFile,' is written' )  

#
#-------------------------------------------------------
  print ('\ntrack_length = %d \n' % track_length) 
#-------------------------------------------------------
#
# Saving of data y-trajectories:
#
  saveTracksFlag = 1
  saveTracksFlag = 0       # no saving
  
  nmbrYtracks = np.size(indx_i)
  
  if (saveTracksFlag == 1):
     for m in range(1,y0_numb):  
#
# Opening the output files for fist trajectories with different initial y0: 
#
        saveDataFile = nameDevice
        saveDataFile += '_first_yTracks_y0-'+'{:d}'.format(int(100.*y0rel[m]))+'_coldBeam.dat'
        print ('\nOpen output file "%s"...' % saveDataFile)
        saveDataFile_flag=0
        try:
           outfile = open(saveDataFile,'w')
           saveDataFile_flag=1
        except:
           print ('Problem to open output file "%s"' % saveDataFile)
#
# Writing the results to output file: 
#
        strOut = '\n         '+nameDevice+':  First y-Tracks (micrometers) for y0_rel = '
        strOut +='{:f}'.format(y0rel[m])+'\n\n' 
        outfile.write (strOut)
        outfile.write ('Transversal Velocity of Electrons: rmsTran = %5.3e m/s (T = %5.3f eV)' % \
                       (Delta_e_tr,eT_tr_eV))
        outfile.write ('\nLongitudinal Velocity of Electrons: rmsLong = %5.3e m/s (T = %5.3f eV)' % \
                       (Delta_e_par,eT_par_eV))
        outfile.write ('\nTime of Cooling = %5.3e s During %d Steps' % (Tsim,Nstep))
        outfile.write ('\nBfield=%4.2f T, Zion = %d, Plasma Density = %5.3e 1/m^3' % (Bz,Z,n_e))
        outfile.write ('\nLarmor Frequency = %e 1/s, Plasma Frequency = %e 1/s' % (w_L,w_p))
        outfile.write ('\nLarmor Period = %e s, Plasma Period = %e s;' % (T_L,T_p))
        outfile.write ('\nTypical Larmor Radius = %e m; Nsplit = %d; nmbrYtracks = %d' % \
                       (r_L,int(Nsplit_arr[0]),nmbrYtracks))  
        strOut = '\n\nVion/VeTrnsv: '
        for ii in range(np.size(indx_i)):
           i = indx_i[ii]
           strOut += '   '+'{:f}'.format(v_i_rel[i]) 
        strOut += '\nPoint of Track              y-Tracks\n'	 
        outfile.write(strOut)
        for j in range(track_length):
           strOut = '\n{:4d}'.format(j)+'          '
	   for ii in range(np.size(indx_i)):
              i = indx_i[ii]
	      strOut += '   '+'{:f}'.format(1.e6*y_traject[j,0,i,m])
#	if (m == 1):   
#	   print (strOut)
           outfile.write(strOut)
        outfile.close()
        print ('Close output file "%s"...' % saveDataFile)
  
     for m in range(1,y0_numb):  
#
# Opening the output files for last trajectories with different initial y0: 
#
        saveDataFile = nameDevice
        saveDataFile += '_last_yTracks_y0-'+'{:d}'.format(int(100.*y0rel[m]))+'_coldBeam.dat'
        print ('\nOpen output file "%s"...' % saveDataFile)
        saveDataFile_flag=0
        try:
           outfile = open(saveDataFile,'w')
           saveDataFile_flag=1
        except:
           print ('Problem to open output file "%s"' % saveDataFile)
#
# Writing the results to output file: 
#
        strOut = '\n         '+nameDevice+':  Last y-Tracks (micrometers) for y0_rel = '
        strOut +='{:f}'.format(y0rel[m])+'\n\n' 
        outfile.write (strOut)
        outfile.write ('Transversal Velocity of Electrons: rmsTran = %5.3e m/s (T = %5.3f eV)' % \
                       (Delta_e_tr,eT_tr_eV))
        outfile.write ('\nLongitudinal Velocity of Electrons: rmsLong = %5.3e m/s (T = %5.3f eV)' % \
                       (Delta_e_par,eT_par_eV))
        outfile.write ('\nTime of Cooling = %5.3e s During %d Steps' % (Tsim,Nstep))
        outfile.write ('\nBfield=%4.2f T, Zion = %d, Plasma Density = %5.3e 1/m^3' % (Bz,Z,n_e))
        outfile.write ('\nLarmor Frequency = %e 1/s, Plasma Frequency = %e 1/s' % (w_L,w_p))
        outfile.write ('\nLarmor Period = %e s, Plasma Period = %e s;' % (T_L,T_p))
        outfile.write ('\nTypical Larmor Radius = %e m; Nsplit = %d; nmbrYtracks = %d' % \
                       (r_L,int(Nsplit_arr[0]),nmbrYtracks))  
        strOut = '\n\nVion/VeTrnsv: '
        for ii in range(np.size(indx_i)):
           i = indx_i[ii]
           strOut += '   '+'{:f}'.format(v_i_rel[i]) 
        strOut += '\nPoint of Track              y-Tracks\n'	 
        outfile.write(strOut)
        for j in range(track_length):
           strOut = '\n{:4d}'.format(j)+'          '
	   for ii in range(np.size(indx_i)):
              i = indx_i[ii]
	      strOut += '   '+'{:f}'.format(1.e6*y_traject[j,1,i,m])
           outfile.write(strOut)
        outfile.close()
        print ('Close output file "%s"...' % saveDataFile)
  
  plt.show()
  sys.exit()

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
    print ("Nsplit = 1? Really?")

#  for k in range(Nsplit):  
#     print ('z_split[%d] = %e' % (k,z_split[k]))   
  return z_split  
  
# For a sphere, dS(z)/dz = 2*pi*r = const; but instead of uniformly sampling the [-r,r], we are
# choosing the z coordinates of the initial conditions so as to match the lowest-order z-moments of
# the spherical shell distribution: for Nsplit odd, z-moments up to order Nsplit are matched, for 
# Nsplit even, z-moments up to order Nsplit+1 are matched 

def aveOverNearestNeibs01(n_e, v0, C0, T, Nstep, k_max, Nsplit,z_split,y0_rel): 

# Output: t_arr,F_ave_arrLong,F_free_arrLong,F_ave_arrTrnsv,F_free_arrTrnsv, \
#         zFirst,zLast,vzFirst,vzLast

  dt = T /np.float64(Nstep) 
  t_arr = np.zeros(Nstep+2, dtype=np.float64) 
  t_arr[-1] = T 
  t_arr[1:Nstep+1] = 0.5*dt +np.linspace(0., T -dt, Nstep) 
  F_cumul_arr_Long = 0.0 *t_arr 
  F_cumul_arr_Trnsv = 0.0 *t_arr 
  
  for ik in np.arange(k_max): 
    k = ik +1 
# expectation value of distance to the k-th nearest neighbor, assuming 
# isotropic const density: 
    rk = np.power(3.*k/(4.*np.pi*n_e), 1./3.)  
#
# Taking into account the initial transversal coordinate y0:
#
    rk = rk*np.sqrt(1.-y0_rel*y0_rel)
    z_k = rk *z_split
    y0 = rk*y0_rel  
    
    for i_sp in np.arange(Nsplit): 
      D0 = np.sqrt(rk*rk -z_k[i_sp]*z_k[i_sp]) 
#
# Return from magDyn1D_leapfr02:
# z, vz, t_arr, z_arr, y_arr, Fz_arr_Long, Fz_arr_Trnsv, t_v_arr, vz_arr, vy_arr   
#
      a1, a2, z_arr, y_arr, y_arr, Fz_arr_Long,  Fz_arr_Trnsv, a5, vz_arr, vy_arr = \
      magDyn1D_leapfr02(z_k[i_sp], y0, v0, C0, D0, T, Nstep)
      F_free_arr_Long,F_free_arr_Trnsv = unpert1D(z_k[i_sp], y0, v0, C0, D0, T, Nstep)
      F_cumul_arr_Long += (Fz_arr_Long - F_free_arr_Long) /np.float64(Nsplit)
      F_cumul_arr_Trnsv += (Fz_arr_Trnsv - F_free_arr_Trnsv) /np.float64(Nsplit)



  
  #F_cumul_arr /= np.float64(k_max) 
    if (k == 1): 
#       print( 'First trajectory (for k=',k,')')
       zFirst = z_arr
       vzFirst = vz_arr 
       yFirst = y_arr
       vyFirst = vy_arr 
    if (k == k_max): 
#       print( 'Last trajectory (for k=',k,')')
       zLast = z_arr
       vzLast = vz_arr 
       yLast = y_arr
       vyLast = vy_arr 
  return t_arr, F_cumul_arr_Long,  F_free_arr_Long, \
                F_cumul_arr_Trnsv, F_free_arr_Trnsv, \
		zFirst, zLast, vzFirst, vzLast, yFirst, yLast, vyFirst, vyLast 

def unpert1D(z_ini, y0, v0, C0, D0, T, Nstep): 
  
  dt = T /np.float64(Nstep) 
  t_arr = np.zeros(Nstep+2, dtype=np.float64) 
  t_arr[-1] = T 
  t_arr[1:Nstep+1] = 0.5*dt +np.linspace(0., T -dt, Nstep) 
  z_free_arr = z_ini+v0*t_arr 
  r = np.sqrt(D0*D0+z_free_arr*z_free_arr+y0*y0) 
  F_free_arr_Long = C0 *z_free_arr /(r*r*r)  # from the electron on the ion, hence the "+" sign 
  F_free_arr_Trnsv = C0 *y0 /(r*r*r)  # from the electron on the ion, hence the "+" sign 
  
  return F_free_arr_Long,F_free_arr_Trnsv 




def magDyn1D_leapfr02(z_ini, y0, v0, C0, D0, T, Nstep): 
  dt = T /np.float64(Nstep) 
  t_arr = np.zeros(Nstep+2, dtype=np.float64) 
  z_arr = 1.0*t_arr  
  y_arr = 1.0*t_arr  
  vz_arr = np.zeros(Nstep+1, dtype=np.float64) 
  vy_arr = np.zeros(Nstep+1, dtype=np.float64) 
  t_v_arr = np.linspace(0.0, T, Nstep+1) 
  Fz_arr_Long = 1.0*t_arr  # from electron on the ion 
  Fz_arr_Trnsv = 1.0*t_arr  # from electron on the ion 
  
  z = z_ini 
  y = y0
  vz = v0
  vy = 0. 
  #print( v0)
  
  z_arr[0] = z 
  vz_arr[0] = vz 
  y_arr[0] = y 
  vy_arr[0] = vy 
  r = np.sqrt(D0*D0+z*z+y*y) 
  r3 = r*r*r
  Fz_arr_Long[0] = C0 *z / r3   #  "+", because the force from electron on the ion 
  Fz_arr_Trnsv[0] = C0 *y / r3  #  "+", because the force from electron on the ion 
  
  z += vz * dt/2.
  y += vy * dt/2.
  
  t_arr[1] = dt/2. 
  z_arr[1] = z 
  y_arr[1] = y 
  r = np.sqrt(D0*D0 +z*z+y*y) 
  r3 = r*r*r
  Fz_arr_Long[1] = C0 *z / r3
  Fz_arr_Trnsv[1] = C0 *y / r3
  
  for istep in np.arange(Nstep-1):
    vz += -dt *C0 *z / r3  # C0 > 0 
    z += vz *dt 
    vy += -dt *C0 *y / r3  # C0 > 0 
    y += vy *dt 
    r = np.sqrt(D0*D0 +z*z+y*y) 
    r3 = r*r*r 
    
    t_arr[2+istep] = (1.5 +istep) *dt 
    z_arr[2+istep] = z 
    y_arr[2+istep] = y 
    Fz_arr_Long[2+istep] = C0 *z / r3 
    Fz_arr_Trnsv[2+istep] = C0 *y / r3 
    vz_arr[1+istep] = vz 
    vy_arr[1+istep] = vy 
  
  vz += -dt *C0 *z / r3  # C0 > 0 
  z += vz *dt/2.
  vy += -dt *C0 *y / r3  # C0 > 0 
  y += vy *dt/2.
  
  t_arr[Nstep+1] = T  
  z_arr[Nstep+1] = z 
  y_arr[Nstep+1] = y 
  r = np.sqrt(D0*D0 +z*z+y*y) 
  r3 = r*r*r
  Fz_arr_Long[Nstep+1] = C0 *z / r3 
  Fz_arr_Trnsv[Nstep+1] = C0 *y / r3 
  vz_arr[Nstep] = vz 
  vy_arr[Nstep] = vy 
  
  return z, vz, t_arr, z_arr, y_arr, Fz_arr_Long, Fz_arr_Trnsv, t_v_arr, vz_arr, vy_arr   



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

