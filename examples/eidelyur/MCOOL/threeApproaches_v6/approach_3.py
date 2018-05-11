#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Approach_3: dragging with averaging over nLarmorAvrgng larmor rotation +
#             "Magnus expansion" method to calculate the transferred momenta
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#------- Start of calculations for approach_3 --------------
         if runFlagApproach_3 == 1:
            timeStart=os.times()
            if trackNumb_3 == 0:
	       rhoFirstTurn=rhoCrrnt[iA,iB]
	       rhoLarmorFirstTurn=rho_larm[iA,iB]
#
# 6-vectors for ion and electron and distance 'b' between them for the first trajectory and 
# (T)rajectory with (M)aximal (T)ransferred (dpx) ( trajectory TMTdpx);
# (for checking only; indices 0-5 for electron, indices 6-11 for ion,
#  index=12 for 'b', index=13 for action and index=14 for dy_gc):
#
               prtclCoorFirst_3=np.zeros((15,timePoints_3[trackNumb_3])) # first trajectory                 
            prtclCoorCrrnt_3=np.zeros((15,timePoints_3[trackNumb_3])) # current trajectory                
            prtclCoorMaxAbsDpx_3=np.zeros((15,timePoints_3[trackNumb_3]))      # trajectory TMTdpx
# Current distance from origin of the coordinate system to electron along the trajectory; cm
            bCrrnt_3=np.zeros(timePoints_3[trackNumb_3])           # cm
# Current log10 of two important ratios:
# First - ratio R_larmor/b; dimensionless
            larmR_bCrrnt_3=np.zeros(timePoints_3[trackNumb_3])       
# Second - ratio potential_energy/kinetic_energy; dimensionless
            uPot_enrgKinCrrnt_3=np.zeros(timePoints_3[trackNumb_3])  
# deltaPapprch_3=dpApprch_3Crrnt
            dpApprch_3Crrnt=np.zeros((3,timePoints_3[trackNumb_3]))
            for m in range(6): 
               z_ionCrrnt_3[m]=0.                     # Current initial zero-vector for ion
               z_elecCrrnt_3[m]=0.                    # Zeroing out of vector for electron
# Current initial vector for electron:
            z_elecCrrnt_3[Ix]=rhoCrrnt[iA,iB]+rho_larm[iA,iB]      # x, cm
            z_elecCrrnt_3[Iz]=-halfLintr[iA,iB]                    # z, cm
            z_elecCrrnt_3[Ipy]=m_elec*evTran[iA,iB]                # py, g*cm/sec
            z_elecCrrnt_3[Ipz]=m_elec*eVrmsLong                    # pz, g*cm/sec
	    z_elecCrrnt_gc=toGuidingCenter(z_elecCrrnt_3)          # transfer to system of guiding center
#	    if iA == 0 and iB == 0:
#	       print 'z_elecCrrnt_3: ', z_elecCrrnt_3 
#	       print 'z_elecCrrnt_gc: ', z_elecCrrnt_gc 
#-----------------------------------------------
# Main action - dragging of the current trajectories (for given i and j)
#
            for k in range(int(timePoints_3[trackNumb_3])):
#	       if k < 100: 
#	          print 'k=%d: x=%e, y=%e, z=%e' % \
#                      (k,1.e+4*z_elecCrrnt_gc[0],1.e+4*z_elecCrrnt_gc[2],1.e+4*z_elecCrrnt_gc[4]) 
#	       if k > timePoints_2[trackNumb_3]-100: 
#	          print 'k=%d: x=%e, y=%e, z=%e' % \
#                       (k,1.e+4*z_elecCrrnt_gc[0],1.e+4*z_elecCrrnt_gc[2],1.e+4*z_elecCrrnt_gc[4]) 
#
# dragging both paticles through first half of trajectory:
#
 	       z_elecCrrnt_gc=np.dot(matr_elec_3,z_elecCrrnt_gc)   # electron
 	       z_ionCrrnt_3=np.dot(matr_ion_3,z_ionCrrnt_3)        # ion
#
# dragging both paticles through interaction point:
#
	       dpIon,dpElec,action,dy_gc=MagnusExpansionCollision(z_elecCrrnt_gc,z_ionCrrnt_3,timeStep_3) 
#	       if trackNumb_3 == 0:
#	          print 'point %d: dpgcElec=%e, dpzElec=%e' % 
#                       (pointTrack_3[0],dpElec[1],dpElec[2])
	       for ic in range(3):
####
#### Taking into account transfer of momentum for both particles:
####
####	          z_ionCrrnt_3[2*ic+1] += dpIon[ic]   
####	          z_elecCrrnt_gc[2*ic+1] += dpElec[ic]
#   
# Current values to calculate deltaPapprch_3:  
#
	          dpApprch_3Crrnt[ic,k]=dpIon[ic]                  # g*cm/sec 
#
# dragging both paticles through second half of trajectory:
#
 	       z_elecCrrnt_gc=np.dot(matr_elec_3,z_elecCrrnt_gc)   # electron
 	       z_ionCrrnt_2=np.dot(matr_ion_3,z_ionCrrnt_3)        # ion
#	       crrntPoint=pointTrack_3[trackNumb_3]
#	       if iA == 0 and iB == 0 and crrntPoint < 10:
#                     print 'k, z_ionCrrnt_3', (k,z_ionCrrnt_3)
	       z_elecCrrnt_3=fromGuidingCenter(z_elecCrrnt_gc)     # transfer from system of guiding center 
# Current distance between ion and electron; cm:
 	       bCrrnt_3[k]=np.sqrt((z_ionCrrnt_3[0]-z_elecCrrnt_3[0])**2+ \
	                           (z_ionCrrnt_3[2]-z_elecCrrnt_3[2])**2+ \
			           (z_ionCrrnt_3[4]-z_elecCrrnt_3[4])**2)
# To be prepared to draw future TMTdpx trajectory (for checking only):
               for ic in range(6):
                  prtclCoorCrrnt_3[ic,pointTrack_3[trackNumb_3]]=z_elecCrrnt_3[ic]  # 6-vector for electron
                  prtclCoorCrrnt_3[6+ic,pointTrack_3[trackNumb_3]]=z_ionCrrnt_3[ic] # 6-vector for ion 
	       prtclCoorCrrnt_3[12,pointTrack_3[trackNumb_3]]=bCrrnt_3[k]           # cm
	       prtclCoorCrrnt_3[13,pointTrack_3[trackNumb_3]]=action                # g*cm^2/sec
	       prtclCoorCrrnt_3[14,pointTrack_3[trackNumb_3]]=dy_gc                 # cm
# To draw only first trajector3 (for checking only):
               if trackNumb_3 == 0:
                  for ic in range(6):
                     prtclCoorFirst_3[ic,pointTrack_3[trackNumb_3]]=z_elecCrrnt_3[ic]      # 6-vector for electron
                     prtclCoorFirst_3[6+ic,pointTrack_3[trackNumb_3]]=z_ionCrrnt_3[ic]     # 6-vector for ion 
	          prtclCoorFirst_3[12,pointTrack_3[trackNumb_3]]=bCrrnt_3[k]               # cm
	          prtclCoorFirst_3[13,pointTrack_3[trackNumb_3]]=action                    # g*cm^2/sec
	          prtclCoorFirst_3[14,pointTrack_3[trackNumb_3]]=dy_gc                     # cm
	          if maxYcoorElec_3 < abs(z_elecCrrnt_3[2]):
	             maxYcoorElec_3=abs(z_elecCrrnt_3[2])
	          if maxYcoorIon_3 < abs(z_ionCrrnt_3[2]):
	             maxYcoorIon_3=abs(z_ionCrrnt_3[2])
               pointTrack_3[trackNumb_3] += 1
# End of dragging of the current trajectory	  
#-----------------------------------------------
# To draw the distribution of transferred dp (g*cm/sec):
               dpxTotal_3[iA,iB] +=-dpApprch_3Crrnt[0,k] 
               dpyTotal_3[iA,iB] +=-dpApprch_3Crrnt[1,k] 
               dpzTotal_3[iA,iB] +=-dpApprch_3Crrnt[2,k] 
#
# Total transferred dpx  for current track:
#
	       totAbsDpxApprch_3 += abs(dpApprch_3Crrnt[0,k])
# 
# End of dragging of the current trajectory	  
#-----------------------------------------------
#
# Accumulate transferred momenta and other data: 
#
            if trackNumb_3 == 0: 
# First definition of the total distance from origin of the
# coordinate system to electron along the trajectory; cm:
 	       b_3=bCrrnt_3                                        # cm
# First definition of the log10 of two important ratios; dimensionless:  
	       uPot_enrgKin_3=uPot_enrgKinCrrnt_3                  # ratio A 
	       if maxUpot_enrgKin_3 < uPot_enrgKinCrrnt_3[k]:
	          maxUpot_enrgKin_3=uPot_enrgKinCrrnt_3[k]
	       if minUpot_enrgKin_3 > uPot_enrgKinCrrnt_3[k]:
	          minUpot_enrgKin_3=uPot_enrgKinCrrnt_3[k]
	       larmR_b_3=larmR_bCrrnt_3                            # ratio B
	       if maxLarmR_b_3 < larmR_bCrrnt_3[k]:
	          maxLarmR_b_3=larmR_bCrrnt_3[k]
	       if minLarmR_b_3 > larmR_bCrrnt_3[k]:
	          minLarmR_b_3=larmR_bCrrnt_3[k]
#### First definition of the values to calculate deltaPapprch_3 (g*cm/sec):
###               dpxApprch_3=dpApprch_3Crrnt[0,:] 
###               dpyApprch_3=dpApprch_3Crrnt[1,:]
###               dpzApprch_3=dpApprch_3Crrnt[2,:] 
            else:  
#### Total distance from origin of the coordinate system to electron along the trajectory:
 	       b_3=np.concatenate((b_3,bCrrnt_3),axis=0)           # cm
# Total log10 of two important ratios; dimensionless:  
	       uPot_enrgKin_3=np.concatenate((uPot_enrgKin_3,uPot_enrgKinCrrnt_3),axis=0)        
	       larmR_b_3=np.concatenate((larmR_b_3,larmR_bCrrnt_3),axis=0)                  
#### Total values to calculate deltaPapprch_3 (g*cm/sec):
###               dpxApprch_3=np.concatenate((dpxApprch_3,dpApprch_3Crrnt[0,:]),axis=0)
###               dpyApprch_3=np.concatenate((dpyApprch_3,dpApprch_3Crrnt[1,:]),axis=0)
###               dpzApprch_3=np.concatenate((dpzApprch_3,dpApprch_3Crrnt[2,:]),axis=0)
#
# To draw TMTdpx trajectory (for checking only):
#
            if maxAbsDpxApprch_3 < totAbsDpxApprch_3:
	       maxAbsDpxApprch_3=totAbsDpxApprch_3
	       indxAmaxAbsDpxApprch_3=iA
	       indxBmaxAbsDpxApprch_3=iB
	       trackNumbMaxAbsDpxApprch_3=trackNumb_3
	       prtclCoorMaxAbsDpx_3=prtclCoorCrrnt_3.copy()   
	       rhoMaxAbsDpxTurn_3=rhoCrrnt[iA,iB]
	       rhoLarmorMaxAbsDpxTurn_3=rho_larm[iA,iB]
###	       print 'iA=%d, iB=%d: TMT-track %d, points %d' % \
###	             (indxAmaxAbsDpxApprch_3,indxBmaxAbsDpxApprch_3, \
###		     trackNumbMaxAbsDpxApprch_3,pointTrack_3[trackNumbMaxAbsDpxApprch_3]) 
###	       print 'timePoints.shape: ', \
###	             (prtclCoorMaxAbsDpx_3.shape[0],prtclCoorMaxAbsDpx_3.shape[1])
# End of all calculations for approach_3
            sumPoints_3 += pointTrack_3[trackNumb_3]
# 
# To draw the distribution of transferred dp:
#
	    if dpxTotalMax_3 < abs(dpxTotal_3[iA,iB]):
	       dpxTotalMax_3=abs(dpxTotal_3[iA,iB])
	    if dpyTotalMax_3 < abs(dpyTotal_3[iA,iB]):
	       dpyTotalMax_3=abs(dpyTotal_3[iA,iB])
	    if dpzTotalMax_3 < abs(dpzTotal_3[iA,iB]):
	       dpzTotalMax_3=abs(dpzTotal_3[iA,iB])
            print 'Track %d: dpxTotalMax_3=%e, dpyTotalMax_3=%e, dpzTotalMax_3=%e' % \
	          (trackNumb_3,dpxTotalMax_3,dpyTotalMax_3,dpzTotalMax_3)
            timeEnd=os.times()
	    cpuTime_3[trackNumb_3]=1.e+6*(float(timeEnd[0])-float(timeStart[0]))  # CPU time , mks
            cpuTimeTotal += cpuTime_3[trackNumb_3] 
#
#------- End of approach_3 --------------
#
