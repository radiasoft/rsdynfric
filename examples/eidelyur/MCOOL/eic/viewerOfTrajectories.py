def main(): 
  
#
# Opening the input file: 
#
   nameDevice = 'EIC'
   inputFile=nameDevice+'_first_yTracks_y0-50_coldBeam.dat'
   print ('Open input file "%s"...' % inputFile)
   inpfileFlag=0
   try:
      inpfile = open(inputFile,'r')
      inpfileFlag=1
   except:
      print ('Problem to open input file "%s"' % inputFile)
   if inpfileFlag == 1:
      print ('No problem to open input file "%s"' % inputFile)

   lines=0                                                            # Number of current line from input file   
   pointOfTrack=0 
   
   y0_relFlag = 0
   rmsLongFlag = 0
   trackLengthFlag = 0
   nmbrYtracksFlag = 0                                                      # Number of current value of any types of Data
# Data arrays (with a margin):	     
   vIon_rel = []
   yTracks = []
   while True:
      lineData=inpfile.readline()
      if not lineData:
         break
      lines += 1
#      print ('line=%d: %s' % (lines,lineData))
      if (lines < 14):  
         words=lineData.split()
	 if (words != ''):
            nWords=len(words)
#            print ('Data from %d: words=%s, number of entries = %d' % (lines,words,nWords))
            for k in range(nWords):
               if (words[k] == 'y0_rel'):
                  y0_rel = float(words[k+2])
                  y0_relFlag = 1
               if (words[k] == 'rmsLong'):
                  rmsLong = float(words[k+2])
                  rmsLongFlag = 1
               if (words[k] == 'Steps'):
                  trackLength = int(words[k-1])+2
                  trackLengthFlag = 1
               if (words[k] == 'nmbrYtracks'):
                  nmbrYtracks = int(words[k+2])
                  nmbrYtracksFlag = 1
               if (words[k] == 'Vion/VeTrnsv:'):
                  for j in range(nmbrYtracks):
                     vIon_rel[j] = float(words[j+1])   		  	       
      if (lines == 11): 
         if (y0_relFlag*rmsLongFlag*trackLengthFlag*nmbrYtracksFlag == 0):
	    print 'Wrong data in the Header:' 
	    print 'y0_relFlag = ', y0_relFlag     
	    print 'rmsLongFlag = ', rmsLongFlag     
	    print 'trackLengthFlag = ', trackLengthFlag     
	    print 'nmbrYtracksFlag = ', nmbrYtracksFlag  
            sys.exit()
         else:
	    print 'y0_rel = ', y0_rel     
	    print 'rmsLong = ', rmsLong     
	    print 'trackLength = ', trackLength     
	    print 'nmbrYtracks = ', nmbrYtracks 
# Data arrays:	     
            vIon_rel = np.zeros(nmbrYtracks)
            yTracks = np.zeros((trackLength,nmbrYtracks))
#      if (lines == 14): 
#         print 'vIon_rel: ', vIon_rel
      if (lines > 14):
         words=lineData.split()
         nWords=len(words)
#         print ('Data from %d: words=%s, number of entries = %d' % (lines,words,nWords))
         for j in range(nmbrYtracks):
            yTracks[pointOfTrack,j]=float(words[j+1])
         pointOfTrack += 1

   inpfile.close()
   print ('Close input file "%s"' % inputFile)

   arg_x = np.linspace(start=0,stop=trackLength,num=trackLength,dtype=int)
   typeLines = ['-r','-b','-m','-k','-g','.r','.b','.m','.k','.g', \
               '-r','-b','-m','-k','-g','.r','.b','.m','.k','.g']

# For checking of the rrading:	       
   fig181 = plt.figure(181)
   plt.xlabel('Time Step',color='m',fontsize=14)
   plt.ylabel('y, $\mu$m',color='m',fontsize=14)
   titleHeader = nameDevice
   titleHeader += (' First Trajectories: $y_{0,rel}$=%5.3f' % y0_rel)
   plt.title(titleHeader,color='m',fontsize=14)
   plt.grid(True)
   for i in range(0,nmbrYtracks,3):
      plt.plot(arg_x,yTracks[0:trackLength,i],typeLines[i])

   plt.show()
   sys.exit()

if __name__=="__main__":
  import os, sys
  import numpy as np
  import math
  import matplotlib.pyplot as plt 
  from matplotlib.legend_handler import HandlerLine2D

  main() 
  
