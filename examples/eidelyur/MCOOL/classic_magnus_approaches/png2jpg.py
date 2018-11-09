# from __future__ import division


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
from scipy.integrate import quad, nquad, dblquad

from scipy.constants import pi
from scipy import optimize
from statistics import mean
from array import array

from scipy.constants import speed_of_light as clight
from scipy.constants import epsilon_0 as eps0
from scipy.constants import mu_0 as mu0
from scipy.constants import elementary_charge as qe
from scipy.constants import electron_mass as me
from scipy.constants import proton_mass as mp
from scipy.constants import Boltzmann as kB


# files = os.popen('ls picturesCMA_v7').readlines()
files = os.popen('ls picturesCMA_v7').readlines()
# print ('files = ', files)

numbFiles = len(files)

print ('numbFiles = %d' % numbFiles)

print ('     Files:\n')
nFileWrttn = 0
for n in range(numbFiles):
   fileName = files[n]
   lenName = len(fileName)
#   print ('file(%d) = %s length = %d' % (n,fileName,lenName))
   tail = fileName[lenName-4:lenName-1]
#   print ('tail = "%s"' % tail)
   if (tail != 'png'):
      print ('Not png file(%d) = %s' % (n,fileName))
   else:
      pngName = 'picturesCMA_v7/' + fileName[0:lenName-1]
      jpgName = 'picturesCMA_v7/' + fileName[0:lenName-5] + '.jpg'
#      print ('file(%d) = "%s" --> "%s"' % (n,pngName,jpgName))
      pngFileFlag=0
      try:
         pngFile = open(pngName,'rb')
         pngFileFlag=1
      except:
         print ('Problem to open input file %s' % pngName)
      if pngFileFlag == 1:
#         print ('No problem to open input file "%s"' % pngName)
         fileData = pngFile.read()
         pngFile.close()
         print ('Close input file "%s"' % pngName)
         try:
            jpgFile = open(jpgName,'wb')
#            print ('No problem to open output file "%s"' % jpgName)
            jpgFile.write (fileData)
            jpgFile.close()
            nFileWrttn = nFileWrttn + 1
            print ('Written and closed output file "%s"' % jpgName)
         except:
            print ('Problem to open output file "%s"' % jpgName)

print ('Total written files = %d' % nFileWrttn)

sys.exit()


