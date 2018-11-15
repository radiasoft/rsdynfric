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

# nbrForRename = ['5020cma','5030cma','6000cma','6100cma','6200cma','6301cma','6302cma', \
#                 '6400cma','6402cma','7000cma','7100cma','7200cma','7300cma']
# nbrForRename = ['1345cma','1445cma','1545cma','1645cma','5201cma','5202cma','5250cma', \
#                 '5251cma','5350cma','5351cma','5401cma','5402cma','5450cma','5451cma', \
#                 '6000cma','6002cma','6100cma','6102cma','6200cma','6202cma','6301cma', \
#                 '6302cma','6400cma','6402cma']
# nbrForRename = ['506cma','606cma','636cma','666cma','706cma','736cma','766cma','806cma', \
#                 '836cma','866cma','906cma','1006cma']
nbrForRename = ['503cma','603cma','633cma','663cma','703cma','733cma','763cma','803cma', \
                '833cma','863cma','903cma','1003cma']

totalRename = len(nbrForRename)

print ('totalRename = %d' % totalRename)

print ('     Files:\n')
nFileWrttn = 0
for n in range(numbFiles):
   fileName = files[n]
   lenName = len(fileName)
#   print ('file(%d) = %s length = %d' % (n,fileName,lenName))
   tail = fileName[lenName-4:lenName-1]
#   print ('tail = "%s"' % tail)
   if (tail != 'png'):
#      print ('Not png file(%d) = %s' % (n,fileName))
      continue
   else:
#      print ('File(%d) = %s has *png type' % (n,fileName))
      for k in range(totalRename):
         flagRename = fileName.find(nbrForRename[k])
#         print ('or k=%d file %s in finding: result = %d' % (k,fileName,flagRename))
         if (flagRename > 0): 
            print ('file "%s" will be renamed' % fileName)
            pngName = 'picturesCMA_v7/' + fileName[0:lenName-1]
            jpgName = 'picturesCMA_v7/' + fileName[0:lenName-5] + '.jpg'
            print ('file(%d) = %s --> %s' % (n,pngName,jpgName))
            pngFileFlag=0
            try:
               pngFile = open(pngName,'rb')
               pngFileFlag=1
            except:
               print ('Problem to open input file %s' % pngName)
            if pngFileFlag == 1:
#               print ('No problem to open input file "%s"' % pngName)
               fileData = pngFile.read()
               pngFile.close()
               print ('Close input file "%s"' % pngName)
               try:
                  jpgFile = open(jpgName,'wb')
#                  print ('No problem to open output file "%s"' % jpgName)
                  jpgFile.write (fileData)
                  jpgFile.close()
                  nFileWrttn = nFileWrttn + 1
                  print ('Written and closed output file %d: %s' % (nFileWrttn,jpgName))
                  break
               except:
                  print ('Problem to open output file "%s"' % jpgName)

print ('Total written files = %d' % nFileWrttn)

sys.exit()


