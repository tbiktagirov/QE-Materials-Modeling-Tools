"""

The script parses .xml files generated by pw_export.x
to extract wave-functions FFT coefficients and the corresponding
k-space coordinates.

"""

import numpy as np
from xml.dom import minidom


class WfcDat():
   def __init__(self):
      self.recip = []
      #self.kpoint = [0,0,0]
      self.G = []
      self.Coef = []
      #self.R = []
      #self.wfc = []


   def parsewfc(self, filewfc, nbnd):
      xmldoc = minidom.parse(filewfc)
      tmp = []
      for i in range(1,nbnd+1):
         nodename = "Wfc."+str(i)
         node = xmldoc.getElementsByTagName(nodename)[0]
         C = getCoef(self, node)
         print(nodename)
         tmp = np.append(tmp, C)
      self.Coef = tmp.reshape((nbnd,-1)) # [band, index]
      return self.Coef
#      wfc = np.append(wfc, xmldoc.getElementsByTagName(nodename))


   def getG(self, filegrid, recip):
      self.G = []
      self.recip = recip
      xmldoc = minidom.parse(filegrid)
      tmp = []
      nodelist = xmldoc.getElementsByTagName("grid")[0].childNodes
      for node in nodelist:
         tmp.append(node.data)
      lines = tmp[0].split('\n')
      for x in lines:
         try:
            int(x.split()[0])
            coef = [int(x.split()[0]),int(x.split()[1]),int(x.split()[2])]
            vect = recip[0,:]*coef[0] + recip[1,:]*coef[1] + recip[2,:]*coef[2]
            self.G = np.append(self.G,vect)
         except:
            pass
      self.G = self.G.reshape((-1,3)) #[index, :xyz]
      return self.G


def getCoef(self, node):
         nodelist = node.childNodes
         tmp = []
         C = []
         for node in nodelist:
            tmp.append(node.data)
         #self.Coef = tmp
         lines = tmp[0].split('\n')
         #print(tmp2)
         for x in lines:
            try:
               float(x.split(',')[0])
               coefval = complex(float(x.split(',')[0]),float(x.split(',')[1]))
               C = np.append(C, coefval)
            except:
               pass
         return C


def getrecip(latvec):
   recip = np.zeros((3,3))
   vol = np.dot(latvec[0][:],np.cross(latvec[1,:],latvec[2,:]))
   recip[0,:] = 2*np.pi*np.cross(latvec[1,:],latvec[2,:])/vol
   recip[1,:] = 2*np.pi*np.cross(latvec[2,:],latvec[0,:])/vol
   recip[2,:] = 2*np.pi*np.cross(latvec[0,:],latvec[1,:])/vol
   return recip


if ( __name__ == "__main__"):
   """
   An example of usage. 
   """
   alat = 10
   latvec=np.array([[1, 0, 0],[0, 1, 0],[0, 0, 1]]) * alat
   recip = getrecip(latvec)   
   filewfc = 'wfc.1'
   filegrid = 'grid.1'
   nbnd = 10

   tmp = WfcDat()
   Coef = tmp.parsewfc(filewfc,nbnd)
   G = tmp.getG(filegrid,recip) 
