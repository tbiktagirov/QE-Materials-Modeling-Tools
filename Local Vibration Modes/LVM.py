
"""
LVM (by T. Biktagirov): 

Extracts local vibration modes from Quantum Espresso ph.x output.

Usage: run LVM cut-off DYNFILE EXCOORFILE(optional)
   cut-off: cut-off value for the localization measure
   DYNFILE: .dyn file containing dynamical matrices
   EXCOORFILE(optional): '.xsf', '.xyz', or geometry optimization output '.out' 
      containing coordinates of the system in an excited state;
      if specified, projections of the displacement vectors onto 
	  the polarization vectors are evaluated (in A^2/amu^1/2)

Returns:
   {LVM number, LVM energy (meV), Projection (angstrom^2/amu^1/2)}
   {Phonon mode energy (meV), Localisation measure (a.u.)}
   
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

class Phonon():
   def __init__(self):
      self.numat = 0
      self.freq = []
      self.locparam = []
      self.proj = []
      self.mass = []
      self.coor = []


   def parse(self, filedyn, fileex):
      f = open(filedyn,'r')
      masstype = []
      strnum = 0
      m = 0
      n = 0
      for i in f:
         if strnum == 2:
           self.numat = int(i.split()[1])
         if ("'    ") in i:
           masstype=np.append(masstype,float(i.split()[3]))
           m = 1
         if (m != 0) and (m <= self.numat) and not (("'    ") in i):
           type = int(i.split()[1])
           self.mass = np.append(self.mass, masstype[(type-1)])
           tmp = [float(i.split()[2]), float(i.split()[3]), float(i.split()[4])]
           self.coor = np.append(self.coor, tmp)
           m += 1
         if ("omega(" in i) or ("freq" in i):
           freq_i = float(i.split()[-2])
           self.freq = np.append(self.freq,freq_i)
           n = 0
           polquad = []
           polsqr = []
         if ("(" in i) and not ("q" in i) and not ("omega" in i) and not ("freq" in i):
	       #"omega" and "freq" keywords are related to the different versions of ph.x
           n += 1
           polar = [float(i.split()[1]), float(i.split()[3]), float(i.split()[5])]
           polquad = np.append(polquad,np.dot(polar,polar)**4)
           polsqr = np.append(polsqr,np.dot(polar,polar)**2)
           if n == self.numat:
              self.locparam = np.append(self.locparam,sum(polquad)/sum(polsqr))
              if fileex:
                 q_i = self.get_proj(polar, fileex)
                 self.proj = np.append(self.proj,q_i)
         strnum += 1
      #self.freq -= self.freq[0] #for gamma point only
      self.freq = self.freq*1.23981e-1
      print(self.proj)
      return self.freq, self.locparam, self.proj


   def get_proj(self, polar, fileex):
      f = open(fileex,'r')
      strnum = 0
      n = 0
      t = 0
      coor_ex = []
      for i in f:
         if ".out" in fileex:
            if "Begin final coordinates" in i:
              strnum = 0
              t = 1
            if (t == 1) and (strnum > 8) and (n < self.numat):
              n += 1
              tmp = [float(i.split()[1]), float(i.split()[2]), float(i.split()[3])]
              coor_ex = np.append(coor_ex, tmp)
         if ".xsf" in fileex:
            if (strnum > 6)  and (n < self.numat):
              n += 1
              tmp = [float(i.split()[1]), float(i.split()[2]), float(i.split()[3])]
              coor_ex = np.append(coor_ex, tmp)
            if self.numat == n:
              break
         if ".xyz" in fileex:
            if (strnum > 1) and (n < self.numat):
              n += 1
              tmp = [float(i.split()[1]), float(i.split()[2]), float(i.split()[3])]
              coor_ex = np.append(coor_ex, tmp)
         strnum += 1
      q_i = 0
      tmp = coor_ex - self.coor
      dR = np.reshape(tmp, (self.numat,3))
      k = 0
      while k < self.numat:
         q_i += np.dot(dR[k,:],polar) / (self.mass[k])**(1/2)
         k += 1
      return q_i


def get_lvm(freq, locparam, cutoff):
      lvmind = np.where(locparam > np.average(locparam)*cutoff)
      lvm = freq[(lvmind)] 
      return lvm


def get_lvm_proj(proj, locparam, cutoff):
      lvmind = np.where(locparam > np.average(locparam)*cutoff)
      lvm_proj = proj[(lvmind)] 
      return lvm_proj


def main(cutoff, filedyn, fileex):
   tmp = Phonon()
   freq, locparam, proj = tmp.parse(filedyn, fileex)
   output1 = np.column_stack((freq, locparam))
   np.savetxt(filedyn+'.locparam.csv',output1,delimiter=",")

   lvm = get_lvm(freq, locparam, cutoff)
   if fileex:
      lvm_proj = get_lvm_proj(proj, locparam, cutoff)
      output2 = np.column_stack((np.arange(1, lvm.shape[0]+1), lvm, lvm_proj))
   else:
      output2 = np.column_stack((np.arange(1, lvm.shape[0]+1), lvm))
   np.savetxt(filedyn+'.lvm.csv',output2,delimiter=",")

   plt.bar(freq, locparam/np.average(locparam), color='r', label='LVM')
   plt.hold('on')
   plt.plot([0,freq[-1]+10],[cutoff,cutoff],'k--')
   plt.text(1,cutoff,'cutoff')
   plt.xlabel('Frequency, meV')
   plt.ylabel('Localization/mean')
   plt.show()


if __name__ == "__main__":
   if (len(sys.argv) == 3):
      filedyn = sys.argv[2]
      fileex = ''
      cutoff = int(sys.argv[1])
   elif (len(sys.argv) == 4):
      filedyn = sys.argv[2]
      fileex = sys.argv[3]
      cutoff = int(sys.argv[1])
   else:
      print("Usage: run LVM cut-off DYNFILE EXCOORFILE(optional)")
      sys.exit()
   main(cutoff, filedyn, fileex)
