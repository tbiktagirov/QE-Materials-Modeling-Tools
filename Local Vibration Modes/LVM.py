"""
LVM (T. Biktagirov): 
Extracts local vibration modes from Quantum Espresso ph.x output.
For detailed theory, see
   1. Alkauskas, A.; Buckley, B. B.; Awschalom, D. D.; Van de Walle, C. G. New J. Phys. 2014, 16, 073026. 
   2. Markham, J. J. Rev. Mod. Phys. 1959, 31, 956.
   3. von Bardeleben, H.J.; Cantin, J.L.; Gerstmann, U.; Schmidt, W.G.; Biktagirov, T. Nano Lett. 2021, 21, 8119-8125.
Usage: run LVM cut-off DYNFILE EXCOORFILE(optional)
   cut-off: cut-off value for the localization measure
   DYNFILE: file containing dynamical matrices (dyn.G)
   EXCOORFILE(optional): '.xsf', '.xyz', or geometry optimization output '.out' 
      containing coordinates of the system in an excited state;
      if specified, projections of the displacement vectors onto 
	  the polarization vectors are evaluated (in A^2/amu^1/2)
Returns:
   File 1 (DYNFILE+'.locparam_hr.dat'):
   {Phonon mode energy (meV), Localisation measure (a.u.), partial Huang-Rhys factors (a.u.)}
   File 2 (DYNFILE+'.spec_dens.dat'):
   {Energy (meV), Spectral density (a.u.)}
   File 3 (DYNFILE+'.spec_func.dat'):
   {Energy (meV), Shape of the phonon sideband (a.u.) with ZPL energy set to zero}   
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.integrate import quad

class Phonon():
   def __init__(self):
      self.numat     = 0
      self.freq      = []
      self.locparam  = []
      self.huangrhys = []
      self.mass      = []
      self.coor      = []

   # Function: Parces DYNFILE (filedyn) and EXCOORFILE (fileex) files
   # Output (tuple): 
   #    self.freq: phonon frequencies, in meV 
   #    self.locparam: phonon localizations (IPR - 1/N), dimensionless, 
   #    self.huangrhys, partial Huang-Rhys factor, dimensionless
   def parse(self, filedyn, fileex):
      angstr2m = 1e-10
      amu2kg   = 1.66054e-27 # kg
      me       = 911.44 # 0.5 a.u.
      hbar     = 1.054571817e-34 # m2 kg / s
      bohr2ang = 0.5291772085876858
      cm2mev   = 0.12398
      cm2hz    = 0.02998*1e12
	   
      f = open(filedyn,'r')
      m = 0
      n = 0
      t = 0
      strnum   = 0
      masstype = []
      for i in f:
         if strnum == 2:
           self.numat = int(i.split()[1])
         if ("'    ") in i:
           #note: masses in dyn.G are in me (electron mass) units
           masstype = np.append(masstype, float(i.split()[3])/me)
           m  = 1
         if (m != 0) and (m <= self.numat) and not (("'    ") in i):
           type      = int(i.split()[1])
           self.mass = np.append(self.mass, masstype[(type-1)])
           tmp       = [float(i.split()[2])*bohr2ang, float(i.split()[3])*bohr2ang, float(i.split()[4])*bohr2ang] 
           self.coor = np.append(self.coor, tmp )
           m += 1
         if ("omega(" in i) or ("freq" in i):
           freq_i    = float(i.split()[-2]) #in cm-1
           self.freq = np.append(self.freq,freq_i)
           polquad = []
           pol     = []
           n  = 0
           if fileex:
                 dR = self.get_dR(fileex)
                 q_i = 0
           t  = 1
           continue
         if (t == 1 ) and (n < self.numat) and not ("omega" in i) and not ("freq" in i):
	   #note: "omega" and "freq" keywords are related to the different versions of ph.x (5 vs 6, 7)
           pol_tmp = [float(i.split()[1]), float(i.split()[3]), float(i.split()[5])]
           polquad = np.append(polquad,np.sqrt(np.dot(pol_tmp,pol_tmp))**4)
           if fileex:
              q_i += np.dot(dR[n,:],pol_tmp)*angstr2m * (self.mass[n]*amu2kg)**(1/2) # in m*kg^1/2
           n += 1
           if n == self.numat:
              self.locparam = np.append(self.locparam,sum(polquad)) #Inverse Participation Ratio
              t = 0
              if fileex:
                 s_i = 0.5/hbar*q_i**2*(freq_i*cm2hz)
                 self.huangrhys = np.append(self.huangrhys,s_i)
         strnum += 1
#      self.freq -= self.freq[0] #for gamma-point calculations only
      self.freq  = self.freq*cm2mev 
      return self.freq, self.locparam, self.huangrhys


   #Function: Calculates atomic displacements (dR) between ground-state and excited state geometries
   def get_dR(self, fileex):
      f = open(fileex,'r')
      n = 0
      t = 0
      strnum  = 0
      coor_ex = []
      for i in f:
         if ".out" in fileex:
            if "Begin final coordinates" in i:
              strnum = 0
              t  = 1
            if (t == 1) and (strnum > 2) and (n < self.numat):
              tmp     = [float(i.split()[1]), float(i.split()[2]), float(i.split()[3])]
              coor_ex = np.append(coor_ex, tmp)
              n += 1
         if ".xsf" in fileex:
            if (strnum > 6)  and (n < self.numat):
              tmp     = [float(i.split()[1]), float(i.split()[2]), float(i.split()[3])]
              coor_ex = np.append(coor_ex, tmp)
              n += 1
            if self.numat == n:
              break
         if ".xyz" in fileex:
            if (strnum > 1) and (n < self.numat):
              tmp = [float(i.split()[1]), float(i.split()[2]), float(i.split()[3])]
              coor_ex = np.append(coor_ex, tmp)
              n += 1
         strnum += 1
      q_i = 0
      tmp = coor_ex - self.coor
      dR  = np.reshape(tmp, (self.numat,3))
      return dR

# Function: Calculate spectral density
# Input:
   #Emin, Emax: energy range, meV
   #SIGMA:      broadening parameter, meV
def get_spectral_density(freq, huangrhys, Emin, Emax, Npoints, SIGMA):
   Enrg = np.linspace(Emin,Emax,Npoints)
   S = np.zeros(Npoints)
   for i in range(0,Npoints):
      s_tmp = 0.0
      for k in range(0,np.size(freq)):
         s_tmp += huangrhys[k]/SIGMA/np.sqrt(np.pi)*np.exp(-0.5*(Enrg[i]-freq[k])**2/SIGMA/SIGMA)
      S[i] = s_tmp
   return Enrg, S


def main(cutoff, filedyn, fileex, SIGMA):
   tmp = Phonon()
   freq, locparam, huangrhys = tmp.parse(filedyn, fileex)

   # Generate output file 1:
   # column1 - Frequencies, meV
   # column2 - Localizations (IPR - 1/N), dimensionless
   # column3 - partial Huang-Rhys factors, dimensionless
   output1 = np.column_stack((freq, locparam - 1/511, huangrhys))
   np.savetxt(filedyn+'.locparam_hr.dat', output1, delimiter=" ")

   # Compute total Huang-Rhys factor
   print("total HR, S = ", sum(huangrhys))

   # Generate output file 2:
   # column1 - Energy, meV
   # column2 - Spectral density
   Enrg, S = get_spectral_density(freq, huangrhys, 0, 300, 1024, SIGMA)
   output2 = np.column_stack((Enrg, S))
   np.savetxt(filedyn+'.spec_dens.dat', output2, delimiter=" ")

   # Generating function:
   S_t     = np.fft.fft(S)
   Freq    = Enrg*241799050402.417 # meV to Hz
   T       = np.fft.fftfreq(Freq.shape[-1])
   S_t_0   = np.sum(huangrhys)
   G       = np.exp(S_t - S_t_0)

   # Generate output file 3:
   # column1 - Energy, eV
   # column2 - Shape of the phonon sideband
   A       = np.fft.ifft(G)
   output3 = np.column_stack((Enrg*1e-3, np.abs(A)))
   np.savetxt(filedyn+'.spec_funct.dat', output3, delimiter=" ")


if __name__ == "__main__":
   SIGMA = 5         # broadening parameter, meV
   if (len(sys.argv) == 3):
      filedyn = sys.argv[2]
      fileex  = ''
      cutoff  = int(sys.argv[1])
   elif (len(sys.argv) == 4):
      filedyn = sys.argv[2]
      fileex  = sys.argv[3]
      cutoff  = int(sys.argv[1])
   else:
      print("Usage: run LVM cut-off DYNFILE EXCOORFILE(optional)")
      sys.exit()
   main(cutoff, filedyn, fileex, SIGMA)
