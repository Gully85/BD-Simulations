# calculates properties of a set of parameters 
# for BD simulation of capillarily interacting particles

import numpy as np

N1 = int(input("N1 = "))
N2 = int(input("N2 = "))

N = N1+N2

rho = float(input("total density rho="))

L = np.sqrt(N/rho)

print("The box length is L=", L)

lam = float(input("capillary length lambda="))
T = float(input("T*="))

kTgamma_f2 = lam**2*T/2

rhoc = 4./np.pi
kJ = np.sqrt(1./kTgamma_f2) * np.sqrt(rho/(1-3*(rho/rhoc)**2 - 8*(rho/rhoc)**3   ))

print("With this, Jeans' wavenumber is kJ=", kJ)

kc = kJ*np.sqrt(1-1./(lam*kJ)**2)
print("and the critical perturbation size is ", 2.*np.pi/kc)
