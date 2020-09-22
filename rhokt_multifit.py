"""
Script file to read gridded densities from 
files localDens/rho<t>.txt, perform a 2-dim FFT,
and for each k (vector) write time evolution to 
files rhok_timeEv/rhok_<kx>_<ky>.txt. Negative signs 
are possible, like in rhok_-3_-2.txt. The kx and ky in 
the filename are in units of dq=2pi/L

TODO now: Perform curve_fit on the time evolutions, 
report on fit quality and how close they are to expectation.
"""

from scipy.optimize import curve_fit

# ParProgress2tmp.py should be executed first. That will write 
# files "tmp" and "tmp_denstest". Usually done by the "make obs" call
with open("tmp", "r") as f:
    line = f.readline()
    words = line.split("\t")
    obs_anzahl = int(words[0]) #number of reached CPs, not necessarily obs_anzahl in parameter.h
    L = float(words[1])
    lambdac = float(words[2])
    snapgrid_num = int(words[4])
    obs_dt = float(words[5])
    N1 = int(words[6])
    N2 = int(words[7])
    f2_epsgam = float(words[8])
    T = float(words[9])

rho = (N1+N2)/L**2

import numpy as np
from numpy import exp
from calc_rhoIso import readGrid, ft_and_shift
# readGrid takes filename and returns two 2dim-arrays of size (snapgrid_num x snapgrid_num), rho1 and rho2
# ft_and_shift takes a 2dim array of that size, returns one of that size

# giant arrays to store all rho(k,t) in. srhokt is (rho1+rho2)(k,t), drhokt is (rho1-rho2)(k,t). "sum" and "difference".
rho1kt = np.zeros((snapgrid_num, snapgrid_num, obs_anzahl))
rho2kt = np.zeros((snapgrid_num, snapgrid_num, obs_anzahl))
srhokt = np.zeros((snapgrid_num, snapgrid_num, obs_anzahl))
drhokt = np.zeros((snapgrid_num, snapgrid_num, obs_anzahl))

# hard disk effect H(rho) on timescales. In units of 1/sigma**2.
rhomax = 4./np.pi
def HDE(rho):
    return -1./(rhomax-rho) -rho/(rhomax-rho)**2 + 2*rho**2/(rhomax-rho)**3

# Greens Function of the capillary interactions. In units of inverse sigma. Parameter is (k squared) = qx**2 + qy**2
def W(k2):
    return f2_epsgam/ (k2 + lambdac**-2)

def exponential(t, N0, lam, sat):
    """Function to be fitted. Assuming the shape f(t) = N0*exp(lam*t) + sat, exponential decay towards saturation value sat"""
    return N0*exp(lam*t) + sat

for obs in range(obs_anzahl):
    filename = "localDens/rho{}.txt".format(obs)
    # print(filename)
    
    rho1x,rho2x = readGrid(filename)
    rho1k = ft_and_shift(rho1x)
    rho2k = ft_and_shift(rho2x)
    srhok = np.absolute(rho1k + rho2k)
    drhok = np.absolute(rho1k - rho2k)
    
    rho1kt[:,:,obs] = np.absolute(rho1k)
    rho2kt[:,:,obs] = np.absolute(rho2k)
    srhokt[:,:,obs] = srhok
    drhokt[:,:,obs] = drhok
    
dq = 2 * np.pi / L

# for testing purposes, now commented out. Testing q = dq*(1,1), ind is the index for that q.
if True:
    ind = 1 + snapgrid_num//2
    #print(rho1kt[ind,ind,:])
    t   = np.linspace(0, obs_dt*(obs_anzahl-1), obs_anzahl)
    r1t = rho1kt[ind,ind,:]
    r2t = rho2kt[ind,ind,:]
    srt = srhokt[ind,ind,:]
    drt = drhokt[ind,ind,:]
    writeout = np.array([t,r1t,r2t,srt,drt]).transpose()
    np.savetxt('rhok_timeEv/rhok_1_1.txt', writeout, delimiter='\t', 
               header="Time evolution of rho(k,t), k=dq*(see filename), dq=2pi/L={}. \nFormat: t rho1 rho2 (rho1+rho2) (rho1-rho2)\n".format(dq))
    qx=dq
    qy=dq
    q2 = qx*qx + qy*qy
    prediction_sum = -q2 * (1. + rho*HDE(rho)) * T
    prediction_dif = -q2 * (T - rho*W(q2))
    
    N0_st = 0.1*rho
    lam_st = prediction_sum
    saturation_st = rho
    
    optimal, covar = curve_fit(exponential, t, srt, p0=(N0_st, lam_st, saturation_st))
    #print(optimal)
    
    import matplotlib.pyplot as plt
    plt.plot(t, srt, 'ro')
    plt.plot(t, exponential(t, *optimal), 'r-')
    plt.plot(t, exponential(t, optimal[0], lam_st, optimal[2]), 'r--')
    print("Total dens Mode: predicted lam={}, fitted lam={}".format(lam_st, optimal[1]))
    
    N0_st = 0.1*rho
    lam_st = prediction_dif
    saturation_st = 0
    
    optimal, covar = curve_fit(exponential, t, drt, p0=(N0_st, lam_st, saturation_st))
    plt.plot(t, drt, 'bo')
    plt.plot(t, exponential(t, *optimal), 'b-')
    plt.plot(t, exponential(t, optimal[0], lam_st, optimal[2]), 'b--')
    print("Diff dens Mode: predicted lam={}, fitted lam={}".format(lam_st, optimal[1]))
    plt.show()
