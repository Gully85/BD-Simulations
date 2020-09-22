"""
Script file to read gridded densities from 
files localDens/rho<t>.txt, perform a 2-dim FFT,
and for each k (vector) write time evolution to 
files rhok_timeEv/rhok_<kx>_<ky>.txt. Negative signs 
are possible, like in rhok_-3_-2.txt. The kx and ky in 
the filename are in units of dq=2pi/L
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

# hold q (abs value) and fitted decay constants (sum and dif mode)
fitresults = np.zeros((snapgrid_num*snapgrid_num, 3))

failed_fits     = 0
successful_fits = 0
# loop 2dim over all k
for i in range(snapgrid_num):
    qx = dq*(i - snapgrid_num//2)
    if i%10 == 0:
        print("Progress: {}/{}. Successful fits: {}/{}".format(i,snapgrid_num, successful_fits, (successful_fits+failed_fits)))
    for j in range(snapgrid_num):
        qy = dq*(j - snapgrid_num//2)
        
        index = i*snapgrid_num + j
        #print(rho1kt[ind,ind,:])
        t   = np.linspace(0, obs_dt*(obs_anzahl-1), obs_anzahl)
        r1t = rho1kt[i,j,:]
        r2t = rho2kt[i,j,:]
        srt = srhokt[i,j,:]
        drt = drhokt[i,j,:]
        writeout = np.array([t,r1t,r2t,srt,drt]).transpose()
        np.savetxt('rhok_timeEv/rhok_{}_{}.txt'.format(i-snapgrid_num//2, j-snapgrid_num//2), writeout, delimiter='\t', 
                header="Time evolution of rho(k,t), k=dq*(see filename), dq=2pi/L={}. \nFormat: t rho1 rho2 (rho1+rho2) (rho1-rho2)\n".format(dq))
        #qx=dq
        #qy=dq
        q2 = qx*qx + qy*qy
        fitresults[index,0] = np.sqrt(q2)
        prediction_sum = -q2 * (1. + rho*HDE(rho)) * T
        prediction_dif = -q2 * (T - rho*W(q2))
        
        N0_st = 0.1*rho
        lam_st = prediction_sum
        saturation_st = rho
        
        try:
            optimal, covar = curve_fit(exponential, t, srt, p0=(N0_st, lam_st, saturation_st))
            fitresults[index,1] = optimal[1]
            successful_fits = successful_fits+1
        except RuntimeError:
            #print("Fit failed: Sum mode, q=({},{})".format(i-snapgrid_num//2, j-snapgrid_num//2))
            fitresults[index,1] = np.nan
            failed_fits = failed_fits+1
        #print(optimal)
        
        
        #import matplotlib.pyplot as plt
        #plt.plot(t, srt, 'ro')
        #plt.plot(t, exponential(t, *optimal), 'r-')
        #plt.plot(t, exponential(t, optimal[0], lam_st, optimal[2]), 'r--')
        #print("Total dens Mode: predicted lam={}, fitted lam={}".format(lam_st, optimal[1]))
        
        N0_st = 0.1*rho
        lam_st = prediction_dif
        saturation_st = 0
        
        try:
            optimal, covar = curve_fit(exponential, t, drt, p0=(N0_st, lam_st, saturation_st))
            fitresults[index,2] = optimal[1]
            successful_fits=successful_fits+1
        except RuntimeError:
            #print("Fit failed: Difference mode, q=({},{})".format(i-snapgrid_num//2, j-snapgrid_num//2))
            fitresults[index,2]
            failed_fits = failed_fits+1
        
        #plt.plot(t, drt, 'bo')
        #plt.plot(t, exponential(t, *optimal), 'b-')
        #plt.plot(t, exponential(t, optimal[0], lam_st, optimal[2]), 'b--')
        #print("Diff dens Mode: predicted lam={}, fitted lam={}".format(lam_st, optimal[1]))
        #plt.show()

#writeout = np.array()
np.savetxt("fitresults.txt", fitresults, delimiter='\t', header="Results of exponential fits to sum (rho1+rho2)(k) and and dif (rho1-rho2)(k) density modes \nFormat: q lam_s lam_d, where lam_s is the decay constant of the sum mode and lam_d the one of the difference mode")

# Produce binned version: Combine q's with similar abs value.
qmax = dq*(snapgrid_num//2)*np.sqrt(2.)
qbin = 0.05
q_numbins = int(qmax/qbin)+1

# 5 infos per qbin: |q|, lam_s avg, lam_s errorbar, lam_d avg, lam_d errorbar
fitresults_binned = np.zeros((q_numbins, 5))
for qi in range(q_numbins):
    print("Binning Progress: {}/{}".format(qi, q_numbins))
    # current qbin: (qi-0.5)*qbin < q < (qi+0.5)*qbin
    fitresults_binned[qi,0] = qi*qbin
    
    # quick and dirty: walk twice through fitresults. Sum up and count contribution on first walkthrough
    # qssum means: For current q-value, sum mode, sum of the decay constants. qscont is for current q-value, sum mode, number of contributions
    qssum=0.0
    qscont=0
    qdsum=0.0
    qdcont=0
    for i in range(snapgrid_num):
        qx = dq*(i-snapgrid_num//2)
        for j in range(snapgrid_num):
            qy = dq*(j-snapgrid_num//2)
            q2 = qx*qx + qy*qy
            if not all([(qi-0.5)**2 * qbin**2 < q2, (qi+0.5)**2 * qbin**2 >= q2]):
                continue
             # index in fitresults[:,]. fitresults[,1] is the decay constant of the sum mode, fitresults[,2] the one of the difference mode
            index = i*snapgrid_num+j
            
            if np.isfinite(fitresults[index, 1]):
                qssum = qssum + fitresults[index,1]
                qscont = qscont+1
            if np.isfinite(fitresults[index,2]):
                qdsum = qdsum + fitresults[index,2]
                qdcont = qdcont+1
    
    # qsavg means: for current q-value, sum mode, average of decay constants
    if qscont > 0:
        qsavg = qssum/qscont
    else:
        qsavg = np.nan
    
    if qdcont > 0:
        qdavg = qdsum/qdcont
    else:
        qdavg = np.nan
    
    # second walkthrough: Sum quadratic deviations from average
    # qserr2 means: For current q, sum of (lam-qsavg)**2 for all contributing lam
    qserr2=0.0
    qderr2=0.0
    for i in range(snapgrid_num):
        qx = dq*(i-snapgrid_num//2)
        for j in range(snapgrid_num):
            qy = dq*(j-snapgrid_num//2)
            q2 = qx*qx + qy*qy
            if not all([(qi-0.5)**2 *qbin**2 < q2, (qi+0.5)**2 *qbin**2 >= q2]):
                continue
            index = i*snapgrid_num+j
            
            if np.isfinite(fitresults[index,1]):
                qserr2 = qserr2 + (fitresults[index,1]-qsavg)**2
            if np.isfinite(fitresults[index,2]):
                qderr2 = qderr2 + (fitresults[index,2]-qdavg)**2
    
    # calculate actual variance: (sum of squares of deviations from avg) / sqrt(cont*(cont-1))
    if qscont <= 1:
        fitresults_binned[qi,1] = np.nan
        fitresults_binned[qi,2] = np.nan
    else:
        fitresults_binned[qi,1] = qsavg
        fitresults_binned[qi,2] = qserr2 / np.sqrt(qscont*(qscont-1))
    
    if qdcont <= 1:
        fitresults_binned[qi,3] = np.nan
        fitresults_binned[qi,4] = np.nan
    else:
        fitresults_binned[qi,3] = qdavg
        fitresults_binned[qi,4] = qderr2 / np.sqrt(qdcont*(qdcont-1))

np.savetxt("fitresults_binned.txt", fitresults_binned, delimiter="\t", header="Results of exponential fits to sum (rho1+rho2)(k) and difference (rho1-rho2)(k) modes.\nSimilar (abs values of) q's have been binned together. \nFormat: q lam_s errorbar lam_d errorbar, where lam is the decay/growth constant of exponential")
