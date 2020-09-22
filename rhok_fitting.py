
"""Provides functions for fitting an exponential decay (or rise) 
to FFT'ed density profiles. 
"""

import rhok_functions as rf

import numpy as np

# read parameters from parameter.h
from readParameters import *
(N1,N2,L,lam,T,obs_dt,obs_anzahl) = readParameters()

def read_and_FFT(filename, N, numt, FFT_samples=1024):
    """Read numt blocks of N positions from file. 
    Return an array rho[qx, qy, t].
    qx and qy are 0...FFT_samples-1, t is 0...numt-1"""
    
    # store all the FT'ed rho(q,t). First two indices are q/dq, 
    # third index is t/(dt*time_bin)
    rhoqt = np.zeros((FFT_samples, FFT_samples, numt), dtype='complex128')

    with open(filename, "r") as file:
        for i in range(numt):
            pos = rf.readPos(N1, file)
            rhok = rf.calc_rhok(pos, L, FFT_samples=FFT_samples)
            #rhok = np.absolute(rhok)
            rhoqt[:,:,i] += rhok
    
    return rhoqt

def multifit(rhoqt, silent=False):
    """Takes an array rhoqt like the one read_and_FFT() returns.
    Performs an exponential fit over every qx,qy.
    Returns arrays ampl, decay, weights and allq
    where allq[index] = [qx,qy] for ampl[index]
    
    If the fit failed, there is no entry in the arrays."""
    
    FFT_samples=rhoqt.shape[0]
    numt = len(rhoqt[0,0,:])
    
    # Store fitting parameters ampl,decay,weights for every q.
    # These are arrays of a single index, matching the qx,qy 
    # in the parallel array allq
    Ampl_arr   = np.zeros(FFT_samples**2, dtype=float)
    Decay_arr  = np.zeros(FFT_samples**2, dtype=float)
    Weights_arr= np.zeros(FFT_samples**2, dtype=float)
    allq       = np.zeros((FFT_samples**2,2), dtype=int)
    failedfits =-np.zeros(FFT_samples**2, dtype=int) # -1's are placeholders
                                                     # so the array doesn't 
    numfails = 0                                     # need to grow in-place
    numsuccess=0

    ts = np.linspace(0, numt, num=numt)
    # function to fit
    exp_func = lambda t,ampl,decay,offset: ampl*np.exp(-decay*t)+offset
    # jacobian: derivatives w.r.t. ampl,decay,offset
    jac_func = lambda t,ampl,decay,offset: np.transpose([np.exp(-decay*t),
                                                         -ampl*t*np.exp(-decay*t),
                                                         np.ones_like(t)
                                                        ])
    from scipy.optimize import curve_fit
    import warnings
    from scipy.optimize import OptimizeWarning
    warnings.filterwarnings("error")
    
    # for each qx,qy perform the fit
    for index in range(FFT_samples**2):
        
        #index = qx*FFT_samples + qy
        qx = index // FFT_samples
        qy = index %  FFT_samples
        vals = rhoqt[qx,qy]
        # curve_fit can throw warnings/errors if the fit fails. In this case, 
        # skip this [qx,qy]. Do not create an entry to the arrays, count fails
        # and afterwards shrink the arrays
        try:
            (ampl,decay,offset),pcov = curve_fit(exp_func, ts, vals, jac=jac_func, maxfev=10000)
            Ampl_arr[numsuccess]    = ampl
            Decay_arr[numsuccess]   = decay
            Weights_arr[numsuccess] = vals[0]
            allq[numsuccess,:]      = [qx,qy]
            numsuccess += 1
        except (RuntimeWarning, OptimizeWarning, RuntimeError):
            numfails   += 1
    
    assert numsuccess + numfails == FFT_samples**2
    
    Ampl_arr  = Ampl_arr   [:numsuccess]
    Decay_arr = Decay_arr  [:numsuccess]
    Weights_arr=Weights_arr[:numsuccess]
    allq      = allq       [:numsuccess]
    
    print("While fitting, ", numfails, "out of", FFT_samples**2, " q didn't work.")
    
    return Ampl_arr,Decay_arr,Weights_arr,allq
def average_by_lengths(arr, allq, q_bin, weights=None):
    """Takes an array arr with one index and a parallel array allq
    which assigns qx,qy to each index.
    Performs an averaging over vectors qx,qy with similar
    length (binning constant q_bin).
    Weights can be provided, they must be the same shape as arr.
    Returns an array with only one index."""
    
    if weights is None:
        weights = np.ones_like(arr)
    elif weights.shape != arr.shape:
        raise ValueError("Mismatch between shapes of arr and weights!"
                         "Shape of weights should be" + arr.shape + ", but it is " + weights.shape)
    
    if allq.shape[0] != arr.shape[0]:
        raise ValueError("Mismatch between arr and allq. They should be same size.",
                         "allq.shape=",allq.shape," and arr.shape=",arr.shape)
    
    numq = arr.shape[0]
    
    # in principle, the highest q could be not in the last entry
    # (happens if the last few fits all fail). Ignore that chance
    # and potentially miss a few q's.
    lastqs = allq[-1] 
    maxqabs = np.linalg.norm(lastqs)
    num_bins = int(np.sqrt(2.)* (maxqabs // q_bin)) + 1
    q_bins = np.linspace(0, maxqabs, num=num_bins)
    
    data_binned = np.zeros(num_bins, dtype=float)
    normalize   = np.zeros_like(data_binned)
    
    import warnings
    for index in range(numq):
        qx,qy = allq[index]
        qabs2 = qx**2 + qy**2
        bin = int(qabs2 // q_bin)
        if bin < 0 or bin >= num_bins:
            #warnings.warn("While binning, the q-vector ({},{}) didn't fit in any bin".format(qx,qy))
            continue
        norm = weights[index]
        data_binned[bin] += norm*arr[index]
        normalize[bin]   += norm

    # check which bins have no entries, delete them from binned arrays
    skip_bin = (normalize==0)
    print(len(skip_bin), "bins were empty and are thus not returned")
    cut_bins = [i for i,v in enumerate(skip_bin) if v]
    
    data_binned = np.delete(data_binned, cut_bins)
    normalize   = np.delete(normalize,   cut_bins)
    q_bins      = np.deleta(q_bins,      cut_bins)
    
    # normalize and return
    data_binned /= normalize
    
    return q_bins,data_binned