
"""Reads File(s) pos1_1.txt and perhaps pos2_1.txt
Calculates rho(k), rho1+rho2(k), rho1-rho2(k) and
compares to exponential prediction with Jeans-Time."""


import numpy as np

from readParameters import *
(N1,N2,L,lam,T,obs_dt,obs_anzahl) = readParameters()

FFT_samples = 1024 #must be even
dq_bin = 0.05
def readPos(N, f):
    """reads N Positions from file f. Returns array [N][dim]"""
    ret = np.zeros((N, 2))
    
    for i in range(N):
        
        line = f.readline()
        words = line.split('\t')
        
        # skip empty and comment lines
        while words[0][0] == '\n' or words[0][0] == "#":
            line = f.readline()
            words = line.split('\t')
        
        ret[i] = [ float(words[1]), float(words[2]) ]
        
    return ret    


def gridTSC(arr, L, FFT_samples):
    """turns an array of positions into an gridded array of density. Returns array of density."""
    
    cellwidth = L / FFT_samples
    drho = 1/cellwidth**2 # contribution of one particle to the density
    
    ret = np.zeros((FFT_samples, FFT_samples))
    
    for [x,y] in arr:
        
        if any(np.array([x,y])<0.) or any(np.array([x,y])>L):
            raise ValueError("Position ({},{}) not within box [0,{}]".format(x,y,L))
        
        # k,l are indices of the cell containing the center of the particle
        k = int(x // cellwidth)
        l = int(y // cellwidth)
        
        # kleft and kright are the neighboring cells. 
        # And lleft, lright for y-direction. "left"/"right" is slightly misleading.
        kleft  = (k-1+FFT_samples)%FFT_samples
        kright = (k+1)%FFT_samples
        lleft  = (l-1+FFT_samples)%FFT_samples
        lright = (l+1)%FFT_samples
        
        # x-direction. x_in is the fraction of the cloud in the center cell k
        d = x - (k+0.5)*cellwidth
        x_in = 0.75 - d**2/cellwidth**2
        # fraction in kleft
        tmp = abs((d+cellwidth)/cellwidth)
        x_left = 0.5*(1.5-tmp)**2
        # fraction in kright
        tmp = abs((d-cellwidth)/cellwidth)
        x_right = 0.5*(1.5-tmp)**2
        
        # same for y-direction
        d = y - (l+0.5)*cellwidth
        y_in = 0.75 - d**2/cellwidth**2
        tmp = abs((d+cellwidth)/cellwidth)
        y_left = 0.5*(1.5-tmp)**2
        tmp = abs((d-cellwidth)/cellwidth)
        y_right = 0.5*(1.5-tmp)**2
        
        ret[kleft, lleft ] += x_left *y_left  *drho
        ret[kleft, l     ] += x_left *y_in    *drho
        ret[kleft, lright] += x_left *y_right *drho
        
        ret[k,     lleft ] += x_in   *y_left  *drho
        ret[k,     l     ] += x_in   *y_in    *drho
        ret[k,     lright] += x_in   *y_right *drho
        
        ret[kright,lleft ] += x_right*y_left  *drho
        ret[kright,l     ] += x_right*y_in    *drho
        ret[kright,lright] += x_right*y_right *drho
    return ret
def calc_rhok(pos, L, FFT_samples=1024):
    """Takes array of positions as produced be readPos().
    Performs a 2-dim FFT, return rho(k) for all 2D vectors k.
    Result has alternating signs, but low k are at the start of the array.
    """
    grid = gridTSC(pos, L, FFT_samples)

    # 2-dim FFT of real-valued arrays
    from numpy.fft import fft2 as fft

    rhok = fft(grid)



    ## alternate signs
    #alt_tile = np.array([[1, -1], [-1, 1]])
    #altsigns = np.tile(alt_tile, (FFT_samples//2, FFT_samples//2))
    #rhok *= altsigns

    # normalization factor for FFT
    rhok *= (L / FFT_samples)**2

    #print(np.shape(rhok))

    # shift FT'ed array so that k=0 is at the start
    rhok_shifted = np.zeros_like(rhok)
    rhok_shifted[FFT_samples//2: , FFT_samples//2: ] = rhok[0:FFT_samples//2, 0:FFT_samples//2]
    rhok_shifted[FFT_samples//2: , 0:FFT_samples//2] = rhok[0:FFT_samples//2, FFT_samples//2: ]
    rhok_shifted[0:FFT_samples//2, FFT_samples//2: ] = rhok[FFT_samples//2: , 0:FFT_samples//2]
    rhok_shifted[0:FFT_samples//2, 0:FFT_samples//2] = rhok[FFT_samples//2: , FFT_samples//2: ]
    
    return rhok_shifted
def rhok_to_iso(rhok):
    """ Takes a (presumably fourier-transformed) array rhok(q),
    Expects a (complex) 2-dim input array, uses the absolute value of input.
    Averages over all q=|q|, performs binning with a binning constant dq_bin.
    Returns an array rho(k) for k from 0 to qmax=FFT_samples*2pi/L.
    Uses absolute values of the complex rho(k).
    """
    dq = 2*np.pi/L
    q_numbins = int( dq*FFT_samples/dq_bin ) +1

    rhok = np.absolute(rhok)
    
    # Strategy: Walk over all rhok_shifted and add them to
    # the corresponding bin. Count number of contributions 
    # in "normalize". At the end, divide each bin by the
    # corresponding number in "normalize".
    
    rhok_final = np.zeros(q_numbins)
    normalize = np.zeros_like(rhok_final)

    for (k,l),rho in np.ndenumerate(rhok):
        qabs2 = dq**2*(k**2 + l**2)
        bin = int(np.sqrt(qabs2 // dq_bin**2))
        if bin < 0 or bin >= q_numbins:
            continue
        rhok_final[bin] += rho
        normalize[bin] += 1

    rhok_final = rhok_final/normalize
    return rhok_final

def draw_rhok(rho):
    """Take 2-dim array rho(k). Flatten it by averaging over vectors q
    with the (almost) same length (binning constant dq_bin). Draw the resulting rho(q)."""
    
    dq = 2*np.pi/L
    q_numbins = int( dq*FFT_samples/dq_bin ) +1
    
    import matplotlib.pyplot as plt
    q_final = np.linspace(0, FFT_samples*2*np.pi/L, num=q_numbins)
    rhok_final = rhok_to_iso(rho)
    plt.plot(q_final[1:], rhok_final[1:])
    #plt.ylim((-500, 500))
    plt.show()