"""Script file to read gridded densities from 
files localDens/rho<t>.txt, perform a 2-dim FFT,
bin the data to become 1-dimensional (|k| instead
of vector k) and write rho1,rho2,rho1+rho2,rho1-rho2
to file rhoIso<t>.txt.

For lowest <t>, the unbinned FTed data are written to file
rhok_start.txt.
"""

# ParProgress2tmp.py should be executed first, so that
# the file "tmp" contains obs_anzahl,L,lam,frames_per_tJ
with open("tmp", "r") as f:
    line = f.readline()
    words = line.split("\t")
    obs_anzahl = int(words[0])
    L = float(words[1])
    snapgrid_num = int(words[4])

import numpy as np

dq_bin = 0.1
dq = 2*np.pi/L
q_numbins = int(dq*snapgrid_num / dq_bin) +1

def readGrid(filename):
    """Take filename. Open file, read gridded
    densities. Return two 2-dim np-arrays. Does not
    perform consistency checks on the 
    coordinates/times in the first columns."""

    rho1 = np.zeros([snapgrid_num,snapgrid_num])
    rho2 = np.zeros([snapgrid_num,snapgrid_num])

    line = 0
    try:
        with open(filename, "r") as f:
            # skip first 3 lines
            f.readline()
            f.readline()
            f.readline()
            
            for i in range(snapgrid_num):
                for j in range(snapgrid_num):
                    words = f.readline().split()
                    r1 = float(words[3])
                    r2 = float(words[4])
                    rho1[i,j] = r1
                    rho2[i,j] = r2
                    line = line+1
    except IndexError: 
        print ("IndexError in Zeile ", line, words)
        
    return rho1,rho2

from numpy.fft import fft2

def ft_and_shift(arr):
    """Take a 2-dimensional np-array, expects snapgrid_num entries per dim.
    Perform a 2-dim FFT and shift (low q to middle of the array)."""
    
    arrq = fft2(arr)
    arrq = np.absolute(arrq)
    ret = np.zeros_like(arr)
    half = snapgrid_num // 2
    
    ret[:half ,:half ] = arrq[ half:, half:]
    ret[:half , half:] = arrq[ half:,:half ]
    ret[ half:,:half ] = arrq[:half , half:]
    ret[ half:, half:] = arrq[:half ,:half ]
    
    return ret

def rhok_to_iso(arr):
    """Takes a (presumable FT'ed) array. Averages
    over all q=|q|, performs binning with a binning-constant dq_bin.
    Returns a 1-dim array"""

    rhok_iso = np.zeros(q_numbins)
    normalize = np.zeros(q_numbins, dtype=int)
    
    for(k,l),rho in np.ndenumerate(arr):
        qx = dq*(k-snapgrid_num//2)
        qy = dq*(l-snapgrid_num//2)
        qabs2 = qx**2 + qy**2
        bin = int(np.sqrt(qabs2//dq_bin**2))
        if bin<0 or bin>= q_numbins:
            continue
        rhok_iso[bin] += rho
        normalize[bin] += 1
    
    return rhok_iso/normalize
#from calc_rhoIso import *



for t in range(1,obs_anzahl):
    rho1,rho2 = readGrid("localDens/rho{}.txt".format(t))
    
    rho1q = ft_and_shift(rho1)
    rho2q = ft_and_shift(rho2)

    rho1q_iso = rhok_to_iso(rho1q)
    rho2q_iso = rhok_to_iso(rho2q)

    outfile = open("rhok_iso/t{}.txt".format(t), "w")
    outfile.write("# Format: q[sigma11] TAB (rho1+rho2)(q) TAB (rho1-rho2)(q) TAB rho1(q) TAB rho2(q) \n\n")

    for i in range(q_numbins):
        q = dq*i
        r1 = rho1q_iso[i]
        r2 = rho2q_iso[i]
        r  = r1+r2
        dr = r1-r2
        outfile.write("{}\t{}\t{}\t{}\t{}\n".format(q, r, dr, r1, r2))

    outfile.close()

if False:    
    # lowest t: Write FTed data to file rhok_start.py
    rho1,rho2 = readGrid("localDens/rho1.txt")

    #rho1q = fft2(rho1)
    #rho2q = fft2(rho2)

    #rho1q = np.absolute(rho1q)
    #rho2q = np.absolute(rho2q)

    rho1q = ft_and_shift(rho1)
    rho2q = ft_and_shift(rho2)

    outfile = open("rhok_test.txt", "w")
    outfile.write("# Format: qx[sigma11] TAB qy[sigma11] TAB rho1(q) TAB rho2(q)\n\n")

    for i in range(snapgrid_num):
        qx = dq*(i-snapgrid_num/2)
        for j in range(snapgrid_num):
            qy = dq*(j-snapgrid_num/2)
            outfile.write("{} \t {} \t {} \t {}\n".format(qx, qy, rho1q[i,j], rho2q[i,j]))
    outfile.close()




### Treatment of the noiseless test density: Gauss function with widths 0.1*L in x-direction and 0.15*L in y-direction. 
# Density in file test_noiseless.txt, Format 0 x y rho1 rho2   (here rho1=rho2)
if False:
    rho1,rho2 = readGrid("localDens/test_noiseless.txt")
    rho1q = ft_and_shift(rho1)
    rho1q = np.absolute(rho1q) # absolute value of complex array. Correct to ignore the phase?

    outfile=open("FTnoiseless.txt", "w")
    outfile.write("# Format: qx[sigma11] qy[sigma11] rho1(q) rho2(q)\n\n")

    for i in range(snapgrid_num):
        qx = dq*(i-snapgrid_num/2)
        for j in range(snapgrid_num):
            qy = dq*(j-snapgrid_num/2)
            outfile.write("{} \t {} \t {} \t {} \n".format(qx, qy, rho1q[i,j], rho1q[i,j]))

    outfile.close()
