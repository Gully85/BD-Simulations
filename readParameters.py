"""Read parameters (particle number, box size etc) from file parameter.h
Return them as a big tuple. Then use

from readParameters import *
(N1,N2,L,lam,T,obs_dt,obs_anzahl) = readParameters()

to get them.
"""

import numpy as np

def readParameters():
    N1,N2,L,lam,T,obs_dt,obs_anzahl = -1,-1,-1.,-1.,-1.,-1.,-1
    
    with open("parameter.h", "r") as file:
        for linenr,line in enumerate(file):
            #print(linenr,line[0:4])
            if len(line) == 1:
                continue
            if line[0:5] == "const":
                halves = line.split("=")
                left,right = halves[0],halves[1]
                var = left.split()[2]
                val = right.split(";")[0]
                #print(var, val)
                if var == "N1":
                    N1 = int(val)
                elif var == "N2":
                    N2 = int(val)
                elif var == "L":
                    L = float(val)
                elif var == "T":
                    T = float(val)
                elif var == "lambda_kapillar":
                    lam = float(val)
                elif var == "obs_dt":
                    obs_dt = float(val)
                elif var == "obs_anzahl":
                    # obs_anzahl is allowed to be composed of numbers
                    # in the form a*b*c. 
                    nums = val.split("*")
                    nums2 = [int(i) for i in nums]
                    #print(nums2)
                    obs_anzahl = np.prod(nums2)
                    #nums2= int(nums)
                    #obs_anzahl = np.prod(nums)
                    
                
    if any(np.array([N1,N2,L,lam,T,obs_dt,obs_anzahl]) < 0):
        print (np.array([N1,N2,L,lam,T,obs_dt,obs_anzahl]) < 0)
        raise ValueError("at least one of the needed variables wasn't found in parameter.h. \n"
                         "N1={}, N2={}, L={}, lambda={}, obs_dt={}, obs_anzahl={}".format(
                        N1,N2,L,lam, obs_dt, obs_anzahl))
    return N1,N2,L,lam,T,obs_dt,obs_anzahl
    