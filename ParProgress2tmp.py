"""script to read parameters (box size etc) from 
file parameter.h and progress from file out.txt

Writes everything the gnuplot-file plot_animation.plt needs 
to file tmp. The file is not commented or anything and should 
be deleted by the Makefile afterwards.
"""

import numpy as np


# initialize variables to -1, so we can check whether all were found
L,obs_dt,N1,N2,f2_epsgam,snapgrid_num = -1,-1,-1,-1,-1,-1

with open("parameter.h" ,"r") as fileP:
    for linenr,line in enumerate(fileP):
        
        # skip empty lines
        if 1 == len(line):
            continue
        
        # find lines starting with "const"
        if "const" == line[0:5]:
            halves = line.split("=")
            left,right = halves[0],halves[1]
            var = left.split()[2]
            val = right.split(";")[0]
            if "L" == var:
                L = float(val)
            elif "lambda_kapillar" == var:
                lam = float(val)
            elif "obs_dt" == var:
                obs_dt = float(val)
            elif "N1" == var:
                N1 = int(val)
            elif "N2" == var:
                N2 = int(val)
            elif "kapillar_vorfaktor" == var:
                f2_epsgam = float(val)
            elif "snapgrid_crude" == var:
                snapgrid_crude = float(val)
                snapgrid_num = 2*(int) (L/(2*snapgrid_crude) + 0.5)

    if any(np.array([L,obs_dt,N1,N2,f2_epsgam,snapgrid_num]) < 0.0):
        raise ValueError("at least one of the needed variables wasn't found \n"
                         "in parameter.h. N1={}, N2={}, L={}, lam={}, obs_dt={}, kap_vorfaktor={}, snapgrid_num={}".format(
                         N1,N2,L,lam,obs_dt,f2_epsgam,snapgrid_num))

# if this point is reached, reading from parameter.h was successful.
rho = (N1+N2)/(L*L)
tJ = 1./(rho * f2_epsgam)
print(rho, tJ)
frames_per_tJ = tJ / obs_dt

# number of checkpoints already reached
progress = -1

with open("out.txt", "r") as fileO:
    for linenr,line in enumerate(fileO):
        if 1 == len(line):
            continue
        words = line.split()
        if len(words) < 3:
            #print(words)
            continue
        if "recording" == words[2]:
            progress = int(words[4])

# write gathered information to file
### TEMPORARY: WRITE 100 INSTEAD OF THE ACTUAL PROGRESS
outfile = open("tmp", "w")
outfile.write("{}\t{}\t{}\t{}\t{}".format(progress, L, lam, frames_per_tJ,snapgrid_num))
#outfile.write("{}\t{}\t{}\t{}\t{}".format(200, L, lam, frames_per_tJ,snapgrid_num))