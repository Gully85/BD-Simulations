# gnuplot script to plot timescales for density modes (as a function of wavenumber k)
# write resulting plot to files timescales_of_k_bigbox.pdf and timescales_of_k_dense.pdf

reset

# only positive k are meaningful. Everything is a function of k**2 anyway.
set xr [0  : 1.5]
set yr [-0.25 : 0.05]

set border 1+2+8
set xtics nomirror
set ytics 0.1 nomirror

set arrow from 0,0 to 1.5,0 nohead
set label 51 "decay (stable)"    at 1.1,-0.005 font ",12"
set label 52 "growth (unstable)"  at 1.1,0.005  font ",12"
set label 53 "slow  fast" at 1.53,0.005 rotate by 90 
set arrow 53 from 1.57,0.005 to 1.57,0.04 head lw 2
set label 54 "fast          slow" at 1.53,-0.13 rotate by 90
set arrow 54 from 1.57,-0.005 to 1.57,-0.15 head lw 2

set lmargin 13
set xlabel "k[1/{/Symbol s}]"
set ylabel "-{/Symbol G}k^2{/Symbol l}[1/T_{BD}]" norotate offset -1,6

set samples 2000

unset key

# bigbox runs have different Teff. 0.03, 0.04, 0.07, 0.1, 0.15, 0.22, 0.3, 0.35
# these are the kT/eps for the runs
kteps_003 = 10.41666
kteps_004 = 13.88888
kteps_007 = 24.30555
kteps_010 = 34.72222
kteps_015 = 52.08333
kteps_022 = 76.38888
kteps_030 =104.16666
kteps_035 =121.52777

# linestyles. 11-19 are for total density mode (just -k^2 decay). 
# Taken from ColorBrewer sequential 9-class Yellow-Green
set style line 11 lw 3 lt 1 lc rgb "#ffffe5"
set style line 12 lw 3 lt 1 lc rgb "#f7fcb9"
set style line 13 lw 3 lt 1 lc rgb "#d9f0a3"
set style line 14 lw 3 lt 1 lc rgb "#addd8e"
set style line 15 lw 3 lt 1 lc rgb "#78c679"
set style line 16 lw 3 lt 1 lc rgb "#41ab5d"
set style line 17 lw 3 lt 1 lc rgb "#238443"
set style line 18 lw 3 lt 1 lc rgb "#006837"
set style line 19 lw 3 lt 1 lc rgb "#004529"

# and 21-29 are for demixing mode (more complicated function lambda(k))
# Taken from ColorBrewer sequential 9-class Red-Purple
set style line 21 lw 3 lt 1 lc rgb "#fff7f3"
set style line 22 lw 3 lt 1 lc rgb "#fde0dd"
set style line 23 lw 3 lt 1 lc rgb "#fcc5c0"
set style line 24 lw 3 lt 1 lc rgb "#fa9fb5"
set style line 25 lw 3 lt 1 lc rgb "#f768a1"
set style line 26 lw 3 lt 1 lc rgb "#dd3497"
set style line 27 lw 3 lt 1 lc rgb "#ae017e"
set style line 28 lw 3 lt 1 lc rgb "#7a0177"
set style line 29 lw 3 lt 1 lc rgb "#49006a"

# Timescale prediction in Linear Stability Analysis was: T_LS = 1./(Gamma k**2 lambda), where 
# lambda=kT-rho*W(k) for density mode and lambda=kT(1+rho*H(rho)) for demixing mode.
# Meaningful is the ratio between this term and the BD timescale T_BD = sigma**2/(eps Gamma)
# Plotted is T_LS/T_BD = 1./(k sigma)**2 1./(kT/eps) 1./(1+rho*H(rho)) [density mode]
# and resp. T_LS/T_BD  = 1./(k sigma)**2 eps/(kT-rho W(k)) [demixing mode]

# capillary length lambda. Not to be confused with inverse timescale 
# lambda, eigenvalue of the matrix in linear stability analysis. 
# Capillary length is called "lam". In units of sigma.
lam_bigbox = 25.

# H(rho) is constant for all the bigbox runs, and also a (different) constant for all the dense runs.
rho_bigbox = 0.1
rmax = 4./pi
kap_factor = 10./9
H(rho) = -1./(rmax-rho) - rho/(rmax-rho)**2 + 2*rmax**2/(rmax-rho)**3
rhoH_bigbox = rho_bigbox * H(rho_bigbox)
rho_dense = 0.4
rhoH_dense  = rho_dense  * H(rho_dense)

# W(k) is the Greens Function for the capillary interaction. With units, it is W(k) = f^2/gamma 1./(k**2+lam_bigbox**-2)
W(k) = kap_factor * 1./(k**2+lam_bigbox**-2)

# defined without the 1./(kT/eps)
lambda_densMode(k) = -k**2 * (1+rhoH_bigbox)

lambda_demixMode(k, kTeps)= -k**2 * (kTeps - rho_bigbox*W(k))


# label temperature and dens-Mode vs demix-Mode
set label 11 "density mode {/Symbol r}_1+{/Symbol r}_2" at 0.8,-0.08 rotate by -40 tc rgb "#006837"
set label 12 "T_{eff}=0.03" at 1.25,-0.18 rotate by -45 tc rgb "#006837" font ",16"
set label 13 "T_{eff}=0.15" at 1.28,-0.03 rotate by -10 tc rgb "#addd8e" front font ",14"
set arrow 11 from 0.6,-0.04 to 0.65,0.01 lc rgb "black" lw 3 head front
set label 14 "T_{eff}" at 0.65,0.01

set label 21 "demixing mode {/Symbol r}_1-{/Symbol r}_2" at 0.35,-0.23 tc rgb "#7a0177"
set arrow 21 from 0.18,-0.2 to 0.02,-0.21 lc rgb "black" lw 3 head front
set label 24 "T_{eff}" at 0.01,-0.22 front


# for zoom-in: Highlight areas. The rectangle [oldleft, oldbot] x [oldright, oldtop] is plotted again in [newleft, newbot] x [newright, newtop]
oldleft  = 0
oldbot   =-0.005
oldright = 0.11
oldtop   = 0.045
set style line 2 lw 2 lc rgb "black"
set style arrow 2 ls 2 nohead
set arrow 61 from oldleft,oldbot to oldleft,oldtop as 2
set arrow 62 from oldleft,oldtop to oldright,oldtop as 2
set arrow 63 from oldright,oldtop to oldright,oldbot as 2
set arrow 64 from oldright,oldbot to oldleft,oldbot as 2

newleft = 0.24
newbot = -0.2
newright = 0.8
newtop = -0.1
set arrow 71 from newleft,newbot to newleft,newtop as 2
set arrow 72 from newleft,newtop to newright,newtop as 2
set arrow 73 from newright,newtop to newright,newbot as 2
set arrow 74 from newright,newbot to newleft,newbot as 2

set arrow 75 from oldright,oldtop to newright,newtop as 2 front
set arrow 76 from oldleft,oldbot  to newleft,newbot as 2 front


outfile = "timescales_of_k_bigbox.eps"
set term postscript enhanced eps color solid font ",20"
set output outfile

set multiplot
plot \
lambda_densMode(x) * 1./kteps_035 ls 11 ,\
lambda_densMode(x) * 1./kteps_030 ls 12 ,\
lambda_densMode(x) * 1./kteps_022 ls 13 ,\
lambda_densMode(x) * 1./kteps_015 ls 14 ,\
lambda_densMode(x) * 1./kteps_010 ls 15 ,\
lambda_densMode(x) * 1./kteps_007 ls 16 ,\
lambda_densMode(x) * 1./kteps_004 ls 17 ,\
lambda_densMode(x) * 1./kteps_003 ls 18 ,\
\
lambda_demixMode(x, kteps_035)    ls 21 ,\
lambda_demixMode(x, kteps_030)    ls 22 ,\
lambda_demixMode(x, kteps_022)    ls 23 ,\
lambda_demixMode(x, kteps_015)    ls 24 ,\
lambda_demixMode(x, kteps_010)    ls 25 ,\
lambda_demixMode(x, kteps_007)    ls 26 ,\
lambda_demixMode(x, kteps_004)    ls 27 ,\
lambda_demixMode(x, kteps_003)    ls 28

set origin 0.335,0.265
set size 0.272,0.275
set xr [oldleft : oldright]
set yr [oldbot  : oldtop]
set bmargin 0; set tmargin 0; set lmargin 0; set rmargin 0
unset border
unset ytics; set xtics out font ",14" offset 0,0.5; set ytics (0.04) font ",14" nomirror offset 0.8,0
set x2tics 2.*pi/350. out font ",14" format ""
unset ylabel; unset xlabel
unset label; unset arrow
set label at 0.005,1.05*oldtop "dk" font ",12"
set arrow from oldleft,0 to oldright,0 nohead
set label "T_{eff}=0.03" at 0.08,0.03 tc rgb "#7a0177" rotate by -60 font ",14"
set label "T_{eff}=0.15" at 0.015,0.004 tc rgb "#fa9fb5" rotate by -15 font ",8"
replot

system "epstopdf ".outfile
system "rm ".outfile





