# gnuplot script file to draw particle positions, 2dim density, and FTed 2dim density.
# The goal is to validate the fourier-transformation, so we can reliably use the rho(k)

reset

# Strategy: Multiplot with several plots. 
# 1: Box with particle positions
# 2: Box with gridded density, colorbar.
# 2b The same along x axis, with "prediction"
# 3: Box with FT'ed gridded density, colorbar.
# 3b The same along kx axis, with "prediction"
# (4 FT'ed density, directions averaged out)


# read parameters from tmp-file. That file is written by python3 ParProgress2tmp.py
set table "/dev/null"
plot 'tmp' u (L=$2, lam=$3, frames_per_tJ=$4, FFT_num=$5):(1)
unset table

# read other parameters from tmp_denstest. That file is also written by ParProgress2tmp.py
set table "/dev/null"
plot 'tmp_denstest' u (gaussdens_x=$1,gaussdens_y=$2):(1)
unset table


# double-ended arrow indicating lambda in the snapshot
set style arrow 1 heads size 2,90 front lw 4

# scale particle radius with box size
particle_radius = 150./L

# colloids (dots) in the snapshot
set style line 1 lc rgb "#ad3000" pt 7 ps particle_radius lt 1 lw 2
set style line 2 lc rgb "#0060ad" pt 7 ps particle_radius lt 1 lw 2

# lines when a function is plotted
set style line 3 lc rgb "#225ea8" lw 3
set style line 4 lc rgb "#4daf4a" lw 3

### more linestyles here


unset xlabel
unset ylabel
unset key

# max expected density in particles/(sigma**2). Manual fine-tuning.
rmax = 0.2

# normalized total rho. The interval [0 : rmax] is mapped to [0 : 1], with clipping
ntr(rho1,rho2) = rho1+rho2 > rmax ? 1 : (rho1+rho2)/rmax

# weaker version of RGBalpha: Use alpha values [ alpmin : alpmax ] instead of 0..255. At rmax or higher, this gives alpmax.
alpmin=0; alpmax=200
#a(rho1,rho2) = int(alpmin + (alpmax-alpmin)*ntr(rho1,rho2))


# max density is #e31a1c. Red part is 0xe3=227, green part is 0x1a=26, blue part is 0x1c=28
# min density is #ffffff, white. 255 in each byte.
# 
# ra,ga,ba are red,green,blue with "pseudo-alpha": shifted towards white to become pale as alpha is low
ra(rho1,rho2) = 255 - (255-227)*ntr(rho1,rho2)
ga(rho1,rho2) = 255 - (255- 26)*ntr(rho1,rho2)
ba(rho1,rho2) = 255 - (255- 28)*ntr(rho1,rho2)


set term pngcairo font ",28" size 2000,1000

set xrange [0:L]
set yrange [0:L]

unset xtics
unset ytics

# files with particle coordinates. Format t x y
posdatei1 = "pos1_1.txt"
posdatei2 = "pos2_1.txt"

# file with gridded densities rho1,rho2 in real space. Format t x y rho1 rho2
rhoxdatei = "localDens/test_noiseless.txt"

outfile= "rhok_testresult_noiseless.png"
set output outfile

set multiplot

# first plot: Particle positions. Top left
set origin 0,0.5
set size 0.33,0.5

#plot posdatei1 u 2:3 ls 1, posdatei2 u 2:3 ls 1
set label 7 "no positions in\nnoiseless test" at 0.1*L,0.3*L
p NaN

# second plot: Particle density in real space. Top center. Draw lines x=L/2, y=L/2 for cuts
set origin 0.33,0.5
set size 0.33,0.5

# linestyle for these lines
set style arrow 21 ls 3 nohead
set style arrow 22 lw 3 lc rgb "black" nohead

set arrow 21 from L/2.,0.1*L  to L/2.,0.9*L  arrowstyle 21 front
set arrow 22 from 0.1*L,L/2.  to 0.9*L,L/2.  arrowstyle 22 front
#set arrow 23 from 0.1*L,0.1*L to 0.9*L,0.9*L arrowstyle 21 front
#set arrow 24 from 0.1*L,0.9*L to 0.9*L,0.1*L arrowstyle 21 front

plot rhoxdatei using 2:3:(ra($4,$5)):(ga($4,$5)):(ba($4,$5)) with rgbimage

#unset arrow 21; unset arrow 22; unset arrow 23

# 2b plot: Particle density in real space, projections along x=L/2, y=L/2, x=y. And expectation. Bottom row, next to the second plot
set origin 0.33,0
set size 0.33,0.5

set xr [0 : L]
set yr [0 : 1.1]

set ytics
set xtics
set xlabel "x"
#set ylabel "{/Symbol r}(x)" rotate

# almost: Returns true if the arguments differ less than 0.01*L
almost(a,b) = abs(a-b)<0.01*L

# linestyle for function that depicts expected density
# #225ea8
set style line 25 lc rgb "#a0aed3" lw 2
set style line 26 lc rgb "#080808" lw 2

plot \
exp(-0.5*(x-L/2.)**2/(gaussdens_y**2)) w l ls 25,\
exp(-0.5*(x-L/2.)**2/(gaussdens_x**2)) w l ls 26,\
rhoxdatei using (almost(L/2., $2) ? $3          : NaN):4 w    p      ls 3 ,\
rhoxdatei using (almost(L/2., $3) ? $2          : NaN):4 with points ls 3 lc rgb "black"#,\
#rhoxdatei using (almost($2,$3)   ? ($2-L/2.)*sqrt(2.)+L/2. : NaN):$4 w p ls 3 ,\
#rhoxdatei using (almost($2,L-$3) ? ($2-L/2.)*sqrt(2.)+L/2. : NaN):$4 w p ls 3


# third plot: Particle density in k-space (absolute value of complex rho(k)). Top right. 
set origin 0.66,0.5
set size 0.33,0.5

qmax = 2.*pi/L * FFT_num

set arrow 21 from -0.4*qmax,0 to 0.4*qmax,0 arrowstyle 22 front
set arrow 22 from 0,-0.4*qmax to 0,0.4*qmax arrowstyle 21 front

unset xlabel
#set ylabel "k_y{/Symbol s}" rotate


set xr [-qmax/2 : qmax/2]
set yr [-qmax/2 : qmax/2]

rhokdatei = "FTnoiseless.txt"
rmax = 80

plot rhokdatei using 1:2:(ra($3,$4)):(ga($3,$4)):(ba($3,$4)) with rgbimage

unset arrow 21; unset arrow 22

# 3b plot: Particle density in k-space, projections along x=L/2, y=L/2, x=y. Bottom right.
set origin 0.66,0
set size 0.33,0.5

set xr [-0.6*qmax : 0.6*qmax]
set yr [0 : *]

set xlabel "k"
#set ylabel "{/Symbol r}(k)" rotate

# redefine the "almost" function to have a reasonable tolerance
dq = 2.*pi/L
almost(a,b) = abs(a-b)<0.3*dq

# read parameters of the test functions from file tmp_denstest
# TODO

print gaussdens_x
print gaussdens_y
set samples 10000

print 2*pi*gaussdens_x * gaussdens_y

# ls 25 ist blau, ls 26 ist grau/schwarz

plot \
(FFT_num/L)**2 * 2*pi*gaussdens_x*gaussdens_y * exp(-0.5 * (x*gaussdens_y)**2) ls 25,\
(FFT_num/L)**2 * 2*pi*gaussdens_x*gaussdens_y * exp(-0.5 * (x*gaussdens_x)**2) ls 26,\
rhokdatei using (almost(0., $2) ? $1 : NaN):3 with points ls 3 lc rgb "black",\
rhokdatei using (almost(0., $1) ? $2 : NaN):3 with points ls 3 #,\
#rhokdatei using (almost($1,$2)      ? ($1-qmax/2.)*sqrt(2.)+qmax/2. : NaN) :($3+$4) w p ls 4 ,\
#rhokdatei using (almost($1,qmax-$2) ? ($1-qmax/2.)*sqrt(2.)+qmax/2. : NaN) :($3+$4) w p ls 4 


# final: add a few arrows and labels to indicate procedure
set origin 0,0
set size 1,1
set xr [0:1]; set yr [0:1]
unset border
unset xtics; unset ytics
unset xlabel; unset ylabel

set style arrow 11 lc rgb "black" lw 3 head 

set arrow from 0.3,0.75 to 0.35,0.75 arrowstyle 11 front
set label "grid" at 0.3,0.78
set arrow from 0.64,0.75 to 0.7,0.75 arrowstyle 11 front
set label "FFT" at 0.64,0.78
set arrow from 0.5,0.55 to 0.5,0.45 arrowstyle 11 front
set label "cuts" at 0.51,0.5

plot 0 lc rgb "white"

unset multiplot
