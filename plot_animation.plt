reset


set xlabel 'x'
set ylabel 'y'



# get parameters from tmp-file. That file is supposed to be written by python3 ParProgress2tmp.py
set terminal unknown
stats 'tmp' u (anzahl_frames=$1, L=$2, lam=$3, frames_per_tJ=$4):1 nooutput

## TEMP
#anzahl_frames = 10

#print anzahl_frames
#print L
#print lam
#print frames_per_tJ

# double-ended arrow indicating lambda in the snapshot
set style arrow 1 heads size 2,90 front lw 5

# scale particle radius with box size
particle_radius = 150./L

# colloids (dots) in the snapshot
set style line 1 lc rgb "#ad3000" pt 7 ps particle_radius lt 1 lw 2
set style line 2 lc rgb "#0060ad" pt 7 ps particle_radius lt 1 lw 2

# green line, inter-species g(r)
set style line 3 lc rgb "#30a000" pt 7 ps 1.5 lt 1 lw 2

# black line, ratio g11(r)/g12(r)
set style line 4 lc rgb 'black'   pt 7 ps 1.5 lt 1 lw 2

# violet line for rho1+rho2
set style line 11 lc rgb "dark-violet" lt 1 lw 2

# light-orange line for rho1-rho2
set style line 12 lc rgb "#fdbf6f" lt 1 lw 2

# red line for species 1
set style line 13 lc rgb "#d53e4f" lt 1 lw 2

# blue line for species 2
set style line 14 lc rgb "#3288bd" lt 1 lw 2


# set xlabel 'x/{/Symbol s}'
# set ylabel 'y/{/Symbol s}'

unset xlabel
unset ylabel

# principial max packing for 2dim hard disks. Only used in axis range/labels of histograms, not in the coloring procedure
rhomax = 2./sqrt(3.)

# density: invisible for rho1+rho2=0, full opacity for rho1+rho2 > maxdens*rho0
# color: #d53e4f for (rho1-rho2)/(rho1+rho2) < -maxdemix, #3288bd for (rho1-rho2)/(rho1+rho2) > maxdemix
maxdemix = 1.0
rho0 = 0.01
# which ranges of rho are expected? e.g. rhorange=5 means a range from rho0/5 to rho0*5. Must be greater than 1.
rhorange = 2

# deltarho, normalized to -1..1
dr(rho1,rho2) = rho1+rho2<0.0001 ? 0 : (rho1-rho2)/(rho1+rho2)
# normalized deltarho, -maxdemix..maxdemix is mapped to -1..+1
min(x,y) = x<y ? x : y
max(x,y) = x>y ? x : y
ndr(rho1,rho2) = max(-maxdemix, min(dr(rho1,rho2), maxdemix))/maxdemix

# normalized totalrho, rho0/2 .. maxdens*rho0 is mapped to 0..1
#ntr(rho1,rho2) = min(rho0*maxdens, rho1+rho2) / (rho0*maxdens)
ntr(rho1,rho2) = rho1+rho2<rho0/rhorange ? 0 : rho1+rho2<rho0*rhorange ? (rho1+rho2-rho0/rhorange)/(rhorange-1./rhorange)/rho0 : 1

# red: 0xd5 = 229 for ndr=-1, 0x32 = 50 for ndr=+1. Linear interpolation
r(rho1,rho2) = int(140 + 90*ndr(rho1,rho2))
# green: 0x3e = 62 for ndr=-1, 0x88=136 for ndr=+1
g(rho1,rho2) = int(99  - 37*ndr(rho1,rho2))
# blue: 0x4f = 79 for ndr=-1, 0xbd = 189 for ndr=+1
b(rho1,rho2) = int(134 - 55*ndr(rho1,rho2))
# alpha: normalized totalrho, normalized to 0..255
# a(rho1,rho2) = int(255*ntr(rho1,rho2))
# "weaker" version of alpha: instead of 0..255, map to alpmin..alpmax
alpmin = 0
alpmax = 200
a(rho1,rho2) = int(alpmin + (alpmax-alpmin)*ntr(rho1,rho2))

# ra, ga, ba are red,green,blue with "pseudo-alpha": shifted towards white to become pale as alpha is low
ra(rho1,rho2) = r(rho1,rho2) + (255-a(rho1,rho2))/255.*(255-r(rho1,rho2))
ga(rho1,rho2) = g(rho1,rho2) + (255-a(rho1,rho2))/255.*(255-g(rho1,rho2))
ba(rho1,rho2) = b(rho1,rho2) + (255-a(rho1,rho2))/255.*(255-b(rho1,rho2))



set term pngcairo font ",28" size 2000,1000

set xrange [0:L]
set yrange [0:L]

unset xtics
unset ytics

datei1 = "pos1_1.txt"
datei2 = "pos2_1.txt"

do for [stepsize in "100 30 5 1"]{ \
system "mkdir -p animation" . stepsize ;\
do for [k=1:anzahl_frames-1:stepsize]{ \
t_inTJ = k / frames_per_tJ ;\
t_str = sprintf("t = %f t_J (%d/%d)", t_inTJ,k,anzahl_frames) ;\
set output "animation".stepsize. "/frame_". sprintf("%03.0f",k) . ".png"; \
set multiplot; \
set label 84 t_str at graph 0.3,1.05 front; \
set origin 0,0; set size 0.48,0.96; \
set arrow 99 from L/10,L/2 to L/10+lam,L/2 arrowstyle 1 ;\
set label 99 at L/10+lam/2,L*0.52 "{/Symbol l}" font ",40" front ;\
unset grid;\
set xr [0:L]; set yr [0:L] ;\
unset xlabel; unset xtics; unset ylabel; unset ytics; unset key; \
plot "localDens/rho".k.".txt" using 2:3:(ra($4,$5)):(ga($4,$5)):(ba($4,$5)) with rgbimage; \
unset arrow 99; unset label 99; unset label 84; \
set size 0.16,0.33; set origin 0.45, 0 ;\
set xr [0:rhomax]; set yr [0: 0.05]; set xtics ("0" 0, "{/Symbol r}_{max}" rhomax, "" 0.5*rhomax) font ",12" offset 0,0.5 nomirror; unset ytics; \
plot "rhoHist.txt" every :::k::k u 2:3 w filledcurves y=0 lc rgb "#d53e4f" ;\
set origin 0.57,0; \
plot "rhoHist.txt" every :::k::k u 2:4 w filledcurves y=0 lc rgb "#3288bd" ;\
set origin 0.68,0;\
plot "rhoHist.txt" every :::k::k u 2:5 w filledcurves y=0 lc rgb "dark-violet" ;\
set origin 0.8,0;\
set xr [-1:1]; set xtics ("blue" -1, "1:1" 0, "red" 1); \
set label 87 "histogram\nof local\ndensities" font ",14" at graph 1.05,0.5; \
plot "drhoHist.txt" every :::k::k using 2:3 w filledcurves y=0 lc rgb "#fdbf6f";\
unset label 87
set origin 0.45,0.3; set size 0.28,0.33;\
set xr [*:*]; set yr [*:*]; set xtics autofreq; set ytics autofreq font ",15"; \
plot "rhok_iso/t".k.".txt" u 1:4 w l ls 13, "" u 1:5 w l ls 14 ;\
unset ytics; \
set origin 0.68,0.3 ; set size 0.28,0.33; \
set label 87 "{/Symbol r}({/Symbol s}k)" at graph 1.05,0.5; \
plot "rhok_iso/t".k.".txt" u 1:2 w l ls 11, "" u 1:3 w l ls 12 ;\
unset label 87; \
set origin 0.43,0.55; set size 0.3,0.35 ;\
set xr [0:L/2]; set xtics format "%g" font ",12"; set yr [-1:5]; set ytics autofreq format "%g" font ",15"; \
set label 84 "red"  tc rgb "#d53e4f" at graph 0.2,1.13 front; \
set label 85 "and"  tc rgb "black"   at graph 0.4,1.13 front; \
set label 86 "blue" tc rgb "#3288bd" at graph 0.6,1.13 front; \
plot "korrfunk_run1.txt" every :::k::k using 2:3 ls 13, "" every :::k::k using 2:5 ls 14;\
set origin 0.68,0.55; set size 0.28,0.35 ;\
set ytics format "";\
set label 84 "sum"  tc rgb "dark-violet" at graph 0.2,1.13 front; \
set label 85 "and"  tc rgb "black"       at graph 0.4,1.13 front; \
set label 86 "diff" tc rgb "#fdbf6f"     at graph 0.6,1.13 front; \
set label 87 "g(r/{/Symbol s})" at graph 1.05,0.5;\
plot "korrfunk_run1.txt" every :::k::k using 2:(0.25*($3+$5+2*$4)) ls 11, "" every :::k::k using 2:(0.25*($3+$5-2*$4)) ls 12;
unset label 84; unset label 85; unset label 86; unset label 87;\
unset multiplot
}}



#set size 0.25,0.5; set origin 0.5, 0.45; \
#unset key ;\
#unset arrow 99; unset label 99; set grid y; set xtics format ""; \
#set xr [*:*]; set yr [0:0.05]; set ytics format "%g" font ",15"; set xtics format "%g" font ",15"; \
#plot "rhoHist.txt" every :::k::k using 2:3 w boxes ;\
#set size 0.25,0.5; set origin 0.75,0.45; \
#set ytics format "";\
#plot "drhoHist.txt" every :::k::k using 2:3 w boxes;\
#set size 0.5,0.5; set origin 0.5,0; set yr [0:5];\
#set key top right ;\
#set xtics format "%g"; set xlabel "r/{/Symbol s}";\
#plot "korrfunk_run1.txt" every :::k::k using 2:3 ls 1 ps 1.5 ti "red-red", "" every :::k::k using 2:5 ls 2 ps 1.5 ti "blue-blue", "" every :::k::k using 2:4 ls 3 ti "red-blue";\



# old plot command with dots (particles) AND densities.
# plot "localDens/rho".k.".txt" using 2:3:(ra($4,$5)):(ga($4,$5)):(ba($4,$5)) with rgbimage, datei1 using 2:3 every :::k::k ls 1, datei2 using 2:3 every :::k::k ls 2; \

# old plot command for top-right, sum and ratio of g(r)'s
plot "korrfunk_run1.txt" every :::k::k using 2:(($3+$5)/(2*$4)) lc rgb "black" lw 3 ti "ratio in-species/mixed", "" every :::k::k using 2:(($3+$5)/2 - $4) ti 'difference in-species/mixed' lc rgb "#ff7f00" lw 3;\

set term qt
set origin 0,0
set size 1,1
set xr [0:L]
set yr [0:L]
unset ytics




#set term postscript enhanced eps color solid font ",28" size 3,3
#
#do for [k=0:(anzahl_frames-1)/frame_skip] { \
#set output 'animation_frame.'.k.'.eps' ; \
#plot datei1 using 2:3 every :::k*frame_skip::k*frame_skip ls 2, datei2 using 2:3 every :::k*frame_skip::k*frame_skip ls 1 ;\
#print 'animation_frame.'.k.'.eps'
#}
#
## alle in pdf convertieren
#! ls *.eps | xargs -n1 epstopdf
#
## alle eps löschen
#! rm *.eps
#
#
#set term qt


#!mkdir animation10 animation3 animation1
#do for [stepsize in "10 3 1"]{ \
#do for [k=0:anzahl_frames-1:stepsize]{ \
#unset key; \
#t_inTJ = k / frames_per_tJ ; \
#ti_str = sprintf("t = %f t_J (%d/%d)", t_inTJ,k,anzahl_frames) ;\
#set output "animation".stepsize. "/frame_". sprintf("%03.0f",k) . ".png"; \
#set multiplot layout 1,2 title ti_str ;\
#set arrow 99 from L/10,L/2 to L/10+lam,L/2 arrowstyle 1
#set label 99 at L/10+lam/2,L*0.52 "{/Symbol l}" font ",40"
#set xr [0:L]; set yr [0:L]; \
#unset xlabel; unset xtics; unset ytics; unset key; \
#plot datei1 using 2:3 every :::k::k ls 1, datei2 using 2:3 every :::k::k ls 2; \
#set key top right ;\
#unset arrow 99; unset label 99; \
#set xr [0.0001 : L*0.5];\
#set xlabel "r/{/Symbol s}" ;\
#set yr [0 : 5] ;\
#set ytics add 1 ;\
#set xtics ;\
#plot 'korrfunk_run1.txt' every :::k::k using 2:($3*$2/32) ls 1 ti 'blue-blue', '' every :::k::k using 2:($5*$2/32) ls 2 ti 'red-red', '' every :::k::k using 2:($4*$2/32) ls 3 ti #'red-blue', '' every :::k::k using 2:(($3+$5)/(2*$4)) lc rgb 'black' lw 3 ti 'order parameter (g_{11}+g_{22}/g_{12})';\
#unset multiplot;
#} }
