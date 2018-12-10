reset


set xlabel 'x'
set ylabel 'y'



# get parameters from tmp-file. That file is supposed to be written by python3 ParProgress2tmp.py
set terminal unknown
stats 'tmp' u (anzahl_frames=$1, L=$2, lam=$3, frames_per_tJ=$4):1 nooutput

print anzahl_frames
print L
print lam
print frames_per_tJ


set style arrow 1 heads size 2,90 front lw 5


set style line 1 lc rgb "#0060ad" pt 7 ps 1.5 lt 1 lw 2
set style line 2 lc rgb "#ad3000" pt 7 ps 1.5 lt 1 lw 2
set style line 3 lc rgb "#30a000" pt 7 ps 1.5 lt 1 lw 2
set style line 4 lc rgb 'black'   pt 7 ps 1.5 lt 1 lw 2

# set xlabel 'x/{/Symbol s}'
# set ylabel 'y/{/Symbol s}'

unset xlabel
unset ylabel




set term png font ",28" size 2000,1000

set xrange [0:L]
set yrange [0:L]

unset xtics
unset ytics

datei1 = "pos1_1.txt"
datei2 = "pos2_1.txt"

!mkdir animation10 animation3 animation1
do for [stepsize in "10 3 1"]{ \
do for [k=1:anzahl_frames-1:stepsize]{ \
print k;\
t_inTJ = k / frames_per_tJ ;\
ti_str = sprintf("                                       t = %f t_J (%d/%d)", t_inTJ,k,anzahl_frames) ;\
set output "animation".stepsize. "/frame_". sprintf("%03.0f",k) . ".png"; \
set multiplot title ti_str;\
set origin 0,0; set size 0.5,1; \
set arrow 99 from L/10,L/2 to L/10+lam,L/2 arrowstyle 1 ;\
set label 99 at L/10+lam/2,L*0.52 "{/Symbol l}" font ",40" ;\
unset grid;\
set xr [0:L]; set yr [0:L] ;\
unset xlabel; unset xtics; unset ylabel; unset ytics; unset key; \
plot datei1 using 2:3 every :::k::k ls 1, datei2 using 2:3 every :::k::k ls 2; \
set size 0.5,0.5; set origin 0.5, 0.45; \
set key top right ;\
unset arrow 99; unset label 99; set grid y; set xtics 10 format ""; \
set xr [0.0001: L/2]; set yr [-1:3]; set ytics; \
plot "korrfunk_run1.txt" every :::k::k using 2:(($3+$5)/(2*$4)) lc rgb "black" lw 3 ti "ratio in-species/mixed", "" every :::k::k using 2:(($3+$5 - 2*$4)*$2/64) ti 'difference in-species/mixed' lc rgb "#fdc086" lw 3;\
set size 0.5,0.5; set origin 0.5,0; set yr [0:5];\
set xtics 10 format "%g"; set xlabel "r/{/Symbol s}";\
plot "korrfunk_run1.txt" every :::k::k using 2:($3*$2/64) ls 1 ti "red-red", "" every :::k::k using 2:($5*$2/64) ls 2 ti "blue-blue", "" every :::k::k using 2:($4*$2/64) ls 3 ti "red-blue";\
unset multiplot
}}







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
