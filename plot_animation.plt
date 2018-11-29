reset

unset key
set xlabel 'x'
set ylabel 'y'
set size square



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




set term png font ",28" size 2400,1200

set xrange [0:L]
set yrange [0:L]

unset xtics
unset ytics

datei1 = "pos1_1.txt"
datei2 = "pos2_1.txt"



!mkdir animation10

do for [k=0:anzahl_frames-1:10]{ \
t_inTJ = k / frames_per_tJ ; \
ti_str = sprintf("t = %f t_J (%d/%d)", t_inTJ,k,anzahl_frames) ;\
set output "animation10/frame_". sprintf("%03.0f",k) . ".png"; \
set multiplot layout 1,2 title ti_str ;\
set arrow 99 from L/10,L/2 to L/10+lam,L/2 arrowstyle 1
set label 99 at L/10+lam/2,L*0.52 "{/Symbol l}" font ",40"
set xr [0:L]; set yr [0:L]; \
unset xlabel; unset xtics; unset ytics; \
plot datei1 using 2:3 every :::k::k ls 1, datei2 using 2:3 every :::k::k ls 2; \
unset arrow 99; unset label 99; \
set xr [0.0001 : L*0.5];\
set xlabel "r/{/Symbol s}" ;\
set yr [0 : 5] ;\
set ytics add 1 ;\
set xtics ;\
plot 'korrfunk_run1.txt' every :::k::k using 2:($3*$2/32) ls 1, '' every :::k::k using 2:($5*$2/32) ls 2, '' every :::k::k using 2:($4*$2/32) ls 3, '' every :::k::k using 2:(($3+$5)/(2*$4)) lc rgb 'black' lw 3; \
unset multiplot;
}

!mkdir animation3
do for [k=0:anzahl_frames-1:3]{ \
t_inTJ = k / frames_per_tJ ; \
ti_str = sprintf("t = %f t_J (%d/%d)", t_inTJ,k,anzahl_frames) ;\
set output "animation3/frame_". sprintf("%03.0f",k) . ".png"; \
set multiplot layout 1,2 title ti_str ;\
set arrow 99 from L/10,L/2 to L/10+lam,L/2 arrowstyle 1
set label 99 at L/10+lam/2,L*0.52 "{/Symbol l}" font ",40"
set xr [0:L]; set yr [0:L]; \
unset xlabel; unset xtics; unset ytics; \
plot datei1 using 2:3 every :::k::k ls 1, datei2 using 2:3 every :::k::k ls 2; \
unset arrow 99; unset label 99; \
set xr [0.0001 : L*0.5];\
set xlabel "r/{/Symbol s}" ;\
set yr [0 : 5] ;\
set ytics add 1 ;\
set xtics ;\
plot 'korrfunk_run1.txt' every :::k::k using 2:($3*$2/32) ls 1, '' every :::k::k using 2:($5*$2/32) ls 2, '' every :::k::k using 2:($4*$2/32) ls 3, '' every :::k::k using 2:(($3+$5)/(2*$4)) lc rgb 'black' lw 3; \
unset multiplot;
}

!mkdir animation1
do for [k=0:anzahl_frames-1:1]{ \
t_inTJ = k / frames_per_tJ ; \
ti_str = sprintf("t = %f t_J (%d/%d)", t_inTJ,k,anzahl_frames) ;\
set output "animation1/frame_". sprintf("%03.0f",k) . ".png"; \
set multiplot layout 1,2 title ti_str ;\
set arrow 99 from L/10,L/2 to L/10+lam,L/2 arrowstyle 1
set label 99 at L/10+lam/2,L*0.52 "{/Symbol l}" font ",40"
set xr [0:L]; set yr [0:L]; \
unset xlabel; unset xtics; unset ytics; \
plot datei1 using 2:3 every :::k::k ls 1, datei2 using 2:3 every :::k::k ls 2; \
unset arrow 99; unset label 99; \
set xr [0.0001 : L*0.5];\
set xlabel "r/{/Symbol s}" ;\
set yr [0 : 5] ;\
set ytics add 1 ;\
set xtics ;\
plot 'korrfunk_run1.txt' every :::k::k using 2:($3*$2/32) ls 1, '' every :::k::k using 2:($5*$2/32) ls 2, '' every :::k::k using 2:($4*$2/32) ls 3, '' every :::k::k using 2:(($3+$5)/(2*$4)) lc rgb 'black' lw 3; \
unset multiplot;
}

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
