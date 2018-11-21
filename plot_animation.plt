reset

unset key
set xlabel 'x'
set ylabel 'y'
set size square

anzahl_frames = 500

# number of checkpoints per Jeans-time
frames_per_tJ = 1.11 / 0.05

# fps = 10


# box length
L = 100.0
# lambda is the capillary length
lam = 60.0

set style arrow 1 heads size 2,90 front lw 5
set arrow from L/10,L/2 to L/10+lam,L/2 arrowstyle 1
set label at L/10+lam/2,L*0.52 "{/Symbol l}" font ",40"


set style line 1 lc rgb "#0060ad" pt 7 ps 1.5 lt 1 lw 2
set style line 2 lc rgb "#ad3000" pt 7 ps 1.5 lt 1 lw 2

# set xlabel 'x/{/Symbol s}'
# set ylabel 'y/{/Symbol s}'

unset xlabel
unset ylabel


set xrange [0:L]
set yrange [0:L]

unset xtics
unset ytics

datei1 = "pos1_1.txt"
datei2 = "pos2_1.txt"

set term png font ",28" size 1200,1200

!mkdir animation10

do for [k=0:anzahl_frames-1:10]{ \
t_inTJ = k / frames_per_tJ ; \
ti_str = sprintf("t = %f t_J (%d/%d)", t_inTJ,k,anzahl_frames)
set title ti_str ;\
set output "animation10/frame_". sprintf("%03.0f",k) . ".png"; \
plot datei1 using 2:3 every :::k::k ls 1, datei2 using 2:3 every :::k::k ls 2; \
}

!mkdir animation3
do for [k=0:anzahl_frames-1:3]{ \
t_inTJ = k / frames_per_tJ ; \
ti_str = sprintf("t = %f t_J (%d/%d)", t_inTJ,k,anzahl_frames)
set title ti_str ;\
set output "animation3/frame_". sprintf("%03.0f",k) . ".png"; \
plot datei1 using 2:3 every :::k::k ls 1, datei2 using 2:3 every :::k::k ls 2; \
}

!mkdir animation1
do for [k=0:anzahl_frames-1:1]{ \
t_inTJ = k / frames_per_tJ ; \
ti_str = sprintf("t = %f t_J (%d/%d)", t_inTJ,k,anzahl_frames)
set title ti_str ;\
set output "animation1/frame_". sprintf("%03.0f",k) . ".png"; \
plot datei1 using 2:3 every :::k::k ls 1, datei2 using 2:3 every :::k::k ls 2; \
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