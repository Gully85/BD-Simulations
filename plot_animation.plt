unset key
set xlabel 'x'
set ylabel 'y'
set size square

anzahl_frames = 2000

#zeichne nur jeden n-ten Frame
frame_skip = 20


fps = 10



L = 141.42

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

set term qt

do for [k=0:(anzahl_frames-1)/frame_skip] { \
plot datei1 using 2:3 every :::k*frame_skip::k*frame_skip ls 1, datei2 using 2:3 every :::k*frame_skip::k*frame_skip ls 2 ; \
pause 1./fps }


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