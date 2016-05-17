# zeichnet Schnappschüsse eines Runs, animiert

reset


set xlabel 'x'
set ylabel 'y'
set size square


anzahl_frames = 100
fps = 3


L = 400

set style line 1 lw 2
set style line 2 lw 2 lc 3

set xrange [0:L]
set yrange [0:L]

#edit, um anderen Run anzuschauen
datei1 = "pos1_5.txt"
datei2 = "pos2_5.txt"

do for [k=0:anzahl_frames]{ \
plot datei1 using 2:3 every :::k::k ls 1, datei2 using 2:3 every :::k::k ls 2 ; \
pause 1./fps \
}