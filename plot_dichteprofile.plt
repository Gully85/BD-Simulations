# zeichnet Dichteprofile für beide Teilchentypen. Animation für zeitlichen Verlauf.

reset

anzahl_frames = 100
fps = 3


set xrange [0:100]

datei = "dichteprofil.txt"

set style line 1 lw 3 
set style line 2 lw 3 lc 3


do for [k=0:anzahl_frames]{ \
plot datei u 1:3:4 every :::k::k ls 1 w yerrorbars title "Typ 1", datei u 1:5:6 every :::k::k ls 2 w yerrorbars title "Typ 2"; \
pause 1./fps
}