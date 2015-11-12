# zeichnet den Verlauf der Staerke der Kapillarkraft und die Vorhersage, fuer verschiedene Gridding-Schemata



# Vorhersage
K1 = 'bessel_aus_NR/K1.txt'


# kapillarlaenge
lam = 20.0


set xrange [0:3]
set yrange [0:0.9]


plot K1 using ($1*lam):($2/2/pi/lam) with lines, 'Fkap_NGP.txt' lw 2 lc rgb "#008000", 'Fkap_CIC.txt' lw 2 lc 3, 'Fkap_TSC.txt' lw 2 lc rgb "#000000"
