# zeichnet den Verlauf der Staerke der Kapillarkraft und die Vorhersage, fuer verschiedene Gridding-Schemata



# Vorhersage
K1 = 'bessel_aus_NR/K1.txt'


# kapillarlaenge
lam = 20.0


set xrange [0:1.9]
set yrange [0:0.6]

set xtics font ",20"
set ytics font ",20"

set key font ",20" spacing 1.5

set xlabel "r/{/Symbol l}" font ",20"
set ylabel "F [a.u.]" rotate

plot K1 using ($1*lam):($2/2/pi/lam) with lines lw 2 title 'analytic solution', 'Fkap_NGP.txt' lw 2 lc rgb "#00A000" title 'Nearest Grid Point', 'Fkap_CIC.txt' lw 2 lc 3 title 'Cloud-In-Cell', 'Fkap_TSC.txt' lw 2 lc rgb "#000000" title 'Triangle-Shaped Cloud'

set term postscript enhanced eps color solid font ",20"
set output "Kraefte_je_Gridschema.eps"
replot
set term wxt

!epstopdf Kraefte_je_Gridschema.eps
!rm Kraefte_je_Gridschema.eps
