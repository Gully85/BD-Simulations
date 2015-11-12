# zeichnet Kapillar- und WCA-Kraefte ueber Abstand. Vorhersage und Ergebnis des Programms



# Vorhersage Kapillar
K1 = 'bessel_aus_NR/K1.txt'
# Vorhersage WCA. Ableitung des Potentials
WCA(x) = x<2.0**(1./6.)? -4.0*6.0*(2.0*x**-13 - x**-7) : 0

# Datei mit Ergebnissen des Programms
dat = 'Fvonr.txt'

## Parameter, muessen mit denen in parameter.h uebereinstimmen

# kapillarlaenge
lam = 15.0

# Vorfaktor der Kapillarkraft
kap_vor = 100.0



## Punkte fuer die berechneten, blasse Linien fuer die Vorhersagen
set style line 1 lc rgb "#FF0000" lw 2
set style line 2 lc rgb "#FFA0A0" lw 2

set style line 3 lc rgb "#00FF00" lw 2
set style line 4 lc rgb "#C0FFC0" lw 2

set style line 5 lc rgb "#000000" lw 2
set style line 6 lc rgb "#B0B0B0" lw 2


set xr [0:2]
set yr [-100:50]


plot dat using 1:2 ls 1 with points title "Kapillar", K1 using ($1*lam):(kap_vor*$2/2/pi/lam) ls 2 with lines noti,\
dat using 1:3 ls 3 with points title "WCA", WCA(x) ls 4 noti,\
dat using 1:4 ls 5 with points title "Summe", K1 using ($1*lam):(kap_vor*$2/2/pi/lam + WCA($1*lam)) ls 6 with lines noti
