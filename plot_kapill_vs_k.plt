# zeichnet Effekte der Kapillarwechselwirkung (Hoehenprofil u(r) und Kraefte F(r)) entlang x-Achse gegen die analytische Vorhersage K0(x) bzw K1(x)

reset


set style line 1 lw 5
set style line 3 lw 2


# K0 und K1 sind modifizierte Besselfunktionen. Aus Numerical Recipes
K0 = 'bessel_aus_NR/K0.txt'
K1 = 'bessel_aus_NR/K1.txt'

# Groesse der Simulationsbox
L = 100

#Kapillarlaenge
lam = 7.0

#logscale?
set xrange [0.01:L]
set yrange [0.001:20]

# set logscale



# Hoehenprofil entlang x-Achse, gegen K0/2pi
# plot '../Fx.txt' using ($1-0.5*L):($2==0.5*L? $3 : NaN) ls 1, K0 using ($1*lam):($2/(2*pi)) ls 3




# y-Komponente der Kraft entlang y-Achse, gegen K1/(2pi lambda)
 plot 'Fy.txt' using ($2-0.5*L):($1 ==0.5*L? -$3 : NaN) ls 1,K1 using ($1*lam):($2/(2*pi*lam)) ls 3
