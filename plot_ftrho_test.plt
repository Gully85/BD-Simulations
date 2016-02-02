# zeichnet die fouriertransformierte Dichte, und die Vorhersage für Dichteverteilung=Kreisscheibe


# Teilchenzahl
N=10000

# Boxgröße, je Richtung
L=600

# Radius der Kreisscheibe
a = 15.0



rho1 = N/(a*a*pi)

plot 'ftrho.txt' using 1:3 title 're', 'ftrho.txt' using 1:5 title 'im', rho1*2*pi*a*a*besj1(x*a)/(x*a)
