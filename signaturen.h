#pragma once


// Index-Wrapping, um 2-dimensionale Felder mit einem Index anzusprechen. FFTW erwartet Felder mit nur einem Index
int iw(int i, int j);

// diese drei Funktionen fuehren density-Gridding durch, dh berechnen aus N Teilchenpositionen im Feld r[][] eine Dichteverteilung rhox[][].
// (der zweite Index ist Imaginaerteil, alle Eintraege stehen in rhox[.][0], alle rhox[.][1] werden nicht veraendert.)
// es wird Index-Wrapping mit int iw(int, int) verwendet. Alle drei schreiben in rhox.
void gridDensity_NGP(fftw_complex* rhox, double** r); //Nearest Grid Point
void gridDensity_CIC(fftw_complex* rhox, double** r); //Cloud In Cell
void gridDensity_TSC(fftw_complex* rhox, double** r); //Triangle-Shaped Cloud


// diese drei Funktionen fuehren inverses density-Gridding durch, dh berechnen aus N Teilchenpositionen und den Kraeften an Gitterpunkten die Kraefte auf Teilchen.
void inv_gridDensity_NGP(double** Fkap, fftw_complex* Fx, fftw_complex* Fy, double** r);
void inv_gridDensity_CIC(double** Fkap, fftw_complex* Fx, fftw_complex* Fy, double** r);
void inv_gridDensity_TSC(double** Fkap, fftw_complex* Fx, fftw_complex* Fy, double** r);


// initialisiert Felder, plant Fouriertrafo
void kapkraefte_init();
//berechnet Kapillarkraefte (mittels Fouriertransformation), schreibt sie in Fkap
void berechne_kapkraefte(double** r, double** Fkap);

