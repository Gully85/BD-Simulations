#pragma once


// Index-Wrapping, um 2-dimensionale Felder mit einem Index anzusprechen. FFTW erwartet Felder mit nur einem Index
int iw(int i, int j);

// diese drei Funktionen fuehren density-Gridding durch, dh berechnen aus N Teilchenpositionen im Feld r[][] eine Dichteverteilung rhox[][].
// (der zweite Index ist Imaginaerteil, alle Eintraege stehen in rhox[.][0], alle rhox[.][1] werden nicht veraendert.)
// es wird Index-Wrapping mit int iw(int, int) verwendet. Alle drei schreiben in rhox.
void gridDensity_NGP(fftw_complex* rhox, int** rr_git, double** rr_rel); //Nearest Grid Point
void gridDensity_CIC(fftw_complex* rhox, int** rr_git, double** rr_rel); //Cloud In Cell
void gridDensity_TSC(fftw_complex* rhox, int** rr_git, double** rr_rel); //Triangle-Shaped Cloud


// diese drei Funktionen fuehren inverses density-Gridding durch, dh berechnen aus N Teilchenpositionen und den Kraeften an Gitterpunkten die Kraefte auf Teilchen.
void inv_gridDensity_NGP(double** Fkap, fftw_complex* Fx, fftw_complex* Fy, int** rr_git, double** rr_rel); 
void inv_gridDensity_CIC(double** Fkap, fftw_complex* Fx, fftw_complex* Fy, int** rr_git, double** rr_rel); 
void inv_gridDensity_TSC(double** Fkap, fftw_complex* Fx, fftw_complex* Fy, int** rr_git, double** rr_rel); 

//reserviere Speicher fuer Felder in der main
void main_init();
// initialisiert Teilchenpositionen auf Zufallspositionen. Schreibt in r_git und r_rel
void init_zufallspos(); 
// initialisiert Felder/Nachbarlisten fuer WCA-Kraefte. Erwartet, dass in r_git schon die Gittervektoren der Teilchenpositionen stehen.
void WCA_init(); 
// berechnet WCA-Kraefte. Schreibt sie in F_WCA[][2]
void berechne_WCAkraefte(double** F_WCA); 
// initialisiert Felder fuer Kapillarkraefte, plant Fouriertrafo
void kapkraefte_init();
//berechnet Kapillarkraefte (mittels Fouriertransformation), schreibt sie in Fkap
void berechne_kapkraefte(int** rr_git, double** rr_rel, double** Fkap); 
//Schreibt 2*anzahl Zufallszahlen mit Varianz 1.0 in dn[anzahl][2]
void berechne_zufallskraefte(int anzahl, double** dn);

//Greensfunktion an der Stelle qjk. Index-Wrapping, Verschiebung.
double G(int j, int k);

//berechnet das Quadrat des Abstands zwischen Teilchen i und Teilchen j (gleichen Typs). Berücksichtigt periodische Randbedingungen.
double abstand2(int i, int j, int** r_git, double** r_rel);