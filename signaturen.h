#pragma once

#include <fftw3.h>
#include <vector>
#include <algorithm>

using std::vector; using std::min;

extern const double dt_max;

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

//reserviere Speicher und setze Teilchen auf Startpositionen
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
//Fuehre einen Zeitschritt durch: Berechne Kraefte, ermittle optimale Dauer, bewege Teilchen, aktualisiere ggf Nachbarlisten. Gibt Dauer zurueck.
double zeitschritt(double dt = dt_max);
//bestimme optimale Dauer des Zeitschritts, so dass max_reisedistanz eingehalten wird
double optimaler_zeitschritt(double** F_WCA, double** Fkap, double** F_noise, double deltat, int N1); 
void erwListe_rem(vector<int>& liste, int i); //streiche Eintrag aus erwListe
//Greensfunktion an der Stelle qjk. Index-Wrapping, Verschiebung.
double G(int j, int k);

//reserviere Speicher für Paarkorrelationsfunktion
void init_korrelationsfunktion();//void init_korrelationsfunktion
//reserviere Speicher und vorab-Berechnungen für fouriertransformierte Dichte
void init_ftrho();
//schreibt aktuelles rho(k) in ftrho_re[ar][t] und ftrho_im
void record_ftrho_unkorrigiert(int ar, int t);
//dividiert Korrekturen aus ftrho raus
void korrigiere_ftrho();
//suche Position im vector, an der die Zahl a steht. Wenn nicht drin, gebe -1 zurück
int suche(int a, vector<int> v);


//Schreibt aktuelle Korrelationsfunktion in g11[ar][t].
void record_korrelationsfunktion(int ar, int t);
void auswerten_korrelationsfunktion();//ruft statistik auf und schreibt Ergebnisse in Datei
void statistik_1(vector<double>*& input, vector<double>& mittelwerte, vector<double>& varianzen, int anzahl);



//berechnet das Quadrat des Abstands zwischen Teilchen i und Teilchen j (gleichen Typs). Berücksichtigt periodische Randbedingungen.
double abstand2(int i, int j);
