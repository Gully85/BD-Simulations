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
void gridDensity_NGP(fftw_complex* rhox, double gewicht1, double gewicht2); //Nearest Grid Point
void gridDensity_CIC(fftw_complex* rhox, double gewicht1, double gewicht2); //Cloud In Cell
void gridDensity_TSC(fftw_complex* rhox, double gewicht1, double gewicht2); //Triangle-Shaped Cloud


// diese drei Funktionen fuehren inverses density-Gridding durch, dh berechnen aus N Teilchenpositionen und den Kraeften an Gitterpunkten die Kraefte auf Teilchen.
void inv_gridDensity_NGP();
void inv_gridDensity_CIC();
void inv_gridDensity_TSC();

//schreibt Jeans-Zeit und Jeans-Länge in cout
void ausgabe_jeansgroessen();

//reserviere Speicher und setze Teilchen auf Startpositionen
void main_init();
// initialisiert Teilchenpositionen auf Zufallspositionen. Schreibt in r_git und r_rel
void init_zufallspos(); 
// initialisiert Teilchenpositionen auf Positionen, die in Datei stehen. Dateiname wird in parameter.h festgelegt. 
void init_pos_aus_datei();
// initialisiert Teilchenpositionen auf ein regelmäßiges quadratisches Gitter. 
void init_gitterstart();
// initialisiert alle Teilchenpositionen auf (L/2, L/2), also alle am selben Ort
void init_allegleich();
// initialisiert alle Teilchenpositionen innerhalb eines Kreises (Radius 20) in der Mitte der Box
void init_kreisscheibe();
// schreibt aktuelle Teilchenpositionen in Datei.
void pos_schreiben_einedatei();
//schreibt aktuelle Teilchenpositionen in Dateien. FILE-Pointer muss übergeben werden.
void pos_schreiben(double t, FILE* datei1, FILE* datei2);

// initialisiert Felder/Nachbarlisten fuer WCA-Kraefte. Erwartet, dass in r_git schon die Gittervektoren der Teilchenpositionen stehen.
void WCA_init(); 
// berechnet WCA-Kraefte. Schreibt sie in F_WCA[][2]
void berechne_WCA11(); 
void berechne_WCA22();
void addiere_WCA12();
// initialisiert Felder fuer Kapillarkraefte, plant Fouriertrafo
void kapkraefte_init();
//berechnet Kapillarkraefte (mittels Fouriertransformation), schreibt sie in Fkap
void berechne_kapkraefte(); 
//Schreibt 2*anzahl Zufallszahlen mit Varianz 1.0 in dn[anzahl][2]
void berechne_zufallskraefte(int anzahl, double** dn);
//Fuehre einen Zeitschritt durch: Berechne Kraefte, ermittle optimale Dauer, bewege Teilchen, aktualisiere ggf Nachbarlisten. Gibt Dauer zurueck.
double zeitschritt(double dt = dt_max);
//bestimme optimale Dauer des Zeitschritts, so dass max_reisedistanz eingehalten wird
double optimaler_zeitschritt(double** F_WCA, double** Fkap, double** F_noise, double deltat, int N); 

//erneuere erwNachbarlisten. Beim Aufruf sind r_rel<0 oder r_rel>nachList_Breite zulässig, wird behoben
void refresh_erwNachbar();

//streiche Eintrag aus erwListe
void erwListe_rem(vector<int>& liste, int i); 
//Greensfunktion an der Stelle qjk. Index-Wrapping, Verschiebung.
double G(int j, int k);

//reserviere Speicher für Paarkorrelationsfunktion
void init_korrelationsfunktion();//void init_korrelationsfunktion
//reserviere Speicher und vorab-Berechnungen für fouriertransformierte Dichte
void init_ftrho();
//reserviere Speicher und plane FFTs für rho(k) via FFTW. Variable vector<double>rhoFFTW[run][t][k/dq].
void init_rhoFFTW(); 
//schreibt aktuelles rho(k) in ftrho_re[ar][t] und ftrho_im
void record_ftrho_unkorrigiert(int ar, int t);
//berechnet aktuelles rho(k) via FFTW, schreibt es in rhoFFTW[ar][t][.]
void record_rhoFFTW(int run, int t);
//dividiert Korrekturen aus ftrho raus
void korrigiere_ftrho();
//ruft korrigiere() und statistik() auf und schreibt Ergebnisse in Datei. 
void auswerten_ftrho();
//berechnet Mittelwerte und Fehler von rhoFFT, schreibt in Datei rhoFFT.txt
void auswerten_rhoFFTW();
//berechnet eine Variante von rhoFFT, bei der zuerst jeder Run mit seinem Startwert normiert wird und dann erst über Runs gemittelt. Schreibt in rhoFFTW_re und _im
void auswerten_rhoFFTW_normjerun();


//suche Position im vector, an der die Zahl a steht. Wenn nicht drin, gebe -1 zurück
int suche(int a, vector<int> v);


//Schreibt aktuelle Korrelationsfunktion in g11[ar][t].
void record_korrelationsfunktion11(int ar, int t);
void record_korrelationsfunktion12(int ar, int t);
void record_korrelationsfunktion22(int ar, int t);
void auswerten_korrelationsfunktion();//ruft statistik auf und schreibt Ergebnisse in Datei. Es wird g11 und g22 ausgewertet.
void auswerten_korrelationsfunktion_mixed(); //das gleiche für g12
void statistik_1(vector<double>**& input, vector<double>*& mittelwerte, vector<double>*& varianzen, int anzahl);



//berechnet das Quadrat des Abstands zwischen Teilchen i(Typ 1) und Teilchen j(Typ 1). Berücksichtigt periodische Randbedingungen.
double abstand2_11(int i, int j);
double abstand2_12(int i, int j); // i ist Typ1, j ist Typ2
double abstand2_22(int i, int j);
