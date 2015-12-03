#include "signaturen.h"
#include "zufall.cpp"
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <math.h>
#include "gridRoutinen.cpp"
#include <fftw3.h>
#include "dynamik_methoden.cpp"
#include <algorithm>

#include "auswertung.cpp"


using std::cout; using std::endl; using std::flush; using std::vector; using std::min;


// importiere Variablen aus parameter.h, siehe dort, was die Variable tut.
extern const int N;
extern const double L;

extern const int densGrid_Schema, densGrid_Zellen;
extern const double densGrid_Breite;

extern const int nachList_Zellen;
extern const double nachList_Breite;

extern const double lambda_kapillar;

extern const int obs_anzahl;
extern const double obs_dt;

extern const int runs;



/// Globale Felder. Mit Null initialisiert.

// Teilchenpositionen, Gittervektor. Erster Index TeilchenNr, zweiter Index 0 fuer x und 1 fuer y-Richtung
int** r_git = NULL;
// dasgleiche, Vektor innerhalb der Zelle
double** r_rel = NULL;



int main(){

	init_rng(); //Random Seed
/*
	main_init();
	double t= zeitschritt();
	
	cout << "t="<<t<<endl;
/*/
for(int run=0; run<runs; run++){
	main_init();
	cout << "run nr " << run << endl;
	
	cout << "record obs-Punkt 0 von " << obs_anzahl << endl;
	record_korrelationsfunktion(run, 0);
	for(int obs_nr=1; obs_nr<obs_anzahl; obs_nr++){
		double t=0.0;
		
		//Zeitschritte bis t=obs_dt
		while(t < obs_dt){
			t += zeitschritt(obs_dt-t);
// 			cout << "t="<<t<<endl;
		}//while obs-Punkt noch nicht erreicht
		
		cout << "record obs-Punkt " << obs_nr << " von " << obs_anzahl << endl;
		record_korrelationsfunktion(run, obs_nr);
	}//for obs_nr

}//for run


//Statistik, und Ergebnis in Datei schreiben
auswerten_korrelationsfunktion();


// */


return 0;

}//int main







//index-Wrapping, ermoeglicht es 1dim-Arrays mit zwei Indices anzusprechen. FFTW erwartet Arrays mit nur einem Index. Keine period. Randbed.
int iw(int i, int j){
	if (i>=densGrid_Zellen || j>= densGrid_Zellen || i<0 || j<0){
		cout << "Index out of Bounds! Aufruf war iw(i=" <<i<<", j=" <<j<<") \n" << endl << flush;
		return -1; //erzeugt absichtlich Segmentation Fault
	}//if index out-of-bounds

	return i*densGrid_Zellen + j;

}//int index

//berechnet das Quadrat des Abstands zwischen Teilchen i und Teilchen j (gleichen Typs). Berücksichtigt periodische Randbedingungen.
double abstand2(int i, int j){
	
	//Absolutpositionen
	double xi_abs = r_git[i][0]*nachList_Breite + r_rel[i][0];
	double xj_abs = r_git[j][0]*nachList_Breite + r_rel[j][0];
	double yi_abs = r_git[i][1]*nachList_Breite + r_rel[i][1];
	double yj_abs = r_git[j][1]*nachList_Breite + r_rel[j][1];
	
	// dx^2, durch die periodischen Randbed gibt es drei Moeglichkeiten
	double tmp1 = xi_abs - xj_abs;
	double tmp2 = xi_abs - xj_abs + L;
	double tmp3 = xi_abs - xj_abs - L;
	double dx2 = min(tmp1*tmp1, min(tmp2*tmp2, tmp3*tmp3));
	
	// das gleiche fuer dy^2
	tmp1 = yi_abs - yj_abs;
	tmp2 = yi_abs - yj_abs + L;
	tmp3 = yi_abs - yj_abs - L;
	double dy2 = min(tmp1*tmp1, min(tmp2*tmp2, tmp3*tmp3));
	
	return dx2 + dy2;
}//double abstand2