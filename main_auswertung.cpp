#include "signaturen_auswertung.h"
#include "parameter.h"
#include <iostream>

using std::cout; using std::endl; using std::flush;

extern const double obs_dt;
extern const int obs_anzahl;
extern const int N1, N2;
extern const int jobs;
//Teilchenpositionen. Erster Index Job (0 bis jobs-1), zweiter Index Zeit (in Einheiten obs_dt, 0 bis obs_anzahl-1), dritter Index Teilchen (bis N1 bzw N2), vierter Index 0/1 für x/y.
double**** r1_abs = NULL;
double**** r2_abs = NULL;



int main(){
	
	//reserviere Speicher
	r1_abs = new double***[jobs];
	if(0 != N2) r2_abs = new double***[jobs];
	
	for(int job=0; job<jobs; job++){
		r1_abs[job] = new double**[obs_anzahl];
		if(0 != N2) r2_abs[job] = new double**[obs_anzahl];
		
		for(int t=0; t<obs_anzahl; t++){
			r1_abs[job][t] = new double*[N1];
			if(0 != N2) r2_abs[job][t] = new double*[N2];
			
			for(int teilchen=0; teilchen<N1; teilchen++)
				r1_abs[job][t][teilchen] = new double[2];
			for(int teilchen=0; teilchen<N2; teilchen++)
				r2_abs[job][t][teilchen] = new double[2];
		}//for i bis obs_anzahl
	}
	cout << "Lese Positionen..." << endl;
	//lese Positionen (zu allen Zeiten) aus Dateien
	rvont_lesen();
	
	//cout << "Berechne Mittelpunkte und Dichteprofile..." << endl;
	
	//berechne_dichteprofile();
	
	cout << "Berechne Korrelationsfunktionen..." << endl;
	berechne_korrelationsfunktionen();
	
	cout << "Berechne Mittelwerte und Fehler..." << endl;
	//schreibe Ergebnisse
	//auswerten_dichteprofile();
	auswerten_korrelationsfunktionen();
	
	return 0;
	
}//int main