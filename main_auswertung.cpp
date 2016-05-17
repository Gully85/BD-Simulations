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
		
		for(int i=0; i<obs_anzahl; i++){
			r1_abs[job][i] = new double*[N1];
			if(0 != N2) r2_abs[job][i] = new double*[N2];
			
			for(int j=0; j<N1; j++)
				r1_abs[job][i][j] = new double[2];
			for(int j=0; j<N2; j++)
				r2_abs[job][i][j] = new double[2];
		}//for i bis obs_anzahl
	}
	
	//lese Positionen (zu allen Zeiten) aus Dateien
	rvont_lesen();
	
	//berechne geforderte Observable
	berechne_dichteprofile();
	
	//schreibe Ergebnisse
	auswerten_dichteprofile();
	
	return 0;
	
}//int main