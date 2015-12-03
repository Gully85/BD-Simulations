//stellt Methoden bereit, die in der Auswertung benutzt werden

#pragma once

#include "parameter.h"
#define _USE_MATH_DEFINES
#include <math.h>


extern const double korr_dr;
extern const int korr_bins;
extern const double obs_dt;
extern const int obs_anzahl;
extern const int runs;

extern int** r_git;
extern double** r_rel;
extern const int N;


using std::vector;


//Paarkorrelationsfunktion. Erster Index Durchlauf, zweiter Index Zeit in Einheiten obs_dt, dritter Index r in Einheiten korr_dr
vector<double>** g11;
// vector<vector<double> > g11;

//Paarkorrelationsfunktion, über Runs gemittelt. Erster Index Zeit in Einheiten obs_dt, zweiter Index r in Einheiten korr_dr
vector<double>* g11_mittel;

//Fehlerbalken Paarkorrelationsfunktion, über Runs gemittelt. Erster Index Zeit in Einheiten obs_dt, zweiter Index r in Einheiten korr_dr
vector<double>* g11_fehler;

void init_korrelationsfunktion(){
	
	if(g11 != NULL) return;
	
	g11 = new vector<double>*[runs];
	
	for(int i=0; i<runs; i++){
		g11[i] = new vector<double>[obs_anzahl];
	}//for i
	
}//void init_korrelationsfunktion

//Schreibt aktuelle Korrelationsfunktion in g11[ar][t].
void record_korrelationsfunktion(int ar, int t){
	
	//Nullen
	g11[ar][t].assign(korr_bins, 0.0); //Nullen an alle Stellen
	
	//konstanter Vorfaktor
	
	const double rho = N/L/L;
	const double vorfaktor = 2.0/((N-1)*M_PI*rho*korr_dr*korr_dr);
	
	//Schleife über alle Teilchenpaare
	for(int i=1; i<N; i++)
	for(int j=0; j<i; j++){
		
		double a2 = abstand2(i,j);
		
		if(a2 > korr_bins*korr_bins*korr_dr*korr_dr) continue;
		
		int index = (int) (sqrt(a2)/korr_dr);
		
		if(index < 0 || index >= korr_bins) continue;
		
		g11[ar][t][index] += vorfaktor /(2*index+1);
	}//for j
	
	
}//void record_korrelationsfunktion


//berechnet Mittelwerte und Varianzen von input[][][]. Mittelung über den ersten Index, der bis runs geht. Der zweite geht bis anzahl.
void statistik_1(vector<double>**& input, vector<double>*& mittelwerte, vector<double>*& varianzen, int anzahl){
	int ar; //aktueller run
	int i; // aktuelle Stelle im vector, vorletzter index
	int j; // aktuelle Stelle im vector, letzter Index
	
	//erster Index Zeit in Einheiten obs_dt, zweiter Index r in Einheiten korr_dr
	if(mittelwerte == NULL){
		mittelwerte = new vector<double>[anzahl];
		varianzen = new vector<double>[anzahl];
	}//if mittelwerte==NULL
	else cout << "Warnung: statistik_1() aufgerufen, Ergebnisvektoren schon benutzt!" << std::endl;
	
	for(i=0; i<anzahl; i++){ //zweiter Index im Input, erster in mittelwerte
		for(j=0; j<input[0][0].size(); j++){ //dritter Index im Input, zweiter in Mittelwerte
			//Mittelwert bestimmen
			mittelwerte[i].push_back(0.0); //an Stelle [i][j] eine Null
			for(ar=0; ar<runs; ar++)
				mittelwerte[i][j] += input[ar][i][j];
			mittelwerte[i][j] /= runs; //Summe durch Anzahl
			
			//Fehler bestimmen
			varianzen[i].push_back(0.0);
			for(ar=0; ar<runs; ar++)
				varianzen[i][j] += (input[ar][i][j]-mittelwerte[i][j])*(input[ar][i][j]-mittelwerte[i][j]);
			varianzen[i][j] /= runs; //Summe der Fehlerquadrate durch Anzahl. Ist Varianz^2.
			varianzen[i][j] = sqrt(varianzen[i][j]); //Varianz
			varianzen[i][j] /= sqrt((double)runs); //Standardfehler
			
		}//for j, dritter Index im input bzw zweiter in mittelwerte/varianzen
		
	}//for i, zweiter Index im Input bzw erster in mittelwerte
	
}//void statistik_1

// Korrelationsfunktion auswerten: Mittelung, Fehlerbalken berechnen und alles in Datei schreiben
void auswerten_korrelationsfunktion(){
	
	//berechne Mittelwerte und Fehler. Schreibe diese in g11_mittel und g11_fehler
	statistik_1(g11, g11_mittel, g11_fehler, obs_anzahl);
	
	//schreibe in Datei
	FILE* out = fopen("g11.txt", "w");
	fprintf(out, "# Format: r TAB t TAB g11(r,t) TAB Fehlerbalken\n");

	for(int j=0; j<obs_anzahl; j++){
		double t = j*obs_dt;
		
		for(int i=0; i<korr_bins; i++){
			double r = i*korr_dr;
			fprintf(out, "%g \t %g \t %g \t %g \n", r, t, g11_mittel[j][i], g11_fehler[j][i]);
		}//for i bis korr_bins
		fprintf(out, "\n");
	}//for j bis obs_anzahl
}//void auswerten_korrelationsfunktion





