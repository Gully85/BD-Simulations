//stellt Methoden bereit, die in der Auswertung benutzt werden

#pragma once

#include "parameter.h"
#define _USE_MATH_DEFINES
#include <math.h>


extern const double korr_dr;
extern const int korr_bins;
extern const int ftrho_qbins;
extern const double ftrho_dq;
extern const double obs_dt;
extern const int obs_anzahl;
extern const int runs;

extern int** r_git;
extern double** r_rel;
extern const double nachList_Breite;
extern const int N;


using std::vector;


//Paarkorrelationsfunktion. Erster Index Durchlauf, zweiter Index Zeit in Einheiten obs_dt, dritter Index r in Einheiten korr_dr
vector<double>** g11 = NULL;
// vector<vector<double> > g11;

//Paarkorrelationsfunktion, über Runs gemittelt. Erster Index Zeit in Einheiten obs_dt, zweiter Index r in Einheiten korr_dr
vector<double>* g11_mittel = NULL;

//Fehlerbalken Paarkorrelationsfunktion, über Runs gemittelt. Erster Index Zeit in Einheiten obs_dt, zweiter Index r in Einheiten korr_dr
vector<double>* g11_fehler = NULL;

//fouriertransformierte Dichte, Realteil und Imaginärteil. Erster Index Durchlauf, zweiter Index Zeit/obs_dt, dritter Index r (nachschlagen in order[])
vector<double>** ftrho_re = NULL;
vector<double>** ftrho_im = NULL;

//Hilfsfelder für fouriertransformierte Dichte
vector<int> ftrho_order;
vector<int> ftrho_qabs2;
vector<int> ftrho_beitraege;



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


void init_ftrho(){
	
	if(ftrho_re != NULL) return;
	
	ftrho_order.clear();
	ftrho_qabs2.clear();
	ftrho_beitraege.clear();
	
	const int qbins2 = ftrho_qbins*ftrho_qbins;
	
	int stelle_in_order=0;
	
	for(int i=0; i<ftrho_qbins; i++)
	for(int j=0; j<ftrho_qbins; j++){
		int ind2 = i*i + j*j;
		if(ind2 > qbins2) break;
		
		int b = suche(ind2, ftrho_qabs2); // -1 wenn nicht drin, sonst Index wo ind2 in qabs steht
		
		if(-1 == b){
			ftrho_qabs2.push_back(ind2);
			ftrho_beitraege.push_back(1);
			ftrho_order.push_back(stelle_in_order++);
		}//if nicht drin
		else{
			ftrho_beitraege[b]++;
			ftrho_order.push_back(b);
		}//else drin		
	}//for i,j
	
	//reserviere Vektoren für alle record()-Aufrufe
	ftrho_re = new vector<double>*[runs];
	ftrho_im = new vector<double>*[runs];
	for(int i=0; i<runs; i++){
		ftrho_re[i] = new vector<double>[obs_anzahl];
		ftrho_im[i] = new vector<double>[obs_anzahl];
	}//for i bis runs
	
}//void init_ftrho

//schreibt aktuelles rho(k) in ftrho_re[ar][t] und ftrho_im
void record_ftrho_unkorrigiert(int ar, int t){
	
	//Vektor auf richtige Länge bringen, Nullen reinschreiben
	ftrho_re[ar][t].assign(ftrho_qabs2.size(), 0.0);
	ftrho_im[ar][t].assign(ftrho_qabs2.size(), 0.0);
	
	const int qbins2 = ftrho_qbins*ftrho_qbins;
	
	int stelle_in_order=0;
	
	for(int i=0; i<ftrho_qbins; i++)
	for(int j=0; j<ftrho_qbins; j++){
		int ind2 = i*i + j*j;
		
		if(ind2 > qbins2) break;
		
		//berechne die Summen cos(k*r_i) und sin(k*r_i)
		double sum_c = 0.0;
		double sum_s = 0.0;
		for(int teilchenNr=0; teilchenNr<N; teilchenNr++){
			double x = r_git[teilchenNr][0]*nachList_Breite + r_rel[teilchenNr][0];
			double y = r_git[teilchenNr][1]*nachList_Breite + r_rel[teilchenNr][0];
			
			double skalarprodukt = ftrho_dq*(x*i + y*j);
			
			sum_c += cos(skalarprodukt);
			sum_s += sin(skalarprodukt);
		}//for teilchen
		
		//addiere an richtige Stelle in S[.] diese Summe
		int stelle = ftrho_order[stelle_in_order++];
		ftrho_re[ar][t][stelle] += sum_c;
		ftrho_im[ar][t][stelle] += sum_s;
		
	}//for i,j
}//void record_ftrho

//dividiert Korrekturen aus ftrho raus. Für alle runs und obs_anzahl, nur einmal aufrufen
void korrigiere_ftrho(){
	for(int ar=0; ar<runs; ar++)
	for(int t=0; t<obs_anzahl; t++){
		double faktor = L*L;
		for(int i=0; i<ftrho_re[0][0].size(); i++){
			ftrho_re[ar][t][i] /= faktor*ftrho_beitraege[i];
			ftrho_im[ar][t][i] /= faktor*ftrho_beitraege[i];
		}//for i bis size
	}//for run, t
}//void korrigiere_ftrho

//suche Position im vector, an der die Zahl a steht. Wenn nicht drin, gebe -1 zurück
int suche(int a, vector<int> v){
	
	for(int i=0; i<v.size(); i++)
		if(v[i] == a)
			return i;
		
	return -1;
		
}//int suche

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




