//stellt Methoden bereit, die in reinen Auswertungs-Runs die Observablen berechnen 

#include "parameter.h"
#include "signaturen_auswertung.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <vector>
using std::vector;

#include <cstdio>



//Erster Index Job (0 bis jobs-1), zweiter Index Zeit (in Einheiten obs_dt, 0 bis obs_anzahl-1), dritter Index Teilchen (bis N1 bzw N2), vierter Index 0/1 für x/y.
extern double**** r1_abs;
extern double**** r2_abs;

extern const double Gamma2_1;
extern const double rhoring_dr;
extern const int rhoring_bins;


extern const int jobs;

extern const int N1, N2;



namespace Observablen{
	//Dichteprofil. Erster Index Run bzw Job, zweiter Index Zeit in Einheiten obs_dt, dritter Index r in Einheiten rhoring_dr
	vector<double>** rho1ring = NULL;
	vector<double>** rho2ring = NULL;
	
	//Mobilitätsschwerpunkte. Erster Index Run bzw Job, zweiter Index Zeit in Einheiten obs_dt, dritter Index 0/1 für x/y
	double*** Sgamma = NULL;

}


//berechnet für alle Zeiten und alle Runs/Jobs die Mobilitätsschwerpunkte (sum(typ1) r + 1/Gamma2_1 sum(typ2) r)/(N1 + N2/Gamma2_1). Schreibt diese in Sgamma
void berechne_mobischwerpunkte(){
	
	using namespace Observablen;
	
	//reserviere Speicher für Sgamma
	Sgamma = new double**[jobs];
	for(int job=0; job<jobs; job++){
		Sgamma[job] = new double*[obs_anzahl];
		for(int zeit=0; zeit<obs_anzahl; zeit++)
			Sgamma[job][zeit] = new double[2];
	}//for job
	
	
	//berechne Mobischwerpunkte
	for(int job=0; job<jobs; job++)
	for(int zeit=0; zeit<obs_anzahl; zeit++){
		double x=0.0;
		double y=0.0;
		
		double nenner = N1 + N2/(Gamma2_1);
		
		for(int teilchen=0; teilchen<N2; teilchen++){
			x += r2_abs[job][zeit][teilchen][0];
			y += r2_abs[job][zeit][teilchen][1];
		}//for teilchen bis N2
		x/=Gamma2_1;
		y/=Gamma2_1;
		for(int teilchen=0; teilchen<N1; teilchen++){
			x += r1_abs[job][zeit][teilchen][0];
			y += r1_abs[job][zeit][teilchen][1];
		}//for teilchen bis N1
		
		Sgamma[job][zeit][0] = x/nenner;
		Sgamma[job][zeit][1] = y/nenner;
		
	}//for job,zeit
	
}//void berechne_mobischwerpunkte





//berechnet für jeden Run und jede Zeit das Dichteprofil rho1ring und rho2ring. Schreibt es in rho1ring[][][]
void berechne_dichteprofile(){
	using namespace Observablen;
	
	berechne_mobischwerpunkte();
	
	rho1ring = new vector<double>*[jobs];
	if(0 != N2) rho2ring = new vector<double>*[jobs];
	
	
	//für jeden Run...
	for(int job=0; job<jobs; job++){
		rho1ring[job] = new vector<double>[obs_anzahl];
		if(0 != N2) rho2ring[job] = new vector<double>[obs_anzahl];
		//für jede Zeit...
		for(int zeit=0; zeit<obs_anzahl; zeit++){
			rho1ring[job][zeit].assign(rhoring_bins, 0.0);
			
			const double vorfaktor1 = L*L/(N1*rhoring_dr*rhoring_dr*M_PI);
			
			//für jedes Teilchen Typ 1...
			for(int teilchen=0; teilchen<N1; teilchen++){
				//Abstand des Teilchens vom Mobilitätsschwerpunkt. KEINE periodischen Randbedingungen! (weil nicht nötig, Teilchen und Mobischwerpunkt befinden sich in der Mitte der Box)
				double x = r1_abs[job][zeit][teilchen][0];
				double y = r1_abs[job][zeit][teilchen][1];
				double dx = x - Sgamma[job][zeit][0];
				double dy = y - Sgamma[job][zeit][1];
				
				double r = sqrt(dx*dx + dy*dy);
				
				int bin = (int) (r/rhoring_dr);
				if(0 > bin || bin >= rhoring_bins) continue;
				
				rho1ring[job][zeit][bin] += vorfaktor1/(2*bin+1);
				
			}//for teilchen Typ 1
			
			
			if(0 == N2) continue;
			
			rho2ring[job][zeit].assign(rhoring_bins, 0.0);
			
			const double vorfaktor2 = L*L/(N2*rhoring_dr*rhoring_dr*M_PI);
			
			//für jedes Teilchen Typ 2...
			for(int teilchen=0; teilchen<N2; teilchen++){
				//Abstand des Teilchens vom Mobilitätsschwerpunkt. KEINE periodischen Randbedingungen! (weil nicht nötig, Teilchen und Mobischwerpunkt befinden sich in der Mitte der Box)
				double x = r2_abs[job][zeit][teilchen][0];
				double y = r2_abs[job][zeit][teilchen][1];
				double dx = x - Sgamma[job][zeit][0];
				double dy = y - Sgamma[job][zeit][1];
				
				double r = sqrt(dx*dx + dy*dy);
				
				int bin = (int) (r/rhoring_dr);
				if(0 > bin || bin >= rhoring_bins) continue;
				
				rho2ring[job][zeit][bin] += vorfaktor2/(2*bin+1);
				
			}//for Teilchen Typ 2
		}//for job
	}//for zeit
	
}//void berechne_dichteprofile


//Dichteprofile auswerten, in Dateien schreiben
void auswerten_dichteprofile(){
	using namespace Observablen;
	
	//lokale Felder: Mittelwerte und Fehler von Typ1- und Typ2-Dichteprofilen
	vector<double>* rho1_mittel = new vector<double>[obs_anzahl];
	vector<double>* rho1_fehler = new vector<double>[obs_anzahl];
	vector<double>* rho2_mittel = new vector<double>[obs_anzahl];
	vector<double>* rho2_fehler = new vector<double>[obs_anzahl];
	
	//berechne Mittelwerte und Fehler von rho1
	statistik_1(rho1ring, rho1_mittel, rho1_fehler, obs_anzahl);
	statistik_1(rho2ring, rho2_mittel, rho2_fehler, obs_anzahl);
	
	//schreibe in Datei
	FILE* out = fopen("dichteprofil.txt", "w");
	fprintf(out, "#Format: r TAB t TAB rho1 TAB fehler TAB rho2 TAB fehler \n\n");
	
	
	for(int zeit=0; zeit<obs_anzahl; zeit++){
		double t = zeit*obs_dt;
		
		for(int r_bin=0; r_bin<rhoring_bins; r_bin++){
			double r = (r_bin+0.5)*rhoring_dr;
	
			double r1m = rho1_mittel[zeit][r_bin]; //rho1, Mittelwert
			double r1f = rho1_fehler[zeit][r_bin]; //rho1, Fehler
			double r2m = rho2_mittel[zeit][r_bin];
			double r2f = rho2_fehler[zeit][r_bin];
			fprintf(out, "%g \t %g \t %g \t %g \t %g \t %g \n", r, t, r1m, r1f, r2m, r2f);
		}//for zeit
		fprintf(out, "\n");
	}//for zeit
	
	fclose(out);
}//void auswerten_dichteprofile


//berechnet Mittelwerte und Varianzen von input[][][]. Mittelung über den ersten Index, der bis runs geht. Der zweite geht bis anzahl. mittelwerte und varianzen müssen schon (leer) initialisiert sein.
void statistik_1(vector<double>**& input, vector<double>*& mittelwerte, vector<double>*& varianzen, int anzahl){
	int ar; //aktueller run
	int i; // aktuelle Stelle im vector, vorletzter index
	int j; // aktuelle Stelle im vector, letzter Index
	
	
	/*
	//erster Index Zeit in Einheiten obs_dt, zweiter Index r in Einheiten korr_dr
	if(mittelwerte == NULL){
		mittelwerte = new vector<double>[anzahl];
		varianzen = new vector<double>[anzahl];
	}//if mittelwerte==NULL
	else cout << "Warnung: statistik_1() aufgerufen, Ergebnisvektoren schon benutzt!" << std::endl;
	// */
	
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
