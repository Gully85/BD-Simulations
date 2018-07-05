#include "parameter.h"
#include "signaturen2.h"


#include <cstdlib>
#define _USE_MATH_DEFINES
#include <math.h>
//#include "gridRoutinen.cpp"
#include <fftw3.h>
//#include "dynamik_methoden.cpp"
#include <algorithm>
#include <iostream>
#include <omp.h>

//#include "auswertung.cpp"
//#include "positionen_speichernladen.cpp"


using std::cout; using std::endl; using std::flush; using std::vector; using std::min; using std::string;


// importiere Variablen aus parameter.h, siehe dort, was die Variable tut.
extern const int N1;
extern const double L;

extern const int densGrid_Schema, densGrid_Zellen;
extern const double densGrid_Breite;

extern const int nachList_Zellen;
extern const double nachList_Breite;

extern const double lambda_kapillar;
extern const double kapillar_vorfaktor;
extern const double T;

extern const double dq_rhoFFTW;
extern const int rhoFFTW_bins;

extern const int obs_anzahl;
extern const double obs_dt;

extern const int runs;
extern const int maxThreads;

extern const bool auswerten_korrfunk;
extern const bool auswerten_korrfunk_mixed;
extern const bool auswerten_rhovonk;
extern const bool auswerten_rhoviaFFTW;
extern const bool auswerten_rhoFT_normjerun;
extern const bool auswerten_animation;
extern const bool debugmode;

extern const double max_reisedistanz;

/// Globale Felder. Mit Null initialisiert.
/*
// Teilchenpositionen Typ 1, Gittervektor. Erster Index TeilchenNr, zweiter Index 0 fuer x und 1 fuer y-Richtung
int** r1_git = NULL;
// dasgleiche, Vektor innerhalb der Zelle
double** r1_rel = NULL;

//dasgleiche für Typ 2
int** r2_git = NULL;
double** r2_rel = NULL;
*/


int main(){

	init_rng(); //Random Seed
	//fftw_init_threads();

	///Felder, in die die Ergebnisse der Runs geschrieben werden
	//Paarkorrelationsfunktion. Erster Index Durchlauf, zweiter Index Zeit in Einheiten obs_dt, dritter Index r in Einheiten korr_dr
	vector<double>** g11 = NULL;
	vector<double>** g12 = NULL;
	vector<double>** g22 = NULL;
	
	//fouriertransformierte Dichte, Realteil und Imaginärteil. Erster Index Durchlauf, zweiter Index Zeit/obs_dt, dritter Index q (nachschlagen in order[] und qabs2[])
	vector<double>** ftrho1_re = NULL;
	vector<double>** ftrho1_im = NULL;
	vector<double>** ftrho2_re = NULL;
	vector<double>** ftrho2_im = NULL;
	
	//Fouriertransformierte Dichte, berechnet via FFTW. Erster Index Durchlauf, zweiter Index obs/dt, dritter Index q/dq_rhoFFTW. Kein Index-Wrapping, Verschiebung oder Vorzeichen. Binning, q/dq_rhoFFTW läuft von 0 bis rhoFFTW_bins-1
	vector<double>** rho1FFTW_re = NULL;
	vector<double>** rho1FFTW_im = NULL;
	vector<double>** rho2FFTW_re = NULL;
	vector<double>** rho2FFTW_im = NULL;
	
	{ //reserviere Speicher für diese Felder
	g11 = new vector<double>*[runs];
	g12 = new vector<double>*[runs];
	g22 = new vector<double>*[runs];
	ftrho1_re = new vector<double>*[runs];
	ftrho1_im = new vector<double>*[runs];
	ftrho2_re = new vector<double>*[runs];
	ftrho2_im = new vector<double>*[runs];
	rho1FFTW_re = new vector<double>*[runs];
	rho1FFTW_im = new vector<double>*[runs];
	rho2FFTW_re = new vector<double>*[runs];
	rho2FFTW_im = new vector<double>*[runs];
	
	for(int run=0; run<runs; run++){
		g11[run] = new vector<double>[obs_anzahl];
		g12[run] = new vector<double>[obs_anzahl];
		g22[run] = new vector<double>[obs_anzahl];
		ftrho1_re[run] = new vector<double>[obs_anzahl];
		ftrho1_im[run] = new vector<double>[obs_anzahl];
		ftrho2_re[run] = new vector<double>[obs_anzahl];
		ftrho2_im[run] = new vector<double>[obs_anzahl];
		rho1FFTW_re[run] = new vector<double>[obs_anzahl];
		rho1FFTW_im[run] = new vector<double>[obs_anzahl];
		rho2FFTW_re[run] = new vector<double>[obs_anzahl];
		rho2FFTW_im[run] = new vector<double>[obs_anzahl];
	}//for runs
	}

//if(maxThreads != 0)
//    omp_set_num_threads(maxThreads);

//#pragma omp parallel for
for(int run=0; run<runs; run++){
	
    RunZustand theRun;
    //#pragma omp critical
    {
	
	
	theRun.ausgabe_jeansgroessen();
	//main_init();
	theRun.init(run);
	
	
	cout << "Run nr " << run << ": initialisiert." << endl;
    }
        
	for(int count=0; count < 1e4; count++)
            theRun.dyn.zeitschritt(1.0);

}//for run

/*
//Statistik, und Ergebnis in Datei schreiben
if(auswerten_korrfunk) auswerten_korrelationsfunktion();
if(auswerten_korrfunk_mixed) auswerten_korrelationsfunktion_mixed(); //HIER WEITER! 
if(auswerten_rhovonk) auswerten_ftrho();
if(auswerten_rhoviaFFTW) auswerten_rhoFFTW();
if(auswerten_rhoFT_normjerun) auswerten_rhoFFTW_normjerun();

pos_schreiben_einedatei(); //schreibe Positionen am Ende des letzten Runs in Datei pos.txt
*/

//mittele alle Observablen, schreibe Mittelwerte und Fehler in Dateien
//alle_auswerten(); //TODO


return 0;

}//int main




