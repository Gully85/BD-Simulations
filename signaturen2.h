#pragma once

#include <vector>
#include <string>
#include <fftw3.h>
#include <iostream>
#include <sstream>



extern const double nachList_Breite, L;
extern const int densGrid_Zellen;

using std::vector; using std::string; using std::min;
using std::cout; using std::endl; using std::flush;

//Zahl zu String, keine führenden Nullen
string int_to_string(int zahl);

//suche Zahl im vector<int>. Falls nicht drin, wird -1 zurückgegeben
int suche(int a, vector<int> v);

//fftw_complex beispiel;

//enthält alle Informationen/Variablen/Felder, um den Zustand eines Runs zu definieren
class RunZustand{

public:	
	//enthält alle Informationen/Variablen/Felder, um Zeitschritte durchzuführen
	class RunDynamik{
	public:
		RunDynamik(){
			r1_git = NULL;
			r2_git = NULL;
			r1_rel = NULL;
			r2_rel = NULL;
			
			rhox = NULL;
			rhok = NULL;
			sinxG = NULL;
			sinyG = NULL;
			
			Fxk = NULL;
			Fyk = NULL;
			Fx = NULL;
			Fy = NULL;
			
			erwNachbarn1 = NULL;
			erwNachbarn2 = NULL;
			
			F1kap = NULL;
			F2kap = NULL;
			F1_WCA = NULL;
			F2_WCA = NULL;
			F1_noise = NULL;
			F2_noise = NULL;
			
			forward_plan = NULL;
			backx_plan = NULL;
			backy_plan = NULL;
		} 
		virtual ~RunDynamik()=default;
		
		//Ein Zeitschritt, Länge so dass max_reisedistanz eingehalten wird, höchstens tmax
		double zeitschritt(double tmax);
		
		
		
		//Methoden während eines Runs	
		//Initialisiert Nachbarlisten mit vorhandenen Positionen, plant Fouriertrafos
		void init();
		
		
		void berechne_kapkraefte();
		void berechne_WCA11();
		void berechne_WCA22();
		void addiere_WCA12();
		
		void refresh_erwNachbar();
		
		void berechne_zufallskraefte();
		double optimaler_zeitschritt();
		
		
		//Zeiger auf die Positionen. Soll auf dasselbe Feld zeigen wie RunZustand::r1_git usw.
		int** r1_git;
		int** r2_git;
		double** r1_rel;
		double** r2_rel;
		
	private:
		
		//Initialisierung eines Runs erfordert diese
		void WCA_init(); //Speicher reservieren, Nachbarlisten bauen
		void kapkraefte_init(); //Speicher reservieren, FFTs planen
		
		//Dichte nach Gridding, Ortsraum. Index-Wrapping
		fftw_complex* rhox;
		
		//dasselbe fouriertransformiert. Index-Wrapping, Verschiebung, alternierende Vorzeichen.
		fftw_complex* rhok;
		
		//Konstantes Feld sin(qx dx)/dx G(q). Index-Wrapping, Verschiebung.
		double* sinxG;
		double* sinyG;
		
		//Produkt rhok mal sinxG bzw mal sinyG. Index-Wrapping, Verschiebung, alternierende Vorzeichen.
		fftw_complex* Fxk;
		fftw_complex* Fyk;
		
		//rücktransformierte Kräfte. Index-Wrapping.
		fftw_complex* Fx;
		fftw_complex* Fy;
		
		//Nachbarlisten, je Typ eine. Erster Index ZellenNr in x-Richtung, zweiter Index in y-Richtung. Im vector stehen die Indices aller Teilchen in gleicher oder benachbarter Zelle.
		vector<int>** erwNachbarn1;
		vector<int>** erwNachbarn2;
		
		//Kapillarkräfte, je Typ ein Feld. Erster Index TeilchenNr, zweiter Raumrichtung.
		double** F1kap;
		double** F2kap;
		
		//WCA-Kräfte auf Typ 1. Erster Index TeilchenNr, zweiter Raumrichtung.
		double** F1_WCA;
		double** F2_WCA;
		
		//Zufallskräfte auf Typ 1. Erster Index TeilchenNr, zweiter Raumrichtung.
		double** F1_noise;
		double** F2_noise;
		
		//FFTW-Pläne für Berechnung der Kapillarkräfte
		fftw_plan forward_plan;
		fftw_plan backx_plan;
		fftw_plan backy_plan;
		
		//Gridding-Methoden innerhalb der Kapillarkraftberechnung
		void gridDensity_NGP();
		void gridDensity_CIC();
		void gridDensity_TSC();
		void inv_gridDensity_NGP();
		void inv_gridDensity_CIC();
		void inv_gridDensity_TSC();
		
		//berechnet das Quadrat des Abstands zwischen Teilchen i(Typ 1) und Teilchen j(Typ 1). Berücksichtigt periodische Randbedingungen.
		double abstand2_11(int i, int j);

		//berechnet das Quadrat des Abstands zwischen Teilchen i(Typ 1) und Teilchen j(Typ 2). Berücksichtigt periodische Randbedingungen.
		double abstand2_12(int i, int j);
			
		//berechnet das Quadrat des Abstands zwischen Teilchen i(Typ 2) und Teilchen j(Typ 2). Berücksichtigt periodische Randbedingungen.
		double abstand2_22(int i, int j);
		
	};//class RunDynamik

	RunDynamik dyn;
	
	
	//alles nötige für Observablen und Schreiben von Ergebnissen, die den einzelnen Run betreffen. 
	class RunObs{
	public:
		RunObs()=default;
		virtual ~RunObs()=default;
		
		int** r1_git;
		int** r2_git;
		double** r1_rel;
		double** r2_rel;
		
		void obs_init(int run); 
		
		//diese Methoden schreiben je eine Datei, mit den Ergebnissen (eine Observable) dieses Runs
		void schreibe_korrfunk(int run); //schreibt nur g11 und g22
		void schreibe_korrfunk_alle(int run); //schreibt g11, g22, g12
		void schreibe_ftrho(int run); //schreibt ftrho, falls N2!=0 für beide Typen, sonst nur für eins
		void schreibe_rhoFFT(int run); 
		void normalize(); //die Observablen, die mittels Binning gewonnen wurden (zB ftrho, rhoFFTW), werden je Bin durch (Anzahl Beiträge zum Bin) geteilt
		void schreibe_pos(int t);
		
		
		void record_korrelationsfunktion11(int t);
		void record_korrelationsfunktion12(int t);
		void record_korrelationsfunktion22(int t);
		void record_ftrho_unkorrigiert(int t);
		void record_rhoFFTW(int t);
		
		
		
		//TODO eine Methode, um am Ende des Runs alle Ergebnisse in den globalen Vektor zu kopieren
		
	private:
		//Positionen der Teilchen. Erster Index Teilchen, zweiter Raumrichtung. Zeigt auf dieselbe Variable wie RunZustand::Positionen

		
		//Paarkorrelationsfunktion. Erster Index Zeit in Einheiten obs_dt, zweiter Index r in Einheiten korr_dr
		vector<double>* g11;
		vector<double>* g12;
		vector<double>* g22;
		
	protected:
		
		void init_korrelationsfunktion();
		void init_ftrho();
		void init_rhoFFTW();
		
		void gridDensity_NGP1();
		void gridDensity_CIC1();
		void gridDensity_TSC1();
		void gridDensity_NGP2();
		void gridDensity_CIC2();
		void gridDensity_TSC2();
		


		//Dichte im k-Raum, Realteil und Imaginärteil. Erster Index Zeit/obs_dt, zweiter Index q (nachschlagen in order[] und qabs2[])
		vector<double>* ftrho1_re;
		vector<double>* ftrho1_im;
		vector<double>* ftrho2_re;
		vector<double>* ftrho2_im;

		//Hilfsfelder für Dichte im k-Raum
		vector<int> ftrho_order;
		vector<int> ftrho_qabs2;
		vector<int> ftrho_beitraege;



		//Dichte nach Density-Gridding für rho(k) via FFTW
		fftw_complex* rhox; //Ortsraum. Index-Wrapping.
		fftw_complex* rhok; //k-Raum. Index-Wrapping, Verschiebung, alternierende Vorzeichen.

		//Fouriertransformierte Dichte, berechnet via FFTW. Erster Index obs/dt, zweiter Index q/dq_rhoFFTW. Kein Index-Wrapping, Verschiebung oder Vorzeichen. Binning.
		vector<double>* rho1FFTW_re;
		vector<double>* rho1FFTW_im;
		vector<double>* rho2FFTW_re;
		vector<double>* rho2FFTW_im;
		
		//zaehlt für jeden Bin mit, wie viele k drinliegen
		int* rhoFFTW_beitraege;
		
		//speichert je run das rhoFFTW zur Zeit t=0.
		fftw_plan rhoFFTW_plan;
		
		
		
		//Dateizeiger, in den Positionen der Teilchen je Obs-Point geschrieben werden
		FILE* pos1;
		FILE* pos2;
		

		//berechnet das Quadrat des Abstands zwischen Teilchen i(Typ 1) und Teilchen j(Typ 1). Berücksichtigt periodische Randbedingungen.
		double abstand2_11(int i, int j);

		//berechnet das Quadrat des Abstands zwischen Teilchen i(Typ 1) und Teilchen j(Typ 2). Berücksichtigt periodische Randbedingungen.
		double abstand2_12(int i, int j);

		//berechnet das Quadrat des Abstands zwischen Teilchen i(Typ 2) und Teilchen j(Typ 2). Berücksichtigt periodische Randbedingungen.
		double abstand2_22(int i, int j);
		
				
	};//class RunObs
	
	RunObs obs;
	
	
public:
	//Konstruktor: init() aufrufen, Ortsvariablen in der RunDynamik und in RunObs auf die Ortsvariablen von RunZustand zeigen lassen
	RunZustand()=default;
	virtual ~RunZustand()=default;
	
	//es werden const int runs viele Runs gestartet. Nummer dieses Runs.
	int nr;
	double t;
	int obs_nr;
	int schritte_seit_obs;
	
	//Positionen der Teilchen
	int** r1_git;
	int** r2_git;
	double** r1_rel;
	double** r2_rel;
	
	//initialisiere den Run
	void init(int nr);
	

	void ausgabe_jeansgroessen();
	
	void zeitschritte_bis_obs();
	void obs_point(int nr_obspoint); //alle Observablen ausrechnen, speichern
	
	//Schreibe Ergebnisse in Dateien
	void schreibe_obs();
	
	
private:
	
	//liest Startpositionen aus Datei. Es muss zuerst Typ 1 und dann Typ 2 stehen.
	void init_pos_aus_datei();
	
	//initialisiert Teilchenpositionen auf Zufallspositionen
	void init_zufallspos();
	
	//bestimmt zufällig EINE Position, setzt ALLE Teilchen dorthin
	void init_allegleich();
	
	//setzt alle Teilchen zufällig gleichverteilt auf eine Kreisscheibe in der Mitte der Box
	void init_kreisscheibe();
	
	//setzt Teilchen auf ein quadratisches Gitter. Nur ein Teilchentyp, sonst Abbruch. Falls N keine Quadratzahl ist, bleiben die hintersten paar Gitterplätze frei.
	void init_gitterstart();
	

	//berechnet das Quadrat des Abstands zwischen Teilchen i(Typ 1) und Teilchen j(Typ 1). Berücksichtigt periodische Randbedingungen.
	double abstand2_11(int i, int j);

	//berechnet das Quadrat des Abstands zwischen Teilchen i(Typ 1) und Teilchen j(Typ 2). Berücksichtigt periodische Randbedingungen.
	double abstand2_12(int i, int j);

	//berechnet das Quadrat des Abstands zwischen Teilchen i(Typ 2) und Teilchen j(Typ 2). Berücksichtigt periodische Randbedingungen.
	double abstand2_22(int i, int j);
	

}; //class RunZustand








//index-Wrapping, ermoeglicht es 1dim-Arrays mit zwei Indices anzusprechen. FFTW erwartet Arrays mit nur einem Index. Keine period. Randbed.
int iw(int i, int j);

//entfernt einen Eintrag aus übergebener Nachbarliste
void erwListe_rem(vector<int>& liste, int k);

//Greensfunktion, Index-Wrapping
double G(int, int);

//Random Seed aus aktueller time()
void init_rng();

//gleichverteilte Zufallszahl aus [0.0, 1.0]
double zufall_gleichverteilt();

//gleichverteilte Zufallszahl aus [a,b]
double zufall_gleichverteilt_vonbis(double min, double max);

string int_to_string(int);