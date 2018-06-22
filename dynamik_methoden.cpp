// stellt Methoden bereit, die die Bewegung der Kolloide behandeln.
// ausserdem Methoden zur Initialisierung der Dynamik



#include "parameter.h"
#include "signaturen2.h"

#include <fftw3.h>
#include <math.h>
#include <cstdio>
#include <iostream>
#include <vector>

using std::cout; using std::endl; using std::flush;
using std::vector;


extern const double densGrid_Breite, dq, lambda_kapillar, L, zweihoch1_6, kapillar_vorfaktor, dt_max, max_reisedistanz, T;
extern const int densGrid_Zellen, densGrid_Schema;
extern const int N1, N2;
extern const double obs_dt;
extern const bool noWCA, noRNG, quickInit, restrictRadial;

/*
extern int** r1_git;
extern double** r1_rel;
extern int** r2_git;
extern double** r2_rel;
*/

extern const double sigma11_22, sigma11_12, Gamma2_1, eps22_11, eps12_11;
extern const double f2_f1;

extern const int startpos_methode; //Zufall=1, aus Datei=2, Gitterstart=3, allegleich=4, Kreisscheibe=5


const double dx = densGrid_Breite;
const int Z = densGrid_Zellen;

/*
/// semi-globale Felder/Variablen, erreichbar fuer die Methoden in dieser Datei. Namespace = Dateiname
namespace dynamik_methoden{
fftw_complex* rhox = NULL; // Dichte vor  Fouriertrafo FORWARD. Index-Wrapping
fftw_complex* rhok = NULL; // Dichte nach Fouriertrafo FORWARD. Index-Wrapping, Verschiebung und alternierende Vorzeichen
double* sinxG = NULL;	   // der Term sin(qx dx)/dx G(q). Index-Wrapping und Verschiebung
double* sinyG = NULL;	   // der Term sin(qy dx)/dx G(q). Index-Wrapping und Verschiebung
fftw_complex* Fxk = NULL;  // rhok mal i sinxG. Index-Wrapping, Verschiebung und alternierende Vorzeichen
fftw_complex* Fyk = NULL;  // rhok mal i sinyG. Index-Wrapping, Verschiebung und alternierende Vorzeichen
fftw_complex* Fx = NULL;   // Kraefte (x-komp) nach Fouriertrafo BACKWARD. Index-Wrapping
fftw_complex* Fy = NULL;   // Kraefte (y-komp) nach Fouriertrafo BACKWARD. Index-Wrapping
vector<int>** erwNachbarn1 = NULL; // Erweiterte Nachbarliste. Erster Index = ZellenNr in x-Richtung, zweiter=y-Richtung. Im vector stehen die Indices aller Teilchen, die in der gleichen oder benachbarten Zellen sind.
vector<int>** erwNachbarn2 = NULL; //hier stehen die Indices der Typ2-Teilchen drin

double** F1kap = NULL; // Kapillarkraefte. Erster Index TeilchenNr, zweiter Index Raumrichtung
double** F2kap = NULL;
double** F1_WCA = NULL; //WCA-Kraefte. Erster Index TeilchenNr, zweiter Index Raumrichtung
double** F2_WCA = NULL;
double** F1_noise = NULL; //Zufallskraefte. Erster Index TeilchenNr, zweiter Index Raumrichtung
double** F2_noise = NULL;

const int Z = densGrid_Zellen;

fftw_plan forward_plan=NULL;
fftw_plan backx_plan=NULL;
fftw_plan backy_plan=NULL;
}//namespace dynMeth
*/

//Fuehre einen Zeitschritt durch: Berechne Kraefte, ermittle optimale Dauer, bewege Teilchen, aktualisiere ggf Nachbarlisten. Gibt Dauer zurueck.
double RunZustand::RunDynamik::zeitschritt(double tmax){ 
	
	//using namespace dynamik_methoden;
	
	double deltat=tmax; //tatsaechliche Dauer. Koennte kleiner werden als tmax
	
	if(tmax <= 0.0){
		cout << "Fehler: Negative Zeitschrittweite " << tmax << " gefordert, breche ab!" << endl << flush;
		return 1.0e10;
	}//if tmax negativ, alles abbrechen
	
	//Berechne alle benoetigten Kraefte
	berechne_WCA11(); //schreibt in F1_WCA
	berechne_WCA22();
	addiere_WCA12();
	berechne_kapkraefte();
        if(!noWCA)
            berechne_zufallskraefte();
        else{
            for(int i=0; i<N1; i++){
                F1_noise[i][0] = 0.0;
                F1_noise[i][1] = 0.0;
            }//for i
            for(int i=0; i<N2; i++){
                F2_noise[i][0] = 0.0;
                F2_noise[i][1] = 0.0;
            }//for i
        }//else ohne Zufallskräfte
	
	/* Test: Gebe Teilchenpositionen und aktuelle Kraefte aus. 
	for(int teilchen=0; teilchen<N; teilchen++){
		cout << teilchen << " an (" << r_git[teilchen][0]*nachList_Breite+r_rel[teilchen][0] << "," << r_git[teilchen][1]*nachList_Breite+r_rel[teilchen][1] << "):\tWCA=(" << F_WCA[teilchen][0] << "," << F_WCA[teilchen][1]<< "),\tkap=(" << Fkap[teilchen][0]<<","<<Fkap[teilchen][1]<<"), noise=("<<F_noise[teilchen][0]<<","<<F_noise[teilchen][1]<<")" << endl;
	}//for teilchen
	// */
	
	//bestimme optimalen Zeitschritt, so dass kein Teilchen weiter als max_reisedistanz bewegt wird
	//deltat = optimaler_zeitschritt(F1_WCA, F1kap, F1_noise, deltat, N1);
	//deltat = optimaler_zeitschritt(F2_WCA, F2kap, F2_noise, deltat, N2);
	deltat = optimaler_zeitschritt();
	
// 	cout << endl << "Zeitschritt: " << deltat << endl;
	
	//alle Teilchen Typ 1 bewegen
	for(int teilchen=0; teilchen<N1; teilchen++){
		double dx = deltat * (F1_WCA[teilchen][0] + F1kap[teilchen][0]) + sqrt(2.0*T*deltat)*F1_noise[teilchen][0];
		double dy = deltat * (F1_WCA[teilchen][1] + F1kap[teilchen][1]) + sqrt(2.0*T*deltat)*F1_noise[teilchen][1];
		
		//* Test: dr ausgeben
// 		cout << "Teilchen " << teilchen << ": dr=("<<dx<<","<<dy<<")"<< endl;
		// */
		
		r1_rel[teilchen][0] += dx;
		r1_rel[teilchen][1] += dy;
		
	}//for teilchen bis N1
	
	//Typ 2 bewegen
	for(int teilchen=0; teilchen<N2; teilchen++){
		double dx = deltat * (F2_WCA[teilchen][0] + F2kap[teilchen][0]) + sqrt(2.0*T*deltat)*F2_noise[teilchen][0];
		double dy = deltat * (F2_WCA[teilchen][1] + F2kap[teilchen][1]) + sqrt(2.0*T*deltat)*F2_noise[teilchen][1];
		
		r2_rel[teilchen][0] += dx;
		r2_rel[teilchen][1] += dy;
                
                if(restrictRadial){
                    //messe abstand zur mitte
                    const double mitte = 0.5*L;
                    double x = r2_git[teilchen][0]*nachList_Breite + r2_rel[teilchen][0];
                    double y = r2_git[teilchen][1]*nachList_Breite + r2_rel[teilchen][1];
                    
                    double dxmitte = mitte - x;
                    double dymitte = mitte - y;
                    double abst = sqrt(dxmitte*dxmitte + dymitte*dymitte);
                    
                    
                    //aus Abstand und teilchenNr ergibt 
                    double x_soll = mitte + abst * cos(2*M_PI*(int)teilchen/N2);
                    double y_soll = mitte + abst * sin(2*M_PI*(int)teilchen/N2);
                    
                    r2_rel[teilchen][0] += (x_soll - x);
                    r2_rel[teilchen][1] += (y_soll - y);
                }//if restrictRadial
                
	}//for teilchen bis N2
	
	//erwListen aktualisieren: Falls ein Teilchen seine Zelle verlassen hat (dh r_rel<0 oder r_rel>nachList_Breite), streichen/hinzufügen
	refresh_erwNachbar();
	
	return deltat;
}//double zeitschritt


// Fuehre einen Zeitschritt durch: Berechnet Kräfte, ermittle optimale Dauer, bewege Teilchen, aktualisiere Nachbarlisten. 
// Schreibe zusätzliche Informationen in TimestepInfo
double RunZustand::RunDynamik::zeitschritt_debug(double tmax, TimestepInfo &info){
    double deltat = tmax; //kann nur kleiner werden
    
    if (tmax <= 0.0){
        cout << "Fehler: Negative Zeitschrittweite "<<tmax<<" gefordert, breche ab!" << endl << flush;
        return 1.0e10;
    }//if tmax negativ
    
    info.reset();
    
    //berechne alle Kräfte
    berechne_WCA11_debug(info);
    berechne_WCA22_debug(info);
    addiere_WCA12_debug(info);
    berechne_kapkraefte_debug(info);
    
    berechne_zufallskraefte_debug(info);
    
    deltat = optimaler_zeitschritt();
    info.dt = deltat;
    
    //alle Teilchen Typ 1 bewegen
    for(int teilchen=0; teilchen<N1; teilchen++){
        double dx = deltat * (F1_WCA[teilchen][0] + F1kap[teilchen][0]) + sqrt(2.0*T*deltat)*F1_noise[teilchen][0];
        double dy = deltat * (F1_WCA[teilchen][1] + F1kap[teilchen][1]) + sqrt(2.0*T*deltat)*F1_noise[teilchen][1];
        
        r1_rel[teilchen][0] += dx;
        r1_rel[teilchen][1] += dy;
    }//for teilchen bis N1
    
    // alle Teilchen Typ 2 bewegen
    for(int teilchen=0; teilchen<N2; teilchen++){
        double dx = deltat * (F2_WCA[teilchen][0] + F2kap[teilchen][0]) + sqrt(2.0*T*deltat)*F2_noise[teilchen][0];
        double dy = deltat * (F2_WCA[teilchen][1] + F2kap[teilchen][1]) + sqrt(2.0*T*deltat)*F2_noise[teilchen][1];
        
        r2_rel[teilchen][0] += dx;
        r2_rel[teilchen][1] += dy;
    }//for teilchen bis N1
    
    refresh_erwNachbar_debug(info);
    
    info.ausgabe();
    
    return deltat;
    
    
    
    
}//double zeitschritt_debug



//erneuere erwNachbarlisten. Beim Aufruf sind r_rel<0 oder r_rel>nachList_Breite zulässig, wird behoben
void RunZustand::RunDynamik::refresh_erwNachbar(){
	// using namespace dynamik_methoden;
	//Typ 1 über Zellgrenze bewegt?
	for(int teilchen=0; teilchen<N1; teilchen++){
		if(r1_rel[teilchen][0] < 0.0 || r1_rel[teilchen][0] > nachList_Breite
		|| r1_rel[teilchen][1] < 0.0 || r1_rel[teilchen][1] > nachList_Breite){
			//Indices der alten Zelle
			int ic = r1_git[teilchen][0];
			int jc = r1_git[teilchen][1];
			
			//Zellen, in denen das Teilchen bisher eingetragen war
			int x[3]; int y[3];
			x[0] = (ic-1 + nachList_Zellen)%nachList_Zellen;
			x[1] = ic;
			x[2] = (ic+1)%nachList_Zellen;
			y[0] = (jc-1 + nachList_Zellen)%nachList_Zellen;
			y[1] = jc;
			y[2] = (jc+1)%nachList_Zellen;
			
			//aus allen entfernen
			for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				erwListe_rem(erwNachbarn1[x[i]][y[j]], teilchen);
			
			//Indices der neuen Zelle. Koennen noch negative oder zu grosse Werte haben. floor() rundet auch negative Zahlen korrekt ab.
			ic = (int) floor((r1_git[teilchen][0]*nachList_Breite + r1_rel[teilchen][0])/nachList_Breite);
			jc = (int) floor((r1_git[teilchen][1]*nachList_Breite + r1_rel[teilchen][1])/nachList_Breite);
			
			// korrigiere negative bzw zu grosse Zellindices, durch period. Randbed.
			ic = (ic+nachList_Zellen)%nachList_Zellen;
			jc = (jc+nachList_Zellen)%nachList_Zellen;
			
			//trage neue Position in Gitter- und Relativvektor ein
			r1_rel[teilchen][0] = fmod(r1_rel[teilchen][0] + nachList_Breite, nachList_Breite); //fmod=modulo
			r1_rel[teilchen][1] = fmod(r1_rel[teilchen][1] + nachList_Breite, nachList_Breite);
			r1_git[teilchen][0] = ic;
			r1_git[teilchen][1] = jc;
			
			//fuege Teilchen den neuen Nachbarlisten hinzu
			//Zellen, in die das Teilchen eingetragen werden soll. ic,jc sind schon aktuell
			x[0] = (ic-1+nachList_Zellen)%nachList_Zellen;
			x[1] = ic;
			x[2] = (ic+1)%nachList_Zellen;
			y[0] = (jc-1+nachList_Zellen)%nachList_Zellen;
			y[1] = jc;
			y[2] = (jc+1)%nachList_Zellen;
			
			for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				erwNachbarn1[x[i]][y[j]].push_back(teilchen);
			
			
		}//if Typ1 aus Zelle bewegt
	}//for teilchen bis N1
	
	
	//Typ 2 über Zellgrenze bewegt?
	for(int teilchen=0; teilchen<N2; teilchen++){
		if(r2_rel[teilchen][0] < 0.0 || r2_rel[teilchen][0] > nachList_Breite
		|| r2_rel[teilchen][1] < 0.0 || r2_rel[teilchen][1] > nachList_Breite){
			//Indices der alten Zelle
			int ic = r2_git[teilchen][0];
			int jc = r2_git[teilchen][1];
			
			//Zellen, in denen das Teilchen bisher eingetragen war
			int x[3]; int y[3];
			x[0] = (ic-1 + nachList_Zellen)%nachList_Zellen;
			x[1] = ic;
			x[2] = (ic+1)%nachList_Zellen;
			y[0] = (jc-1 + nachList_Zellen)%nachList_Zellen;
			y[1] = jc;
			y[2] = (jc+1)%nachList_Zellen;
			
			//aus allen entfernen
			for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				erwListe_rem(erwNachbarn2[x[i]][y[j]], teilchen);
			
			//Indices der neuen Zelle. Koennen noch negative oder zu grosse Werte haben. floor() rundet auch negative Zahlen korrekt ab.
			ic = (int) floor((r2_git[teilchen][0]*nachList_Breite + r2_rel[teilchen][0])/nachList_Breite);
			jc = (int) floor((r2_git[teilchen][1]*nachList_Breite + r2_rel[teilchen][1])/nachList_Breite);
			
			// korrigiere negative bzw zu grosse Zellindices, durch period. Randbed.
			ic = (ic+nachList_Zellen)%nachList_Zellen;
			jc = (jc+nachList_Zellen)%nachList_Zellen;
			
			//trage neue Position in Gitter- und Relativvektor ein
			r2_rel[teilchen][0] = fmod(r2_rel[teilchen][0] + nachList_Breite, nachList_Breite); //fmod=modulo
			r2_rel[teilchen][1] = fmod(r2_rel[teilchen][1] + nachList_Breite, nachList_Breite);
			r2_git[teilchen][0] = ic;
			r2_git[teilchen][1] = jc;
			
			//fuege Teilchen den neuen Nachbarlisten hinzu
			//Zellen, in die das Teilchen eingetragen werden soll. ic,jc sind schon aktuell
			x[0] = (ic-1+nachList_Zellen)%nachList_Zellen;
			x[1] = ic;
			x[2] = (ic+1)%nachList_Zellen;
			y[0] = (jc-1+nachList_Zellen)%nachList_Zellen;
			y[1] = jc;
			y[2] = (jc+1)%nachList_Zellen;
			
			for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				erwNachbarn2[x[i]][y[j]].push_back(teilchen);
			
			
		}//if aus Zelle bewegt
	}//for teilchen bis N2
	
}//void refresh_erwNachbar

void RunZustand::RunDynamik::refresh_erwNachbar_debug(TimestepInfo &info){
    // using namespace dynamik_methoden;
	//Typ 1 über Zellgrenze bewegt?
	for(int teilchen=0; teilchen<N1; teilchen++){
		if(r1_rel[teilchen][0] < 0.0 || r1_rel[teilchen][0] > nachList_Breite
		|| r1_rel[teilchen][1] < 0.0 || r1_rel[teilchen][1] > nachList_Breite){
                    
                    info.num_Zellwechsel1++;
                    
			//Indices der alten Zelle
			int ic = r1_git[teilchen][0];
			int jc = r1_git[teilchen][1];
			
			//Zellen, in denen das Teilchen bisher eingetragen war
			int x[3]; int y[3];
			x[0] = (ic-1 + nachList_Zellen)%nachList_Zellen;
			x[1] = ic;
			x[2] = (ic+1)%nachList_Zellen;
			y[0] = (jc-1 + nachList_Zellen)%nachList_Zellen;
			y[1] = jc;
			y[2] = (jc+1)%nachList_Zellen;
			
			//aus allen entfernen
			for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				erwListe_rem(erwNachbarn1[x[i]][y[j]], teilchen);
			
			//Indices der neuen Zelle. Koennen noch negative oder zu grosse Werte haben. floor() rundet auch negative Zahlen korrekt ab.
			ic = (int) floor((r1_git[teilchen][0]*nachList_Breite + r1_rel[teilchen][0])/nachList_Breite);
			jc = (int) floor((r1_git[teilchen][1]*nachList_Breite + r1_rel[teilchen][1])/nachList_Breite);
			
			// korrigiere negative bzw zu grosse Zellindices, durch period. Randbed.
			ic = (ic+nachList_Zellen)%nachList_Zellen;
			jc = (jc+nachList_Zellen)%nachList_Zellen;
			
			//trage neue Position in Gitter- und Relativvektor ein
			r1_rel[teilchen][0] = fmod(r1_rel[teilchen][0] + nachList_Breite, nachList_Breite); //fmod=modulo
			r1_rel[teilchen][1] = fmod(r1_rel[teilchen][1] + nachList_Breite, nachList_Breite);
			r1_git[teilchen][0] = ic;
			r1_git[teilchen][1] = jc;
			
			//fuege Teilchen den neuen Nachbarlisten hinzu
			//Zellen, in die das Teilchen eingetragen werden soll. ic,jc sind schon aktuell
			x[0] = (ic-1+nachList_Zellen)%nachList_Zellen;
			x[1] = ic;
			x[2] = (ic+1)%nachList_Zellen;
			y[0] = (jc-1+nachList_Zellen)%nachList_Zellen;
			y[1] = jc;
			y[2] = (jc+1)%nachList_Zellen;
			
			for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				erwNachbarn1[x[i]][y[j]].push_back(teilchen);
			
			
		}//if Typ1 aus Zelle bewegt
	}//for teilchen bis N1
	
	
	//Typ 2 über Zellgrenze bewegt?
	for(int teilchen=0; teilchen<N2; teilchen++){
		if(r2_rel[teilchen][0] < 0.0 || r2_rel[teilchen][0] > nachList_Breite
		|| r2_rel[teilchen][1] < 0.0 || r2_rel[teilchen][1] > nachList_Breite){
                    
                    info.num_Zellwechsel2++;
                    
			//Indices der alten Zelle
			int ic = r2_git[teilchen][0];
			int jc = r2_git[teilchen][1];
			
			//Zellen, in denen das Teilchen bisher eingetragen war
			int x[3]; int y[3];
			x[0] = (ic-1 + nachList_Zellen)%nachList_Zellen;
			x[1] = ic;
			x[2] = (ic+1)%nachList_Zellen;
			y[0] = (jc-1 + nachList_Zellen)%nachList_Zellen;
			y[1] = jc;
			y[2] = (jc+1)%nachList_Zellen;
			
			//aus allen entfernen
			for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				erwListe_rem(erwNachbarn2[x[i]][y[j]], teilchen);
			
			//Indices der neuen Zelle. Koennen noch negative oder zu grosse Werte haben. floor() rundet auch negative Zahlen korrekt ab.
			ic = (int) floor((r2_git[teilchen][0]*nachList_Breite + r2_rel[teilchen][0])/nachList_Breite);
			jc = (int) floor((r2_git[teilchen][1]*nachList_Breite + r2_rel[teilchen][1])/nachList_Breite);
			
			// korrigiere negative bzw zu grosse Zellindices, durch period. Randbed.
			ic = (ic+nachList_Zellen)%nachList_Zellen;
			jc = (jc+nachList_Zellen)%nachList_Zellen;
			
			//trage neue Position in Gitter- und Relativvektor ein
			r2_rel[teilchen][0] = fmod(r2_rel[teilchen][0] + nachList_Breite, nachList_Breite); //fmod=modulo
			r2_rel[teilchen][1] = fmod(r2_rel[teilchen][1] + nachList_Breite, nachList_Breite);
			r2_git[teilchen][0] = ic;
			r2_git[teilchen][1] = jc;
			
			//fuege Teilchen den neuen Nachbarlisten hinzu
			//Zellen, in die das Teilchen eingetragen werden soll. ic,jc sind schon aktuell
			x[0] = (ic-1+nachList_Zellen)%nachList_Zellen;
			x[1] = ic;
			x[2] = (ic+1)%nachList_Zellen;
			y[0] = (jc-1+nachList_Zellen)%nachList_Zellen;
			y[1] = jc;
			y[2] = (jc+1)%nachList_Zellen;
			
			for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				erwNachbarn2[x[i]][y[j]].push_back(teilchen);
			
			
		}//if aus Zelle bewegt
	}//for teilchen bis N2
	
}//void refresh_erwNachbar_debug

void RunZustand::zeitschritte_bis_obs(){
	while(t < obs_dt){
		t += dyn.zeitschritt(obs_dt-t);
		schritte_seit_obs++;
	}//while
}//void RunDynamik::zeitschritte_bis_obs


// initialisiert Felder/Nachbarlisten fuer WCA-Kraefte. Erwartet, dass in r_git schon die Gittervektoren der Teilchenpositionen stehen.
void RunZustand::RunDynamik::WCA_init(){
	//using namespace dynamik_methoden;
	//vector<int>** erwNachbarn ist die erweiterte Nachbarliste. Erster Index = ZellenNr in x-Richtung, zweiter=y-Richtung. Im vector stehen die Indices aller Teilchen, die in der gleichen oder benachbarten Zellen sind.
	
	//WCA-Kraefte. Erster Index Teilchen, zweiter Index Raumrichtung
	F1_WCA = new double*[N1];
	F2_WCA = new double*[N2];
	
	//reserviere Speicher
	if(erwNachbarn1 == NULL){
		erwNachbarn1 = new vector<int>*[nachList_Zellen];
		erwNachbarn2 = new vector<int>*[nachList_Zellen];
		for(int i=0; i<nachList_Zellen; i++){
			erwNachbarn1[i] = new vector<int>[nachList_Zellen];
			erwNachbarn2[i] = new vector<int>[nachList_Zellen];
		}//for i
	}//if erwNachbarn==NULL
	else
		for(int i=0; i<nachList_Zellen; i++)
		for(int j=0; j<nachList_Zellen; j++){
			erwNachbarn1[i][j].clear();
			erwNachbarn2[i][j].clear();
		}//for i,j
	int* x = new int[3]; //links davon, exakt, rechts davon
	int* y = new int[3];

	//trage jedes Typ1-Teilchen in seine Zelle, und die 8 Nachbarzellen ein
	for(int teilchen=0; teilchen<N1; teilchen++){
		int k = r1_git[teilchen][0];
		int l = r1_git[teilchen][1];
		
		//ermittle beteiligte Zellen
		x[0] = (k-1+nachList_Zellen)%nachList_Zellen;
		x[1] = k;
		x[2] = (k+1)%nachList_Zellen;
		y[0] = (l-1+nachList_Zellen)%nachList_Zellen;
		y[1] = l;
		y[2] = (l+1)%nachList_Zellen;
		
		//fuege Teilchen diesen Zellen hinzu
		for(int i=0;i<3;i++)
		  for(int j=0;j<3;j++)
		    erwNachbarn1[x[i]][y[j]].push_back(teilchen);
		
	}//for teilchen
	
	//trage jedes Typ2-Teilchen in seine Zelle und Nachbarzellen ein
	for(int teilchen=0; teilchen<N2; teilchen++){
		int k = r2_git[teilchen][0];
		int l = r2_git[teilchen][1];
		
		//ermittle beteiligte Zellen
		x[0] = (k-1+nachList_Zellen)%nachList_Zellen;
		x[1] = k;
		x[2] = (k+1)%nachList_Zellen;
		y[0] = (l-1+nachList_Zellen)%nachList_Zellen;
		y[1] = l;
		y[2] = (l+1)%nachList_Zellen;
		
		//fuege Teilchen diesen Zellen hinzu
		for(int i=0;i<3;i++)
		  for(int j=0;j<3;j++)
		    erwNachbarn2[x[i]][y[j]].push_back(teilchen);

	}//for teilchen
	

} //void WCA_init

//initialisiert Felder, plant Fouriertrafos
void RunZustand::RunDynamik::kapkraefte_init(){
	//using namespace dynamik_methoden;
/// alloziere Felder
	
	//Kapillarkraefte. Erster Index TeilchenNr, zweiter Index Raumrichtung
	F1kap = new double*[N1];
	F2kap = new double*[N2];

	//Zufallskraefte, gleiche Indices
	F1_noise = new double*[N1];
	F2_noise = new double*[N2];
	for(int i=0; i<N1; i++){
		F1kap[i] = new double[2];
		F1_WCA[i]= new double[2];
		F1_noise[i] = new double[2];
	}//for i
	for(int i=0; i<N2; i++){
		F2kap[i] = new double[2];
		F2_WCA[i]= new double[2];
		F2_noise[i] = new double[2];
	}//for i

	
	//Dichte vor FFT. Index-Wrapping. rhox[iw(j,k)][0] ist die Dichte in Zelle (j,k). rhox[][1] ist der Imaginaerteil, Null.
	if(rhox == NULL){
		rhox = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Z*Z);

	//Imaginaerteil auf Null setzen. Der wird sich nie aendern.
	for(int j=0; j<Z; j++)
	for(int k=0; k<Z; k++)
		rhox[iw(j,k)][1]=0.0;

	}//if rhox==NULL

	// Dichte nach FFT. Index-Wrapping, Verschiebung, Vorzeichen. rhok[iw(j,k)][0] ist der Realteil von (-1)^(j+k) mal rhoschlange(qjk), wobei qjk = dq* ((j+Zellen/2)%Zellen - Zellen/2 ex) und analog k ey
	if(rhok == NULL)
		rhok = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Z*Z);

	if(Fxk == NULL){
		// Kraefte vor FFT-1. Index-Wrapping, Verschiebung, Vorzeichen.
		Fxk = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Z*Z);
		Fyk = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Z*Z);

		// Kraefte nach FFT-1. Index-Wrapping.
		Fx  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Z*Z);
		Fy  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Z*Z);
	}//if Fxk==NULL

	// der Term sin(qx dx)/dx G(q). Index-Wrapping, Verschiebung. Analog sinyG
	if(sinxG == NULL){
		sinxG = new double[Z*Z];
		sinyG = new double[Z*Z];

/// schreibe Werte (zB G) in die Felder, die konstant sind
		int i,j,k;
		double q;
		for(j=0; j<Z; j++)
		for(k=0; k<Z; k++){
			q = dq* ((j+(int)(Z*0.5))%Z - 0.5*Z );
			sinxG[iw(j,k)] = kapillar_vorfaktor * sin(q*dx)/dx * G(j,k)/Z/Z;
			//sinxG[iw(j,k)] = G(j,k); //liefert spaeter das Hoehenprofil u(r)

			q = dq* ((k+(int)(Z*0.5))%Z - 0.5*Z );
			sinyG[iw(j,k)] = kapillar_vorfaktor * sin(q*dx)/dx * G(j,k)/Z/Z; //liefert spaeter die y-Komponente der Kraft. Entlang der y-Achse sollte das K1 entsprechen.
		}//for j,k
	}//if sinxG==NULL
	
/// plane FFTs
	//#pragma omp critical
	{
            if(quickInit){
                forward_plan = fftw_plan_dft_2d(Z, Z, rhox, rhok, FFTW_FORWARD,  FFTW_ESTIMATE);
                backx_plan   = fftw_plan_dft_2d(Z, Z, Fxk,  Fx,   FFTW_BACKWARD, FFTW_ESTIMATE);
                backy_plan   = fftw_plan_dft_2d(Z, Z, Fyk,  Fy,   FFTW_BACKWARD, FFTW_ESTIMATE);
            }//if quickInit
            else{
                forward_plan = fftw_plan_dft_2d(Z, Z, rhox, rhok, FFTW_FORWARD,  FFTW_PATIENT);
                backx_plan   = fftw_plan_dft_2d(Z, Z, Fxk,  Fx,   FFTW_BACKWARD, FFTW_PATIENT);
                backy_plan   = fftw_plan_dft_2d(Z, Z, Fyk,  Fy,   FFTW_BACKWARD, FFTW_PATIENT);
                
            }//else nicht quickInit
	
	}
	
}//void kapkraefte_init


// berechnet WCA-Kraefte. Schreibt sie in F1_WCA[N1][2]
void RunZustand::RunDynamik::berechne_WCA11(){
	//using namespace dynamik_methoden;
	int i; //Teilchen
	vector<int>::iterator j; //wechselwirkendes Teilchen. Iteriert ueber die Nachbarliste von i's Zelle
	double a2; //Abstandsquadrat der beiden wechselwirkenden Teilchen

	if(F1_WCA == NULL){
		cout << "Fehler: WCA-Kraefteberechnung aufgerufen, aber nicht initialisiert!" << endl;
		return;
	}//if nicht initialisiert

	//setze F_WCA auf Null
	for(i=0; i<N1; i++){
		F1_WCA[i][0]=0.0;
		F1_WCA[i][1]=0.0;
	}//for
	
	//Schleife ueber Teilchen. Fuer jedes, iteriere durch seine Nachbarn und addiere Kraefte.
	for(i=0; i<N1; i++){
		int ic = r1_git[i][0];
		int jc = r1_git[i][1];
		
		for(j = erwNachbarn1[ic][jc].begin(); j!=erwNachbarn1[ic][jc].end(); j++){
			if(i == *j) continue; //ueberspringe WW mit sich selbst
			
			//Abstand berechnen, mit periodischen Randbedingungen
			a2 = abstand2_11(i, *j);
			
			//falls zu weit entfernt, WW ueberspringen
			if(a2 > zweihoch1_6*zweihoch1_6) continue;
			
			//Absolutkoordinaten
			double x1 = r1_git[ i][0]*nachList_Breite + r1_rel[ i][0];
			double y1 = r1_git[ i][1]*nachList_Breite + r1_rel[ i][1];
			double x2 = r1_git[*j][0]*nachList_Breite + r1_rel[*j][0];
			double y2 = r1_git[*j][1]*nachList_Breite + r1_rel[*j][1];
			
			double dx = x1-x2;
			double dy = y1-y2;
			
			//periodische Randbedingungen. Falls dx^2 + dy^2 == a2 ist, ist dies nicht noetig.
			if(dx*dx + dy*dy > a2){
				//x-Richtung nahe am Rand?
				if((dx-L)*(dx-L) < dx*dx)
					x2 += L;
				else if((dx+L)*(dx+L) < dx*dx)
					x2 -= L;
				dx = x1-x2;
				
				//y-Richtung nahe am Rand?
				if((dy-L)*(dy-L) < dy*dy)
					y2 += L;
				else if((dy+L)*(dy+L) < dy*dy)
					y2 -= L;
				dy = y1-y2;
				
			}//if periodische Randbedingungen
			
			// a^(-8) und a^(-14)
			double a_8 = 1.0/(a2*a2*a2*a2);
			double a_14= 1.0/(a2*a2*a2*a2*a2*a2*a2);
			
			//addiere Kraefte zu F_WCA[i]. Die 4 kommt aus dem Lennard-Jones-Potential, die 6 aus d/dr r^(-6)
			F1_WCA[i][0] += 4.0*6.0*(2.0*a_14 - a_8) * dx;
			F1_WCA[i][1] += 4.0*6.0*(2.0*a_14 - a_8) * dy;
			
			
			
		}//for j
	}//for i
}//void berechne_WCA11

void RunZustand::RunDynamik::berechne_WCA11_debug(TimestepInfo &info){
    int i; //Aufteilchen
    vector<int>::iterator j; //wechselwirkendes Teilchen
    double a2; //Abstandsquadrat
    
    if(F1_WCA == NULL){
        cout << "Fehler: WCA-Kraftberechnung aufgerufen, aber nicht initialisiert!" << endl;
        return;
    }//if nicht initialisiert
    
    //Kräfte auf Null setzen
    for(i=0; i<N1; i++){
        F1_WCA[i][0]=0.0;
        F1_WCA[i][1]=0.0;
    }//for
    
    //Schleife über Teilchen
    for(i=0; i<N1; i++){
        int ic = r1_git[i][0];
        int jc = r1_git[i][0];
        
        //Schleife über wechselwirkende Teilchen
        for(j=erwNachbarn1[ic][jc].begin(); j != erwNachbarn1[ic][jc].end(); j++){
            if(i == *j) continue; //überspringe WW mit sich selbst
            
            info.num_WCA_berechnet++;
            a2 = abstand2_11(i, *j); //Abstandsquadrat, berücksichtigt periodische Randbedingungen
            
            //falls zu weit entfernt, WW ueberspringen
            if(a2 > zweihoch1_6*zweihoch1_6) continue;
            
            //Absolutkoordinaten
            double x1 = r1_git[ i][0]*nachList_Breite + r1_rel[ i][0];
            double y1 = r1_git[ i][1]*nachList_Breite + r1_rel[ i][1];
            double x2 = r1_git[*j][0]*nachList_Breite + r1_rel[*j][0];
            double y2 = r1_git[*j][1]*nachList_Breite + r1_rel[*j][1];
            
            double dx = x1-x2;
            double dy = y1-y2;
            
            //periodische Randbedingungen. Falls dx^2 + dy^2 == a2 ist, ist dies nicht noetig.
            if(dx*dx + dy*dy > a2){
                //x-Richtung nahe am Rand?
                if((dx-L)*(dx-L) < dx*dx)
                        x2 += L;
                else if((dx+L)*(dx+L) < dx*dx)
                        x2 -= L;
                dx = x1-x2;
                
                //y-Richtung nahe am Rand?
                if((dy-L)*(dy-L) < dy*dy)
                        y2 += L;
                else if((dy+L)*(dy+L) < dy*dy)
                        y2 -= L;
                dy = y1-y2;
                    
            }//if periodische Randbedingungen
            
            // a^(-8) und a^(-14)
            double a_8 = 1.0/(a2*a2*a2*a2);
            double a_14= 1.0/(a2*a2*a2*a2*a2*a2*a2);
            
            //addiere Kraefte zu F_WCA[i]. Die 4 kommt aus dem Lennard-Jones-Potential, die 6 aus d/dr r^(-6)
            F1_WCA[i][0] += 4.0*6.0*(2.0*a_14 - a_8) * dx;
            F1_WCA[i][1] += 4.0*6.0*(2.0*a_14 - a_8) * dy;
        }//for j
    }//for i
}//void berechne_WCA11_debug


//berechne WCA-Kräfte aus 22-Wechselwirkung, schreibe sie in F2_WCA[N2][2]
void RunZustand::RunDynamik::berechne_WCA22(){
	
	//using namespace dynamik_methoden;
	int i; //Teilchen
	vector<int>::iterator j; //wechselwirkendes Teilchen. Iteriert ueber die Nachbarliste von i's Zelle
	double a2; //Abstandsquadrat der beiden wechselwirkenden Teilchen

	if(F2_WCA == NULL){
		cout << "Fehler: WCA-Kraefteberechnung aufgerufen, aber nicht initialisiert!" << endl;
		return;
	}//if nicht initialisiert

	//setze F_WCA auf Null
	for(i=0; i<N2; i++){
		F2_WCA[i][0]=0.0;
		F2_WCA[i][1]=0.0;
	}//for
	
	//Schleife ueber Teilchen. Fuer jedes, iteriere durch seine Nachbarn und addiere Kraefte.
	for(i=0; i<N2; i++){
		int ic = r2_git[i][0];
		int jc = r2_git[i][1];
		
		for(j = erwNachbarn2[ic][jc].begin(); j!=erwNachbarn2[ic][jc].end(); j++){
			if(i == *j) continue; //ueberspringe WW mit sich selbst
			
			//Abstand berechnen, mit periodischen Randbedingungen
			a2 = abstand2_22(i, *j);
			
			//falls zu weit entfernt, WW ueberspringen
			if(a2 * sigma11_22*sigma11_22 > zweihoch1_6*zweihoch1_6 ) continue;
			
			//Absolutkoordinaten
			double x1 = r2_git[ i][0]*nachList_Breite + r2_rel[ i][0];
			double y1 = r2_git[ i][1]*nachList_Breite + r2_rel[ i][1];
			double x2 = r2_git[*j][0]*nachList_Breite + r2_rel[*j][0];
			double y2 = r2_git[*j][1]*nachList_Breite + r2_rel[*j][1];
			
			double dx = x1-x2;
			double dy = y1-y2;
			
			//periodische Randbedingungen. Falls dx^2 + dy^2 == a2 ist, ist dies nicht noetig.
			if(dx*dx + dy*dy > a2){
				//x-Richtung nahe am Rand?
				if((dx-L)*(dx-L) < dx*dx)
					x2 += L;
				else if((dx+L)*(dx+L) < dx*dx)
					x2 -= L;
				dx = x1-x2;
				
				//y-Richtung nahe am Rand?
				if((dy-L)*(dy-L) < dy*dy)
					y2 += L;
				else if((dy+L)*(dy+L) < dy*dy)
					y2 -= L;
				dy = y1-y2;
				
			}//if periodische Randbedingungen
			
			a2 *= sigma11_22*sigma11_22;
			dx *= sigma11_22;
			dy *= sigma11_22;
			
			
			// a^(-8) und a^(-14)
			double a_8 = 1.0/(a2*a2*a2*a2);
			double a_14= 1.0/(a2*a2*a2*a2*a2*a2*a2);
			
			//addiere Kraefte zu F_WCA[i]. Die 4 kommt aus dem Lennard-Jones-Potential, die 6 aus d/dr r^(-6)
			F2_WCA[i][0] += eps22_11*sigma11_22*4.0*6.0*(2*a_14 - a_8)*dx;
			F2_WCA[i][1] += eps22_11*sigma11_22*4.0*6.0*(2*a_14 - a_8)*dy;
			
		}//for j
	}//for i
	
}//void berechne_WCA22

void RunZustand::RunDynamik::berechne_WCA22_debug(TimestepInfo &info){
    //using namespace dynamik_methoden;
    int i; //Teilchen
    vector<int>::iterator j; //wechselwirkendes Teilchen. Iteriert ueber die Nachbarliste von i's Zelle
    double a2; //Abstandsquadrat der beiden wechselwirkenden Teilchen

    if(F2_WCA == NULL){
            cout << "Fehler: WCA-Kraefteberechnung aufgerufen, aber nicht initialisiert!" << endl;
            return;
    }//if nicht initialisiert

    //setze F_WCA auf Null
    for(i=0; i<N2; i++){
            F2_WCA[i][0]=0.0;
            F2_WCA[i][1]=0.0;
    }//for
    
    //Schleife ueber Teilchen. Fuer jedes, iteriere durch seine Nachbarn und addiere Kraefte.
    for(i=0; i<N2; i++){
            int ic = r2_git[i][0];
            int jc = r2_git[i][1];
            
            for(j = erwNachbarn2[ic][jc].begin(); j!=erwNachbarn2[ic][jc].end(); j++){
                    if(i == *j) continue; //ueberspringe WW mit sich selbst
                    
                    info.num_WCA_berechnet++;
                    //Abstand berechnen, mit periodischen Randbedingungen
                    a2 = abstand2_22(i, *j);
                    
                    //falls zu weit entfernt, WW ueberspringen
                    if(a2 * sigma11_22*sigma11_22 > zweihoch1_6*zweihoch1_6 ) continue;
                    
                    //Absolutkoordinaten
                    double x1 = r2_git[ i][0]*nachList_Breite + r2_rel[ i][0];
                    double y1 = r2_git[ i][1]*nachList_Breite + r2_rel[ i][1];
                    double x2 = r2_git[*j][0]*nachList_Breite + r2_rel[*j][0];
                    double y2 = r2_git[*j][1]*nachList_Breite + r2_rel[*j][1];
                    
                    double dx = x1-x2;
                    double dy = y1-y2;
                    
                    //periodische Randbedingungen. Falls dx^2 + dy^2 == a2 ist, ist dies nicht noetig.
                    if(dx*dx + dy*dy > a2){
                            //x-Richtung nahe am Rand?
                            if((dx-L)*(dx-L) < dx*dx)
                                    x2 += L;
                            else if((dx+L)*(dx+L) < dx*dx)
                                    x2 -= L;
                            dx = x1-x2;
                            
                            //y-Richtung nahe am Rand?
                            if((dy-L)*(dy-L) < dy*dy)
                                    y2 += L;
                            else if((dy+L)*(dy+L) < dy*dy)
                                    y2 -= L;
                            dy = y1-y2;
                            
                    }//if periodische Randbedingungen
                    
                    a2 *= sigma11_22*sigma11_22;
                    dx *= sigma11_22;
                    dy *= sigma11_22;
                    
                    
                    // a^(-8) und a^(-14)
                    double a_8 = 1.0/(a2*a2*a2*a2);
                    double a_14= 1.0/(a2*a2*a2*a2*a2*a2*a2);
                    
                    //addiere Kraefte zu F_WCA[i]. Die 4 kommt aus dem Lennard-Jones-Potential, die 6 aus d/dr r^(-6)
                    F2_WCA[i][0] += eps22_11*sigma11_22*4.0*6.0*(2*a_14 - a_8)*dx;
                    F2_WCA[i][1] += eps22_11*sigma11_22*4.0*6.0*(2*a_14 - a_8)*dy;
                    
            }//for j
    }//for i
    
}//void berechne_WCA22_debug

//berechne WCA-Kräfte aus 12-Wechselwirkung, addiere sie in F1_WCA[N1][2] und F2_WCA[N2][2] 
void RunZustand::RunDynamik::addiere_WCA12(){
	//using namespace dynamik_methoden;
	
	int i,j2; // i läuft über Typ 1, j läuft über Typ 2
	vector<int>::iterator j; //läuft über die Typ2-Nachbarliste
	double a2; //abstandsquadrat
	int ic, jc; //Indices der Zelle, wo das Typ1-Teilchen drin ist
	
	for(i=0; i<N1; i++){
		ic = r1_git[i][0];
		jc = r1_git[i][1];
		
		for(j=erwNachbarn2[ic][jc].begin(); j!=erwNachbarn2[ic][jc].end(); j++){
			j2 = *j;
			a2 = abstand2_12(i,j2);
			
			//Falls Abstand zu groß, überspringen
			if(a2 * sigma11_12*sigma11_12> zweihoch1_6*zweihoch1_6) continue;
			
			//Absolutkoordinaten
			double x1 = r1_git[i ][0]*nachList_Breite + r1_rel[i ][0];
			double y1 = r1_git[i ][1]*nachList_Breite + r1_rel[i ][1];
			double x2 = r2_git[j2][0]*nachList_Breite + r2_rel[j2][0];
			double y2 = r2_git[j2][1]*nachList_Breite + r2_rel[j2][1];
			
			double dx = x1-x2;
			double dy = y1-y2;
			
			//Periodische Randbedingungen: Falls nahe am Rand, x2 y2 verschieben
			if(dx*dx+dy*dy > a2){
				//x-Richtung nahe am Rand?
				//dx = x1 - x2. dx-L ist also x1-(x2+L), und wenn dx-L kürzer ist als dx, setze x2 auf x2+L.
				if((dx-L)*(dx-L) < dx*dx)
					x2 += L;
				else if((dx+L)*(dx+L) < dx*dx)
					x2 -= L;
				
				//y-Richtung nahe am Rand?
				if((dy-L)*(dy-L) < dy*dy)
					y2 += L;
				else if((dy+L)*(dy+L) < dy*dy)
					y2 -= L;
				
			}//if nahe am Rand
			
			//skalieren mit sigma11_12
			a2 *= sigma11_12*sigma11_12;
			dx *= sigma11_12;
			dy *= sigma11_12;
			
			// a hoch -8 und hoch -14
			double a_8 = 1.0/(a2*a2*a2*a2);
			double a_14= 1.0/(a2*a2*a2*a2*a2*a2*a2);
			
			
			F1_WCA[i ][0] += eps12_11*sigma11_12* 4.0*6.0*(2*a_14 - a_8)*dx;
			F1_WCA[i ][1] += eps12_11*sigma11_12* 4.0*6.0*(2*a_14-a_8)*dy;
			F2_WCA[j2][0] -= eps12_11*sigma11_12* 4.0*6.0*(2*a_14-a_8)*dx;
			F2_WCA[j2][1] -= eps12_11*sigma11_12* 4.0*6.0*(2*a_14-a_8)*dy;
		}//for j durch die Nachbarliste
	}//for i, Typ 1
	
}//void addiere_WCA12

void RunZustand::RunDynamik::addiere_WCA12_debug(TimestepInfo &info){
    int i,j2; // i läuft über Typ 1, j läuft über Typ 2
    vector<int>::iterator j; //läuft über die Typ2-Nachbarliste
    double a2; //abstandsquadrat
    int ic, jc; //Indices der Zelle, wo das Typ1-Teilchen drin ist
    
    for(i=0; i<N1; i++){
            ic = r1_git[i][0];
            jc = r1_git[i][1];
            
            for(j=erwNachbarn2[ic][jc].begin(); j!=erwNachbarn2[ic][jc].end(); j++){
                    j2 = *j;
                    a2 = abstand2_12(i,j2);
                    info.num_WCA_berechnet++;
                    
                    //Falls Abstand zu groß, überspringen
                    if(a2 * sigma11_12*sigma11_12> zweihoch1_6*zweihoch1_6) continue;
                    
                    //Absolutkoordinaten
                    double x1 = r1_git[i ][0]*nachList_Breite + r1_rel[i ][0];
                    double y1 = r1_git[i ][1]*nachList_Breite + r1_rel[i ][1];
                    double x2 = r2_git[j2][0]*nachList_Breite + r2_rel[j2][0];
                    double y2 = r2_git[j2][1]*nachList_Breite + r2_rel[j2][1];
                    
                    double dx = x1-x2;
                    double dy = y1-y2;
                    
                    //Periodische Randbedingungen: Falls nahe am Rand, x2 y2 verschieben
                    if(dx*dx+dy*dy > a2){
                            //x-Richtung nahe am Rand?
                            //dx = x1 - x2. dx-L ist also x1-(x2+L), und wenn dx-L kürzer ist als dx, setze x2 auf x2+L.
                            if((dx-L)*(dx-L) < dx*dx)
                                    x2 += L;
                            else if((dx+L)*(dx+L) < dx*dx)
                                    x2 -= L;
                            
                            //y-Richtung nahe am Rand?
                            if((dy-L)*(dy-L) < dy*dy)
                                    y2 += L;
                            else if((dy+L)*(dy+L) < dy*dy)
                                    y2 -= L;
                            
                    }//if nahe am Rand
                    
                    //skalieren mit sigma11_12
                    a2 *= sigma11_12*sigma11_12;
                    dx *= sigma11_12;
                    dy *= sigma11_12;
                    
                    // a hoch -8 und hoch -14
                    double a_8 = 1.0/(a2*a2*a2*a2);
                    double a_14= 1.0/(a2*a2*a2*a2*a2*a2*a2);
                    
                    
                    F1_WCA[i ][0] += eps12_11*sigma11_12* 4.0*6.0*(2*a_14 - a_8)*dx;
                    F1_WCA[i ][1] += eps12_11*sigma11_12* 4.0*6.0*(2*a_14-a_8)*dy;
                    F2_WCA[j2][0] -= eps12_11*sigma11_12* 4.0*6.0*(2*a_14-a_8)*dx;
                    F2_WCA[j2][1] -= eps12_11*sigma11_12* 4.0*6.0*(2*a_14-a_8)*dy;
            }//for j durch die Nachbarliste
    }//for i, Typ 1
    
    //// Debug-spezifisch: Suche größte Kraft, für info.maxWCAx etc
    // Typ 1
    int maxTyp=1;
    int maxIndex=-1;
    double maxAbs2 = 0.0;
    for(i=0; i<N1; i++){
        double x = F1_WCA[i][0];
        double y = F1_WCA[i][1];
        if(x*x + y*y > maxAbs2){
            maxIndex=i; maxAbs2=x*x+y*y;
        }
    }//for i bis N1
    for(i=0; i<N2; i++){
        double x = F2_WCA[i][0];
        double y = F2_WCA[i][0];
        if(x*x + y*y > maxAbs2){
            maxTyp=2; maxIndex=i; maxAbs2=x*x+y*y;
        }
    }//for i bis N2
    info.maxWCAtype=maxTyp;
    info.maxWCAindex=maxIndex;
    if(maxTyp == 1){ info.maxWCAx = F1_WCA[maxIndex][0]; info.maxWCAy = F1_WCA[maxIndex][1]; }
    else           { info.maxWCAx = F2_WCA[maxIndex][0]; info.maxWCAy = F2_WCA[maxIndex][1]; }
    
    int imax, jmax;
    if(maxTyp == 1){ imax = r1_git[maxIndex][0]; jmax = r1_git[maxIndex][1]; }
    else           { imax = r2_git[maxIndex][0]; jmax = r2_git[maxIndex][1]; }
    info.maxWCA_numWW1 = erwNachbarn1[imax][jmax].size();
    info.maxWCA_numWW2 = erwNachbarn2[imax][jmax].size();
    
}//void addiere_WCA12_debug


//berechnet Kapillarkraefte (mittels Fouriertransformation), schreibt sie in Fkap. 
void RunZustand::RunDynamik::berechne_kapkraefte(){
	
	//using namespace dynamik_methoden;
	int j,l;
	double x,y;

/// Density-Gridding: Schreibe Dichte in rhox[][0]
	switch(densGrid_Schema){
		case 0: gridDensity_NGP(); break;
		case 1: gridDensity_CIC(); break;
		case 2: gridDensity_TSC(); break;
		default: cout << "Density-Gridding-Schema (NearestGridPoint/CloudInCell/TSC) nicht erkannt! densGrid_Schema="<<densGrid_Schema<<", zulaessig sind nur 0,1,2." << endl; 
	}//switch

/// FFT: Schreibe transformierte Dichte in rhok[][]
	fftw_execute(forward_plan);

/// Multiplikation mit i sin()/dx G: Schreibe Ergebnis in Fxk[][] und Fyk[][]. Komplexe Multiplikation, Gsinx rein imaginaer

	for(j=0; j<Z*Z; j++){
		Fxk[j][0] = - sinxG[j] * rhok[j][1]; // -im*im
		Fxk[j][1] =   sinxG[j] * rhok[j][0]; // im*re

		Fyk[j][0] = - sinyG[j] * rhok[j][1]; 
		Fyk[j][1] =   sinyG[j] * rhok[j][0];
	}//for j



/// FFT-1: Schreibe ruecktransformierte Kraefte in Fx[][0] und Fy[][0]
	fftw_execute(backx_plan);
	fftw_execute(backy_plan);


/// inverse Density-Gridding: Schreibe Kraefte in Fkap
	switch(densGrid_Schema){
		case 0: inv_gridDensity_NGP(); break; 
		case 1: inv_gridDensity_CIC(); break;
		case 2: inv_gridDensity_TSC(); break;
		default: cout << "Density-Gridding-Schema (NearestGridPoint/CloudInCell/TSC) nicht erkannt! densGrid_Schema="<<densGrid_Schema<<", zulaessig sind nur 0,1,2." << endl; 
	}//switch

}//void berechne_kapkraefte

//berechne Kapillarkraefte (via FFT), schreibe sie in Fkap. Informationen in die Felder von info
void RunZustand::RunDynamik::berechne_kapkraefte_debug(TimestepInfo &info){
    int j,l;
    double x,y;
    
/// Density-Gridding
    switch(densGrid_Schema){
        case 0: gridDensity_NGP(); break;
        case 1: gridDensity_CIC(); break;
        case 2: gridDensity_TSC(); break;
	default: cout << "Density-Gridding-Schema (NearestGridPoint/CloudInCell/TSC) nicht erkannt! densGrid_Schema="<<densGrid_Schema<<", zulaessig sind nur 0,1,2." << endl; 
    }//switch
    
/// FFT: Transformierte Dichte in rhok[][]
    fftw_execute(forward_plan);
    
/// Multiplikation mit i sin()/dx G: Schreibe Ergebnis in Fxk[][] und Fyk[][]. Komplexe Multiplikation. Bedenke, dass G*sin() rein imaginär ist.
    for(j=0; j<Z*Z; j++){
        Fxk[j][0] = -sinxG[j] * rhok[j][1]; // -im*im
        Fxk[j][1] =  sinxG[j] * rhok[j][0]; // im*re
        
        Fyk[j][0] = -sinyG[j] * rhok[j][1];
        Fyk[j][1] =  sinyG[j] * rhok[j][0];
    }//for j
    
/// FFT-1: Schreibe rücktransformierte Kräfte in Fx[][0] und Fy[][0]
    fftw_execute(backx_plan);
    fftw_execute(backy_plan);
    
/// inverse Density_Gridding: Kräfte in Fkap
    switch(densGrid_Schema){
		case 0: inv_gridDensity_NGP(); break; 
		case 1: inv_gridDensity_CIC(); break;
		case 2: inv_gridDensity_TSC(); break;
		default: cout << "Density-Gridding-Schema (NearestGridPoint/CloudInCell/TSC) nicht erkannt! densGrid_Schema="<<densGrid_Schema<<", zulaessig sind nur 0,1,2." << endl; 
    }//switch
    
/// suche betrags-größte Kraft, für info
    double maxBetrag2 = 0.0;
    int maxBetragIndex=-1;
    int maxBetragTyp=1;
    for(l=0; l<N1; l++){
        double betrag2 = F1kap[l][0]*F1kap[l][0] + F1kap[l][1]*F1kap[l][1];
        if (betrag2 > maxBetrag2){maxBetragIndex=l; maxBetrag2=betrag2;}
    }//for l bis N1
    
    for(l=0; l<N2; l++){
        double betrag2 = F2kap[l][0]*F2kap[l][0] + F2kap[l][1]*F2kap[l][1];
        if (betrag2 > maxBetrag2){maxBetragTyp=2; maxBetragIndex=l; maxBetrag2=betrag2;}
    }//for l bis N2
    
    info.maxKapindex=maxBetragIndex;

    if(maxBetragTyp == 1){
        info.maxKaptype = 1;
        info.maxKapx = F1kap[maxBetragIndex][0];
        info.maxKapy = F1kap[maxBetragIndex][1];
    }//if Typ 1
    else {
        info.maxKaptype = 2;
        info.maxKapx = F2kap[maxBetragIndex][0];
        info.maxKapy = F2kap[maxBetragIndex][1];
    }//else Typ 2
    
    
}//void berechne_kapkraefte_debug



//Schreibt Zufallszahlen mit Varianz 1.0 in F1_noise[N1][2] und F2_noise[N2][2]. Nicht gaussverteilt, sondern gleichverteilt.
void RunZustand::RunDynamik::berechne_zufallskraefte(){
	const double wurzelfuenf = sqrt(5.0);
	double x,y;
	for(int i=0; i<N1; i++){
		do{
			x = wurzelfuenf*zufall_gleichverteilt_vonbis(-1.0,1.0);
			y = wurzelfuenf*zufall_gleichverteilt_vonbis(-1.0,1.0);
		} while(x*x + y*y > 5.0);
		
		F1_noise[i][0] = x;
		F1_noise[i][1] = y;
		
	}//for i bis N1
	
	for(int i=0; i<N2; i++){
		do{
			x = wurzelfuenf*zufall_gleichverteilt_vonbis(-1.0,1.0);
			y = wurzelfuenf*zufall_gleichverteilt_vonbis(-1.0,1.0);
		} while(x*x + y*y > 5.0);
		F2_noise[i][0] = x;
		F2_noise[i][1] = y;
	}//for i bis N2
}//void berechne_zufallskraefte

void RunZustand::RunDynamik::berechne_zufallskraefte_debug(TimestepInfo &info){
    
    const double wurzelfuenf = sqrt(5.0);
    double x,y;
    int maxIndex = -1;
    double maxKraft2 = 0.0;
    int maxTyp = 1;
    
    for(int i=0; i<N1; i++){
        do{
            x = wurzelfuenf*zufall_gleichverteilt_vonbis(-1.0, 1.0);
            y = wurzelfuenf*zufall_gleichverteilt_vonbis(-1.0, 1.0);
        } while(x*x + y*y > 5.0);
        F1_noise[i][0] = x;
        F1_noise[i][1] = y;
        if(x*x + y*y > maxKraft2){maxIndex = i; maxKraft2 = x*x+y*y;}
        
    }//for i bis N1
    for(int i=0; i<N2; i++){
        do{
            x = wurzelfuenf*zufall_gleichverteilt_vonbis(-1.0, 1.0);
            y = wurzelfuenf*zufall_gleichverteilt_vonbis(-1.0, 1.0);
        } while(x*x + y*y > 5.0);
        F2_noise[i][0] = x;
        F2_noise[i][1] = y;
        if(x*x + y*y > maxKraft2){maxTyp=2; maxIndex=i; maxKraft2 = x*x + y*y;}
    }//for i bis N2
    
    info.maxRNDtype = maxTyp;
    info.maxRNDindex= maxIndex;
    if (maxTyp == 1){
        info.maxRNDx = F1_noise[maxIndex][0];
        info.maxRNDy = F1_noise[maxIndex][1];
    } // if
    else{
        info.maxRNDx = F2_noise[maxIndex][0];
        info.maxRNDy = F2_noise[maxIndex][1];
    }//else
    
}//void berechne_zufallskraefte_debug


//Greensfunktion an der Stelle qjk. Index-Wrapping, Verschiebung.
double G(int j, int k){
	//using namespace dynamik_methoden;
	
	double qx = dq*(((int)(j+Z*0.5))%Z - 0.5*Z );
	double qy = dq*(((int)(k+Z*0.5))%Z - 0.5*Z );
	if(j==0 && k==0 && lambda_kapillar>0.2*L) return 0.0;
	return 1.0/( 4*sin(0.5*qx*dx)*sin(0.5*qx*dx)/dx/dx + 4*sin(0.5*qy*dx)*sin(0.5*qy*dx)/dx/dx + 1.0/(lambda_kapillar*lambda_kapillar) );
	
	//return - qx*qx - qy*qy;
	//return - 4*sin(0.5*qx*dx)*sin(0.5*qx*dx)/dx/dx - 4*sin(0.5*qy*dx)*sin(0.5*qy*dx)/dx/dx;
}//double G


//streiche Eintrag aus erwListe
void erwListe_rem(vector<int>& liste, int k){
	vector<int>::iterator i;
	
	//springe zum richtigen Eintrag
	for(i=liste.begin(); i!=liste.end(); ++i)
		if(*i == k) break;
	
	/* Test: War k ueberhaupt in der Liste?
	if(*i != k){
		cout << "Fehler in erwListe_rem(): Teilchen " << k << " nicht gefunden! Liste: ";
		for(i=liste.begin(); i<liste.end(); i++){
			cout << *i << " ";
		cout << endl;
		return;
		}//for i
	}//if k nicht in Liste
	// */
	
	liste.erase(i);
}//void erwListe_rem

/*
//bestimme optimale Dauer des Zeitschritts, so dass max_reisedistanz eingehalten wird
double optimaler_zeitschritt(double** F_WCA, double** Fkap, double** F_noise, double deltat, int NN){
	if (0==NN) return deltat;
	
	double ret = deltat; //spaeter return-Wert, kann nur kleiner werden
	
	const double mrd = max_reisedistanz; //kuerzerer Name
	
	//Loesen quadratischer Gleichungen:
	// Aus m < |Fdt + sqrt(2T dt)dn| erhaelt man
	// Falls F==0: dt < m^2/(2T (dn^2))
	// Falls F!=0: dt < (w-|u|)^2 mit u=dn*sqrt(2T)/(8F) und w=sqrt(u^2 + m/(4|F|))
	
	// hier:
	// m = mrd
	// F = F_WCA+Fkap
	// dn = F_noise
	// T=T und dt = ret
	
	
	for(int teilchen=0; teilchen<NN; teilchen++){
		
		// x-Komponente
		double F = F_WCA[teilchen][0] + Fkap[teilchen][0];
		if(fabs(F) < 1.0e-5) //eigentlich F==0.0. Konstante 10^-5 auslagern oder durch was anderes ausdruecken, zB F << F_noise ?
			ret = min(ret, mrd*mrd/(2.0*T*F_noise[teilchen][0]*F_noise[teilchen][0]));
		else{
			double u = F_noise[teilchen][0]*sqrt(0.5*T)/F;
			double w = sqrt(u*u + mrd/fabs(F));
			ret = min(ret, (w-fabs(u))*(w-fabs(u)));
		}//else F ungleich Null
		
		
		
		
		// y-Komponente
		F = F_WCA[teilchen][1] + Fkap[teilchen][1];
		if(fabs(F) < 1.0e-5) //eigentlich F==0.0. Konstante 10^-5 auslagern oder durch was anderes ausdruecken, zB F << F_noise ?
			ret = min(ret, mrd*mrd/(2.0*T*F_noise[teilchen][1]*F_noise[teilchen][1]));
		else{
			double u = F_noise[teilchen][1]*sqrt(0.5*T)/F;
			double w = sqrt(u*u + mrd/fabs(F));
			ret = min(ret, (w-fabs(u))*(w-fabs(u)));
		}//else F ungleich Null
		
		
		
	}//for teilchen
	return ret;
}//double optimaler_zeitschritt 
// */

double RunZustand::RunDynamik::optimaler_zeitschritt(){
	double ret = dt_max; //späterer Return-Wert, kann nur kleiner werden

	const double mrd = max_reisedistanz;
	
	////Loesen quadratischer Gleichungen:
	// Aus m < |Fdt + sqrt(2T dt)dn| erhaelt man
	// Falls F==0: dt < m^2/(2T (dn^2))
	// Falls F!=0: dt < (w-|u|)^2 mit u=dn*sqrt(2T)/(8F) und w=sqrt(u^2 + m/(4|F|))
	
	// hier:
	// m = mrd
	// F = F_WCA+Fkap
	// dn = F_noise
	// T=T und dt = ret
	
	//Erst Typ 1
	for(int teilchen=0; teilchen<N1; teilchen++){
                if(noWCA){
                    if(F1_WCA[teilchen][0] != 0.0 || F1_WCA[teilchen][1] != 0.0)
                        cout << "Warnung, Teilchen "<<teilchen<<" (Typ 1) hat WCA-Kraft." <<endl;
                }//if noWCA, warnen wenn doch WCA-Kraft auftritt
		//x-Richtung
		double F = F1_WCA[teilchen][0] + F1kap[teilchen][0];
		if(fabs(F) < 1.0e-5)
			ret = min(ret, mrd*mrd/(2.0*T*F1_noise[teilchen][0]*F1_noise[teilchen][0]));
		else{
			double u = F1_noise[teilchen][0]*sqrt(0.5*T)/F;
			double w = sqrt(u*u + mrd/fabs(F));
			ret = min(ret, (w-fabs(u))*(w-fabs(u)));
		}//else
		
		//y-Richtung
		F = F1_WCA[teilchen][1] + F1kap[teilchen][1];
		if(fabs(F) < 1.0e-5)
			ret = min(ret, mrd*mrd/(2.0*T*F1_noise[teilchen][1]*F1_noise[teilchen][1]));
		else{
			double u = F1_noise[teilchen][1]*sqrt(0.5*T)/F;
			double w = sqrt(u*u + mrd/fabs(F));
			ret = min(ret, (w-fabs(u))*(w-fabs(u)));
		}//else
	}//for teilchen
	
	//Dann Typ 2
	for(int teilchen=0; teilchen<N2; teilchen++){
                if(noWCA){
                    if(F2_WCA[teilchen][0] != 0.0 || F2_WCA[teilchen][1] != 0.0)
                        cout << "Warnung, Teilchen "<<teilchen<<" (Typ 2) hat WCA-Kraft." <<endl;
                }//if noWCA, warnen wenn doch WCA-Kraft auftritt
		//x-Richtung
		double F = F2_WCA[teilchen][0] + F2kap[teilchen][0];
		if(fabs(F) < 1.0e-5)
			ret = min(ret, mrd*mrd/(2.0*T*F2_noise[teilchen][0]*F2_noise[teilchen][0]));
		else{
			double u = F2_noise[teilchen][0]*sqrt(0.5*T)/F;
			double w = sqrt(u*u + mrd/fabs(F));
			ret = min(ret, (w-fabs(u))*(w-fabs(u)));
		}//else
		
		//y-Richtung
		F = F2_WCA[teilchen][1] + F2kap[teilchen][1];
		if(fabs(F) < 1.0e-5)
			ret = min(ret, mrd*mrd/(2.0*T*F2_noise[teilchen][1]*F2_noise[teilchen][1]));
		else{
			double u = F2_noise[teilchen][1]*sqrt(0.5*T)/F;
			double w = sqrt(u*u + mrd/fabs(F));
			ret = min(ret, (w-fabs(u))*(w-fabs(u)));
		}//else
	}//for teilchen
	
	return ret;
}//double RunZustand::RunDynamik::optimaler_zeitschritt














//liest Startpositionen aus Datei. Es muss zuerst Typ 1 und dann Typ 2 stehen.
void RunZustand::init_pos_aus_datei(){
	const string startpos_dateiname="startpos.txt";
	
	FILE* in = fopen(startpos_dateiname.c_str(), "r");
	
	//erste zwei Zeilen überspringen
	char buffer[100];
	//zwei Zeilen lesen. Wenn beim Lesen von mindestens einer ein Fehler kommt, abbrechen
	if(NULL == fgets(buffer, 100, in) || NULL == fgets(buffer, 100, in)){
		cout << "Fehler beim Lesen der Inputdatei, Zeile 1 oder 2" << endl << flush;
		return;
	}//if Fehler beim Lesen
	
	// lese zeilenweise
	int teilchenNr=0;
	int typ;
	double x,y;
	float xf, yf;
	bool fehler=false;
	for(int zeile=0; zeile<N1; zeile++){
		if(4 != fscanf(in, "%d \t %d \t %g \t %g", &teilchenNr, &typ, &xf, &yf)) fehler=true;
		if(teilchenNr != zeile || typ != 1) fehler=true;
		
		if(fehler){
			cout << "Fehler beim Lesen der Inputdatei, Zeile " << zeile << endl << flush;
			return;
		}//if fehler
		
		x = (double) xf;
		y = (double) yf;
		if(x<0.0 || x>L || y<0.0 || y>L){
			cout << "Fehler beim Laden: Koordinaten von Teilchen "<<teilchenNr<<" (Typ 1), x="<<x<<", y="<<y<<" out of range!" << endl << flush;
			return;
		}//if x oder y out of range
		r1_git[teilchenNr][0] = (int) (x/nachList_Breite);
		r1_git[teilchenNr][1] = (int) (y/nachList_Breite);
		r1_rel[teilchenNr][0] = x - r1_rel[teilchenNr][0]*nachList_Breite;
		r1_rel[teilchenNr][1] = y - r1_rel[teilchenNr][1]*nachList_Breite;
		
	}//for int zeile
	for(int zeile=0; zeile<N2; zeile++){
		if(4 != fscanf(in, "%d \t %d \t %g \t %g", &teilchenNr, &typ, &xf, &yf)) fehler=true;
		if(teilchenNr != zeile || typ != 2) fehler=true; 
		
		if(fehler){
			cout << "Fehler beim Lesen der Inputdatei, Zeile " << zeile+N1<< ". Möglicherweise sind nicht genau " <<N1<<"+"<<N2<<"Teilchen in startpos.txt?" << endl<<flush;
			return;
		}//if fehler
		
		x = (double) xf;
		y = (double) yf;
		if(x<0.0 || x>L || y<0.0 || y>L){
			cout << "Fehler beim Laden: Koordinaten von Teilchen "<<teilchenNr<<" (Typ 2), x="<<x<<", y="<<y<<" out of range!" << endl << flush;
			return;
		}//if x oder y out of range
		
		r2_git[teilchenNr][0] = (int) (x/nachList_Breite);
		r2_git[teilchenNr][1] = (int) (y/nachList_Breite);
		r2_rel[teilchenNr][0] = x - r2_rel[teilchenNr][0]*nachList_Breite;
		r2_rel[teilchenNr][1] = y - r2_rel[teilchenNr][1]*nachList_Breite;

	}//for zeile
	
	fclose(in);
	
}//void init_pos_aus_datei


// initialisiert Teilchenpositionen auf Zufallspositionen. Schreibt sie in r_git und r_rel
void RunZustand::init_zufallspos(){

	
//bestimme zufaellige Positionen in [0,L], setze die Teilchen dort hin
	for(int i=0; i<N1; i++){
		double x = zufall_gleichverteilt_vonbis(0.0, L);
		double y = zufall_gleichverteilt_vonbis(0.0, L);

		int jc = (int) x/nachList_Breite;
		int kc = (int) y/nachList_Breite;

		r1_git[i][0] = jc;
		r1_git[i][1] = kc;
		r1_rel[i][0] = x - jc*nachList_Breite;
		r1_rel[i][1] = y - kc*nachList_Breite;
		
	}//for i
	
	for(int i=0; i<N2; i++){
		double x = zufall_gleichverteilt_vonbis(0.0, L);
		double y = zufall_gleichverteilt_vonbis(0.0, L);
		
		int jc = (int) x/nachList_Breite;
		int kc = (int) y/nachList_Breite;
		
		r2_git[i][0] = jc;
		r2_git[i][1] = kc;
		r2_rel[i][0] = x - jc*nachList_Breite;
		r2_rel[i][1] = y - kc*nachList_Breite;
	}//for i
	
} //void init_zufallspos


// bestimme EINE Position in [0,L], setze alle Teilchen dorthin
void RunZustand::init_allegleich(){
	double x = zufall_gleichverteilt_vonbis(0.0,L);
	double y = zufall_gleichverteilt_vonbis(0.0,L);
	
	// x = 2.0;
	// y = 0.0;
	
	int jc = (int) (x/nachList_Breite);
	int kc = (int) (y/nachList_Breite);
	for(int i=0; i<N1; i++){
		r1_git[i][0] = jc;
		r1_git[i][1] = kc;
		r1_rel[i][0] = x - jc*nachList_Breite;
		r1_rel[i][1] = y - kc*nachList_Breite;
	}//for i
	
	
	//Typ 2 analog
	x = zufall_gleichverteilt_vonbis(0.0, L);
	y = zufall_gleichverteilt_vonbis(0.0, L);
	jc = (int) (x/nachList_Breite);
	kc = (int) (y/nachList_Breite);
	
	for(int i=0; i<N2; i++){
		r2_git[i][0] = jc;
		r2_git[i][1] = kc;
		r2_rel[i][0] = x - jc*nachList_Breite;
		r2_rel[i][1] = y - kc*nachList_Breite;
	}//for i
}//void init_allegleich

// initialisiert alle Teilchenpositionen zufällig gleichverteilt in einem Kreis in der Mitte der Box
void RunZustand::init_kreisscheibe(){
	const double rad = startpos_kreisradius;
	const double rad2 = rad*rad;
// 	double x=2.0;
// 	double y=2.0;
	for(int i=0; i<N1; i++){
		double x,y;
		do{
			x = zufall_gleichverteilt_vonbis(-rad,rad);
			y = zufall_gleichverteilt_vonbis(-rad,rad);
		} while (x*x + y*y > rad2);
		x+= 0.5*L;
		y+= 0.5*L;

		
		int jc = (int) (x/nachList_Breite);
		int kc = (int) (y/nachList_Breite);
		
		r1_git[i][0] = jc;
		r1_git[i][1] = kc;
		r1_rel[i][0] = x - jc*nachList_Breite + 0.03;
		r1_rel[i][1] = y - kc*nachList_Breite + 0.03;
		
	}//for i
	
	for(int i=0; i<N2; i++){
		double x,y;
		do{
			x = zufall_gleichverteilt_vonbis(-rad, rad);
			y = zufall_gleichverteilt_vonbis(-rad, rad);
		} while(x*x + y*y > rad2);
		x+= 0.5*L;
		y+= 0.5*L;
		
		int jc = (int)(x/nachList_Breite);
		int kc = (int)(y/nachList_Breite);
		
		r2_git[i][0] = jc;
		r2_git[i][1] = kc;
		r2_rel[i][0] = x - jc*nachList_Breite + 0.03;
		r2_rel[i][1] = y - kc*nachList_Breite + 0.03;
		
	}//for i
	
}//void init_kreisscheibe


// initialisiert Teilchenpositionen gleichverteilt in der Box. Baut ein quadratisches Gitter auf. Wenn N keine Quadratzahl ist, bleiben dabei Gitterplätze frei.
void RunZustand::init_gitterstart(){
	
	if(0 != N2){
		cout << "Fehler: Gitterstart mit zwei Teilchentypen noch nicht implementiert!" << endl << flush;
		return;
	}//if N2 ungleich Null
	// es sollen N Teilchen auf ein 2dim-Gitter verteilt werden. In jeder Richtung also sqrt(N) viele Teilchen auf die Länge L. Falls N keine Quadratzahl ist, sogar ceil(sqrt(N)).
	double wurzelN = ceil(sqrt(N1));
	double schrittweite = L/wurzelN;
	int index =0;
	for(int i=0; i<wurzelN; i++)
	for(int j=0; j<wurzelN; j++){
		double x = i*schrittweite;
		double y = j*schrittweite;
		
		int ic = (int) (x/nachList_Breite);
		int jc = (int) (y/nachList_Breite);
		
		r1_git[index][0] = ic;
		r1_git[index][1] = jc;
		r1_rel[index][0] = x - ic*nachList_Breite;
		r1_rel[index][1] = y - jc*nachList_Breite;
		
		if(++index == N1) return;
		
	}//for i,j
	
}//void init_gitterstart


void RunZustand::init_kernundring(){
    const double mitte = 0.5*L;
    if(N1 != 1 || N2 < 3){
        cout << "Fehler: Start ist Kern und Ring, aber die Teilchenanzahlen sind nicht (1,viele)!" << endl << flush;
        return;
    }//if
    
    double x, y;
    int ic, jc;
    
    {// Typ 1, gross und in der Mitte
        x = mitte;
        y = mitte;
        
        ic = (int)(x/nachList_Breite);
        jc = (int)(y/nachList_Breite);
        r1_git[0][0] = ic;
        r1_git[0][1] = jc;
        r1_rel[0][0] = x - ic*nachList_Breite;
        r1_rel[0][1] = y - jc*nachList_Breite;
    }
    
    for(int teil=0; teil<N2; teil++){// Typ 2, Ring mit Radius startpos_kreisradius drumherum
        x = mitte + startpos_kreisradius*cos(2*M_PI*(double)teil/N2);
        y = mitte + startpos_kreisradius*sin(2*M_PI*(double)teil/N2);
        
        ic = (int)(x/nachList_Breite);
        jc = (int)(y/nachList_Breite);
        r2_git[teil][0] = ic;
        r2_git[teil][1] = jc;
        r2_rel[teil][0] = x - ic*nachList_Breite;
        r2_rel[teil][1] = y - jc*nachList_Breite;

    }//for teil
    
    
}//void init_kreisundring