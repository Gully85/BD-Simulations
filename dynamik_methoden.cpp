// stellt Methoden bereit, die die Bewegung der Kolloide behandeln.
// ausserdem Methoden zur Initialisierung

#pragma once

#include "parameter.h"
#include "signaturen.h"
#include "zufall.cpp"
#include "gridRoutinen.cpp"
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
extern int** r1_git;
extern double** r1_rel;
extern int** r2_git;
extern double** r2_rel;

extern const double sigma11_22, sigma11_12, Gamma2_1, eps22_11, eps12_11;

extern const double f2_f1;

extern const int startpos_methode; //Zufall=1, aus Datei=2, Gitterstart=3, allegleich=4, Kreisscheibe=5


const double dx = densGrid_Breite;

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

//TODO: Das Aktualisieren der Nachbarlisten auslagern!
//Fuehre einen Zeitschritt durch: Berechne Kraefte, ermittle optimale Dauer, bewege Teilchen, aktualisiere ggf Nachbarlisten. Gibt Dauer zurueck.
double zeitschritt(double tmax){ 
	
	using namespace dynamik_methoden;
	
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
	berechne_zufallskraefte(N1, F1_noise);
	berechne_zufallskraefte(N2, F2_noise);
	
	/* Test: Gebe Teilchenpositionen und aktuelle Kraefte aus. 
	for(int teilchen=0; teilchen<N; teilchen++){
		cout << teilchen << " an (" << r_git[teilchen][0]*nachList_Breite+r_rel[teilchen][0] << "," << r_git[teilchen][1]*nachList_Breite+r_rel[teilchen][1] << "):\tWCA=(" << F_WCA[teilchen][0] << "," << F_WCA[teilchen][1]<< "),\tkap=(" << Fkap[teilchen][0]<<","<<Fkap[teilchen][1]<<"), noise=("<<F_noise[teilchen][0]<<","<<F_noise[teilchen][1]<<")" << endl;
	}//for teilchen
	// */
	
	//bestimme optimalen Zeitschritt, so dass kein Teilchen weiter als max_reisedistanz bewegt wird
	deltat = optimaler_zeitschritt(F1_WCA, F1kap, F1_noise, deltat, N1);
	
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
	}//for teilchen bis N2
	
	//erwListen aktualisieren: Falls ein Teilchen seine Zelle verlassen hat (dh r_rel<0 oder r_rel>nachList_Breite), streichen/hinzufügen
	refresh_erwNachbar();
	
	return deltat;
}//double zeitschritt


//erneuere erwNachbarlisten. Beim Aufruf sind r_rel<0 oder r_rel>nachList_Breite zulässig, wird behoben
void refresh_erwNachbar(){
	using namespace dynamik_methoden;
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




// initialisiert Felder/Nachbarlisten fuer WCA-Kraefte. Erwartet, dass in r_git schon die Gittervektoren der Teilchenpositionen stehen.
void WCA_init(){
	using namespace dynamik_methoden;
	//vector<int>** erwNachbarn ist die erweiterte Nachbarliste. Erster Index = ZellenNr in x-Richtung, zweiter=y-Richtung. Im vector stehen die Indices aller Teilchen, die in der gleichen oder benachbarten Zellen sind.

	
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
void kapkraefte_init(){
	using namespace dynamik_methoden;
/// alloziere Felder
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
	if(forward_plan == NULL){
		forward_plan = fftw_plan_dft_2d(Z, Z, rhox, rhok, FFTW_FORWARD,  FFTW_PATIENT);
		backx_plan   = fftw_plan_dft_2d(Z, Z, Fxk,  Fx,   FFTW_BACKWARD, FFTW_PATIENT);
		backy_plan   = fftw_plan_dft_2d(Z, Z, Fyk,  Fy,   FFTW_BACKWARD, FFTW_PATIENT);
	}//if forward_plan==NULL
}//void kapkraefte_init


// berechnet WCA-Kraefte. Schreibt sie in F1_WCA[N1][2]
void berechne_WCA11(){
	using namespace dynamik_methoden;
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

//berechne WCA-Kräfte aus 22-Wechselwirkung, schreibe sie in F2_WCA[N2][2]
void berechne_WCA22(){
	
	using namespace dynamik_methoden;
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

//berechne WCA-Kräfte aus 12-Wechselwirkung, addiere sie in F1_WCA[N1][2] und F2_WCA[N2][2] 
void addiere_WCA12(){
	using namespace dynamik_methoden;
	
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


//berechnet Kapillarkraefte (mittels Fouriertransformation), schreibt sie in Fkap. 
void berechne_kapkraefte(){
	
	using namespace dynamik_methoden;
	int j,l;
	double x,y;

/// Density-Gridding: Schreibe Dichte in rhox[][0]
	switch(densGrid_Schema){
		case 0: gridDensity_NGP(rhox, 1.0, f2_f1); break; //TODO gridDensity auf global umstellen, namespace nutzen
		case 1: gridDensity_CIC(rhox, 1.0, f2_f1); break;
		case 2: gridDensity_TSC(rhox, 1.0, f2_f1); break;
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
	}//for i



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

//Schreibt 2*anzahl Zufallszahlen mit Varianz 1.0 in dn[anzahl][2]. Nicht gaussverteilt, sondern gleichverteilt.
void berechne_zufallskraefte(int anzahl, double** dn){
	const double wurzelfuenf = sqrt(5.0);
	double x,y;
	for(int i=0; i<anzahl; i++){
		do{
			x = wurzelfuenf*zufall_gleichverteilt_vonbis(-1.0,1.0);
			y = wurzelfuenf*zufall_gleichverteilt_vonbis(-1.0,1.0);
		} while(x*x + y*y > 5.0);
		
		dn[i][0] = x;
		dn[i][1] = y;
		
	}//for i
}//void berechne_zufallskraefte


//reserviere Speicher, setze Teilchen auf Startpositionen, rufe andere init() auf.
void main_init(){
	
	using namespace dynamik_methoden;
	
	ausgabe_jeansgroessen();

	init_korrelationsfunktion();
	init_ftrho();  //Observable rho(k), berechnet via Summe über Teilchen und über Gitter im k-Raum
	init_rhoFFTW();//Observable rho(k), berechnet via density-Gridding und FFTW 
	
	// Teilchenpositionen, Gittervektor. Erster Index TeilchenNr, zweiter Index 0 fuer x und 1 fuer y-Richtung
	if(r1_git == NULL){
		r1_git = new int*[N1];
		// dasgleiche, Vektor innerhalb der Zelle
		r1_rel = new double*[N1];
		r2_git = new int*[N2];
		r2_rel = new double*[N2];
		for(int i=0; i<N1; i++){
			r1_git[i] = new int[2];
			r1_rel[i] = new double[2];
		}//for i
		for(int i=0; i<N2; i++){
			r2_git[i] = new int[2];
			r2_rel[i] = new double[2];
		}//for i
	}//if r_git == NULL
	
	//Kapillarkraefte. Erster Index TeilchenNr, zweiter Index Raumrichtung
	if(F1kap == NULL){
		F1kap = new double*[N1];
		F2kap = new double*[N2];
		//WCA-Kraefte, gleiche Indices
		F1_WCA = new double*[N1];
		F2_WCA = new double*[N2];
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
	}//if F1kap==NULL
	
	switch(startpos_methode){
		case 1:
			init_zufallspos();
			break;
		case 2:
			init_pos_aus_datei();
			break;
		case 3:
			init_gitterstart();
			break;
		case 4:
			init_allegleich();
			break;
		case 5:
			init_kreisscheibe();
			break;
		default:
			cout << "Fehler: startpos_methode="<<startpos_methode<<", erlaubt sind 1,2,3,4,5" << endl;
			return;
	}//switch startpos_methode
	
	
	//init Kraftberechnungen
	WCA_init();
	kapkraefte_init();
	
}//void main_init


//Greensfunktion an der Stelle qjk. Index-Wrapping, Verschiebung.
double G(int j, int k){
	using namespace dynamik_methoden;
	
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













