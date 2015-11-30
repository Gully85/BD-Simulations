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

using std::cout; using std::endl;
using std::vector;


extern const double densGrid_Breite, dq, lambda_kapillar, L, zweihoch1_6, kapillar_vorfaktor, dt_max, max_reisedistanz, T;
extern const int densGrid_Zellen, densGrid_Schema;
extern int** r_git;
extern double** r_rel;

const int Z = densGrid_Zellen;
const double dx = densGrid_Breite;

/// semi-globale Felder/Variablen, erreichbar fuer die Methoden in dieser Datei
fftw_complex* rhox = NULL; // Dichte vor  Fouriertrafo FORWARD. Index-Wrapping
fftw_complex* rhok = NULL; // Dichte nach Fouriertrafo FORWARD. Index-Wrapping, Verschiebung und alternierende Vorzeichen
double* sinxG = NULL;	   // der Term sin(qx dx)/dx G(q). Index-Wrapping und Verschiebung
double* sinyG = NULL;	   // der Term sin(qy dx)/dx G(q). Index-Wrapping und Verschiebung
fftw_complex* Fxk = NULL;  // rhok mal i sinxG. Index-Wrapping, Verschiebung und alternierende Vorzeichen
fftw_complex* Fyk = NULL;  // rhok mal i sinyG. Index-Wrapping, Verschiebung und alternierende Vorzeichen
fftw_complex* Fx = NULL;   // Kraefte (x-komp) nach Fouriertrafo BACKWARD. Index-Wrapping
fftw_complex* Fy = NULL;   // Kraefte (y-komp) nach Fouriertrafo BACKWARD. Index-Wrapping
vector<int>** erwNachbarn = NULL; // Erweiterte Nachbarliste. Erster Index = ZellenNr in x-Richtung, zweiter=y-Richtung. Im vector stehen die Indices aller Teilchen, die in der gleichen oder benachbarten Zellen sind.
double** Fkap = NULL; // Kapillarkraefte. Erster Index TeilchenNr, zweiter Index Raumrichtung
double** F_WCA = NULL; //WCA-Kraefte. Erster Index TeilchenNr, zweiter Index Raumrichtung
double** F_noise = NULL; //Zufallskraefte. Erster Index TeilchenNr, zweiter Index Raumrichtung



fftw_plan forward_plan, backx_plan, backy_plan;



//Fuehre einen Zeitschritt durch: Berechne Kraefte, ermittle optimale Dauer, bewege Teilchen, aktualisiere ggf Nachbarlisten. Gibt Dauer zurueck.
double zeitschritt(double tmax){ 
	double deltat=tmax; //tatsaechliche Dauer. Koennte kleiner werden als tmax
	
	if(tmax <= 0.0){
		cout << "Fehler: Negative Zeitschrittweite " << tmax << " gefordert, breche ab!" << endl;
		return 1.0e10;
	}//if tmax negativ, alles abbrechen
	
	//Berechne alle benoetigten Kraefte
	berechne_WCAkraefte(F_WCA); //schreibt in F_WCA
	berechne_kapkraefte(r_git, r_rel, Fkap);
	berechne_zufallskraefte(N, F_noise);
	
	//* Test: Gebe Teilchenpositionen und aktuelle Kraefte aus. 
	for(int teilchen=0; teilchen<N; teilchen++){
		cout << teilchen << " an (" << r_git[teilchen][0]*nachList_Breite+r_rel[teilchen][0] << "," << r_git[teilchen][1]*nachList_Breite+r_rel[teilchen][1] << "):\tWCA=(" << F_WCA[teilchen][0] << "," << F_WCA[teilchen][1]<< "),\tkap=(" << Fkap[teilchen][0]<<","<<Fkap[teilchen][1]<<"), noise=("<<F_noise[teilchen][0]<<","<<F_noise[teilchen][1]<<")" << endl;
	}//for teilchen
	// */
	
	//bestimme optimalen Zeitschritt, so dass kein Teilchen weiter als max_reisedistanz bewegt wird
	deltat = optimaler_zeitschritt(F_WCA, Fkap, F_noise, deltat, N);
	
	cout << endl << "Zeitschritt: " << deltat << endl;
	
	//alle Teilchen bewegen
	for(int teilchen=0; teilchen<N; teilchen++){
		double dx = deltat * (F_WCA[teilchen][0] + Fkap[teilchen][0]) + sqrt(2.0*T*deltat)*F_noise[teilchen][0];
		double dy = deltat * (F_WCA[teilchen][1] + Fkap[teilchen][1]) + sqrt(2.0*T*deltat)*F_noise[teilchen][1];
		
		//*Test: dr ausgeben
		cout << "Teilchen " << teilchen << ": dr=("<<dx<<","<<dy<<")"<< endl;
		// */
		
		r_rel[teilchen][0] += dx;
		r_rel[teilchen][1] += dy;
		
		
		//ueber Zellgrenze bewegt?
		if(r_rel[teilchen][0] < 0.0 || r_rel[teilchen][0] > nachList_Breite
		|| r_rel[teilchen][1] < 0.0 || r_rel[teilchen][1] > nachList_Breite){
			//Indices der alten Zelle
			int ic = r_git[teilchen][0];
			int jc = r_git[teilchen][1];
			
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
				erwListe_rem(erwNachbarn[x[i]][y[j]], teilchen);
			
			//Indices der neuen Zelle. Koennen noch negative oder zu grosse Werte haben. floor() rundet auch negative Zahlen korrekt ab.
			ic = (int) floor((r_git[teilchen][0]*nachList_Breite + r_rel[teilchen][0])/nachList_Breite);
			jc = (int) floor((r_git[teilchen][1]*nachList_Breite + r_rel[teilchen][1])/nachList_Breite);
			
			// korrigiere negative bzw zu grosse Zellindices, durch period. Randbed.
			ic = (ic+nachList_Zellen)%nachList_Zellen;
			jc = (jc+nachList_Zellen)%nachList_Zellen;
			
			//trage neue Position in Gitter- und Relativvektor ein
			r_rel[teilchen][0] = fmod(r_rel[teilchen][0] + nachList_Breite, nachList_Breite); //fmod=modulo
			r_rel[teilchen][1] = fmod(r_rel[teilchen][1] + nachList_Breite, nachList_Breite);
			r_git[teilchen][0] = ic;
			r_git[teilchen][1] = jc;
			
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
				erwNachbarn[x[i]][y[j]].push_back(teilchen);
			
			
		}//if aus Zelle bewegt
		
	}//for i
	
}//double zeitschritt



// initialisiert Teilchenpositionen auf Zufallspositionen. Schreibt sie in r_git und r_rel
void init_zufallspos(){

//reserviere r_git und r_rel
	r_git = new int*[N];
	r_rel = new double*[N];
	for(int i=0; i<N; i++){
		r_git[i] = new int[2];
		r_rel[i] = new double[2];
	}//for i

//bestimme zufaellige Positionen in [0,L], setze die Teilchen dort hin
	for(int i=0; i<N; i++){
		double x = zufall_gleichverteilt_vonbis(0.0, L);
		double y = zufall_gleichverteilt_vonbis(0.0, L);

		int jc = (int) x/nachList_Breite;
		int kc = (int) y/nachList_Breite;

		r_git[i][0] = jc;
		r_git[i][1] = kc;
		r_rel[i][0] = x - jc*nachList_Breite;
		r_rel[i][1] = y - kc*nachList_Breite;
		
	}//for i
	
} //void init_zufallspos

// initialisiert Felder/Nachbarlisten fuer WCA-Kraefte. Erwartet, dass in r_git schon die Gittervektoren der Teilchenpositionen stehen.
void WCA_init(){

	//vector<int>** erwNachbarn ist die erweiterte Nachbarliste. Erster Index = ZellenNr in x-Richtung, zweiter=y-Richtung. Im vector stehen die Indices aller Teilchen, die in der gleichen oder benachbarten Zellen sind.

	//reserviere Speicher
	erwNachbarn = new vector<int>*[nachList_Zellen];
	for(int i=0; i<nachList_Zellen; i++){
		erwNachbarn[i] = new vector<int>[nachList_Zellen];
	}//for i

	int* x = new int[3]; //links davon, exakt, rechts davon
	int* y = new int[3];

	//trage jedes Teilchen in seine Zelle, und die 8 Nachbarzellen ein
	for(int teilchen=0; teilchen<N; teilchen++){
		int k = r_git[teilchen][0];
		int l = r_git[teilchen][1];
		
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
		    erwNachbarn[x[i]][y[j]].push_back(teilchen);
		
	}//for teilchen

} //void WCA_init

//initialisiert Felder, plant Fouriertrafos
void kapkraefte_init(){

/// alloziere Felder
	//Dichte vor FFT. Index-Wrapping. rhox[iw(j,k)][0] ist die Dichte in Zelle (j,k). rhox[][1] ist der Imaginaerteil, Null.
	rhox = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Z*Z);

	//Imaginaerteil auf Null setzen. Der wird sich nie aendern.
	for(int j=0; j<Z; j++)
		for(int k=0; k<Z; k++)
			rhox[iw(j,k)][1]=0.0;

	// Dichte nach FFT. Index-Wrapping, Verschiebung, Vorzeichen. rhok[iw(j,k)][0] ist der Realteil von (-1)^(j+k) mal rhoschlange(qjk), wobei qjk = dq* ((j+Zellen/2)%Zellen - Zellen/2 ex) und analog k ey
	rhok = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Z*Z);

	// der Term sin(qx dx)/dx G(q). Index-Wrapping, Verschiebung. Analog sinyG
	sinxG = new double[Z*Z];
	sinyG = new double[Z*Z];

	// Kraefte vor FFT-1. Index-Wrapping, Verschiebung, Vorzeichen.
	Fxk = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Z*Z);
	Fyk = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Z*Z);

	// Kraefte nach FFT-1. Index-Wrapping.
	Fx  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Z*Z);
	Fy  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Z*Z);



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

/// plane FFTs
	forward_plan = fftw_plan_dft_2d(Z, Z, rhox, rhok, FFTW_FORWARD,  FFTW_MEASURE);
	backx_plan   = fftw_plan_dft_2d(Z, Z, Fxk,  Fx,   FFTW_BACKWARD, FFTW_MEASURE);
	backy_plan   = fftw_plan_dft_2d(Z, Z, Fyk,  Fy,   FFTW_BACKWARD, FFTW_MEASURE);
}//void kapkraefte_init


// berechnet WCA-Kraefte. Schreibt sie in F_WCA[][2]
void berechne_WCAkraefte(double** F_WCA){
	int i; //Teilchen
	vector<int>::iterator j; //wechselwirkendes Teilchen. Iteriert ueber die Nachbarliste von i's Zelle
	double a2; //Abstandsquadrat der beiden wechselwirkenden Teilchen

	if(F_WCA == NULL){
		cout << "Fehler: WCA-Kraefteberechnung aufgerufen, aber nicht initialisiert!" << endl;
		return;
	}//if nicht initialisiert

	//setze F_WCA auf Null
	for(i=0; i<N; i++){
		F_WCA[i][0]=0.0;
		F_WCA[i][1]=0.0;
	}//for
	
	//Schleife ueber Teilchen. Fuer jedes, iteriere durch seine Nachbarn und addiere Kraefte.
	for(i=0; i<N; i++){
		int ic = r_git[i][0];
		int jc = r_git[i][1];
		
		for(j = erwNachbarn[ic][jc].begin(); j!=erwNachbarn[ic][jc].end(); j++){
			if(i == *j) continue; //ueberspringe WW mit sich selbst
			
			//Abstand berechnen, mit periodischen Randbedingungen
			a2 = abstand2(i, *j, r_git, r_rel);
			
			//falls zu weit entfernt, WW ueberspringen
			if(a2 > zweihoch1_6*zweihoch1_6) continue;
			
			//Absolutkoordinaten
			double x1 = r_git[ i][0]*nachList_Breite + r_rel[ i][0];
			double y1 = r_git[ i][1]*nachList_Breite + r_rel[ i][1];
			double x2 = r_git[*j][0]*nachList_Breite + r_rel[*j][0];
			double y2 = r_git[*j][1]*nachList_Breite + r_rel[*j][1];
			
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
			F_WCA[i][0] += 4.0*6.0*(2.0*a_14 - a_8) * dx;
			F_WCA[i][1] += 4.0*6.0*(2.0*a_14 - a_8) * dy;
			
			
			
		}//for j
	}//for i
}//void berechne_WCAkraefte



//berechnet Kapillarkraefte (mittels Fouriertransformation), schreibt sie in Fkap
void berechne_kapkraefte(int** rr_git, double** rr_rel, double** Fkap){
	int j,l;
	double x,y;

/// Density-Gridding: Schreibe Dichte in rhox[][0]
	switch(densGrid_Schema){
		case 0: gridDensity_NGP(rhox, rr_git, rr_rel); break;
		case 1: gridDensity_CIC(rhox, rr_git, rr_rel); break;
		case 2: gridDensity_TSC(rhox, rr_git, rr_rel); break;
		default: cout << "Density-Gridding-Schema (NearestGridPoint/CloudInCell/TSC) nicht erkannt! densGrid_Schema="<<densGrid_Schema<<", zulaessig sind nur 0,1,2." << endl; 
	}//switch

/* Test: Schreibe Dichte in Datei rho.txt
	out = fopen("rho.txt", "w");
	for(j=0; j<Z; j++){
		for(l=0; l<Z; l++)
			fprintf(out, "%g \t %g \t %g \n", j*dx, l*dx, rhox[iw(j,l)][0]);
		fprintf(out, "\n");
	}//for j
	fclose(out);	
// */

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
		case 0: inv_gridDensity_NGP(Fkap, Fx, Fy, rr_git, rr_rel); break;
		case 1: inv_gridDensity_CIC(Fkap, Fx, Fy, rr_git, rr_rel); break;
		case 2: inv_gridDensity_TSC(Fkap, Fx, Fy, rr_git, rr_rel); break;
		default: cout << "Density-Gridding-Schema (NearestGridPoint/CloudInCell/TSC) nicht erkannt! densGrid_Schema="<<densGrid_Schema<<", zulaessig sind nur 0,1,2." << endl; 
	}//switch
	
	/* Test: Schreibe Kapillarkraefte ins Terminal
	for(i=0; i<N; i++){
		cout << "Kapillarkraft auf Teilchen " << i << ": (" << Fkap[i][0] << ","<< Fkap[i][1] << ")" << endl;
	}//for i
	// */

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
	
	//Random Seed
	init_random();
	
	// Teilchenpositionen, Gittervektor. Erster Index TeilchenNr, zweiter Index 0 fuer x und 1 fuer y-Richtung
	r_git = new int*[N];
	// dasgleiche, Vektor innerhalb der Zelle
	r_rel = new double*[N];
	for(int i=0; i<N; i++){
		r_git[i] = new int[2];
		r_rel[i] = new double[2];
	}//for i
	
	//Kapillarkraefte. Erster Index TeilchenNr, zweiter Index Raumrichtung
	Fkap = new double*[N];
	//WCA-Kraefte, gleiche Indices
	F_WCA= new double*[N];
	//Zufallskraefte, gleiche Indices
	F_noise = new double*[N];
	for(int i=0; i<N; i++){
		Fkap[i] = new double[2];
		F_WCA[i]= new double[2];
		F_noise[i] = new double[2];
	}//for i
	
	// Teilchen auf Startpositionen setzen. Zunaechst sind das nur Zufallspositionen
	init_zufallspos();
	
	//init Kraftberechnungen
	WCA_init();
	kapkraefte_init();
	
}//void main_init





//Greensfunktion an der Stelle qjk. Index-Wrapping, Verschiebung.
double G(int j, int k){
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
		if(F < 1.0e-5) //eigentlich F==0.0. Konstante 10^-5 auslagern oder durch was anderes ausdruecken, zB F << F_noise ?
			ret = min(ret, mrd*mrd/(2.0*T*F_noise[teilchen][0]*F_noise[teilchen][0]));
		else{
			double u = F_noise[teilchen][0]*sqrt(0.5*T)/F;
			double w = sqrt(u*u + mrd/fabs(F));
			ret = min(ret, (w-fabs(u))*(w-fabs(u)));
		}//else F ungleich Null
		
		
		
		
		// y-Komponente
		F = F_WCA[teilchen][1] + Fkap[teilchen][1];
		if(F < 1.0e-5) //eigentlich F==0.0. Konstante 10^-5 auslagern oder durch was anderes ausdruecken, zB F << F_noise ?
			ret = min(ret, mrd*mrd/(2.0*T*F_noise[teilchen][1]*F_noise[teilchen][1]));
		else{
			double u = F_noise[teilchen][1]*sqrt(0.5*T)/F;
			double w = sqrt(u*u + mrd/fabs(F));
			ret = min(ret, (w-fabs(u))*(w-fabs(u)));
		}//else F ungleich Null
		
		
		
	}//for teilchen
	return ret;
}//double optimaler_zeitschritt 













