#pragma once

#include <cstdio>
#include <cmath>
#include "parameter.h"
#include <iostream>
#include <fftw3.h>
#include "signaturen.h"
#include "dynamik_methoden.cpp"

using std::cout; using std::endl;

extern const int N,densGrid_Zellen;
extern const double L,densGrid_Breite;
extern const double nachList_Breite;

extern int** r1_git;
extern double** r1_rel;
extern int** r2_git;
extern double** r2_rel;

extern const double f2_f1;

//Beitrag eines Kolloids zur Dichte
const double drho = 1.0/(densGrid_Breite*densGrid_Breite);

// berechnet aus Teilchenpositionen (Gewichtung je Typ als Parameter) die Dichte auf Gitterpunkten. Schreibt in Rho, verwendet Index-Wrapping. Nearest Grid Point.
void gridDensity_NGP(fftw_complex* rho, double gewicht1, double gewicht2){ 

	int i;
	double x,y;
	int k,l;

	//rho auf Null setzen
	for(k=0; k<densGrid_Zellen; k++)
		for(l=0; l<densGrid_Zellen; l++)
			rho[iw(k,l)][0] = 0.0;

	//Schleife ueber Teilchen Typ 1
	for(i=0; i<N1; i++){
		//x-Richtung
		x = r1_git[i][0]*nachList_Breite + r1_rel[i][0];
		// Gitterpunkt k ist bei x=(k+0.5)*breite
		// die x, die am dichtesten bei k sind, gehen von k*breite bis (k+1)*breite
		k = (int) (x/densGrid_Breite);



		//y-Richtung
		y = r1_git[i][1]*nachList_Breite + r1_rel[i][1];
		// die y, die am dichtesten bei l liegen, gegen von l*breite bis (l+1)*breite
		l = (int) (y/densGrid_Breite);

		//kompletten Anteil auf rho[k][l] addieren
		//rho[k][l] += 1.0;
		rho[iw(k,l)][0] += drho*gewicht1;
	}//for i
	
	//Schleife über Teilchen Typ 2
	if(0 == N2) return;
	for(i=0; i<N2; i++){
		//x-Richtung
		x = r2_git[i][0]*nachList_Breite + r2_rel[i][0];
		k = (int) (x/densGrid_Breite);
		
		//y-Richtung
		y = r2_git[i][1]*nachList_Breite + r2_rel[i][1];
		l = (int) (y/densGrid_Breite);
		
		rho[iw(k,l)][0] += drho*gewicht2;
	}//for i bis N2

}//void gridDens_NGP


// berechnet aus Teilchenpos und Kraeften an Gitterpunkten die Kraefte auf Teilchen. Schreibt in Fkap. Nearest Grid Point.
//void inv_gridDensity_NGP(double** Fkap, fftw_complex* Fx, fftw_complex* Fy, int** rr_git, double** rr_rel){
void inv_gridDensity_NGP(){
	
	using namespace dynamik_methoden;
	
	int i,k,l;
	double x,y;

/// setze Kraefte auf Null
	for(i=0; i<N1; i++){
		F1kap[i][0] = 0.0;
		F1kap[i][1] = 0.0;
	}//for i bis N1
	for(i=0; i<N2; i++){
		F2kap[i][0] = 0.0;
		F2kap[i][1] = 0.0;
	}//for i bis N2

/// Schleife ueber Teilchen Typ 1
	for(i=0; i<N1; i++){
		x = r1_git[i][0]*nachList_Breite + r1_rel[i][0];
		//dichtester Gitterpunkt
		k = (int)(x/densGrid_Breite);

		y = r1_git[i][1]*nachList_Breite + r1_rel[i][1];
		//dichtester Gitterpunkt
		l = (int)(y/densGrid_Breite);

		//kompletter Anteil ist Fx[k] und Fy[l]
		F1kap[i][0] = Fx[iw(k,l)][0];
		F1kap[i][1] = Fy[iw(k,l)][0];

	}//for i bis N1
	
	if(0 == N2) return;
	for(i=0; i<N2; i++){
		x = r2_git[i][0]*nachList_Breite + r2_rel[i][0];
		k = (int)(x/densGrid_Breite);
		
		y = r2_git[i][1]*nachList_Breite + r2_rel[i][1];
		l = (int)(y/densGrid_Breite);
		
		F2kap[i][0] = f2_f1 * Fx[iw(k,l)][0];
		F2kap[i][1] = f2_f1 * Fx[iw(k,l)][1];
	}//for i bis N2

}//void inv_gridDensity_NGP



// berechnet aus Teilchenpositionen die Dichte auf Gitterpunkten. Schreibt in Rho. Cloud-In-Cell.
void gridDensity_CIC(fftw_complex* rho, double gewicht1, double gewicht2){

	int i, k, l; //i Teilchen, k und l Gitter
	int kl, kr, ll, lr; //Indices der benachbarten Zellen (kl heisst "von k eins nach links"). Bedenke periodische Randbedingungen.
	double x,y;
	double xexakt, xlinks, xrechts; //Anteil des Teilchens in der linken, rechten und mittleren Zelle
	double yexakt, ylinks, yrechts; //dasselbe fuer die andere Richtung. Die Begriffe "links" und "rechts" sind irrefuehrend.

	double tmp;//Abstand des Teilchens vom Gittermittelpunkt. Im Aufschrieb ist dies d.

	//rho auf Null setzen
	for(k=0; k<densGrid_Zellen; k++)
		for(l=0; l<densGrid_Zellen; l++)
			rho[iw(k,l)][0] = 0.0;

	//Schleife ueber Teilchen Typ 1
	for(i=0; i<N1; i++){
		x = r1_git[i][0]*nachList_Breite + r1_rel[i][0];
		y = r1_git[i][1]*nachList_Breite + r1_rel[i][1];

		/// Indices der Zellen, die Anteile des Teilchens enthalten (koennen)
		//x-Richtung
		k = (int) (x/densGrid_Breite);
		kl = (k-1+densGrid_Zellen)%densGrid_Zellen;
		kr= (k+1)%densGrid_Zellen;

		//y-Richtung
		l = (int)(y/densGrid_Breite);
		ll = (l-1+densGrid_Zellen)%densGrid_Zellen;
		lr= (l+1)%densGrid_Zellen;


		///x-Richtung: Anteil des Teilchens in Zelle k
		tmp = x/densGrid_Breite - (k+0.5);
		xexakt = 1.0 - fabs(tmp);
		//falls tmp<0, ist das Teilchen teilw in linker Zelle. Sonst teilw in rechter Zelle.
		if(tmp<0.0){
			xlinks = -tmp;
			xrechts= 0.0;
		}//if tmp<0, dh Teilchen in der linken Haelfte der Zelle
		else{
			xlinks = 0.0;
			xrechts= tmp;
		}//else, dh Teilchen in der rechten Haelfte der Zelle

		///gleiches fuer y-Richtung
		tmp = y/densGrid_Breite - (l+0.5);
		yexakt = 1.0 - fabs(tmp);
		//Unterscheidung, ob Teilchen in der oberen oder unteren Haelfte der Zelle ist
		if(tmp<0.0){
			ylinks = -tmp;
			yrechts= 0.0;
		}//if tmp<0, dh Teilchen in der unteren Haelfte der Zelle
		else{
			ylinks = 0.0;
			yrechts= tmp;
		}//else, dh Teilchen in der oberen Haelfte der Zelle


		/// nun sind alle Informationen beisammen. Erhoehe die entsprechenden neun Zellen von rho[][]
		rho[iw(kl,ll)][0] += xlinks*ylinks  *drho*gewicht1;
		rho[iw(kl,l )][0] += xlinks*yexakt  *drho*gewicht1;
		rho[iw(kl,lr)][0] += xlinks*yrechts *drho*gewicht1;
		
		rho[iw(k ,ll)][0] += xexakt*ylinks  *drho*gewicht1;
		rho[iw(k ,l )][0] += xexakt*yexakt  *drho*gewicht1;
		rho[iw(k ,lr)][0] += xexakt*yrechts *drho*gewicht1;

		rho[iw(kr,ll)][0] += xrechts*ylinks *drho*gewicht1;
		rho[iw(kr,l )][0] += xrechts*yexakt *drho*gewicht1;
		rho[iw(kr,lr)][0] += xrechts*yrechts*drho*gewicht1;


		/* Test
			cout << "Teilchen "<<i<<" am Ort ("<<x<<","<<y<<"): x-Aufteilung ("<<xlinks<<" "<<xexakt<<" "<<xrechts<<"), y-Aufteilung ("<<ylinks<<" "<<yexakt<<" "<<yrechts<<")." << endl;
			cout << "Beteiligte Zellen: k=("<<kl<<" "<<k<<" "<<kr<<"), l=("<<ll<<" "<<l<<" "<<lr<<")." << endl << endl;
		// */


	}//for i
	
	if(0 == N2) return;
	//Schleife ueber Teilchen Typ 2
	for(i=0; i<N2; i++){
		x = r2_git[i][0]*nachList_Breite + r2_rel[i][0];
		y = r2_git[i][1]*nachList_Breite + r2_rel[i][1];

		/// Indices der Zellen, die Anteile des Teilchens enthalten (koennen)
		//x-Richtung
		k = (int) (x/densGrid_Breite);
		kl = (k-1+densGrid_Zellen)%densGrid_Zellen;
		kr= (k+1)%densGrid_Zellen;

		//y-Richtung
		l = (int)(y/densGrid_Breite);
		ll = (l-1+densGrid_Zellen)%densGrid_Zellen;
		lr= (l+1)%densGrid_Zellen;


		///x-Richtung: Anteil des Teilchens in Zelle k
		tmp = x/densGrid_Breite - (k+0.5);
		xexakt = 1.0 - fabs(tmp);
		//falls tmp<0, ist das Teilchen teilw in linker Zelle. Sonst teilw in rechter Zelle.
		if(tmp<0.0){
			xlinks = -tmp;
			xrechts= 0.0;
		}//if tmp<0, dh Teilchen in der linken Haelfte der Zelle
		else{
			xlinks = 0.0;
			xrechts= tmp;
		}//else, dh Teilchen in der rechten Haelfte der Zelle

		///gleiches fuer y-Richtung
		tmp = y/densGrid_Breite - (l+0.5);
		yexakt = 1.0 - fabs(tmp);
		//Unterscheidung, ob Teilchen in der oberen oder unteren Haelfte der Zelle ist
		if(tmp<0.0){
			ylinks = -tmp;
			yrechts= 0.0;
		}//if tmp<0, dh Teilchen in der unteren Haelfte der Zelle
		else{
			ylinks = 0.0;
			yrechts= tmp;
		}//else, dh Teilchen in der oberen Haelfte der Zelle


		/// nun sind alle Informationen beisammen. Erhoehe die entsprechenden neun Zellen von rho[][]
		rho[iw(kl,ll)][0] += xlinks*ylinks  *drho*gewicht2;
		rho[iw(kl,l )][0] += xlinks*yexakt  *drho*gewicht2;
		rho[iw(kl,lr)][0] += xlinks*yrechts *drho*gewicht2;
		
		rho[iw(k ,ll)][0] += xexakt*ylinks  *drho*gewicht2;
		rho[iw(k ,l )][0] += xexakt*yexakt  *drho*gewicht2;
		rho[iw(k ,lr)][0] += xexakt*yrechts *drho*gewicht2;

		rho[iw(kr,ll)][0] += xrechts*ylinks *drho*gewicht2;
		rho[iw(kr,l )][0] += xrechts*yexakt *drho*gewicht2;
		rho[iw(kr,lr)][0] += xrechts*yrechts*drho*gewicht2;


		/* Test
			cout << "Typ2-Teilchen "<<i<<" am Ort ("<<x<<","<<y<<"): x-Aufteilung ("<<xlinks<<" "<<xexakt<<" "<<xrechts<<"), y-Aufteilung ("<<ylinks<<" "<<yexakt<<" "<<yrechts<<")." << endl;
			cout << "Beteiligte Zellen: k=("<<kl<<" "<<k<<" "<<kr<<"), l=("<<ll<<" "<<l<<" "<<lr<<")." << endl << endl;
		// */


	}//for i


}//void gridDens_CIC


// berechnet aus Teilchenpos und Kraeften an Gitterpunkten die Kraefte auf Teilchen. Schreibt in Fkap. Cloud In Cell.
void inv_gridDensity_CIC(){

	using namespace dynamik_methoden;
// 	using dynamik_methoden::Fy;
// 	using dynamik_methoden::F1kap;
	
	int i, k, l, kl, kr, ll, lr;
	double x,y;

	double xexakt, xlinks, xrechts;
	double yexakt, ylinks, yrechts;
	double tmp;

	//alle Kraefte auf Null setzen
	for(i=0; i<N1; i++){
		F1kap[i][0] = 0.0;
		F1kap[i][1] = 0.0;
	}//for i bis N

	//Schleife ueber Teilchen Typ 1
	for(i=0; i<N1; i++){
		x = r1_git[i][0]*nachList_Breite + r1_rel[i][0];
		k = (int)(x/densGrid_Breite);
		kl = (k-1+densGrid_Zellen)%densGrid_Zellen;
		kr = (k+1)%densGrid_Zellen;

		y = r1_git[i][1]*nachList_Breite + r1_rel[i][1];
		l = (int)(y/densGrid_Breite);
		ll = (l-1+densGrid_Zellen)%densGrid_Zellen;
		lr = (l+1)%densGrid_Zellen;

		//x-Richtung: Anteil des Teilchens in Zelle k
		tmp = x/densGrid_Breite - (k+0.5);
		xexakt = 1.0 - fabs(tmp);
		//falls tmp<0, ist das Teilchen teilw in linker Zelle. Sonst teilw in rechter Zelle.
		if(tmp<0.0){
			xlinks = -tmp;
			xrechts= 0.0;
		}//if tmp<0, dh Teilchen in der linken Haelfte der Zelle
		else{
			xlinks = 0.0;
			xrechts= tmp;
		}//else, dh Teilchen in der rechten Haelfte der Zelle

		///gleiches fuer y-Richtung
		tmp = y/densGrid_Breite - (l+0.5);
		yexakt = 1.0 - fabs(tmp);
		//Unterscheidung, ob Teilchen in der oberen oder unteren Haelfte der Zelle ist
		if(tmp<0.0){
			ylinks = -tmp;
			yrechts= 0.0;
		}//if tmp<0, dh Teilchen in der unteren Haelfte der Zelle
		else{
			ylinks = 0.0;
			yrechts= tmp;
		}//else, dh Teilchen in der oberen Haelfte der Zelle

		//in xexakt,xlinks,xrechts,yexakt,ylinks,yrechts steckt alle benoetigte Information
		F1kap[i][0] =  xlinks*(ylinks*Fx[iw(kl,ll)][0]  +  yexakt*Fx[iw(kl,l)][0]  +  yrechts*Fx[iw(kl,lr)][0]);
		F1kap[i][1] =  xlinks*(ylinks*Fy[iw(kl,ll)][0]  +  yexakt*Fy[iw(kl,l)][0]  +  yrechts*Fy[iw(kl,lr)][0]);

		F1kap[i][0] += xexakt*(ylinks*Fx[iw(k ,ll)][0]  +  yexakt*Fx[iw(k ,l)][0]  +  yrechts*Fx[iw(k ,lr)][0]);
		F1kap[i][1] += xexakt*(ylinks*Fy[iw(k ,ll)][0]  +  yexakt*Fy[iw(k ,l)][0]  +  yrechts*Fy[iw(k ,lr)][0]);

		F1kap[i][0] +=xrechts*(ylinks*Fx[iw(kr,ll)][0]  +  yexakt*Fx[iw(kr,l)][0]  +  yrechts*Fx[iw(kr,lr)][0]);
		F1kap[i][1] +=xrechts*(ylinks*Fy[iw(kr,ll)][0]  +  yexakt*Fy[iw(kr,l)][0]  +  yrechts*Fy[iw(kr,lr)][0]);

	}//for i

	if(0 == N2) return;
	//Schleife ueber Teilchen Typ 2
	for(i=0; i<N2; i++){
		x = r2_git[i][0]*nachList_Breite + r2_rel[i][0];
		k = (int)(x/densGrid_Breite);
		kl = (k-1+densGrid_Zellen)%densGrid_Zellen;
		kr = (k+1)%densGrid_Zellen;

		y = r2_git[i][1]*nachList_Breite + r2_rel[i][1];
		l = (int)(y/densGrid_Breite);
		ll = (l-1+densGrid_Zellen)%densGrid_Zellen;
		lr = (l+1)%densGrid_Zellen;

		//x-Richtung: Anteil des Teilchens in Zelle k
		tmp = x/densGrid_Breite - (k+0.5);
		xexakt = 1.0 - fabs(tmp);
		//falls tmp<0, ist das Teilchen teilw in linker Zelle. Sonst teilw in rechter Zelle.
		if(tmp<0.0){
			xlinks = -tmp;
			xrechts= 0.0;
		}//if tmp<0, dh Teilchen in der linken Haelfte der Zelle
		else{
			xlinks = 0.0;
			xrechts= tmp;
		}//else, dh Teilchen in der rechten Haelfte der Zelle

		///gleiches fuer y-Richtung
		tmp = y/densGrid_Breite - (l+0.5);
		yexakt = 1.0 - fabs(tmp);
		//Unterscheidung, ob Teilchen in der oberen oder unteren Haelfte der Zelle ist
		if(tmp<0.0){
			ylinks = -tmp;
			yrechts= 0.0;
		}//if tmp<0, dh Teilchen in der unteren Haelfte der Zelle
		else{
			ylinks = 0.0;
			yrechts= tmp;
		}//else, dh Teilchen in der oberen Haelfte der Zelle

		//in xexakt,xlinks,xrechts,yexakt,ylinks,yrechts steckt alle benoetigte Information
		F2kap[i][0] =  xlinks*(ylinks*Fx[iw(kl,ll)][0]  +  yexakt*Fx[iw(kl,l)][0]  +  yrechts*Fx[iw(kl,lr)][0]);
		F2kap[i][1] =  xlinks*(ylinks*Fy[iw(kl,ll)][0]  +  yexakt*Fy[iw(kl,l)][0]  +  yrechts*Fy[iw(kl,lr)][0]);

		F2kap[i][0] += xexakt*(ylinks*Fx[iw(k ,ll)][0]  +  yexakt*Fx[iw(k ,l)][0]  +  yrechts*Fx[iw(k ,lr)][0]);
		F2kap[i][1] += xexakt*(ylinks*Fy[iw(k ,ll)][0]  +  yexakt*Fy[iw(k ,l)][0]  +  yrechts*Fy[iw(k ,lr)][0]);

		F2kap[i][0] +=xrechts*(ylinks*Fx[iw(kr,ll)][0]  +  yexakt*Fx[iw(kr,l)][0]  +  yrechts*Fx[iw(kr,lr)][0]);
		F2kap[i][1] +=xrechts*(ylinks*Fy[iw(kr,ll)][0]  +  yexakt*Fy[iw(kr,l)][0]  +  yrechts*Fy[iw(kr,lr)][0]);
		
		F2kap[i][0] *= f2_f1;
		F2kap[i][1] *= f2_f1;

	}//for i



}//void inv_gridDensity_CIC



// berechnet aus Teilchenpositionen die Dichte auf Gitterpunkten. Schreibt in Rho. Triangle-Shaped Cloud.
void gridDensity_TSC(fftw_complex* rho, double gewicht1, double gewicht2){
	
	int i, k, l; //i Teilchen, k und l Gitter
	int kl, kr, ll, lr; //Indices der benachbarten Zellen (kl heisst "von k eins nach links"). Bedenke periodische Randbedingungen.
	double x,y;
	double d; //Abstand des Teilchens vom Gittermittelpunkt
	double xexakt, xlinks, xrechts; //Anteil des Teilchens in der linken, rechten und mittleren Zelle
	double yexakt, ylinks, yrechts; //dasselbe fuer die andere Richtung. Die Begriffe "links" und "rechts" sind irrefuehrend.

	double tmp;//Abstand des Teilchens vom Mittelpunkt der Nachbarzelle. Im Aufschrieb ist dies dtilde oder dquer.

	//rho auf Null setzen
	for(k=0; k<densGrid_Zellen; k++)
		for(l=0; l<densGrid_Zellen; l++)
			rho[iw(k,l)][0] = 0.0;

	//Schleife ueber Teilchen Typ 1
	for(i=0; i<N1; i++){
		x = r1_git[i][0]*nachList_Breite + r1_rel[i][0];
		y = r1_git[i][1]*nachList_Breite + r1_rel[i][1];

		/// Indices der Zellen, die Anteile des Teilchens enthalten
		//x-Richtung
		k = (int) (x/densGrid_Breite);
		kl = (k-1+densGrid_Zellen)%densGrid_Zellen;
		kr= (k+1)%densGrid_Zellen;

		//y-Richtung
		l = (int)(y/densGrid_Breite);
		ll = (l-1+densGrid_Zellen)%densGrid_Zellen;
		lr= (l+1)%densGrid_Zellen;


		///x-Richtung: Anteil des Teilchens in Zelle k
		d = x - (k+0.5)*densGrid_Breite;
		xexakt = 0.75 - d*d/densGrid_Breite/densGrid_Breite;
		 // Anteil in linker Zelle
		tmp=fabs((d+densGrid_Breite)/densGrid_Breite);
		xlinks = 0.5*(1.5-tmp)*(1.5-tmp);
		 // Anteil in rechter Zelle
		tmp=fabs((d-densGrid_Breite)/densGrid_Breite);
		xrechts = 0.5*(1.5-tmp)*(1.5-tmp);

		///gleiches fuer y-Richtung
		d = y - (l+0.5)*densGrid_Breite;
		yexakt = 0.75 - d*d/densGrid_Breite/densGrid_Breite;
		//untere Zelle
		tmp=fabs((d+densGrid_Breite)/densGrid_Breite);
		ylinks = 0.5*(1.5-tmp)*(1.5-tmp);
		 // Anteil in rechter Zelle
		tmp=fabs((d-densGrid_Breite)/densGrid_Breite);
		yrechts = 0.5*(1.5-tmp)*(1.5-tmp);



		/// nun sind alle Informationen beisammen. Erhoehe die entsprechenden neun Zellen von rho[][]
		rho[iw(kl,ll)][0] += xlinks*ylinks  *drho*gewicht1;
		rho[iw(kl,l )][0] += xlinks*yexakt  *drho*gewicht1;
		rho[iw(kl,lr)][0] += xlinks*yrechts *drho*gewicht1;
		
		rho[iw(k ,ll)][0] += xexakt*ylinks  *drho*gewicht1;
		rho[iw(k ,l )][0] += xexakt*yexakt  *drho*gewicht1;
		rho[iw(k ,lr)][0] += xexakt*yrechts *drho*gewicht1;

		rho[iw(kr,ll)][0] += xrechts*ylinks *drho*gewicht1;
		rho[iw(kr,l )][0] += xrechts*yexakt *drho*gewicht1;
		rho[iw(kr,lr)][0] += xrechts*yrechts*drho*gewicht1;
		/* Test
			cout << "Teilchen "<<i<<" am Ort ("<<x<<","<<y<<"): x-Aufteilung ("<<xlinks<<" "<<xexakt<<" "<<xrechts<<"), y-Aufteilung ("<<ylinks<<" "<<yexakt<<" "<<yrechts<<")." << endl;
			cout << "Beteiligte Zellen: k=("<<kl<<" "<<k<<" "<<kr<<"), l=("<<ll<<" "<<l<<" "<<lr<<")." << endl << endl;
		// */


	}//for i
	
	if(0 == N2) return;
	//Schleife ueber Teilchen Typ 2
	for(i=0; i<N2; i++){
		x = r2_git[i][0]*nachList_Breite + r2_rel[i][0];
		y = r2_git[i][1]*nachList_Breite + r2_rel[i][1];

		/// Indices der Zellen, die Anteile des Teilchens enthalten
		//x-Richtung
		k = (int) (x/densGrid_Breite);
		kl = (k-1+densGrid_Zellen)%densGrid_Zellen;
		kr= (k+1)%densGrid_Zellen;

		//y-Richtung
		l = (int)(y/densGrid_Breite);
		ll = (l-1+densGrid_Zellen)%densGrid_Zellen;
		lr= (l+1)%densGrid_Zellen;


		///x-Richtung: Anteil des Teilchens in Zelle k
		d = x - (k+0.5)*densGrid_Breite;
		xexakt = 0.75 - d*d/densGrid_Breite/densGrid_Breite;
		 // Anteil in linker Zelle
		tmp=fabs((d+densGrid_Breite)/densGrid_Breite);
		xlinks = 0.5*(1.5-tmp)*(1.5-tmp);
		 // Anteil in rechter Zelle
		tmp=fabs((d-densGrid_Breite)/densGrid_Breite);
		xrechts = 0.5*(1.5-tmp)*(1.5-tmp);

		///gleiches fuer y-Richtung
		d = y - (l+0.5)*densGrid_Breite;
		yexakt = 0.75 - d*d/densGrid_Breite/densGrid_Breite;
		//untere Zelle
		tmp=fabs((d+densGrid_Breite)/densGrid_Breite);
		ylinks = 0.5*(1.5-tmp)*(1.5-tmp);
		 // Anteil in rechter Zelle
		tmp=fabs((d-densGrid_Breite)/densGrid_Breite);
		yrechts = 0.5*(1.5-tmp)*(1.5-tmp);



		/// nun sind alle Informationen beisammen. Erhoehe die entsprechenden neun Zellen von rho[][]
		rho[iw(kl,ll)][0] += xlinks*ylinks  *drho*gewicht2;
		rho[iw(kl,l )][0] += xlinks*yexakt  *drho*gewicht2;
		rho[iw(kl,lr)][0] += xlinks*yrechts *drho*gewicht2;
		
		rho[iw(k ,ll)][0] += xexakt*ylinks  *drho*gewicht2;
		rho[iw(k ,l )][0] += xexakt*yexakt  *drho*gewicht2;
		rho[iw(k ,lr)][0] += xexakt*yrechts *drho*gewicht2;

		rho[iw(kr,ll)][0] += xrechts*ylinks *drho*gewicht2;
		rho[iw(kr,l )][0] += xrechts*yexakt *drho*gewicht2;
		rho[iw(kr,lr)][0] += xrechts*yrechts*drho*gewicht2;
		/* Test
			cout << "Teilchen "<<i<<" am Ort ("<<x<<","<<y<<"): x-Aufteilung ("<<xlinks<<" "<<xexakt<<" "<<xrechts<<"), y-Aufteilung ("<<ylinks<<" "<<yexakt<<" "<<yrechts<<")." << endl;
			cout << "Beteiligte Zellen: k=("<<kl<<" "<<k<<" "<<kr<<"), l=("<<ll<<" "<<l<<" "<<lr<<")." << endl << endl;
		// */


	}//for i
	
}//void gridDens_TSC


// berechnet aus Teilchenpos und Kraeften an Gitterpunkten die Kraefte auf Teilchen. Schreibt in Fkap. Triangle Shaped Cloud.
void inv_gridDensity_TSC(){
	using namespace dynamik_methoden;

	int i, k, l, kl, kr, ll, lr;
	double x,y;

	double xexakt, xlinks, xrechts;
	double yexakt, ylinks, yrechts;
	double tmp,d;

	//alle Kraefte auf Null setzen
	for(i=0; i<N1; i++){
		F1kap[i][0] = 0.0;
		F1kap[i][1] = 0.0;
	}//for i bis N1
	for(i=0; i<N2; i++){
		F2kap[i][0] = 0.0;
		F2kap[i][1] = 0.0;
	}//for i bis N2

	//Schleife ueber Teilchen
	for(i=0; i<N1; i++){
		x = r1_git[i][0]*nachList_Breite + r1_rel[i][0];
		k = (int)(x/densGrid_Breite);
		kl = (k-1+densGrid_Zellen)%densGrid_Zellen;
		kr = (k+1)%densGrid_Zellen;

		y = r1_git[i][1]*nachList_Breite + r1_rel[i][1];
		l = (int)(y/densGrid_Breite);
		ll = (l-1+densGrid_Zellen)%densGrid_Zellen;
		lr = (l+1)%densGrid_Zellen;

		///x-Richtung: Anteil des Teilchens in Zelle k
		d = x - (k+0.5)*densGrid_Breite;
		xexakt = 0.75 - d*d/densGrid_Breite/densGrid_Breite;
		 // Anteil in linker Zelle
		tmp=fabs((d+densGrid_Breite)/densGrid_Breite);
		xlinks = 0.5*(1.5-tmp)*(1.5-tmp);
		 // Anteil in rechter Zelle
		tmp=fabs((d-densGrid_Breite)/densGrid_Breite);
		xrechts = 0.5*(1.5-tmp)*(1.5-tmp);

		///gleiches fuer y-Richtung
		d = y - (l+0.5)*densGrid_Breite;
		yexakt = 0.75 - d*d/densGrid_Breite/densGrid_Breite;
		//untere Zelle
		tmp=fabs((d+densGrid_Breite)/densGrid_Breite);
		ylinks = 0.5*(1.5-tmp)*(1.5-tmp);
		 // Anteil in rechter Zelle
		tmp=fabs((d-densGrid_Breite)/densGrid_Breite);
		yrechts = 0.5*(1.5-tmp)*(1.5-tmp);


		
		//in xexakt,xlinks,xrechts,yexakt,ylinks,yrechts steckt alle benoetigte Information
		F1kap[i][0] =  xlinks*(ylinks*Fx[iw(kl,ll)][0]  +  yexakt*Fx[iw(kl,l)][0]  +  yrechts*Fx[iw(kl,lr)][0]);
		F1kap[i][1] =  xlinks*(ylinks*Fy[iw(kl,ll)][0]  +  yexakt*Fy[iw(kl,l)][0]  +  yrechts*Fy[iw(kl,lr)][0]);

		F1kap[i][0] += xexakt*(ylinks*Fx[iw(k ,ll)][0]  +  yexakt*Fx[iw(k ,l)][0]  +  yrechts*Fx[iw(k ,lr)][0]);
		F1kap[i][1] += xexakt*(ylinks*Fy[iw(k ,ll)][0]  +  yexakt*Fy[iw(k ,l)][0]  +  yrechts*Fy[iw(k ,lr)][0]);

		F1kap[i][0] +=xrechts*(ylinks*Fx[iw(kr,ll)][0]  +  yexakt*Fx[iw(kr,l)][0]  +  yrechts*Fx[iw(kr,lr)][0]);
		F1kap[i][1] +=xrechts*(ylinks*Fy[iw(kr,ll)][0]  +  yexakt*Fy[iw(kr,l)][0]  +  yrechts*Fy[iw(kr,lr)][0]);

	}//for i

	if(0 == N2) return;
	//Schleife ueber Teilchen
	for(i=0; i<N2; i++){
		x = r2_git[i][0]*nachList_Breite + r2_rel[i][0];
		k = (int)(x/densGrid_Breite);
		kl = (k-1+densGrid_Zellen)%densGrid_Zellen;
		kr = (k+1)%densGrid_Zellen;

		y = r2_git[i][1]*nachList_Breite + r2_rel[i][1];
		l = (int)(y/densGrid_Breite);
		ll = (l-1+densGrid_Zellen)%densGrid_Zellen;
		lr = (l+1)%densGrid_Zellen;

		///x-Richtung: Anteil des Teilchens in Zelle k
		d = x - (k+0.5)*densGrid_Breite;
		xexakt = 0.75 - d*d/densGrid_Breite/densGrid_Breite;
		 // Anteil in linker Zelle
		tmp=fabs((d+densGrid_Breite)/densGrid_Breite);
		xlinks = 0.5*(1.5-tmp)*(1.5-tmp);
		 // Anteil in rechter Zelle
		tmp=fabs((d-densGrid_Breite)/densGrid_Breite);
		xrechts = 0.5*(1.5-tmp)*(1.5-tmp);

		///gleiches fuer y-Richtung
		d = y - (l+0.5)*densGrid_Breite;
		yexakt = 0.75 - d*d/densGrid_Breite/densGrid_Breite;
		//untere Zelle
		tmp=fabs((d+densGrid_Breite)/densGrid_Breite);
		ylinks = 0.5*(1.5-tmp)*(1.5-tmp);
		 // Anteil in rechter Zelle
		tmp=fabs((d-densGrid_Breite)/densGrid_Breite);
		yrechts = 0.5*(1.5-tmp)*(1.5-tmp);


		
		//in xexakt,xlinks,xrechts,yexakt,ylinks,yrechts steckt alle benoetigte Information
		F2kap[i][0] =  xlinks*(ylinks*Fx[iw(kl,ll)][0]  +  yexakt*Fx[iw(kl,l)][0]  +  yrechts*Fx[iw(kl,lr)][0]);
		F2kap[i][1] =  xlinks*(ylinks*Fy[iw(kl,ll)][0]  +  yexakt*Fy[iw(kl,l)][0]  +  yrechts*Fy[iw(kl,lr)][0]);

		F2kap[i][0] += xexakt*(ylinks*Fx[iw(k ,ll)][0]  +  yexakt*Fx[iw(k ,l)][0]  +  yrechts*Fx[iw(k ,lr)][0]);
		F2kap[i][1] += xexakt*(ylinks*Fy[iw(k ,ll)][0]  +  yexakt*Fy[iw(k ,l)][0]  +  yrechts*Fy[iw(k ,lr)][0]);

		F2kap[i][0] +=xrechts*(ylinks*Fx[iw(kr,ll)][0]  +  yexakt*Fx[iw(kr,l)][0]  +  yrechts*Fx[iw(kr,lr)][0]);
		F2kap[i][1] +=xrechts*(ylinks*Fy[iw(kr,ll)][0]  +  yexakt*Fy[iw(kr,l)][0]  +  yrechts*Fy[iw(kr,lr)][0]);
		
		F2kap[i][0] *= f2_f1;
		F2kap[i][1] *= f2_f1;

	}//for i

	
}//void inv_gridDensity_TSC
