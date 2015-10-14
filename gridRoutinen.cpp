#pragma once

#include <cstdio>
#include <cmath>
#include "parameter.h"
#include <iostream>

using std::cout; using std::endl;

extern const int N,densGrid_Zellen;
extern const double L,densGrid_Breite;

// berechnet aus Teilchenpositionen die Dichte auf Gitterpunkten. Schreibt in Rho. Nearest Grid Point.
void gridDensity_NGP(double** rho, double** r){

	int i;
	double x,y;
	int k,l;

	//rho auf Null setzen
	for(k=0; k<densGrid_Zellen; k++)
		for(l=0; l<densGrid_Zellen; l++)
			rho[k][l] = 0.0;

	//Schleife ueber Teilchen
	for(i=0; i<N; i++){
		//x-Richtung
		x = r[i][0];
		// Gitterpunkt k ist bei x=(k+0.5)*breite
		// die x, die am dichtesten bei k sind, gehen von k*breite bis (k+1)*breite
		k = (int) (x/densGrid_Breite);



		//y-Richtung
		y = r[i][1];
		// die y, die am dichtesten bei l liegen, gegen von l*breite bis (l+1)*breite
		l = (int) (y/densGrid_Breite);

		//kompletten Anteil auf rho[k][l] addieren
		rho[k][l] += 1.0;
	}//for i

}//void gridDens_NGP

// berechnet aus Teilchenpositionen die Dichte auf Gitterpunkten. Schreibt in Rho. Cloud-In-Cell.
void gridDensity_CIC(double** rho, double** r){

	int i, k, l; //i Teilchen, k und l Gitter
	int kl, kr, ll, lr; //Indices der benachbarten Zellen (kl heisst "von k eins nach links"). Bedenke periodische Randbedingungen.
	double x,y;
	double xexakt, xlinks, xrechts; //Anteil des Teilchens in der linken, rechten und mittleren Zelle
	double yexakt, ylinks, yrechts; //dasselbe fuer die andere Richtung. Die Begriffe "links" und "rechts" sind irrefuehrend.

	double tmp;//Abstand des Teilchens vom Gittermittelpunkt. Im Aufschrieb ist dies d.

	//rho auf Null setzen
	for(k=0; k<densGrid_Zellen; k++)
		for(l=0; l<densGrid_Zellen; l++)
			rho[k][l] = 0.0;

	//Schleife ueber Teilchen
	for(i=0; i<N; i++){
		x = r[i][0];
		y = r[i][1];

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
		rho[kl][ll] += xlinks*ylinks;
		rho[kl][l]  += xlinks*yexakt;
		rho[kl][lr] += xlinks*yrechts;

		rho[k ][ll] += xexakt*ylinks;
		rho[k ][l]  += xexakt*yexakt;
		rho[k ][lr] += xexakt*yrechts;

		rho[kr][ll] += xrechts*ylinks;
		rho[kr][l]  += xrechts*yexakt;
		rho[kr][lr] += xrechts*yrechts;

		/* Test
			cout << "Teilchen "<<i<<" am Ort ("<<x<<","<<y<<"): x-Aufteilung ("<<xlinks<<" "<<xexakt<<" "<<xrechts<<"), y-Aufteilung ("<<ylinks<<" "<<yexakt<<" "<<yrechts<<")." << endl;
			cout << "Beteiligte Zellen: k=("<<kl<<" "<<k<<" "<<kr<<"), l=("<<ll<<" "<<l<<" "<<lr<<")." << endl << endl;
		// */


	}//for i

}//void gridDens_CIC

// berechnet aus Teilchenpositionen die Dichte auf Gitterpunkten. Schreibt in Rho. Triangle-Shaped Cloud.
void gridDensity_TSC(double** rho, double** r){
	
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
			rho[k][l] = 0.0;

	//Schleife ueber Teilchen
	for(i=0; i<N; i++){
		x = r[i][0];
		y = r[i][1];

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
		rho[kl][ll] += xlinks*ylinks;
		rho[kl][l]  += xlinks*yexakt;
		rho[kl][lr] += xlinks*yrechts;

		rho[k ][ll] += xexakt*ylinks;
		rho[k ][l]  += xexakt*yexakt;
		rho[k ][lr] += xexakt*yrechts;

		rho[kr][ll] += xrechts*ylinks;
		rho[kr][l]  += xrechts*yexakt;
		rho[kr][lr] += xrechts*yrechts;

		/* Test
			cout << "Teilchen "<<i<<" am Ort ("<<x<<","<<y<<"): x-Aufteilung ("<<xlinks<<" "<<xexakt<<" "<<xrechts<<"), y-Aufteilung ("<<ylinks<<" "<<yexakt<<" "<<yrechts<<")." << endl;
			cout << "Beteiligte Zellen: k=("<<kl<<" "<<k<<" "<<kr<<"), l=("<<ll<<" "<<l<<" "<<lr<<")." << endl << endl;
		// */


	}//for i
}//void gridDens_TSC
