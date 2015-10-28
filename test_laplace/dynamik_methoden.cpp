// stellt Methoden bereit, die die Bewegung der Kolloide behandeln.

#pragma once

#include "parameter.h"
#include "signaturen.h"
#include "zufall.cpp"
#include "gridRoutinen.cpp"
#include <fftw3.h>
#include <math.h>
#include <cstdio>
#include <iostream>

using std::cout; using std::endl;


extern const double densGrid_Breite, dq, lambda_kapillar;
extern const int densGrid_Zellen, densGrid_Schema;

const int Z = densGrid_Zellen;
const double dx = densGrid_Breite;

//semi-globale Felder/Variablen, erreichbar fuer die Methoden in dieser Datei
fftw_complex* rhox = NULL; // Dichte vor  Fouriertrafo FORWARD. Index-Wrapping
fftw_complex* rhok = NULL; // Dichte nach Fouriertrafo FORWARD. Index-Wrapping, Verschiebung und alternierende Vorzeichen
double* sinxG = NULL;	   // der Term sin(qx dx)/dx G(q). Index-Wrapping und Verschiebung
double* sinyG = NULL;	   // der Term sin(qy dx)/dx G(q). Index-Wrapping und Verschiebung
fftw_complex* Fxk = NULL;  // rhok mal i sinxG. Index-Wrapping, Verschiebung und alternierende Vorzeichen
fftw_complex* Fyk = NULL;  // rhok mal i sinyG. Index-Wrapping, Verschiebung und alternierende Vorzeichen
fftw_complex* Fx = NULL;   // Kraefte (x-komp) nach Fouriertrafo BACKWARD. Index-Wrapping
fftw_complex* Fy = NULL;   // Kraefte (y-komp) nach Fouriertrafo BACKWARD. Index-Wrapping

fftw_plan forward_plan, backx_plan, backy_plan;

//Greensfunktion an der Stelle qjk. Index-Wrapping, Verschiebung.
double G(int j, int k);



//berechnet Kapillarkraefte (mittels Fouriertransformation), schreibt sie in Fkap
void berechne_kapkraefte(double** r, double** Fkap){
	int j,l;
	double x,y;

/// Density-Gridding: Schreibe Dichte in rhox[][0]
	switch(densGrid_Schema){
		case 0: gridDensity_NGP(rhox, r); break;
		case 1: gridDensity_CIC(rhox, r); break;
		case 2: gridDensity_TSC(rhox, r); break;
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

	//* TESTE LAPLACE: rhox = irgendeine Funktion. laprho = laplace von der Funktion
	double* laprho = new double[Z*Z];
	for(j=0; j<Z; j++)
	for(l=0; l<Z; l++){
		x=j*dx;
		y=l*dx;
		//rhox[iw(j,l)][0] = exp(1.5*cos(4*M_PI*x/L) + 2*cos(2*M_PI*y/L));
		//laprho[iw(j,l)] = -4*M_PI*M_PI/L/L * rhox[iw(j,l)][0] * (6*cos(4*M_PI*x/L) - 9*sin(4*M_PI*x/L)*sin(4*M_PI*x/L) + 2*cos(2*M_PI*y/L) - 4*sin(2*M_PI*y/L)*sin(2*M_PI*y/L));
		rhox[iw(j,l)][0] = cos(2*M_PI/L * (x+y));
		laprho[iw(j,l)] = -8*M_PI*M_PI/L/L * rhox[iw(j,l)][0];
	
		
		rhox[iw(j,l)][1] = 0.0;
	}//for j,l
	// */

/// FFT: Schreibe transformierte Dichte in rhok[][]
	fftw_execute(forward_plan);

/// Multiplikation mit i sin()/dx G: Schreibe Ergebnis in Fxk[][] und Fyk[][]. Komplexe Multiplikation, Gsinx rein imaginaer

	for(j=0; j<Z*Z; j++){
		//Fxk[j][0] = - sinxG[j] * rhok[j][1]; // -im*im
		//Fxk[j][1] =   sinxG[j] * rhok[j][0]; // im*re
		Fxk[j][0] =   sinxG[j] * rhok[j][0];
		Fxk[j][1] =   sinxG[j] * rhok[j][1];

		Fyk[j][0] = - sinyG[j] * rhok[j][1]; 
		Fyk[j][1] =   sinyG[j] * rhok[j][0];
	}//for i



/// FFT-1: Schreibe ruecktransformierte Kraefte in Fx[][0] und Fy[][0]
	fftw_execute(backx_plan);
	fftw_execute(backy_plan);

//* Test: Schreibe Kraefte (auf Gitterzellen) in Dateien Fx.txt und Fy.txt. Und die imaginaerteile (sollten Null sein, weil rho reell und sinG symmetrisch).
	FILE* out  = fopen("Fx.txt", "w");
	FILE* out2 = fopen("Fy.txt", "w");
	FILE* outrk= fopen("Grk.txt", "w");


	//FILE* outu =fopen("u_entlang_xAchse.txt", "w");
	//FILE* outF =fopen("F_entlang_xAchse.txt", "w");
	for(j=0; j<Z; j++){
		for(l=0; l<Z; l++){
			//					              x     y    ruecktransformiert    f                 laplace f
			fprintf(out, "%g \t %g \t %g \t %g \t %g \t %g \n", j*dx, l*dx, Fx[iw(j,l)][0]/Z/Z, rhox[iw(j,l)][0], laprho[iw(j,l)], Fx[iw(j,l)][0]/Z/Z - laprho[iw(j,l)]);
			fprintf(out2, "%g \t %g \t %g \n", j*dx, l*dx, Fy[iw(j,l)][0]/Z/Z);
			//fprintf(imagy, "%g \t %g \t %g \n", j*dx, l*dx, Fy[iw(j,l)][1]);
			double qx = dq*((j+(int)(0.5*Z))%Z - 0.5*Z);
			double qy = dq*((l+(int)(0.5*Z))%Z - 0.5*Z);
			fprintf(outrk,"%g \t %g \t %g \t %g \n", qx, qy, Fxk[iw(j,l)][0]/Z/Z, Fxk[iw(j,l)][1]/Z/Z);
		}//for l
		fprintf(out, "\n");
		fprintf(out2, "\n");
		fprintf(outrk, "\n");
		//fprintf(imagy, "\n");

		//fprintf(outu, "%g \t %g \n", j*dx, Fx[iw(j,Z/2-1)][0]);
		//fprintf(outF, "%g \t %g \n", j*dx, Fy[iw(j,Z/2-1)][0]);
	}//for j
	fclose(out);
	fclose(out2);
	//fclose(imagx);
	//fclose(imagy);
// */

/// inverse Density-Gridding: Schreibe Kraefte in Fkap
	switch(densGrid_Schema){
		case 0: inv_gridDensity_NGP(Fkap, Fx, Fy, r); break;
		case 1: inv_gridDensity_CIC(Fkap, Fx, Fy, r); break;
		case 2: inv_gridDensity_TSC(Fkap, Fx, Fy, r); break;
		default: cout << "Density-Gridding-Schema (NearestGridPoint/CloudInCell/TSC) nicht erkannt! densGrid_Schema="<<densGrid_Schema<<", zulaessig sind nur 0,1,2." << endl; 
	}//switch
	

}//void berechne_kapkraefte











//initialisiert Felder, plant Fouriertrafos
void kapkraefte_init(){

/// alloziere Felder
	//Dichte vor FFT. Index-Wrapping. rhox[iw(j,k)][0] ist die Dichte in Zelle (j,k). rhox[][1] ist der Imaginaerteil, Null.
	rhox = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Z*Z);

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
		//q = dq* ((j+(int)(Z*0.5))%Z - 0.5*Z );
		//sinxG[iw(j,k)] = sin(q*dx)/dx * G(j,k);
		sinxG[iw(j,k)] = G(j,k);

		q = dq* ((k+(int)(Z*0.5))%Z - 0.5*Z );
		sinyG[iw(j,k)] = sin(q*dx)/dx * G(j,k);
	}//for j,k

/// plane FFTs
	forward_plan = fftw_plan_dft_2d(Z, Z, rhox, rhok, FFTW_FORWARD,  FFTW_ESTIMATE);
	backx_plan   = fftw_plan_dft_2d(Z, Z, Fxk,  Fx,   FFTW_BACKWARD, FFTW_ESTIMATE);
	backy_plan   = fftw_plan_dft_2d(Z, Z, Fyk,  Fy,   FFTW_BACKWARD, FFTW_ESTIMATE);
}//void kapkraefte_init



//Greensfunktion an der Stelle qjk. Index-Wrapping, Verschiebung.
double G(int j, int k){
	double qx = dq*(((int)(j+Z*0.5))%Z - 0.5*Z );
	double qy = dq*(((int)(k+Z*0.5))%Z - 0.5*Z );
	if(j==0 && k==0) return 0.0;
	//return 1.0/( sin(0.5*qx*dx)*sin(0.5*qx*dx)/dx/dx + sin(0.5*qy*dx)*sin(0.5*qy*dx)/dx/dx + 1.0/(lambda_kapillar*lambda_kapillar) );
	
	return - qx*qx - qy*qy;
	//return - 4*sin(0.5*qx*dx)*sin(0.5*qx*dx)/dx/dx - 4*sin(0.5*qy*dx)*sin(0.5*qy*dx)/dx/dx;
}//double G

















