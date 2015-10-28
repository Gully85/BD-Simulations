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

/// Density-Gridding: Schreibe Dichte in rhox[][0]
	switch(densGrid_Schema){
		case 0: gridDensity_NGP(rhox, r); break;
		case 1: gridDensity_CIC(rhox, r); break;
		case 2: gridDensity_TSC(rhox, r); break;
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
		q = dq* ((j+(int)(Z*0.5))%Z - 0.5*Z );
		//sinxG[iw(j,k)] = sin(q*dx)/dx * G(j,k);
		sinxG[iw(j,k)] = G(j,k);

		q = dq* ((k+(int)(Z*0.5))%Z - 0.5*Z );
		sinyG[iw(j,k)] = sin(q*dx)/dx * G(j,k);
	}//for j,k

/// plane FFTs
	forward_plan = fftw_plan_dft_2d(Z, Z, rhox, rhok, FFTW_FORWARD,  FFTW_MEASURE);
	backx_plan   = fftw_plan_dft_2d(Z, Z, Fxk,  Fx,   FFTW_BACKWARD, FFTW_MEASURE);
	backy_plan   = fftw_plan_dft_2d(Z, Z, Fyk,  Fy,   FFTW_BACKWARD, FFTW_MEASURE);
}//void kapkraefte_init



//Greensfunktion an der Stelle qjk. Index-Wrapping, Verschiebung.
double G(int j, int k){
	double qx = dq*(((int)(j+Z*0.5))%Z - 0.5*Z );
	double qy = dq*(((int)(k+Z*0.5))%Z - 0.5*Z );
	if(lambda_kapillar < 0.01 && j==0 && k==0) return 0.0;
	return 1.0/( 4*sin(0.5*qx*dx)*sin(0.5*qx*dx)/dx/dx + 4*sin(0.5*qy*dx)*sin(0.5*qy*dx)/dx/dx + 1.0/(lambda_kapillar*lambda_kapillar) );
}//double G

















