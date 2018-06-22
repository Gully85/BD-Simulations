// enthaelt alle veraenderbaren Parameter

#pragma once

#define _USE_MATH_DEFINES
#include <math.h>
#include <string>

using std::string;

/// Parameter, die geaendert werden duerfen

// Teilchenzahl
const int N1 = 1000; //darf nicht Null sein!
const int N2 = 1000;

// Groesse der Simulationsbox (je Raumrichtung) in Einheiten von sigma
const double L = 20.0*sqrt(10.0); // etwa 63

// 0 fuer Nearest Grid Point, 1 fuer Cloud In Cell, 2 fuer Triangular Shaped Cloud
const int densGrid_Schema = 2; 

// Anzahl Samples in der Fouriertransformation. Nur gerade Zahlen erlaubt. Zweierpotenzen für gute Performance.
const int FFT_samples = 256; 


// kapillarlaenge in Einheiten von sigma
const double lambda_kapillar = 10.0; 

// Vorfaktor der Kapillarkraft. Ist f^2/(eps gamma).
const double kapillar_vorfaktor = 0.1; //0.89 ist V0/kT im EPJE-Paper. Es gilt kap_vorfaktor=2pi*V0*T.


// Temperatur, in Einheiten kT/eps
const double T = 20.0;


// maximaler Zeitschritt
const double dt_max = 0.05;

//maximale Reisedistanz in einem Zeitschritt
const double max_reisedistanz = 0.002;

//Verhältnis der Radien Typ1 und Typ2. Sigma11/Sigma22
const double sigma11_22 = 1.0;

//Verhältnis der Radien Sigma11/Sigma12. additive Mischung erfüllt s12 = 0.5(s11+s22)
const double sigma11_12 = 0.5*(1.0 + sigma11_22);

// Verhältnis der Mobilitäten Gamma2/Gamma1. Falls Reibung ~ Radius ist, ist Mobilität ~ 1/Radius
const double Gamma2_1 = sigma11_22;

// Verhältnis der Potentialstärken eps22/eps11
const double eps22_11 = 1.0;

//analog eps12/eps11
const double eps12_11 = 1.0;


//Verhältnis der Kräfte im Referenzzustand
//const double f2_f1 = 1.0/sigma11_22; //gleiche Dichten und Oberflächen
const double f2_f1 = - 1.0; //insgesamt gleichviel Attraktion wie Repulsion


//Warnung ausgeben, wenn WCA-Kräfte auftreten
const bool noWCA = false;
//keine Zufallskräfte, dh keine Diffusion
const bool noRNG = false;
//FFTW nicht so aufwändig initialisieren
const bool quickInit = false;
//Bewegung der Teilchen einschränken: Nur radial, Winkel zur x-Achse nach jedem Zeitschritt korrigieren
const bool restrictRadial = false;



//mitteln über wie viele runs?
const int runs = 1;

//wie viele Jobs soll es geben?
const int jobs = 1;

//wie viele Threads sollen gestartet werden, dh wie viele CPUs verwendet? 0 für unbegrenzt
const int maxThreads = 1;

//Observable aufnehmen in welchem Zeitabstand?
//const double obs_dt = 0.05; //alter Wert
const double obs_dt = 0.05; //

// Obervable aufnehmen wie oft?
const int obs_anzahl = 20*20*10; //


// Binbreite Paarkorrelationsfunktion
const double korr_dr = 0.5;

// größtes r der Paarkorrelationsfunktion
const double korr_rmax = 200.0;


// Binbreite Dichteprofil
const double rhoring_dr = 0.5;

// Anzahl Bins im Dichteprofil
const int rhoring_bins = 500;


//größtes q in rhothilde, in Einheiten dq=4pi/L
const int ftrho_qbins = 30;


// Start: Zufall=1, aus Datei=2, Gitterstart=3, allegleich=4, Kreisscheibe=5, Kern+Ring=6
const int startpos_methode=1;
//falls aus Datei: Name der Datei
const string startpos_dateiname="startpos.txt";
//falls Kreisscheibe oder Kern+Ring: Radius des Kreisss
//const double startpos_kreisradius = 91.6;
const double startpos_kreisradius = 4.5 * lambda_kapillar; //vorhergesagtes Gleichgewicht ist bei 5.31379 lambda

//falls Positionen in EINE Datei geschrieben werden: Name der Datei
const string pos_output_dateiname="pos.txt";


// Welche Observablen sollen aufgenommen werden?
const bool auswerten_korrfunk = false;
const bool auswerten_korrfunk_mixed = false;
const bool auswerten_rhovonk= false; //via Gitter im k-Raum
const bool auswerten_rhoviaFFTW = false;
const bool auswerten_rhoFT_normjerun=false; //nur falls auswerten_rhoviaFFTW gesetzt ist.
const bool auswerten_animation = true; //schreibt im ersten Durchgang Schnappschüsse in Dateien pos1.txt und pos2.txt
const bool auswerten_abstand = false; //Mittelwert des Abstands Typ2-Teilchen von der Boxmitte


/////////////// AB HIER: AENDERN VERBOTEN! /////////////////////////


// die Konstante 2^(1/6)
const double zweihoch1_6 = pow(2.0, 1.0/6.0);

// Anzahl Zellen beim Nachbarlisting
const int nachList_Zellen = (int) floor(L/zweihoch1_6);

// Breite jeder Zelle beim Nachbarlisting
const double nachList_Breite = L/((double)nachList_Zellen);

// Breite jeder Gitterzelle beim density-Gridding
const double densGrid_Breite = L/FFT_samples;

// Anzahl Gitterzellen (je Raumrichtung) fuer Density Gridding (Teilchenpositionen |-> Dichte auf Gitterpunkten)
const int densGrid_Zellen = FFT_samples; //vorher "G"

// Breite jeder Gitterzelle im reziproken Raum
const double dq = 2.0*M_PI/L;

// Binbreite für rho via FFTW. Sollte nicht kleiner sein als dq.
const double dq_rhoFFTW = 2.0*dq;

//Verhältnis dieser beiden dq, sollte kleiner als Eins sein
const double dq_dqFT = dq / dq_rhoFFTW;

//Anzahl Bins für rho(k) via FFTW
const int rhoFFTW_bins = 0.5*FFT_samples*dq/dq_rhoFFTW;

// Abtastrate bzw Bandbreite im reziproken Raum. q's gehen (komponentenweise) von -Omega bis +Omega
const double Omega = M_PI/densGrid_Breite;

// Anzahl Bins Paarkorrelationsfunktion
const int korr_bins = (int) (korr_rmax/korr_dr);

//Breite dq für rhothilde
const double ftrho_dq = 4.0*M_PI/L;
