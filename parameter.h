// enthaelt alle veraenderbaren Parameter

#pragma once

#define _USE_MATH_DEFINES
#include <math.h>
#include <string>

using std::string;

/// Parameter, die geaendert werden duerfen

// Teilchenzahl
const int N1 = 2000; //darf nicht Null sein!
const int N2 = 2000;

// Groesse der Simulationsbox (je Raumrichtung) in Einheiten von sigma
const double L = 300.0; 

// 0 fuer Nearest Grid Point, 1 fuer Cloud In Cell, 2 fuer Triangular Shaped Cloud
const int densGrid_Schema = 2; 

// Anzahl Samples in der Fouriertransformation. Nur gerade Zahlen erlaubt. Zweierpotenzen für gute Performance.
const int FFT_samples = 384; 


// kapillarlaenge in Einheiten von sigma
const double lambda_kapillar = 25.0; 

// Vorfaktor der Kapillarkraft. Ist f^2/(eps gamma).
const double kapillar_vorfaktor = 1.111111111; //0.89 ist V0/kT im EPJE-Paper. Es gilt kap_vorfaktor=2pi*V0*T.
//const double kapillar_vorfaktor = 0.0;


// Temperatur, in Einheiten kT/eps
const double T = 24.3055555;


// maximaler Zeitschritt
const double dt_max = 0.005;

//maximale Reisedistanz in einem Zeitschritt
const double max_reisedistanz = 0.24;

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

// nur ideales Gas. Weder Kapillar- noch WCA-Kräfte
const bool idgas_only = false;
//Warnung ausgeben, wenn WCA-Kräfte auftreten
const bool noWCA = false;
//keine Zufallskräfte, dh keine Diffusion
const bool noRNG = false;
//FFTW nicht so aufwändig initialisieren
const bool quickInit = true;
//Bewegung der Teilchen einschränken: Nur radial, Winkel zur x-Achse nach jedem Zeitschritt korrigieren
const bool restrictRadial = false;
//macht nur wenige Zeitschritte, schreibt zusätzliche Infos, vor allem über Kraftberechnung
const bool debugmode = false;
//macht keine Zeitschritte. Zwei Teilchen, werden auf Positionen gesetzt und Kraft ausgerechnet. Es muss N1=2, N2=0 sein.
const bool krafttestmode = false;

// softening WCA interaction: no two-particle-force can be larger than this
const double maxPairWCA = 1.0e5;

//have every run write its progress at least after this many seconds
const int max_write_interval = 900;

//mitteln über wie viele runs?
const int runs = 1;

//wie viele Jobs soll es geben?
const int jobs = 1;

//wie viele Threads sollen gestartet werden, dh wie viele CPUs verwendet? 0 für unbegrenzt
const int maxThreads = 1;

//Observable aufnehmen in welchem Zeitabstand?
const double obs_dt = 0.2; 

// Obervable aufnehmen wie oft?
const int obs_anzahl = 12; 


// Binbreite Paarkorrelationsfunktion
const double korr_dr = 0.2;

// größtes r der Paarkorrelationsfunktion
const double korr_rmax = L*0.5;

// in Histogram of local densities: number of bins and width of a bin
const int rhohist_bins = 100;
const double rhohist_drho = 0.01;
// bin assignment in rhohist:  x = rhohist_drho*(bin+0.5)
// bin assignment in drhohist: x = rhohist_drho*(bin+0.5-rhohist_bins/2)
// solved for bin: in rhohist bin=(int)(x/rhohist_drho - 0.5)
// and in drhohist:           bin=(int)(x/rhohist_drho - 0.5 + rhohist_bin/2)

// thresholds for rhohist. When a cell has more than X density, it is counted.
const double rho_thresh1 = 0.1;
const double rho_thresh2 = 0.2;
const double rho_thresh3 = 0.5;
// same for drhohist. When in a cell |rho1-rho2|/(rho1+rho2) > X, it is counted.
const double drho_thresh1= 0.1;
const double drho_thresh2= 0.2;
const double drho_thresh3= 0.5;
const double drho_thresh4= 0.99;


// Binbreite Dichteprofil
const double rhoring_dr = 0.5;

// Anzahl Bins im Dichteprofil
const int rhoring_bins = 500;


//größtes q in rhothilde, in Einheiten dq=4pi/L
const int ftrho_qbins = 30;


// Start: Zufall=1, aus Datei=2, Gitterstart=3, allegleich=4, Kreisscheibe=5, Kern+Ring=6, Krafttest=7
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

// snapshot-grid: Crude value of lattice constant. The actual lattice constant will be chosen 
// close to this value such that L is a multiple of this. In units of sigma11
const double snapgrid_crude = 0.5;

// for rhok test: Width of the gaussian distribution that write_gaussdens should produce
const double gaussdens_width_x = 0.1*L;
const double gaussdens_width_y = 0.02*L;

/////////////// AB HIER: AENDERN VERBOTEN! /////////////////////////


// die Konstante 2^(1/6)
const double zweihoch1_6 = pow(2.0, 1.0/6.0);

// snapshot-grid: number of lattice cells (per dimension)
const int snapgrid_num = 2 * (int) (L/(2*snapgrid_crude) + 0.5);
//snapshot-grid: size of each lattice cell
const double snapgrid_width = L / snapgrid_num;

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
const int korr_bins = (int) (korr_rmax/korr_dr)+1;

//Breite dq für rhothilde
const double ftrho_dq = 4.0*M_PI/L;
