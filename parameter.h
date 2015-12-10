// enthaelt alle veraenderbaren Parameter

#pragma once

#define _USE_MATH_DEFINES
#include <math.h>
#include <string>

using std::string;

/// Parameter, die geaendert werden duerfen

// Teilchenzahl
const int N = 8000;

// Groesse der Simulationsbox (je Raumrichtung) in Einheiten von sigma
const double L = 600.0;

// 0 fuer Nearest Grid Point, 1 fuer Cloud In Cell, 2 fuer Triangular Shaped Cloud
const int densGrid_Schema = 2;

// Anzahl Samples in der Fouriertransformation. Nur gerade Zahlen erlaubt.
const int FFT_samples = 512;


// kapillarlaenge in Einheiten von sigma
const double lambda_kapillar = 15.0;

// Vorfaktor der Kapillarkraft. Ist f^2/(eps gamma).
const double kapillar_vorfaktor = 100.0;

// Temperatur, in Einheiten kT/eps
const double T = 0.1;


// maximaler Zeitschritt
const double dt_max = 0.05;

//maximale Reisedistanz in einem Zeitschritt
const double max_reisedistanz = 0.1;



//mitteln über wie viele runs?
const int runs = 1;

//Observable aufnehmen in welchem Zeitabstand?
const double obs_dt = 0.1;

// Obervable aufnehmen wie oft?
const int obs_anzahl = 1;


// Binbreite Paarkorrelationsfunktion
const double korr_dr = 0.05;

// größtes r der Paarkorrelationsfunktion
const double korr_rmax = 5.0;


//größtes q in rhothilde, in Einheiten dq=4pi/L
const int ftrho_qbins = 150;


// Start: Zufall=1, aus Datei=2, Gitterstart=3, allegleich=4, Kreisscheibe=5
const int startpos_methode=5;
//falls aus Datei: Name der Datei
const string startpos_dateiname="startpos.txt";

//falls Positionen in Datei geschrieben werden: Name der Datei
const string pos_output_dateiname="pos.txt";



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

// Abtastrate bzw Bandbreite im reziproken Raum. q's gehen (komponentenweise) von -Omega bis +Omega
const double Omega = M_PI/densGrid_Breite;

// Anzahl Bins Paarkorrelationsfunktion
const int korr_bins = (int) (korr_rmax/korr_dr);

//Breite dq für rhothilde
const double ftrho_dq = 4.0*M_PI/L;