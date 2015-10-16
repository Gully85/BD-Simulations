// enthaelt alle veraenderbaren Parameter

#pragma once

#define _USE_MATH_DEFINES
#include <math.h>

/// Parameter, die geaendert werden duerfen

// Teilchenzahl
const int N = 1;

// Groesse der Simulationsbox (je Raumrichtung) in Einheiten von sigma
const double L = 30.0;

// 0 fuer Nearest Grid Point, 1 fuer Cloud In Cell, 2 fuer Triangular Shaped Cloud
const int densGrid_Schema = 1; //vorher "density_gridding_schema"

// Anzahl Samples in der Fouriertransformation. Nur gerade Zahlen erlaubt.
const int FFT_samples = 256;




// kapillarlaenge in Einheiten von sigma
const double lambda_kapillar = 40.0;


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


