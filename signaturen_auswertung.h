#pragma once

#include <iostream>
#include <string>
#include <vector>

using std::vector;

//lese Teilchenpositionen (mit zugehöriger Zeit) aus Dateien
bool rvont_lesen();

//berechnet Dichteprofile, schreibt je Run in eine Datei
void berechne_dichteprofile();
void berechne_korrelationsfunktionen();


//mittelt Dichteprofile über Jobs, schreibt Mittelwerte und Fehler in Datei
void auswerten_dichteprofile();
void auswerten_korrelationsfunktionen();


//Zahl zu String, ohne führende Nullen
std::string int_to_string(int);

//berechnet Mittelwerte und Varianzen von input[][][]. Mittelung über den ersten Index, der bis runs geht. Der zweite geht bis anzahl. mittelwerte und varianzen müssen schon (leer) initialisiert sein.
void statistik_1(vector<double>**& input, vector<double>*& mittelwerte, vector<double>*& varianzen, int anzahl);


//abstandsquadrat von Teilchen i (Typ 1) und Teilchen j (Typ 1). Berücksichtigt periodische Randbedingungen. Analog andere Typen
double abstand2(double*, double*);


// gleichverteilte Zufallszahl
double zufall_gleichverteilt_vonbis(double min, double max);