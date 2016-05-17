#pragma once

#include <iostream>
#include <string>
#include <vector>

using std::vector;

//lese Teilchenpositionen (mit zugehöriger Zeit) aus Dateien
bool rvont_lesen();

//berechnet Dichteprofile
void berechne_dichteprofile();

//mittelt Dichteprofile über Jobs, schreibt Mittelwerte und Fehler in Datei
void auswerten_dichteprofile();


//Zahl zu String, ohne führende Nullen
std::string int_to_string(int);

//berechnet Mittelwerte und Varianzen von input[][][]. Mittelung über den ersten Index, der bis runs geht. Der zweite geht bis anzahl. mittelwerte und varianzen müssen schon (leer) initialisiert sein.
void statistik_1(vector<double>**& input, vector<double>*& mittelwerte, vector<double>*& varianzen, int anzahl);