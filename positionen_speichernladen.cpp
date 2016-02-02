// stellt Methoden bereit, die zu Beginn der Simulation Teilchen auf die richtigen Positionen setzen.

#pragma once

#include "signaturen.h"
#include "zufall.cpp"


extern int** r_git;
extern double** r_rel;

extern const int N;
extern const double L;
extern const double nachList_Breite;

extern const string pos_output_dateiname;
extern const string startpos_dateiname;
extern const double startpos_kreisradius;




// schreibt aktuelle Teilchenpositionen in Datei
void pos_schreiben(){
	FILE* out = fopen(pos_output_dateiname.c_str(), "w");
	
	fprintf(out, "# Format: teilchenNr TAB x TAB y \n\n");
	
	for(int i=0; i<N; i++){
		double x=r_git[i][0]*nachList_Breite + r_rel[i][0];
		double y=r_git[i][1]*nachList_Breite + r_rel[i][1];
		
		fprintf(out, "%d \t %g \t %g \n", i, x, y);
	}//for i
	
	
	
	
	
	fclose(out);
}//void pos_schreiben

//liest Startpositionen aus Datei
void init_pos_aus_datei(){
	const string startpos_dateiname="startpos.txt";
	
	FILE* in = fopen(startpos_dateiname.c_str(), "r");
	
	//erste zwei Zeilen überspringen
	char buffer[100];
	//zwei Zeilen lesen. Wenn beim Lesen von mindestens einer ein Fehler kommt, abbrechen
	if(NULL == fgets(buffer, 100, in) || NULL == fgets(buffer, 100, in)){
		cout << "Fehler beim Lesen der Inputdatei, Zeile 1 oder 2" << endl;
		return;
	}//if Fehler beim Lesen
	
	// lese zeilenweise
	int teilchenNr=0;
	double x,y;
	float xf, yf;
	bool fehler=false;
	for(int zeile=0; zeile<N; zeile++){
		if(3 != fscanf(in, "%d \t %g \t %g", &teilchenNr, &xf, &yf)) fehler=true;
		if(teilchenNr != zeile) fehler=true;
		
		if(fehler){
			cout << "Fehler beim Lesen der Inputdatei, Zeile " << zeile << endl;
			return;
		}//if fehler
		
		x = (double) xf;
		y = (double) yf;
		r_git[teilchenNr][0] = (int) (x/nachList_Breite);
		r_git[teilchenNr][1] = (int) (y/nachList_Breite);
		r_rel[teilchenNr][0] = x - r_rel[teilchenNr][0]*nachList_Breite;
		r_rel[teilchenNr][1] = y - r_rel[teilchenNr][1]*nachList_Breite;
		
	}//for int zeile
	
	fclose(in);
	
}//void init_pos_aus_datei


// initialisiert Teilchenpositionen auf Zufallspositionen. Schreibt sie in r_git und r_rel
void init_zufallspos(){

	
//bestimme zufaellige Positionen in [0,L], setze die Teilchen dort hin
	for(int i=0; i<N; i++){
		double x = zufall_gleichverteilt_vonbis(0.0, L);
		double y = zufall_gleichverteilt_vonbis(0.0, L);

		int jc = (int) x/nachList_Breite;
		int kc = (int) y/nachList_Breite;

		r_git[i][0] = jc;
		r_git[i][1] = kc;
		r_rel[i][0] = x - jc*nachList_Breite;
		r_rel[i][1] = y - kc*nachList_Breite;
		
	}//for i
	
} //void init_zufallspos


// bestimme EINE Position in [0,L], setze alle Teilchen dorthin
void init_allegleich(){
	double x = zufall_gleichverteilt_vonbis(0.0,L);
	double y = zufall_gleichverteilt_vonbis(0.0,L);
	
	x = 2.0;
	y = 0.0;
	
	int jc = (int) x/nachList_Breite;
	int kc = (int) y/nachList_Breite;
	for(int i=0; i<N; i++){
		r_git[i][0] = jc;
		r_git[i][1] = kc;
		r_rel[i][0] = x - jc*nachList_Breite;
		r_rel[i][1] = y - kc*nachList_Breite;
	}//for i
}//void init_allegleich

// initialisiert alle Teilchenpositionen gleichverteilt in einem Kreis in der Mitte der Box
void init_kreisscheibe(){
	const double rad = startpos_kreisradius;
	const double rad2 = rad*rad;
// 	double x=2.0;
// 	double y=2.0;
	for(int i=0; i<N; i++){
		double x,y;
		do{
			x = zufall_gleichverteilt_vonbis(-rad,rad);
			y = zufall_gleichverteilt_vonbis(-rad,rad);
		} while (x*x + y*y > rad2);
		x+= 0.5*L;
		y+= 0.5*L;

		
		int jc = (int) (x/nachList_Breite);
		int kc = (int) (y/nachList_Breite);
		
		r_git[i][0] = jc;
		r_git[i][1] = kc;
		r_rel[i][0] = x - jc*nachList_Breite + 0.03;
		r_rel[i][1] = y - kc*nachList_Breite + 0.03;
		
	}//for i
	
}//void init_kreisscheibe


// initialisiert Teilchenpositionen gleichverteilt in der Box. Baut ein quadratisches Gitter auf. Wenn N keine Quadratzahl ist, bleiben dabei Gitterplätze frei.
void init_gitterstart(){
	// es sollen N Teilchen auf ein 2dim-Gitter verteilt werden. In jeder Richtung also sqrt(N) viele Teilchen auf die Länge L. Falls N keine Quadratzahl ist, sogar ceil(sqrt(N)).
	double wurzelN = ceil(sqrt(N));
	double schrittweite = L/wurzelN;
	int index =0;
	for(int i=0; i<wurzelN; i++)
	for(int j=0; j<wurzelN; j++){
		double x = i*schrittweite;
		double y = j*schrittweite;
		
		int ic = (int) (x/nachList_Breite);
		int jc = (int) (y/nachList_Breite);
		
		r_git[index][0] = ic;
		r_git[index][1] = jc;
		r_rel[index][0] = x - ic*nachList_Breite;
		r_rel[index][1] = y - jc*nachList_Breite;
		
		if(++index == N) return;
		
	}//for i,j
	
}//void init_gitterstart