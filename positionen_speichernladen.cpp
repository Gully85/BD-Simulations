// stellt Methoden bereit, die zu Beginn der Simulation Teilchen auf die richtigen Positionen setzen.

#pragma once

#include "signaturen.h"
#include <string>


extern int** r1_git;
extern int** r2_git;
extern double** r1_rel;
extern double** r2_rel;

extern const int N1, N2;
extern const double L;
extern const double nachList_Breite;

extern const int jobs;
extern double**** r1_abs;
extern double**** r2_abs;

extern const string startpos_dateiname;
extern const double startpos_kreisradius;




// schreibt aktuelle Teilchenpositionen in Datei
void pos_schreiben_einedatei(){
	FILE* out = fopen(pos_output_dateiname.c_str(), "w");
	
	fprintf(out, "# Es sind %d Typ1- und %d Typ2-Teilchen vorhanden. \n", N1, N2);
	fprintf(out, "# Format: teilchenNr TAB Typ TAB x TAB y \n\n");
	
	for(int i=0; i<N1; i++){
		double x = r1_git[i][0]*nachList_Breite + r1_rel[i][0];
		double y = r1_git[i][1]*nachList_Breite + r1_rel[i][1];
		
		fprintf(out, "%d \t %d \t %g \t %g \n", i, 1, x, y);
	}//for i
	fprintf(out, "\n");
	for(int i=0; i<N2; i++){
		double x = r2_git[i][0]*nachList_Breite + r2_rel[i][0];
		double y = r2_git[i][1]*nachList_Breite + r2_rel[i][1];
		
		fprintf(out, "%d \t %d \t %g \t %g \n", i, 2, x, y);
	}//for i bis N2
	fclose(out);
}//void pos_schreiben

//schreibt aktuelle Teilchenpositionen in zwei Dateien. FILE-Pointer müssen übergeben werden.
void pos_schreiben(double t, FILE* datei1, FILE* datei2){
	
	for(int i=0; i<N1; i++){
		double x = r1_git[i][0]*nachList_Breite + r1_rel[i][0];
		double y = r1_git[i][1]*nachList_Breite + r1_rel[i][1];
		
		//Format: t TAB x TAB y
		fprintf(datei1, "%9.5f \t %g \t %g \n", t, x, y);
	}//for i
	fprintf(datei1, "\n");
	
	if(0 == N2) return;
	
	for(int i=0; i<N2; i++){
		double x = r2_git[i][0]*nachList_Breite + r2_rel[i][0];
		double y = r2_git[i][1]*nachList_Breite + r2_rel[i][1];
		
		//Format: t TAB x TAB y
		fprintf(datei2, "%9.5f \t %g \t %g \n", t, x, y);
	}//for i bis N2
	fprintf(datei2, "\n");
	
}//void pos_schreiben



//liest Startpositionen aus Datei. Es muss zuerst Typ 1 und dann Typ 2 stehen.
void init_pos_aus_datei(){
	const string startpos_dateiname="startpos.txt";
	
	FILE* in = fopen(startpos_dateiname.c_str(), "r");
	
	//erste zwei Zeilen überspringen
	char buffer[100];
	//zwei Zeilen lesen. Wenn beim Lesen von mindestens einer ein Fehler kommt, abbrechen
	if(NULL == fgets(buffer, 100, in) || NULL == fgets(buffer, 100, in)){
		cout << "Fehler beim Lesen der Inputdatei, Zeile 1 oder 2" << endl << flush;
		return;
	}//if Fehler beim Lesen
	
	// lese zeilenweise
	int teilchenNr=0;
	int typ;
	double x,y;
	float xf, yf;
	bool fehler=false;
	for(int zeile=0; zeile<N1; zeile++){
		if(4 != fscanf(in, "%d \t %d \t %g \t %g", &teilchenNr, &typ, &xf, &yf)) fehler=true;
		if(teilchenNr != zeile || typ != 1) fehler=true;
		
		if(fehler){
			cout << "Fehler beim Lesen der Inputdatei, Zeile " << zeile << endl << flush;
			return;
		}//if fehler
		
		x = (double) xf;
		y = (double) yf;
		if(x<0.0 || x>L || y<0.0 || y>L){
			cout << "Fehler beim Laden: Koordinaten von Teilchen "<<teilchenNr<<" (Typ 1), x="<<x<<", y="<<y<<" out of range!" << endl << flush;
			return;
		}//if x oder y out of range
		r1_git[teilchenNr][0] = (int) (x/nachList_Breite);
		r1_git[teilchenNr][1] = (int) (y/nachList_Breite);
		r1_rel[teilchenNr][0] = x - r1_rel[teilchenNr][0]*nachList_Breite;
		r1_rel[teilchenNr][1] = y - r1_rel[teilchenNr][1]*nachList_Breite;
		
	}//for int zeile
	for(int zeile=0; zeile<N2; zeile++){
		if(4 != fscanf(in, "%d \t %d \t %g \t %g", &teilchenNr, &typ, &xf, &yf)) fehler=true;
		if(teilchenNr != zeile || typ != 2) fehler=true; 
		
		if(fehler){
			cout << "Fehler beim Lesen der Inputdatei, Zeile " << zeile+N1<< ". Möglicherweise sind nicht genau " <<N1<<"+"<<N2<<"Teilchen in startpos.txt?" << endl<<flush;
			return;
		}//if fehler
		
		x = (double) xf;
		y = (double) yf;
		if(x<0.0 || x>L || y<0.0 || y>L){
			cout << "Fehler beim Laden: Koordinaten von Teilchen "<<teilchenNr<<" (Typ 2), x="<<x<<", y="<<y<<" out of range!" << endl << flush;
			return;
		}//if x oder y out of range
		
		r2_git[teilchenNr][0] = (int) (x/nachList_Breite);
		r2_git[teilchenNr][1] = (int) (y/nachList_Breite);
		r2_rel[teilchenNr][0] = x - r2_rel[teilchenNr][0]*nachList_Breite;
		r2_rel[teilchenNr][1] = y - r2_rel[teilchenNr][1]*nachList_Breite;

	}//for zeile
	
	fclose(in);
	
}//void init_pos_aus_datei


// initialisiert Teilchenpositionen auf Zufallspositionen. Schreibt sie in r_git und r_rel
void init_zufallspos(){

	
//bestimme zufaellige Positionen in [0,L], setze die Teilchen dort hin
	for(int i=0; i<N1; i++){
		double x = zufall_gleichverteilt_vonbis(0.0, L);
		double y = zufall_gleichverteilt_vonbis(0.0, L);

		int jc = (int) x/nachList_Breite;
		int kc = (int) y/nachList_Breite;

		r1_git[i][0] = jc;
		r1_git[i][1] = kc;
		r1_rel[i][0] = x - jc*nachList_Breite;
		r1_rel[i][1] = y - kc*nachList_Breite;
		
	}//for i
	
	for(int i=0; i<N2; i++){
		double x = zufall_gleichverteilt_vonbis(0.0, L);
		double y = zufall_gleichverteilt_vonbis(0.0, L);
		
		int jc = (int) x/nachList_Breite;
		int kc = (int) y/nachList_Breite;
		
		r2_git[i][0] = jc;
		r2_git[i][1] = kc;
		r2_rel[i][0] = x - jc*nachList_Breite;
		r2_rel[i][1] = y - kc*nachList_Breite;
	}//for i
	
} //void init_zufallspos


// bestimme EINE Position in [0,L], setze alle Teilchen dorthin
void init_allegleich(){
	double x = zufall_gleichverteilt_vonbis(0.0,L);
	double y = zufall_gleichverteilt_vonbis(0.0,L);
	
	// x = 2.0;
	// y = 0.0;
	
	int jc = (int) (x/nachList_Breite);
	int kc = (int) (y/nachList_Breite);
	for(int i=0; i<N1; i++){
		r1_git[i][0] = jc;
		r1_git[i][1] = kc;
		r1_rel[i][0] = x - jc*nachList_Breite;
		r1_rel[i][1] = y - kc*nachList_Breite;
	}//for i
	
	
	//Typ 2 analog
	x = zufall_gleichverteilt_vonbis(0.0, L);
	y = zufall_gleichverteilt_vonbis(0.0, L);
	jc = (int) (x/nachList_Breite);
	kc = (int) (y/nachList_Breite);
	
	for(int i=0; i<N2; i++){
		r2_git[i][0] = jc;
		r2_git[i][1] = kc;
		r2_rel[i][0] = x - jc*nachList_Breite;
		r2_rel[i][1] = y - kc*nachList_Breite;
	}//for i
}//void init_allegleich

// initialisiert alle Teilchenpositionen zufällig gleichverteilt in einem Kreis in der Mitte der Box
void init_kreisscheibe(){
	const double rad = startpos_kreisradius;
	const double rad2 = rad*rad;
// 	double x=2.0;
// 	double y=2.0;
	for(int i=0; i<N1; i++){
		double x,y;
		do{
			x = zufall_gleichverteilt_vonbis(-rad,rad);
			y = zufall_gleichverteilt_vonbis(-rad,rad);
		} while (x*x + y*y > rad2);
		x+= 0.5*L;
		y+= 0.5*L;

		
		int jc = (int) (x/nachList_Breite);
		int kc = (int) (y/nachList_Breite);
		
		r1_git[i][0] = jc;
		r1_git[i][1] = kc;
		r1_rel[i][0] = x - jc*nachList_Breite + 0.03;
		r1_rel[i][1] = y - kc*nachList_Breite + 0.03;
		
	}//for i
	
	for(int i=0; i<N2; i++){
		double x,y;
		do{
			x = zufall_gleichverteilt_vonbis(-rad, rad);
			y = zufall_gleichverteilt_vonbis(-rad, rad);
		} while(x*x + y*y > rad2);
		x+= 0.5*L;
		y+= 0.5*L;
		
		int jc = (int)(x/nachList_Breite);
		int kc = (int)(y/nachList_Breite);
		
		r2_git[i][0] = jc;
		r2_git[i][1] = kc;
		r2_rel[i][0] = x - jc*nachList_Breite + 0.03;
		r2_rel[i][1] = y - kc*nachList_Breite + 0.03;
		
	}//for i
	
}//void init_kreisscheibe


// initialisiert Teilchenpositionen gleichverteilt in der Box. Baut ein quadratisches Gitter auf. Wenn N keine Quadratzahl ist, bleiben dabei Gitterplätze frei.
void init_gitterstart(){
	
	if(0 != N2){
		cout << "Fehler: Gitterstart mit zwei Teilchentypen noch nicht implementiert!" << endl << flush;
		return;
	}//if N2 ungleich Null
	// es sollen N Teilchen auf ein 2dim-Gitter verteilt werden. In jeder Richtung also sqrt(N) viele Teilchen auf die Länge L. Falls N keine Quadratzahl ist, sogar ceil(sqrt(N)).
	double wurzelN = ceil(sqrt(N1));
	double schrittweite = L/wurzelN;
	int index =0;
	for(int i=0; i<wurzelN; i++)
	for(int j=0; j<wurzelN; j++){
		double x = i*schrittweite;
		double y = j*schrittweite;
		
		int ic = (int) (x/nachList_Breite);
		int jc = (int) (y/nachList_Breite);
		
		r1_git[index][0] = ic;
		r1_git[index][1] = jc;
		r1_rel[index][0] = x - ic*nachList_Breite;
		r1_rel[index][1] = y - jc*nachList_Breite;
		
		if(++index == N1) return;
		
	}//for i,j
	
}//void init_gitterstart