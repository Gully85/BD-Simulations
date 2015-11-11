#include "zufall.cpp"
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <math.h>
#include "gridRoutinen.cpp"
#include "signaturen.h"
#include <fftw3.h>
#include "dynamik_methoden.cpp"
#include <algorithm>


using std::cout; using std::endl; using std::flush; using std::vector; using std::min;


// importiere Variablen aus parameter.h, siehe dort, was die Variable tut.
extern const int N;
extern const double L;

extern const int densGrid_Schema, densGrid_Zellen;
extern const double densGrid_Breite;

extern const int nachList_Zellen;
extern const double nachList_Breite;

extern const double lambda_kapillar;



/// Globale Felder. Mit Null initialisiert.

// Teilchenpositionen, x- und y-Koordinate. r[i][0] ist die x-Koordinate der Position von Teilchen i.
double** r = NULL; //wird bald verschwinden, weil durch r_git und r_rel ausdrueckbar

// Teilchenpositionen, Gittervektor. Erster Index TeilchenNr, zweiter Index 0 fuer x und 1 fuer y-Richtung
int** r_git = NULL;
// dasgleiche, Vektor innerhalb der Zelle
double** r_rel = NULL;


// Kapillarkraefte, die gerade auf Teilchen wirken. Erster Index TeilchenNr, zweiter Index Raumrichtung
double** Fkap = NULL;

//WCA-Kraefte, die gerade auf Teilchen wirken. Erster Index TeilchenNr, zweiter Index Raumrichtung
double** F_WCA = NULL;

//Zufallskraefte, gleiche Indices
double** F_noise = NULL;


int main(){

srand(time(NULL));

int i,j,k;
FILE* out;

///reserviere Felder

// Teilchenpositionen, x- und y-Koordinate. r[i][0] ist die 0-Koordinate (x) der Position von Teilchen i.
r = new double*[N];
for(i=0; i<N; i++)
	r[i] = new double[2];

// Kapillarkraefte, die gerade auf Teilchen wirken. Erster Index TeilchenNr, zweiter Index Raumrichtung
Fkap = new double*[N];
for(i=0; i<N; i++)
	Fkap[i] = new double[2];

//WCA-Kraefte, die gerade auf Teilchen wirken. Erster Index TeilchenNr, zweiter Index Raumrichtung
F_WCA = new double*[N];
for(i=0; i<N; i++)
	F_WCA[i] = new double[2];

//Zufallskraefte, gleiche Indices
F_noise = new double*[N];
for(i=0; i<N; i++)
	F_noise[i] = new double[2];


///setze Teilchen auf Positionen

//* Test: Zunaechst nur ein Teilchen an Position x=L/2 + 0.5*densGrid_Breite, y=x. Also genau in der Mitte einer densGrid-Zelle. 
r[0][0] = 0.5*densGrid_Breite + 0.5*L;
r[0][1] = 0.5*L + 0.5*densGrid_Breite;

//r[1][0] = 0.5*L + 2.5*densGrid_Breite;
//r[1][1] = r[0][1];
// */


// bereite Kraftberechnung vor, dh plane Fouriertrafos und reserviere Speicher
kapkraefte_init();

//* Test: Platziere zweites Teilchen im Abstand t (y-Richtung), messe Kraft Fy. Variiere t, trage ueber t auf.
r[1][0] = r[0][0];
FILE* outF = fopen("Fkap.txt", "w");
for(double t=0.1; t<2.0; t+= 0.05){

	r[1][1] = r[0][1] + t;
	//berechne Kapillarkraefte schreibe sie in Fkap
	berechne_kapkraefte(r_git, r_rel, Fkap);
	
	//berechnet Kapillarkraefte (mittels Fouriertransformation), schreibt sie in Fkap
void berechne_kapkraefte(int** rr_git, double** rr_rel, double** Fkap); 

	fprintf(outF, "%g \t %g \n", t, Fkap[0][1]);
}//for t

//* Test: Schreibe Kapillarkraefte ins Terminal
for(i=0; i<N; i++){
	cout << "Kapillarkraft auf Teilchen " << i << ": (" << Fkap[i][0] << ","<< Fkap[i][1] << ")" << endl;
}//for i

return 0;

}//int main




//index-Wrapping, ermoeglicht es 1dim-Arrays mit zwei Indices anzusprechen. FFTW erwartet Arrays mit nur einem Index. NOCH keine period. Randbed, aber die kann man hier einbauen
int iw(int i, int j){
	if (i>=densGrid_Zellen || j>= densGrid_Zellen || i<0 || j<0){
		cout << "Index out of Bounds! Aufruf war iw(i=" <<i<<", j=" <<j<<") \n" << endl << flush;
		return -1; //erzeugt absichtlich Segmentation Fault
	}//if index out-of-bounds

	return i*densGrid_Zellen + j;

}//int index

//berechnet das Quadrat des Abstands zwischen Teilchen i und Teilchen j (gleichen Typs). Berücksichtigt periodische Randbedingungen.
double abstand2(int i, int j, int** rr_git, double** rr_rel){
	
	//Absolutpositionen
	double xi_abs = rr_git[i][0]*nachList_Breite + rr_rel[i][0];
	double xj_abs = rr_git[j][0]*nachList_Breite + rr_rel[j][0];
	double yi_abs = rr_git[i][1]*nachList_Breite + rr_rel[i][1];
	double yj_abs = rr_git[j][1]*nachList_Breite + rr_rel[j][1];
	
	// dx^2, durch die periodischen Randbed gibt es drei Moeglichkeiten
	double tmp1 = xi_abs - xj_abs;
	double tmp2 = xi_abs - xj_abs + L;
	double tmp3 = xi_abs - xj_abs - L;
	double dx2 = min(tmp1*tmp1, min(tmp2*tmp2, tmp3*tmp3));
	
	// das gleiche fuer dy^2
	tmp1 = yi_abs - yj_abs;
	tmp2 = yi_abs - yj_abs + L;
	tmp3 = yi_abs - yj_abs - L;
	double dy2 = min(tmp1*tmp1, min(tmp2*tmp2, tmp3*tmp3));
	
	return dx2 + dy2;
}//double abstand2




