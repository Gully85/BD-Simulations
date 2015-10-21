#include "zufall.cpp"
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <math.h>
#include "gridRoutinen.cpp"
#include "signaturen.h"
#include <fftw3.h>
#include "dynamik_methoden.cpp"

using std::cout; using std::endl; using std::flush;


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
double** r = NULL;

// Kapillarkraefte, die gerade auf Teilchen wirken. Erster Index TeilchenNr, zweiter Index Raumrichtung
double** Fkap = NULL;


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


///setze Teilchen auf Positionen

//* Test: Zunaechst nur zwei Teilchen, eins an Position x=L/2 + 0.5*densGrid_Breite, y=x. Also genau in der Mitte einer densGrid-Zelle. Das andere bei x=L/2 + 2.5*densGrid_Breite, y wie vorher. Damit sind die beiden Teilchen genau  zwei Zellen voneinander entfernt.
r[0][0] = 0.5*L + 0.5*densGrid_Breite;
r[0][1] = 0.5*L + 0.5*densGrid_Breite;

r[1][0] = 0.5*L + 2.5*densGrid_Breite;
r[1][1] = r[0][1];
// */


// bereite Kraftberechnung vor, dh plane Fouriertrafos und reserviere Speicher
kapkraefte_init();

//berechne Kapillarkraefte schreibe sie in Fkap
berechne_kapkraefte(r, Fkap);

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
