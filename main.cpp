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
//double** r = NULL; //wird bald verschwinden, weil durch r_git und r_rel ausdrueckbar

// Teilchenpositionen, Gittervektor. Erster Index TeilchenNr, zweiter Index 0 fuer x und 1 fuer y-Richtung
int** r_git = NULL;
// dasgleiche, Vektor innerhalb der Zelle
double** r_rel = NULL;



int main(){

srand(time(NULL));

int i,j,k;


/// initialisiere Felder, Positionen, Kraftberechnungen
main_init();


/* Test: Platziere zweites Teilchen im Abstand t (y-Richtung), messe Kraft Fy. Variiere t, trage ueber t auf.
r_git[1][0] = r_git[0][0]; //x-Komponente wie Teilchen 0
r_rel[1][0] = r_rel[0][0];

FILE* outF = fopen("Fvonr.txt", "w");
fprintf(outF, "# Format: r TAB Fkap TAB F_WCA TAB summe \n\n");
for(double t=0.1; t<2.0; t+= 0.01){

	x = 0.5*densGrid_Breite + 0.5*L + t; //eigentlich y
	r_git[1][1] = (int) (x/nachList_Breite);
	r_rel[1][1] = x - nachList_Breite*r_git[1][1];
	//berechne Kapillarkraefte, schreibe sie in Fkap
	berechne_kapkraefte(r_git, r_rel, Fkap);
	//berechne WCA-Kraefte, schreibe sie in F_WCA
	berechne_WCAkraefte(F_WCA);
	
	fprintf(outF, "%g \t %g \t %g \t %g\n", t, Fkap[0][1], F_WCA[0][1], Fkap[0][1]+F_WCA[0][1]);
}//for t
//*/



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


