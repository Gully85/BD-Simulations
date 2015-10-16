#include "zufall.cpp"
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <math.h>
#include "gridRoutinen.cpp"
#include "signaturen.h"
#include <fftw3.h>


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
// Dichte auf Gitter im x-Raum. Muss den FFTW-Typ haben, weil dieses Feld fouriertransformiert wird.
fftw_complex* rhox = NULL;
// Dichte auf Gitter im k-Raum. Muss den FFTW-Typ haben, weil dieses Feld fouriertransformiert wird.
fftw_complex* rhok = NULL;

// Greensfunktion, gleiche Indizes wie rhok nach der Fouriertrafo hat. G ist reell.
double* G = NULL;

int main(){

srand(time(NULL));

int i,j,k;
FILE* out;

///reserviere Felder

// Teilchenpositionen, x- und y-Koordinate. r[i][0] ist die 0-Koordinate (x) der Position von Teilchen i.
r = new double*[N];
for(i=0; i<N; i++)
	r[i] = new double[2];

// Dichte auf Gitter im x-Raum. Muss den FFTW-Typ haben, weil dieses Feld fouriertransformiert wird.
rhox = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * densGrid_Zellen*densGrid_Zellen);
// Dichte auf Gitter im k-Raum. Muss den FFTW-Typ haben, weil dieses Feld fouriertransformiert wird.
rhok = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * densGrid_Zellen*densGrid_Zellen);

// es muss Index-Wrapping gemacht werden, mit der Funktion int iw(int i, int j)
// rhox[iw(a,b)][0] ist die Anzahl der Teilchen in der Zelle, deren Mittelpunkt bei x=(a+0.5)*densGrid_Breite, y=(b+0.5)*densGrid_Breite liegt. densGrid_Breite ist also das, was im fftw-Test dx war.
// rhox[.......][1] ist der Imaginaerteil, also Null.

// rhok[iw(c,d)][0] ist die Fouriertransformierte von rho, an der Stelle qx=dq*((c+(int)(0.5*densGrid_Zellen))%densGrid_Zellen - 0.5*densGrid_Zellen), qy entsprechend.
// rhok[.......][1] der Imaginaerteil davon, muss nicht Null sein.


//Greensfunktion, nur Realteil noetig
G = new double[densGrid_Zellen*densGrid_Zellen];



///plane FFTs
fftw_plan ft_voran = fftw_plan_dft_2d(densGrid_Zellen, densGrid_Zellen, rhox, rhok, FFTW_FORWARD, FFTW_MEASURE);
fftw_plan ft_rueck = fftw_plan_dft_2d(densGrid_Zellen, densGrid_Zellen, rhok, rhox, FFTW_BACKWARD,FFTW_MEASURE);


///setze Teilchen auf Positionen

//* Test: Zunaechst nur ein Teilchen, an Position x=L/2 + 0.5*densGrid_Breite, y=x. Also genau in der Mitte einer densGrid-Zelle
r[0][0] = 0.5*(L + densGrid_Breite);
r[0][1] = 0.5*(L + densGrid_Breite);
// */




switch(densGrid_Schema){
	case 0: gridDensity_NGP(rhox, r); break;
	case 1: gridDensity_CIC(rhox, r); break;
	case 2: gridDensity_TSC(rhox, r); break;
	default: cout << "Density-Gridding-Schema (NearestGridPoint/CloudInCell/TSC) nicht erkannt! densGrid_Schema="<<densGrid_Schema<<", zulaessig sind nur 0,1,2." << endl; 
		 return 1;
}//switch


//* Test: Gebe Feld rhox aus
out = fopen("test.txt", "w");
fprintf(out, "#Format: x TAB y TAB rho");
for(j=0; j<densGrid_Zellen; j++){
	for(k=0; k<densGrid_Zellen; k++){
		fprintf(out, "%g \t %g \t %g \n", (j+0.5)*densGrid_Breite, (k+0.5)*densGrid_Breite, rhox[iw(j,k)][0]);
	}//for k
	fprintf(out, "\n");
}//for j
// */


//transformiere rho in den k-Raum. Danach steht die Dichte in rhok[.][0] und rhok[.][1], real- und Imaginaerteil, mit alternierenden Vorzeichen. Bedenke beim Index-Wrapping die Verschiebung um N/2
fftw_execute(ft_voran);

//erzeuge Greensfunktion G(qx, qy) = -1/(sin2(qx dx)/dx2 + sin2(qy dy)/dy2 + 1/lambda2)
//G ist reell. Erzeuge Feld G[] mit demselben Index-Wrapping und Verschiebung
double qx, qy;
double s2x, s2y;
for(j=0; j<densGrid_Zellen; j++){
	qx = dq*((j+(int)(0.5*densGrid_Zellen))%densGrid_Zellen - 0.5*densGrid_Zellen);
	s2x = sin(0.5*qx * densGrid_Breite)*sin(0.5*qx * densGrid_Breite)/(densGrid_Breite*densGrid_Breite);
	for(k=0; k<densGrid_Zellen; k++){
		qy = dq*((k+(int)(0.5*densGrid_Zellen))%densGrid_Zellen - 0.5*densGrid_Zellen);
		s2y = sin(0.5*qy * densGrid_Breite)*sin(0.5*qy * densGrid_Breite)/(densGrid_Breite*densGrid_Breite);
		
		G[iw(j,k)] = -1.0/(s2x + s2y + 1.0/(lambda_kapillar*lambda_kapillar));
	}//for k
}//for j


//multipliziere fouriertransformierte Dichte mit Greensfunktion. Keine Doppelschleife noetig, weil rhok und G dasselbe Index-Wrapping verwenden.
for(j=0; j<densGrid_Zellen*densGrid_Zellen; j++){
	rhok[j][0] *= G[j];
	rhok[j][1] *= G[j];
}//for j bis densGrid_Zellen^2


// transformiere rhok in den x-Raum zurueck
fftw_execute(ft_rueck);
for(j=0; j<densGrid_Zellen*densGrid_Zellen; j++){
	rhox[j][0] /= densGrid_Zellen*densGrid_Zellen;
}//for j

//* Test: Gebe Feld rhox aus
out = fopen("test2.txt", "w");
fprintf(out, "#Format: x TAB y TAB rho");
for(j=0; j<densGrid_Zellen; j++){
	for(k=0; k<densGrid_Zellen; k++){
		fprintf(out, "%g \t %g \t %g \n", (j+0.5)*densGrid_Breite, (k+0.5)*densGrid_Breite, rhox[iw(j,k)][0]);
	}//for k
	fprintf(out, "\n");
}//for j
// */




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
