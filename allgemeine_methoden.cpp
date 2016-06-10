//stellt grundlegende Methoden bereit, die an vielen Stellen benötigt werden. Zum Beispiel Abstandsberechnung oder Index-Wrapping


#include "parameter.h"
#include "signaturen2.h"

#include <iostream>


extern const double nachList_Breite;
extern const double L;
extern const int N1;
extern const double rho;
extern const double kapillar_vorfaktor;
extern const double T;


using std::cout;


//berechnet das Quadrat des Abstands zwischen Teilchen i(Typ 1) und Teilchen j(Typ 1). Berücksichtigt periodische Randbedingungen.
double RunZustand::abstand2_11(int i, int j){
	
	//Absolutpositionen
	double xi_abs = r1_git[i][0]*nachList_Breite + r1_rel[i][0];
	double xj_abs = r1_git[j][0]*nachList_Breite + r1_rel[j][0];
	double yi_abs = r1_git[i][1]*nachList_Breite + r1_rel[i][1];
	double yj_abs = r1_git[j][1]*nachList_Breite + r1_rel[j][1];
	
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
}//double RunZustand::abstand2


//berechnet das Quadrat des Abstands zwischen Teilchen i(Typ 1) und Teilchen j(Typ 2). Berücksichtigt periodische Randbedingungen.
double RunZustand::abstand2_12(int i, int j){
	//Absolutpositionen
	double xi_abs = r1_git[i][0]*nachList_Breite + r1_rel[i][0];
	double xj_abs = r2_git[j][0]*nachList_Breite + r2_rel[j][0];
	double yi_abs = r1_git[i][1]*nachList_Breite + r1_rel[i][1];
	double yj_abs = r2_git[j][1]*nachList_Breite + r2_rel[j][1];
	
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
}//double RunZustand::abstand2_12



//berechnet das Quadrat des Abstands zwischen Teilchen i(Typ 2) und Teilchen j(Typ 2). Berücksichtigt periodische Randbedingungen.
double RunZustand::abstand2_22(int i, int j){
	//Absolutpositionen
	double xi_abs = r2_git[i][0]*nachList_Breite + r2_rel[i][0];
	double xj_abs = r2_git[j][0]*nachList_Breite + r2_rel[j][0];
	double yi_abs = r2_git[i][1]*nachList_Breite + r2_rel[i][1];
	double yj_abs = r2_git[j][1]*nachList_Breite + r2_rel[j][1];
	
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
	
}//double RunZustand::abstand2_22





// Dieselben drei Methoden in den Klassen RunZustand::RunDynamik und RunZustand::RunObs
double RunZustand::RunDynamik::abstand2_11(int i, int j){
	
	//Absolutpositionen
	double xi_abs = r1_git[i][0]*nachList_Breite + r1_rel[i][0];
	double xj_abs = r1_git[j][0]*nachList_Breite + r1_rel[j][0];
	double yi_abs = r1_git[i][1]*nachList_Breite + r1_rel[i][1];
	double yj_abs = r1_git[j][1]*nachList_Breite + r1_rel[j][1];
	
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
}//double RunDynamik::abstand2

double RunZustand::RunDynamik::abstand2_12(int i, int j){
	//Absolutpositionen
	double xi_abs = r1_git[i][0]*nachList_Breite + r1_rel[i][0];
	double xj_abs = r2_git[j][0]*nachList_Breite + r2_rel[j][0];
	double yi_abs = r1_git[i][1]*nachList_Breite + r1_rel[i][1];
	double yj_abs = r2_git[j][1]*nachList_Breite + r2_rel[j][1];
	
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
}//double RunDynamik::abstand2_12

double RunZustand::RunDynamik::abstand2_22(int i, int j){
	//Absolutpositionen
	double xi_abs = r2_git[i][0]*nachList_Breite + r2_rel[i][0];
	double xj_abs = r2_git[j][0]*nachList_Breite + r2_rel[j][0];
	double yi_abs = r2_git[i][1]*nachList_Breite + r2_rel[i][1];
	double yj_abs = r2_git[j][1]*nachList_Breite + r2_rel[j][1];
	
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
	
}//double RunDynamik::abstand2_22


double RunZustand::RunObs::abstand2_11(int i, int j){
	
	//Absolutpositionen
	double xi_abs = r1_git[i][0]*nachList_Breite + r1_rel[i][0];
	double xj_abs = r1_git[j][0]*nachList_Breite + r1_rel[j][0];
	double yi_abs = r1_git[i][1]*nachList_Breite + r1_rel[i][1];
	double yj_abs = r1_git[j][1]*nachList_Breite + r1_rel[j][1];
	
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
}//double RunObs::abstand2

double RunZustand::RunObs::abstand2_12(int i, int j){
	//Absolutpositionen
	double xi_abs = r1_git[i][0]*nachList_Breite + r1_rel[i][0];
	double xj_abs = r2_git[j][0]*nachList_Breite + r2_rel[j][0];
	double yi_abs = r1_git[i][1]*nachList_Breite + r1_rel[i][1];
	double yj_abs = r2_git[j][1]*nachList_Breite + r2_rel[j][1];
	
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
}//double RunZustand::abstand2_12

double RunZustand::RunObs::abstand2_22(int i, int j){
	//Absolutpositionen
	double xi_abs = r2_git[i][0]*nachList_Breite + r2_rel[i][0];
	double xj_abs = r2_git[j][0]*nachList_Breite + r2_rel[j][0];
	double yi_abs = r2_git[i][1]*nachList_Breite + r2_rel[i][1];
	double yj_abs = r2_git[j][1]*nachList_Breite + r2_rel[j][1];
	
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
	
}//double RunObs::abstand2_22





//index-Wrapping, ermoeglicht es 1dim-Arrays mit zwei Indices anzusprechen. FFTW erwartet Arrays mit nur einem Index. Keine period. Randbed.
int iw(int i, int j){
	if (i>=densGrid_Zellen || j>= densGrid_Zellen || i<0 || j<0){
		cout << "Index out of Bounds! Aufruf war iw(i=" <<i<<", j=" <<j<<") \n" << endl << flush;
		return -1; //erzeugt absichtlich Segmentation Fault
	}//if index out-of-bounds

	return i*densGrid_Zellen + j;

}//int IndexWrapper (iw)




//initialisiert einen Run
void RunZustand::init(int nr){
	
	this->nr = nr;
	this->obs_nr = 0;
	
	// Teilchenpositionen. Erster Index TeilchenNr, zweiter Index 0 fuer x und 1 fuer y-Richtung
	r1_git = new int*[N1];
	r1_rel = new double*[N1];
	r2_git = new int*[N2];
	r2_rel = new double*[N2];
	for(int i=0; i<N1; i++){
		r1_git[i] = new int[2];
		r1_rel[i] = new double[2];
	}//for i
	for(int i=0; i<N2; i++){
		r2_git[i] = new int[2];
		r2_rel[i] = new double[2];
	}//for i

	
	switch(startpos_methode){
		case 1:
			init_zufallspos();
			break;
		case 2:
			init_pos_aus_datei();
			break;
		case 3:
			init_gitterstart();
			break;
		case 4:
			init_allegleich();
			break;
		case 5:
			init_kreisscheibe();
			break;
		default:
			cout << "Fehler: startpos_methode="<<startpos_methode<<", erlaubt sind 1,2,3,4,5" << endl;
			return;
	}//switch startpos_methode
	
	
	//in RunDynamik und RunObs sollen r1_git usw auf die Felder aus RunZustand zeigen
	dyn.r1_git = r1_git;
	dyn.r2_git = r2_git;
	dyn.r1_rel = r1_rel;
	dyn.r2_rel = r2_rel;
	
	obs.r1_git = r1_git;
	obs.r2_git = r2_git;
	obs.r1_rel = r1_rel;
	obs.r2_rel = r2_rel;
	
	//init Kraftberechnungen
	dyn.init();
	//init Observablen
	obs.obs_init(nr);

}//void RunZustand::init()


//schreibt Jeans-Zeit und Jeans-Länge in cout
void RunZustand::ausgabe_jeansgroessen(){
	// tJ = 1/rho 1/(f^2/eps gamma)
	
	//Dichte
	const double rho = N1/(L*L);
	//Jeans-Zeit
	const double tJ = 1.0/(rho * kapillar_vorfaktor);
	
	//kritische Dichte, close packed HS in 2-dim
	const double rho_c = 2.0/sqrt(3.0);
	
	
	//der Ausdruck 1/wurzel(2rho_c^2/(rho-rho_c)^2 - 1), kommt in der Jeanslänge vor
	const double lJ_faktor = 1.0/sqrt(2.0*rho_c/((rho-rho_c)*(rho-rho_c)) -1.0);
	
	const double kJ = lJ_faktor*sqrt(rho*kapillar_vorfaktor/T);

	cout << "Jeans-Zeit: \t\t" << tJ << endl;
	cout << "Jeans-Wellenzahl:\t" << kJ << endl;
	cout << "Faktor in kJ: \t\t" << lJ_faktor << endl << endl;
	
	
}//void ausgabe_jeansgroessen



//Zahl zu String, keine führenden Nullen
string int_to_string(int zahl){
	std::stringstream ss;
	ss << zahl;
	//ss << setw(n) << setfill('0') << zahl;
	return ss.str();
}//string int_to_string


void RunZustand::RunDynamik::init(){
	WCA_init(); 
	kapkraefte_init();
}//void RunDynamik::init