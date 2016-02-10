			
//stellt Methoden bereit, die in der Auswertung benutzt werden

#pragma once

#include "parameter.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <fftw3.h>


extern const double korr_dr;
extern const int korr_bins;
extern const int ftrho_qbins;
extern const double ftrho_dq;
extern const double obs_dt;
extern const int obs_anzahl;
extern const int runs;

extern const double dq;
extern const double dq_rhoFFTW;
extern const double dq_dqFT;
extern const int rhoFFTW_bins;
extern const int densGrid_Zellen;
const int Z = densGrid_Zellen;
extern const int densGrid_Schema;

extern int** r1_git;
extern double** r1_rel;

extern int** r2_git;
extern double** r2_rel;

extern const double nachList_Breite;
extern const int N1;
extern const int N2;

extern const bool auswerten_rho1FT_normjerun;
extern const bool auswerten_rho2FT_normjerun;

//TODO durchgehen, überspringen wenn N2==0

using std::vector; using std::cout; using std::endl;

//semi-globale Variablen, sichtbar für die Methoden in dieser Datei. Namespace=Dateiname.
namespace auswertung{
	//Paarkorrelationsfunktion. Erster Index Durchlauf, zweiter Index Zeit in Einheiten obs_dt, dritter Index r in Einheiten korr_dr
	vector<double>** g11 = NULL;
	vector<double>** g12 = NULL;
	vector<double>** g22 = NULL;
	


	//fouriertransformierte Dichte, Realteil und Imaginärteil. Erster Index Durchlauf, zweiter Index Zeit/obs_dt, dritter Index q (nachschlagen in order[] und qabs2[])
	vector<double>** ftrho1_re = NULL;
	vector<double>** ftrho1_im = NULL;
	vector<double>** ftrho2_re = NULL;
	vector<double>** ftrho2_im = NULL;

	//Hilfsfelder für fouriertransformierte Dichte
	vector<int> ftrho_order;
	vector<int> ftrho_qabs2;
	vector<int> ftrho_beitraege;



	//Dichte nach Density-Gridding für rho(k) via FFTW. Muss (datei-)global sein, weil FFT vorab geplant wird.
	fftw_complex* rhox = NULL; //Ortsraum. Index-Wrapping.
	fftw_complex* rhok = NULL; //k-Raum. Index-Wrapping, Verschiebung, alternierende Vorzeichen.

	//Fouriertransformierte Dichte, berechnet via FFTW. Erster Index Durchlauf, zweiter Index obs/dt, dritter Index q/dq_rhoFFTW. Kein Index-Wrapping, Verschiebung oder Vorzeichen. Binning.
	vector<double>** rho1FFTW_re = NULL;
	vector<double>** rho1FFTW_im = NULL;
	vector<double>** rho2FFTW_re = NULL;
	vector<double>** rho2FFTW_im = NULL;
	
	//zaehlt für jeden Bin mit, wie viele k drinliegen
	int* rhoFFTW_beitraege;
	
	//speichert je run das rhoFFTW zur Zeit t=0.

	fftw_plan rhoFFTW_plan = NULL;
	
	//für Auswertung mit Johannes' Skript
	FILE* out_rhojerun = NULL;
	
}//namespace auswertung

void init_korrelationsfunktion(){
	using namespace auswertung;
	
	if(g11 != NULL) return;
	
	g11 = new vector<double>*[runs];
	if(0 != N2){
		g12 = new vector<double>*[runs];
		g22 = new vector<double>*[runs];
	}//if N2
	
	for(int i=0; i<runs; i++){
		g11[i] = new vector<double>[obs_anzahl];
		if(0 != N2){
			g12[i] = new vector<double>[obs_anzahl];
			g22[i] = new vector<double>[obs_anzahl];
		}//if N2
	}//for i
	
}//void init_korrelationsfunktion

//Schreibt aktuelle Korrelationsfunktion in g11[ar][t].
void record_korrelationsfunktion11(int ar, int t){
	using namespace auswertung;
	
	//Nullen
	g11[ar][t].assign(korr_bins, 0.0); //Nullen an alle Stellen
	
	//konstanter Vorfaktor
	const double rho1 = N1/L/L;
	const double vorfaktor = 2.0/((N1-1)*M_PI*rho1*korr_dr*korr_dr);
	
	//Schleife über alle Teilchenpaare
	for(int i=1; i<N1; i++)
	for(int j=0; j<i; j++){
		
		double a2 = abstand2_11(i,j);
		
		// if(a2 > korr_bins*korr_bins*korr_dr*korr_dr) continue;
		
		int index = (int) (sqrt(a2)/korr_dr);
		
		if(index < 0 || index >= korr_bins) continue;
		
		g11[ar][t][index] += vorfaktor /(2*index+1);
	}//for j
	
	
}//void record_korrelationsfunktion11

//schreibt aktuelle Korrelationsfunktion in g12[ar][t].
void record_korrelationsfunktion12(int ar, int t){
	if(0 == N2) return;
	
	using namespace auswertung;
	
	g12[ar][t].assign(korr_bins, 0.0);
	
	//konstanter Vorfaktor
	const double rho2 = N2/L/L;
	const double vorfaktor = 1.0/((N1-1)*M_PI*rho2*korr_dr*korr_dr);
	
	//Schleife über Teilchenpaare
	for(int i=0; i<N1; i++)
	for(int j=0; j<N2; j++){
		double a2 = abstand2_12(i,j); 
		int index = (int) (sqrt(a2)/korr_dr);
		
		if(index < 0 || index >= korr_bins) continue;
		
		g12[ar][t][index] += vorfaktor/(2*index+1);
	}//for i,j
	
	
}//void record_korrelationsfunktion12

void record_korrelationsfunktion22(int ar, int t){
	using namespace auswertung;
	if(0 == N2) return;
	
	//Nullen
	g22[ar][t].assign(korr_bins, 0.0);
	
	//konstanter Vorfaktor
	const double rho2 = N2/L/L;
	const double vorfaktor = 2.0/((N2-1)*M_PI*rho2*korr_dr*korr_dr);
	
	//Schleife über alle Teilchenpaare
	for(int i=1; i<N2; i++)
	for(int j=0; j<i; j++){
		double a2 = abstand2_22(i,j);
		int index = (int) (sqrt(a2)/korr_dr);
		
		if(index<0 || index>=korr_bins) continue;
		
		g22[ar][t][index] += vorfaktor/(2*index+1);
	}//for i,j
}//void record_korrelationsfunktion22

// Korrelationsfunktion auswerten: Mittelung, Fehlerbalken berechnen und alles in Datei schreiben
void auswerten_korrelationsfunktion(){
	using namespace auswertung;
	
	vector<double>* g_mittel = new vector<double>[obs_anzahl];
	vector<double>* g_fehler = new vector<double>[obs_anzahl];
	
	//berechne Mittelwerte und Fehler. Schreibe diese in g11_mittel und g11_fehler
	statistik_1(g11, g_mittel, g_fehler, obs_anzahl);
	
	//schreibe in Datei
	FILE* out = fopen("g11.txt", "w");
	fprintf(out, "# Format: r TAB t TAB g11(r,t) TAB Fehlerbalken\n\n");

	for(int j=0; j<obs_anzahl; j++){
		double t = j*obs_dt;
		
		for(int i=0; i<korr_bins; i++){
			double r = i*korr_dr;
			fprintf(out, "%g \t %g \t %g \t %g \n", r, t, g_mittel[j][i], g_fehler[j][i]);
		}//for i bis korr_bins
		fprintf(out, "\n");
	}//for j bis obs_anzahl
	fclose(out);
	
	if(N2 == 0) return;
	
	statistik_1(g22, g_mittel, g_fehler, obs_anzahl);
	
	out = fopen("g22.txt", "w");
	fprintf(out, "#Format: r TAB t TAB g22(r,t) TAB Fehler \n\n");
	
	for(int j=0; j<obs_anzahl; j++){
		double t = j*obs_dt;
		
		for(int i=0; i<korr_bins; i++){
			double r = i*korr_dr;
			fprintf(out, "%g \t %g \t %g \t %g \n", r, t, g_mittel[j][i], g_fehler[j][i]);
		}//for i
	}//for j bis obs_anzahl
	
	fclose(out);
}//void auswerten_korrelationsfunktion


void auswerten_korrelationsfunktion_mixed(){
	if(N2 == 0) return;
	
	using namespace auswertung;
	
	vector<double>* g_mittel = new vector<double>[obs_anzahl];
	vector<double>* g_fehler = new vector<double>[obs_anzahl];
	
	//berechne Mittelwerte und Fehler, schreibe in g_mittel und g_fehler
	statistik_1(g12, g_mittel, g_fehler, obs_anzahl);
	
	//schreibe in Datei
	FILE* out = fopen("g12.txt", "w");
	fprintf(out, "# Format: r TAB t TAB g12(r,t) TAB Fehler \n\n");
	
	for(int t_int=0; t_int<obs_anzahl; t_int++){
		double t = t_int*obs_dt;
		for(int i=0; i<korr_bins; i++){
			double r=i*korr_dr;
			fprintf(out, "%g \t %g \t %g \t %g \n", r, t, g_mittel[t_int][i], g_fehler[t_int][i]);
		}//for i bis korr_bins
	}//for t_int bis obs_anzahl	
	fclose(out);
}//void auswerten_korrelationsfunktion_mixed


void init_ftrho(){
	using namespace auswertung;
	
	if(ftrho1_re != NULL) return;
	
	
	//Vorbereitung: Berechne die Reihenfolge und die anzahl_beiträge, mit der am Ende normiert werden muss.
	ftrho_order.clear();
	ftrho_qabs2.clear();
	ftrho_beitraege.clear();
	
	const int qbins2 = ftrho_qbins*ftrho_qbins;
	
	int stelle_in_order=0;
	
	for(int i=0; i<ftrho_qbins; i++)
	for(int j=0; j<ftrho_qbins; j++){
		int ind2 = i*i + j*j;
		if(ind2 > qbins2) break;
		
		int b = suche(ind2, ftrho_qabs2); // -1 wenn nicht drin, sonst Index wo ind2 in qabs steht
		
		if(-1 == b){
			ftrho_qabs2.push_back(ind2);
			ftrho_beitraege.push_back(1);
			ftrho_order.push_back(stelle_in_order++);
		}//if nicht drin
		else{
			ftrho_beitraege[b]++;
			ftrho_order.push_back(b);
		}//else drin		
	}//for i,j
	
	//reserviere Vektoren für alle record()-Aufrufe
	ftrho1_re = new vector<double>*[runs];
	ftrho1_im = new vector<double>*[runs];
	if(0 != N2){
		ftrho2_re = new vector<double>*[runs];
		ftrho2_im = new vector<double>*[runs];
	}//if N2
	for(int i=0; i<runs; i++){
		ftrho1_re[i] = new vector<double>[obs_anzahl];
		ftrho1_im[i] = new vector<double>[obs_anzahl];
		if(0 != N2){
			ftrho2_re[i] = new vector<double>[obs_anzahl];
			ftrho2_im[i] = new vector<double>[obs_anzahl];
		}//if N2
	}//for i bis runs
	
}//void init_ftrho

//schreibt aktuelles rho(k) in ftrho_re[ar][t] und ftrho_im
void record_ftrho_unkorrigiert(int ar, int t){
	using namespace auswertung;
	
	//Vektor auf richtige Länge bringen, Nullen reinschreiben
	ftrho1_re[ar][t].assign(ftrho_qabs2.size(), 0.0);
	ftrho1_im[ar][t].assign(ftrho_qabs2.size(), 0.0);
	if(N2 != 0){
		ftrho2_re[ar][t].assign(ftrho_qabs2.size(), 0.0);
		ftrho2_im[ar][t].assign(ftrho_qabs2.size(), 0.0);
	}//if N2
	
	const int qbins2 = ftrho_qbins*ftrho_qbins;
	
	int stelle_in_order=0;
	
	for(int i=0; i<ftrho_qbins; i++)
	for(int j=0; j<ftrho_qbins; j++){
		int ind2 = i*i + j*j;
		
		if(ind2 > qbins2) break;
		
		//berechne die Summen cos(k*r_i) und sin(k*r_i)
		double sum_c = 0.0;
		double sum_s = 0.0;
		for(int teilchenNr=0; teilchenNr<N1; teilchenNr++){
			double x = r1_git[teilchenNr][0]*nachList_Breite + r1_rel[teilchenNr][0];
			double y = r1_git[teilchenNr][1]*nachList_Breite + r1_rel[teilchenNr][0];
			
			double skalarprodukt = ftrho_dq*(x*i + y*j);
			
			sum_c += cos(skalarprodukt);
			sum_s += sin(skalarprodukt);
		}//for teilchen
		
		//addiere an richtige Stelle in S[.] diese Summe
		int stelle = ftrho_order[stelle_in_order++];
		ftrho1_re[ar][t][stelle] += sum_c;
		ftrho1_im[ar][t][stelle] += sum_s;
		
		if(N2 != 0){
			sum_c = 0.0;
			sum_s = 0.0;
			for(int teilchenNr=0; teilchenNr<N2; teilchenNr++){
				double x = r2_git[teilchenNr][0]*nachList_Breite + r2_rel[teilchenNr][0];
				double y = r2_git[teilchenNr][1]*nachList_Breite + r2_rel[teilchenNr][1];
				
				double skalarprodukt = ftrho_dq*(x*i + y*j);
				
				sum_c += cos(skalarprodukt);
				sum_s += sin(skalarprodukt);
			}//for teilchenNr bis N2
			ftrho2_re[ar][t][stelle] += sum_c;
			ftrho2_im[ar][t][stelle] += sum_s;
		}//if N2
		
	}//for i,j
}//void record_ftrho

// FTrho auswerten: ruft korrigiere() und statistik() auf und schreibt Mittelwerte/Fehlerbalken in Datei. 
void auswerten_ftrho(){
	
	using namespace auswertung;
	
	korrigiere_ftrho();
	
	//erster Index Zeit, zweiter Index q-Werte nachschlagen in qabs2
	vector<double>* ftrho_re_mittel = new vector<double>[obs_anzahl];
	vector<double>* ftrho_re_fehler = new vector<double>[obs_anzahl];
	vector<double>* ftrho_im_mittel = new vector<double>[obs_anzahl];
	vector<double>* ftrho_im_fehler = new vector<double>[obs_anzahl];
	
	
	statistik_1(ftrho1_re, ftrho_re_mittel, ftrho_re_fehler, obs_anzahl);
	statistik_1(ftrho1_im, ftrho_im_mittel, ftrho_im_fehler, obs_anzahl);
	
	
	FILE* out = fopen("ftrho1.txt", "w");
	fprintf(out, "# Format: q TAB t TAB re(rhotilde) TAB Fehler davon TAB im(rhotilde) TAB Fehler davon TAB abs TAB fehler davon \n\n");
	for(int t_int=0; t_int<obs_anzahl; t_int++){ // t in Einheiten obs_dt
		double t = t_int*obs_dt;
		for(int i=0; i<ftrho1_re[0][0].size(); i++){ // i läuft durch die q-Werte
		
			double q = ftrho_dq*sqrt(ftrho_qabs2[i]);
			
			double re = ftrho_re_mittel[t_int][i];
			double im = ftrho_im_mittel[t_int][i];
			double refehler = ftrho_re_fehler[t_int][i];
			double imfehler = ftrho_im_fehler[t_int][i];
			
			fprintf(out, "%g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \n", q, t, re, refehler, im, imfehler, sqrt(re*re+im*im), sqrt(refehler*refehler+imfehler*imfehler) );
		}// for i, q-Werte
	}//for t_int, Zeit
	fclose(out);
	
	if(N2 == 0) return;
	
	statistik_1(ftrho2_re, ftrho_re_mittel, ftrho_re_fehler, obs_anzahl);
	statistik_1(ftrho2_im, ftrho_im_mittel, ftrho_im_fehler, obs_anzahl);
	
	out = fopen("ftrho2.txt", "w");
	fprintf(out, "# Format: q TAB t TAB re TAB Fehler TAB im TAB Fehler TAB abs TAB Fehler \n\n");
	for(int t_int=0; t_int<obs_anzahl; t_int++){
		double t = t_int*obs_dt;
		for(int i=0; i<ftrho2_re[0][0].size(); i++){
			double q = ftrho_dq*ftrho_qabs2[i];
			
			double re = ftrho_re_mittel[t_int][i];
			double im = ftrho_im_mittel[t_int][i];
			double refehler = ftrho_re_fehler[t_int][i];
			double imfehler = ftrho_im_fehler[t_int][i];
			
			fprintf(out, "%g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \n", q, t, re, refehler, im, imfehler, sqrt(re*re+im*im), sqrt(refehler*refehler+imfehler*imfehler));
		}//for i, q-Werte
		
	}//for t_int bis obs_anzahl
	
	
	delete[] ftrho_re_mittel;
	delete[] ftrho_re_fehler;
	delete[] ftrho_im_mittel;
	delete[] ftrho_im_fehler;
	
}//void auswerten_ftrho

//dividiert Korrekturen aus ftrho raus. Für alle runs und obs_anzahl, nur einmal aufrufen
void korrigiere_ftrho(){
	using namespace auswertung;
	for(int ar=0; ar<runs; ar++)
	for(int t=0; t<obs_anzahl; t++){
		double faktor = 1.0;
		for(int i=0; i<ftrho1_re[0][0].size(); i++){
			ftrho1_re[ar][t][i] /= faktor*ftrho_beitraege[i];
			ftrho1_im[ar][t][i] /= faktor*ftrho_beitraege[i];
			if(N2 != 0){
				ftrho2_re[ar][t][i] /= faktor*ftrho_beitraege[i];
				ftrho2_im[ar][t][i] /= faktor*ftrho_beitraege[i];
			}
		}//for i bis size
	}//for run, t
}//void korrigiere_ftrho

//suche Position im vector, an der die Zahl a steht. Wenn nicht drin, gebe -1 zurück
int suche(int a, vector<int> v){
	
	for(int i=0; i<v.size(); i++)
		if(v[i] == a)
			return i;
		
	return -1;
		
}//int suche


//reserviere Speicher und plane FFTs für rho(k) via FFTW. Variable vector<double>rhoFFTW[run][t][k/dq_rhoFFTW]
void init_rhoFFTW(){
	using namespace auswertung;
	
	if(NULL != rho1FFTW_re) return;
	
	//// reserviere Speicher
	
	// Dichte nach Density-Gridding
	rhox = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Z*Z);
	rhok = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Z*Z);
		
	// Hier steht das berechnete rho(k) drin. Erster Index Durchlauf, zweiter Index t/obs_dt, dritter Index q/dq_rhoFFTW
	rho1FFTW_re = new vector<double>*[runs];
	rho1FFTW_im = new vector<double>*[runs];
	rho2FFTW_re = new vector<double>*[runs];
	rho2FFTW_im = new vector<double>*[runs];
	for(int i=0; i<runs; i++){
		rho1FFTW_re[i] = new vector<double>[obs_anzahl];
		rho1FFTW_im[i] = new vector<double>[obs_anzahl];
		rho2FFTW_re[i] = new vector<double>[obs_anzahl];
		rho2FFTW_im[i] = new vector<double>[obs_anzahl];
	}//for i bis runs
	
	
	
	//Imaginärteil von rhox auf Null setzen, der wird sich nie ändern
	for(int i=0; i<Z; i++)
	for(int j=0; j<Z; j++)
		rhox[iw(i,j)][1] = 0.0;
	
	rhoFFTW_plan = fftw_plan_dft_2d(Z, Z, rhox, rhok, FFTW_FORWARD, FFTW_PATIENT);
	
	const int rhoFFTW_bins = 0.5*Z*dq_dqFT;
	rhoFFTW_beitraege = new int[rhoFFTW_bins];
	for(int i=0; i<rhoFFTW_bins; i++){
		rhoFFTW_beitraege[i]=0;
	}//for i
	
	//herausfinden, wie viele |k| in jedem Bin landen werden.
	for(int i=0; i<Z; i++){
		double x_dq = ((i+Z/2)%Z - Z/2);
		for(int j=0; j<Z; j++){
			double y_dq = ((j+Z/2)%Z - Z/2);
			
			int bin = (int)(dq_dqFT * sqrt(x_dq*x_dq + y_dq*y_dq) + 0.5);
			if(bin < 0 || bin > rhoFFTW_bins){
				//cout << "Fehler: Index " << bin << " out of bounds! Nur 0 bis " << rhoFFTW_bins-1 << "zulässig!" << endl;
				continue;
			}//if out of bounds
			rhoFFTW_beitraege[bin]++;
		}//for j
	}//for i
	
	/* Test: Feldgröße und beitraege auf Konsole ausgeben
	cout << "rhoFFTW wird in  " << rhoFFTW_bins << " Bins einsortiert." << endl;
	cout << "dq="<<dq<< ", dq_FFT="<<dq_rhoFFTW<< endl;
	
	cout << "Anzahl Beiträge je Bin:" << endl;
	for(int i=0; i<rhoFFTW_bins; i++){
		cout << rhoFFTW_beitraege[i] << "  ";
		
	}//for i bis rhoFFTW_bins
	cout << endl;
	// */
	
	out_rhojerun = fopen("rhok_je_run.txt", "w");
	
}//void init_rhoFFTW

//berechnet aktuelles rho(k) via FFTW, schreibt es in rhoFFTW[ar][t][.]
void record_rhoFFTW(int run, int t){
	
	using namespace auswertung;
	
	//Density-Gridding, dh berechne aus Teilchenpositionen eine kontinuierliche Dichte
	switch(densGrid_Schema){
		case 0: gridDensity_NGP(rhox, 1.0, 0.0); break;
		case 1: gridDensity_CIC(rhox, 1.0, 0.0); break;
		case 2: gridDensity_TSC(rhox, 1.0, 0.0); break;
		default: cout << "Density-Gridding-Schema (NearestGridPoint/CloudInCell/TSC) nicht erkannt! densGrid_Schema="<<densGrid_Schema<<", zulaessig sind nur 0,1,2." << endl; 
	}//switch
	
	//FFT ausführen
	fftw_execute(rhoFFTW_plan);
	//jetzt steht in fftw_complex* rhok die fouriertransformierte Dichte. Index-Wrapping, Verschiebung, alternierende Vorzeichen.
	
	//Vorzeichen korrigieren
	for(int i=0; i<Z; i++)
	for(int j=0; j<Z; j++)
		if((i+j) % 2 == 1){ //falls i+j ungerade ist, Vorzeichen wechseln
			rhok[iw(i,j)][0] *= -1.0; //re
			rhok[iw(i,j)][1] *= -1.0; //im
		}//if i+j ungerade
	
	//Binning, dabei Verschiebung und Index-Wrapping berücksichtigen. 
	rho1FFTW_re[run][t].assign(rhoFFTW_bins, 0.0);
	rho1FFTW_im[run][t].assign(rhoFFTW_bins, 0.0);
	
	for(int i=0; i<Z; i++){
		double x_dq = ((i+Z/2)%Z - Z/2);
		for(int j=0; j<Z; j++){
			double y_dq = ((j+Z/2)%Z - Z/2);
			
			int bin = (int)(dq_dqFT * sqrt(x_dq*x_dq + y_dq*y_dq) + 0.5);
			if(bin < 0 || bin > rhoFFTW_bins){
				//cout << "Fehler: Index " << bin << " out of bounds! Nur 0 bis " << rhoFFTW_bins-1 << "zulässig!" << endl;
				continue;
			}//if out of bounds
			
			rho1FFTW_re[run][t][bin] += rhok[iw(i,j)][0];
			rho1FFTW_im[run][t][bin] += rhok[iw(i,j)][1];
		}//for j
	}//for i
	
	//Skalierung richtig machen: Faktor Z^2/L^2 rausdividieren
	for(int i=0; i<rhoFFTW_bins; i++){
		rho1FFTW_re[run][t][i] *= L*L/(Z*Z)/rhoFFTW_beitraege[i];
		rho1FFTW_im[run][t][i] *= L*L/(Z*Z)/rhoFFTW_beitraege[i];

	}//for i
	
	//Jetzt dasselbe für Typ 2
	if(0 == N2) return;
	
	//Density-Gridding
	switch(densGrid_Schema){
		case 0: gridDensity_NGP(rhox, 0.0, 1.0); break;
		case 1: gridDensity_CIC(rhox, 0.0, 1.0); break;
		case 2: gridDensity_TSC(rhox, 0.0, 1.0); break;
		default: cout << "Density-Gridding-Schema (NearestGridPoint/CloudInCell/TSC) nicht erkannt! densGrid_Schema="<<densGrid_Schema<<", zulaessig sind nur 0,1,2." << endl; 
	}//switch
	
	//FFT ausführen
	fftw_execute(rhoFFTW_plan);
	
	//Vorzeichen korrigieren
	for(int i=0; i<Z; i++)
	for(int j=0; j<Z; j++)
		if((i+j)%2 == 1){
			rhok[iw(i,j)][0] *= -1.0;
			rhok[iw(i,j)][1] *= -1.0;
		}//if i+j ungerade
	
	//Binning, dabei Verschiebung und Index-Wrapping berücksichtigen
	rho2FFTW_re[run][t].assign(rhoFFTW_bins, 0.0);
	rho2FFTW_im[run][t].assign(rhoFFTW_bins, 0.0);
	
	for(int i=0; i<Z; i++){
		double x_dq = ((i+Z/2)%Z - Z/2);
		for(int j=0; j<Z; j++){
			double y_dq = ((j+Z/2)%Z - Z/2);
			
			int bin = (int)(dq_dqFT * sqrt(x_dq*x_dq + y_dq*y_dq) + 0.5);
			if(bin < 0 || bin > rhoFFTW_bins){
				//cout << "Fehler: Index " << bin << " out of bounds! Nur 0 bis " << rhoFFTW_bins-1 << "zulässig!" << endl;
				continue;
			}//if out of bounds
			
			rho2FFTW_re[run][t][bin] += rhok[iw(i,j)][0];
			rho2FFTW_im[run][t][bin] += rhok[iw(i,j)][1];
		}//for j
	}//for i
	
	//Skalierung: Faktor Z^2/L^2 rausdividieren
	for(int i=0; i<rhoFFTW_bins; i++){
		rho2FFTW_re[run][t][i] *= L*L/(Z*Z)/rhoFFTW_beitraege[i];
		rho2FFTW_im[run][t][i] *= L*L/(Z*Z)/rhoFFTW_beitraege[i];
	}//for i
	
}//void record_rhoFFTW

//berechnet Mittelwerte und Fehler von rhoFFT, schreibt in Datei rhoFFT.txt
void auswerten_rhoFFTW(){
	
	using namespace auswertung;
	
	
	//* Für Johannes' Auswertung
	fprintf(out_rhojerun, "schreibe rho je Run einzeln in diese Datei \n");
	
	for(int run=0; run<runs; run++){
		
		for(int t_int=0; t_int<obs_anzahl; t_int++){
			for(int i=0; i<rhoFFTW_bins; i++){
				double renorm = rho1FFTW_re[run][0][i];
				double imnorm = rho1FFTW_im[run][0][i];
				double t = t_int*obs_dt;
				double k = i*dq_rhoFFTW;
				double re = rho1FFTW_re[run][t_int][i];
				double im = rho1FFTW_im[run][t_int][i];
				
				fprintf(out_rhojerun, "%g \t %g \t %d \t %g 0.0 %g 0.0 %g \t %g \n", t, k, i, re, im, renorm, imnorm);
			}//for i bis rhoFFTW_bins
		}//for t_int
		
		fprintf(out_rhojerun, "\n");
	}//for run
	fclose(out_rhojerun);
	// */
	
	//erster Index: t/obs_dt, zweiter Index: q/dq
	vector<double>* rhoFFT_re_mittel = new vector<double>[obs_anzahl];
	vector<double>* rhoFFT_re_fehler = new vector<double>[obs_anzahl];
	vector<double>* rhoFFT_im_mittel = new vector<double>[obs_anzahl];
	vector<double>* rhoFFT_im_fehler = new vector<double>[obs_anzahl];
	
	//berechne Mittelwerte und Fehler, Real- und Imaginärteil einzeln
	statistik_1(rho1FFTW_re, rhoFFT_re_mittel, rhoFFT_re_fehler, obs_anzahl);
	statistik_1(rho1FFTW_im, rhoFFT_im_mittel, rhoFFT_im_fehler, obs_anzahl);
	
	FILE* out = fopen("rho1FFT.txt", "w");
	fprintf(out, "#Format: q TAB t TAB rho_re TAB Fehler davon TAB rho_im TAB Fehler davon TAB rho_abs TAB Fehler davon\n");
	
	for(int t_int=0; t_int<obs_anzahl; t_int++){
		double t = t_int*obs_dt;
		for(int i=0; i<rhoFFTW_bins; i++){
		
			double q = dq_rhoFFTW*i;
			
			double re = rhoFFT_re_mittel[t_int][i];
			double refehler = rhoFFT_re_fehler[t_int][i];
			double im = rhoFFT_im_mittel[t_int][i];
			double imfehler = rhoFFT_im_fehler[t_int][i];
			
			fprintf(out, "%g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \n", q, t, re, refehler, im, imfehler, sqrt(re*re+im*im), sqrt(refehler*refehler + imfehler*imfehler));
			
		}//for i
		fprintf(out, "\n");
	}//for t_int
	fclose(out);
	
	//Jetzt dasselbe für Typ 2
	if(0 == N2) return;
	
	statistik_1(rho2FFTW_re, rhoFFT_re_mittel, rhoFFT_re_fehler, obs_anzahl);
	statistik_1(rho2FFTW_im, rhoFFT_im_mittel, rhoFFT_im_fehler, obs_anzahl);
	
	out = fopen("rho2FFT.txt", "w");
	fprintf(out, "# Format: q TAB t TAB re TAB fehler TAB im TAB fehler TAB abs TAB fehler \n\n");
	
	for(int t_int=0; t_int<obs_anzahl; t_int++){
		double t = t_int*obs_dt;
		for(int i=0; i<rhoFFTW_bins; i++){
			double q = dq_rhoFFTW*i;
			
			double re = rhoFFT_re_mittel[t_int][i];
			double refehler = rhoFFT_re_fehler[t_int][i];
			double im = rhoFFT_im_mittel[t_int][i];
			double imfehler = rhoFFT_im_fehler[t_int][i];
			
			fprintf(out, "%g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \n", q, t, re, refehler, im, imfehler, sqrt(re*re+im*im), sqrt(refehler*refehler + imfehler*imfehler));
		}//for i
		fprintf(out, "\n");
	}//for t_int
	fclose(out);
	
	
	
	delete[] rhoFFT_re_mittel;
	delete[] rhoFFT_re_fehler;
	delete[] rhoFFT_im_mittel;
	delete[] rhoFFT_im_fehler;
	
	if(auswerten_rhoFT_normjerun) auswerten_rhoFFTW_normjerun();
	
}//void auswerten_rhoFFTW

//berechnet eine Variante von rhoFFT, bei der zuerst jeder Run mit seinem Startwert normiert wird und dann erst über Runs gemittelt. Schreibt in rhoFFTW_re und _im
void auswerten_rhoFFTW_normjerun(){
	
	using namespace auswertung;
	
	/// EVTL: Erkenne, welche Werte ignoriert werden sollen, weil r(k,0) zu klein
	
	///normalisiere jeden Run mit seinem Startwert
	//vector<double>rhoFFTWre[run][t][k/dq_rhoFFTW]
	for(int run=0; run<runs; run++){
		for(int i=0; i<rho1FFTW_re[0][0].size(); i++){
			
			double re_norm = 1.0/(rho1FFTW_re[run][0][i]); //eins durch rho(k,0) muss überall dranmultipliziert werden
			double im_norm = 1.0/(rho1FFTW_im[run][0][i]);
			
			for(int t_int=0; t_int<obs_anzahl; t_int++){
				rho1FFTW_re[run][t_int][i] *= re_norm;
				rho1FFTW_im[run][t_int][i] *= im_norm;
			}//for i
		}//for t_int
		
	}//for run
	
	// werte normalisierte Runs aus wie oben
	vector<double>* rhoFFT_re_mittel = new vector<double>[obs_anzahl];
	vector<double>* rhoFFT_im_mittel = new vector<double>[obs_anzahl];
	vector<double>* rhoFFT_re_fehler = new vector<double>[obs_anzahl];
	vector<double>* rhoFFT_im_fehler = new vector<double>[obs_anzahl];	
	
	//berechne Mittelwerte und Fehler, Real- und Imaginärteil einzeln
	statistik_1(rho1FFTW_re, rhoFFT_re_mittel, rhoFFT_re_fehler, obs_anzahl);
	statistik_1(rho1FFTW_im, rhoFFT_im_mittel, rhoFFT_im_fehler, obs_anzahl);
	
	FILE* out = fopen("rho1FFT_njr.txt", "w");
	
	fprintf(out, "#Format: q TAB t TAB re TAB fehler TAB im TAB fehler TAB abs TAB fehler \n \n");
	
	for(int t_int=0; t_int<obs_anzahl; t_int++){
		double t = t_int*obs_dt;
		
		for(int k=0; k<rhoFFTW_bins; k++){
			double q = k*dq_rhoFFTW;
			
			double re = rhoFFT_re_mittel[t_int][k];
			double im = rhoFFT_im_mittel[t_int][k];
			double refehler = rhoFFT_re_fehler[t_int][k];
			double imfehler = rhoFFT_im_fehler[t_int][k];
			
			double abs = sqrt(re*re + im*im);
			double absfehler = sqrt(refehler*refehler + imfehler*imfehler);
			
			fprintf(out, "%g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \n", q, t, re, refehler, im, imfehler, abs, absfehler);
			
		}//for k
		fprintf(out, "\n");
	}//for t_int
	
	fclose(out);
	
	
	//Jetzt dasselbe für Typ 2
	if (0 == N2) return;
	
	//normalisiere jeden Run mit seinem Startwert
	for(int run=0; run<runs; run++){
		for(int i=0; i<rho2FFTW_re[0][0].size(); i++){
			double re_norm = 1.0/rho2FFTW_re[run][0][i];
			double im_norm = 1.0/rho2FFTW_im[run][0][i];
			
			for(int t_int=0; t_int<obs_anzahl; t_int++){
				rho2FFTW_re[run][t_int][i] *= re_norm;
				rho2FFTW_im[run][t_int][i] *= im_norm;
			}//for t_int
		}//for i
	}//for run
	
	// werte normalisierte Runs aus wie üblich. Nutze die bereits vorhandenen rhoFFT_re_mittel etc
	statistik_1(rho2FFTW_re, rhoFFT_re_mittel, rhoFFT_re_fehler, obs_anzahl);
	statistik_1(rho2FFTW_im, rhoFFT_im_mittel, rhoFFT_im_fehler, obs_anzahl);
	
	out = fopen("rho2FFT_njr.txt", "w");
	
	fprintf(out, "# Format: q TAB t TAB re TAB fehler TAB im TAB fehler TAB abs TAB fehler \n\n");
	for(int t_int=0; t_int<obs_anzahl; t_int++){
		double t = t_int*obs_dt;
		
		for(int k=0; k<rhoFFTW_bins; k++){
			double q = k*dq_rhoFFTW;
			
			double re = rhoFFT_re_mittel[t_int][k];
			double im = rhoFFT_im_mittel[t_int][k];
			double refehler = rhoFFT_re_fehler[t_int][k];
			double imfehler = rhoFFT_im_fehler[t_int][k];
			
			double abs = sqrt(re*re + im*im);
			double absfehler = sqrt(refehler*refehler + imfehler*imfehler);
			
			fprintf(out, "%g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \n", q, t, re, refehler, im, imfehler, abs, absfehler);
			
		}//for k
		fprintf(out, "\n");
	}//for t_int
	
	fclose(out);
	
	
	
}//void auswerten_rhoFFTW_normjerun




//berechnet Mittelwerte und Varianzen von input[][][]. Mittelung über den ersten Index, der bis runs geht. Der zweite geht bis anzahl. mittelwerte und varianzen müssen schon (leer) initialisiert sein.
void statistik_1(vector<double>**& input, vector<double>*& mittelwerte, vector<double>*& varianzen, int anzahl){
	int ar; //aktueller run
	int i; // aktuelle Stelle im vector, vorletzter index
	int j; // aktuelle Stelle im vector, letzter Index
	
	
	/*
	//erster Index Zeit in Einheiten obs_dt, zweiter Index r in Einheiten korr_dr
	if(mittelwerte == NULL){
		mittelwerte = new vector<double>[anzahl];
		varianzen = new vector<double>[anzahl];
	}//if mittelwerte==NULL
	else cout << "Warnung: statistik_1() aufgerufen, Ergebnisvektoren schon benutzt!" << std::endl;
	// */
	
	for(i=0; i<anzahl; i++){ //zweiter Index im Input, erster in mittelwerte
		for(j=0; j<input[0][0].size(); j++){ //dritter Index im Input, zweiter in Mittelwerte
			//Mittelwert bestimmen
			mittelwerte[i].push_back(0.0); //an Stelle [i][j] eine Null
			for(ar=0; ar<runs; ar++)
				mittelwerte[i][j] += input[ar][i][j];
			mittelwerte[i][j] /= runs; //Summe durch Anzahl
			
			//Fehler bestimmen
			varianzen[i].push_back(0.0);
			for(ar=0; ar<runs; ar++)
				varianzen[i][j] += (input[ar][i][j]-mittelwerte[i][j])*(input[ar][i][j]-mittelwerte[i][j]);
			varianzen[i][j] /= runs; //Summe der Fehlerquadrate durch Anzahl. Ist Varianz^2.
			varianzen[i][j] = sqrt(varianzen[i][j]); //Varianz
			varianzen[i][j] /= sqrt((double)runs); //Standardfehler
			
		}//for j, dritter Index im input bzw zweiter in mittelwerte/varianzen
		
	}//for i, zweiter Index im Input bzw erster in mittelwerte
	
}//void statistik_1




