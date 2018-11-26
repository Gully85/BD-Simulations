//stellt Methoden bereit, die in reinen Auswertungs-Runs die Observablen berechnen 

#include "parameter.h"
#include "signaturen_auswertung.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <vector>
using std::vector; using std::min;

#include <cstdio>



//Erster Index Job (0 bis jobs-1), zweiter Index Zeit (in Einheiten obs_dt, 0 bis obs_anzahl-1), dritter Index Teilchen (bis N1 bzw N2), vierter Index 0/1 für x/y.
extern double**** r1_abs;
extern double**** r2_abs;

extern const double Gamma2_1;

extern const double rhoring_dr;
extern const int rhoring_bins;
extern const double korr_dr;
extern const int korr_bins;
extern const double obs_dt;
extern const int obs_anzahl;


extern const int jobs;

extern const int N1, N2;

extern const double L;



namespace Observablen{
	//Mobilitätsschwerpunkte. Erster Index Run bzw Job, zweiter Index Zeit in Einheiten obs_dt, dritter Index 0/1 für x/y
	double*** Sgamma = NULL;
	double*** S1 = NULL;
	double*** S2 = NULL;
	
	//Dichteprofil. Erster Index Run bzw Job, zweiter Index Zeit in Einheiten obs_dt, dritter Index r in Einheiten rhoring_dr. Im Abstand vom Mobischwerpunkt, dh von Sgamma[][][]
	vector<double>** rho1ring = NULL;
	vector<double>** rho2ring = NULL;
	
	//rho1zen2 ist die Dichte Typ1-Teilchen im Abstand zum Mittelpunkt der Typ2-Teilchen. Erster Index Job/Run, zweiter Index t/obs_dt, dritter Index r/rhoring_dr
	vector<double>** rho1zen1 = NULL;
	vector<double>** rho1zen2 = NULL;
	vector<double>** rho2zen1 = NULL;
	vector<double>** rho2zen2 = NULL;
	
	//g11, g12, g22 sind Korrelationsfunktionen. Erster Index Job/Run, zweiter Index t/obs_dt, dritter Index r/korr_dr
	vector<double>** g11 = NULL;
	vector<double>** g12 = NULL;
	vector<double>** g22 = NULL;

}

/*
//VERMUTLICH KAPUTT! INDICES VON S und Sgamma!
//berechnet für alle Zeiten und alle Runs/Jobs die Mobilitätsschwerpunkte (sum(typ1) r + 1/Gamma2_1 sum(typ2) r)/(N1 + N2/Gamma2_1). Schreibt diese in Sgamma
void berechne_mobischwerpunkte(){
	
	using namespace Observablen;
	
	//Unterscheidung: Sgamma[][][] ist der Mobilitätsschwerpunkt, dh Schwerpunkt aller Teilchen gewichtet mit 1/Gamma statt mit Masse. S1[][][] ist der Mittelpunkt aller Typ1-Teilchen, S2[][][] der Mittelpunkt aller Typ2-Teilchen
	// Für alle dieser Felder gilt, erster Index Run bzw Job (bis jobs-1), zweiter Index t/obs_dt (bis obs_anzahl-1), dritter Index r/rhoring_dr (bis rhoring_bins-1).
	
	//reserviere Speicher für Sgamma
	Sgamma = new double**[jobs];
	S1 = new double**[jobs];
	S2 = new double**[jobs];
	for(int job=0; job<jobs; job++){
		Sgamma[job] = new double*[obs_anzahl];
		S1[job] = new double*[obs_anzahl];
		S2[job] = new double*[obs_anzahl];
		for(int zeit=0; zeit<obs_anzahl; zeit++){
			Sgamma[job][zeit] = new double[2];
			S1[job][zeit] = new double[2];
			S2[job][zeit] = new double[2];
		}//for zeit
	}//for job
	
	
	//berechne Mobischwerpunkte. Schreibe sie je Job in eine Datei
	for(int job=0; job<jobs; job++){
		string dateiname = "Sgamma_run.txt";
		dateiname.insert(10, int_to_string(job+1));
		FILE* out = fopen(dateiname.c_str(), "w");
		fprintf(out, "# Format: t x y x1 y1 x2 y2 \n# Erklärung: (x,y) ist der Mobilitätsschwerpunkt, (x1,y2) der Mittelpunkt Typ1, (x2,y2) der Mittelpunkt Typ2\n\n");
		
		for(int zeit=0; zeit<obs_anzahl; zeit++){
			double x=0.0;
			double y=0.0;
			S1[job][zeit][0] = 0.0;
			S1[job][zeit][1] = 0.0;
			S2[job][zeit][0] = 0.0;
			S2[job][zeit][1] = 0.0;
			
			const double nenner = N1 + N2/(Gamma2_1);
			
			
			
			for(int teilchen=0; teilchen<N2; teilchen++){
				x += r2_abs[job][zeit][teilchen][0];
				y += r2_abs[job][zeit][teilchen][1];
				S2[job][zeit][0] += r2_abs[job][zeit][teilchen][0];
				S2[job][zeit][1] += r2_abs[job][zeit][teilchen][1];
			}//for teilchen bis N2
			x/=Gamma2_1;
			y/=Gamma2_1;
			S2[job][zeit][0] /= N2;
			S2[job][zeit][1] /= N2;
			
			for(int teilchen=0; teilchen<N1; teilchen++){
				x += r1_abs[job][zeit][teilchen][0];
				y += r1_abs[job][zeit][teilchen][1];
				S1[job][zeit][0] += r1_abs[job][zeit][teilchen][0];
				S1[job][zeit][1] += r1_abs[job][zeit][teilchen][1];
			}//for teilchen bis N1
			S1[job][zeit][0] /= N1;
			S1[job][zeit][1] /= N1;
			
			Sgamma[job][zeit][0] = x/nenner;
			Sgamma[job][zeit][1] = y/nenner;
			//Sgamma[job][zeit][0] = 0.5*L;
			//Sgamma[job][zeit][1] = 0.5*L;
			
			fprintf(out, "%g \t %g \t %g \t %g \t %g \t %g \t %g \n\n", zeit*obs_dt, Sgamma[job][zeit][0], Sgamma[job][zeit][1], S1[job][zeit][0], S1[job][zeit][1], S2[job][zeit][0], S2[job][zeit][1]);
			
		}//for zeit
	}//for job
	
}//void berechne_mobischwerpunkte
// */

//berechnet für alle Runs und alle Zeiten Sgamma, S1 und S2
void berechne_mittelpunkte(){
	using namespace Observablen;
	
	// Null setzen, Speicher reservieren
	S1 = new double**[jobs];
	S2 = new double**[jobs];
	Sgamma = new double**[jobs];
	for(int job=0; job<jobs; job++){
		S1[job] = new double*[obs_anzahl];
		S2[job] = new double*[obs_anzahl];
		Sgamma[job] = new double*[obs_anzahl];
		for(int t=0; t<obs_anzahl; t++){
			S1[job][t] = new double[2];
			S2[job][t] = new double[2];
			Sgamma[job][t] = new double[2];
			
			S1[job][t][0] = 0.0;
			S1[job][t][1] = 0.0;
			S2[job][t][0] = 0.0;
			S2[job][t][1] = 0.0;
			Sgamma[job][t][0] = 0.0;
			Sgamma[job][t][1] = 0.0;
		}//for t bis obs_anzahl
	}//for job
	
	//für jeden Job und Zeitpunkt: für beide Teilchentypen Mittelpunkte berechnen, das ist je Raumrichtung eine Mittelung über cos/sin(2pi/L x)
	for(int job=0; job<jobs; job++){
		//schreibe je Job die Mittelpunkte in eine Datei
		string dateiname = "Sgamma_run.txt";
		dateiname.insert(10, int_to_string(job+1));
		FILE* out = fopen(dateiname.c_str(), "w");
		fprintf(out, "# Format: t x y x1 y1 x2 y2 \n# Erklärung: (x,y) ist der Mobilitätsschwerpunkt, (x1,y2) der Mittelpunkt Typ1, (x2,y2) der Mittelpunkt Typ2\n\n");
		
		for(int t=0; t<obs_anzahl; t++){
			double x1quer_re = 0.0; //Typ 1, x-Richtung aufgerollt, x'-Richtung. xquer in den Notizen
			double x1quer_im = 0.0; //Typ 1, x-Richtung aufgerollt, y'-Richtung. yquer in den Notizen
			double y1quer_re = 0.0; //Typ 1, in y-Richtung aufgerollt
			double y1quer_im = 0.0;
			
			double x2quer_re = 0.0; //Typ 2, x-Richtung aufgerollt, x'-Richtung
			double x2quer_im = 0.0;
			double y2quer_re = 0.0;
			double y2quer_im = 0.0;
			
			const double zweipi_L = 2.0*M_PI/L;
			for(int teilchen=0; teilchen<N1; teilchen++){
				double x = r1_abs[job][t][teilchen][0];
				double y = r1_abs[job][t][teilchen][1];
				
				x1quer_re += cos(x*zweipi_L);
				x1quer_im += sin(x*zweipi_L);
				y1quer_re += cos(y*zweipi_L);
				y1quer_im += sin(y*zweipi_L);
			}//for teilchen Typ 1
			for(int teilchen=0; teilchen<N2; teilchen++){
				double x = r2_abs[job][t][teilchen][0];
				double y = r2_abs[job][t][teilchen][1];
				
				x2quer_re += cos(x*zweipi_L);
				x2quer_im += sin(x*zweipi_L);
				y2quer_re += cos(y*zweipi_L);
				y2quer_im += sin(y*zweipi_L);
			}//for teilchen Typ 2
			double theta;
			
			theta = atan2(-x1quer_im, -x1quer_re) + M_PI;
			S1[job][t][0] = theta / zweipi_L;
			
			theta = atan2(-y1quer_im, -y1quer_re) + M_PI;
			S1[job][t][1] = theta / zweipi_L;
			
			theta = atan2(-x2quer_im, -x2quer_re) + M_PI;
			S2[job][t][0] = theta / zweipi_L;
			
			theta = atan2(-y2quer_im, -y2quer_re) + M_PI;
			S2[job][t][1] = theta / zweipi_L;
			
			//quick and quite dirty: Sgamma = (S1 + S2)/2
			Sgamma[job][t][0] = 0.5*(S1[job][t][0] + S2[job][t][0]);
			Sgamma[job][t][1] = 0.5*(S1[job][t][1] + S2[job][t][1]);
			
			fprintf(out, "%g \t %g \t %g \t %g \t %g \t %g \t %g \n\n", t*obs_dt, Sgamma[job][t][0], Sgamma[job][t][1], S1[job][t][0], S1[job][t][1], S2[job][t][0], S2[job][t][1]);
		}//for t
	}//for int job
	
}//void berechne_mittelpunkte



//berechnet für jeden Run und jede Zeit das Dichteprofil rho1ring und rho2ring. Schreibt es in rho1ring[][][]
void berechne_dichteprofile(){
	using namespace Observablen;
	
	//berechne_mobischwerpunkte();
	berechne_mittelpunkte();
	
	rho1ring = new vector<double>*[jobs];
	rho2ring = new vector<double>*[jobs];
	rho1zen1 = new vector<double>*[jobs];
	rho1zen2 = new vector<double>*[jobs];
	rho2zen1 = new vector<double>*[jobs];
	rho2zen2 = new vector<double>*[jobs];
	
	
	//für jeden Run...
	for(int job=0; job<jobs; job++){
		rho1ring[job] = new vector<double>[obs_anzahl];
		rho2ring[job] = new vector<double>[obs_anzahl];
		rho1zen1[job] = new vector<double>[obs_anzahl];
		rho1zen2[job] = new vector<double>[obs_anzahl];
		rho2zen1[job] = new vector<double>[obs_anzahl];
		rho2zen2[job] = new vector<double>[obs_anzahl];
		
		string dateiname = "rhozen_run.txt";
		dateiname.insert(10, int_to_string(job+1));
		FILE* out = fopen(dateiname.c_str(), "w");
		fprintf(out, "#Format: t r rho1 rho2 rho1zen1 rho1zen2 rho2zen1 rho2zen2 \n #Erklärung: rho1 und rho2(r) sind über Abstand vom Mobischwerpunkt aufgetragen, rho1zen2 ist die Dichte von Typ1-Teilchen aufgetragen über Abstand vom Mittelpunkt Typ2.\n\n ");
		
		//für jede Zeit...
		for(int zeit=0; zeit<obs_anzahl; zeit++){
			
			rho1ring[job][zeit].assign(rhoring_bins, 0.0);
			rho1zen1[job][zeit].assign(rhoring_bins, 0.0);
			rho1zen2[job][zeit].assign(rhoring_bins, 0.0);
			
			const double vorfaktor1 = L*L/(N1*rhoring_dr*rhoring_dr*M_PI);
			
			//für jedes Teilchen Typ 1...
			for(int teilchen=0; teilchen<N1; teilchen++){
				//Abstand des Teilchens vom Mobilitätsschwerpunkt. KEINE periodischen Randbedingungen! (weil nicht nötig, Teilchen und Mobischwerpunkt befinden sich in der Mitte der Box)
				double x = r1_abs[job][zeit][teilchen][0];
				double y = r1_abs[job][zeit][teilchen][1];
				
				double dx = x - Sgamma[job][zeit][0];
				double dy = y - Sgamma[job][zeit][1];
				double r = sqrt(dx*dx + dy*dy);
				int bin = (int) (r/rhoring_dr);
				if(0 > bin || bin >= rhoring_bins) bin=rhoring_bins-1;
				rho1ring[job][zeit][bin] += vorfaktor1/(2*bin+1);
				
				dx = x - S1[job][zeit][0];
				dy = y - S1[job][zeit][1];
				r = sqrt(dx*dx + dy*dy);
				bin = (int) (r/rhoring_dr);
				if(0 > bin || bin >= rhoring_bins) bin=rhoring_bins-1;
				rho1zen1[job][zeit][bin] += vorfaktor1/(2*bin+1);
				
				dx = x - S2[job][zeit][0];
				dy = y - S2[job][zeit][1];
				r = sqrt(dx*dx + dy*dy);
				bin = (int) (r/rhoring_dr);
				if(0 > bin || bin >= rhoring_bins) bin=rhoring_bins-1;
				rho1zen2[job][zeit][bin] += vorfaktor1/(2*bin+1);
				
			}//for teilchen Typ 1
			
			
			rho2ring[job][zeit].assign(rhoring_bins, 0.0);
			rho2zen1[job][zeit].assign(rhoring_bins, 0.0);
			rho2zen2[job][zeit].assign(rhoring_bins, 0.0);
			
			const double vorfaktor2 = L*L/(N2*rhoring_dr*rhoring_dr*M_PI);
			
			//für jedes Teilchen Typ 2...
			for(int teilchen=0; teilchen<N2; teilchen++){
				//Abstand des Teilchens vom Mobilitätsschwerpunkt. KEINE periodischen Randbedingungen! (weil nicht nötig, Teilchen und Mobischwerpunkt befinden sich in der Mitte der Box)
				double x = r2_abs[job][zeit][teilchen][0];
				double y = r2_abs[job][zeit][teilchen][1];
				
				double dx = x - S2[job][zeit][0];
				double dy = y - S2[job][zeit][1];
				double r = sqrt(dx*dx + dy*dy);
				int bin = (int) (r/rhoring_dr);
				if(0 > bin || bin >= rhoring_bins) bin=rhoring_bins-1;
				rho2ring[job][zeit][bin] += vorfaktor2/(2*bin+1);
				
				dx = x - S1[job][zeit][0];
				dy = y - S1[job][zeit][1];
				r = sqrt(dx*dx + dy*dy);
				bin = (int) (r/rhoring_dr);
				if(0 > bin || bin >= rhoring_bins) bin=rhoring_bins-1;
				rho2zen1[job][zeit][bin] += vorfaktor2/(2*bin+1);
				
				dx = x - S2[job][zeit][0];
				dy = y - S2[job][zeit][1];
				r = sqrt(dx*dx + dy*dy);
				bin = (int) (r/rhoring_dr);
				if(0 > bin || bin >= rhoring_bins) bin=rhoring_bins-1;
				rho2zen2[job][zeit][bin] += vorfaktor2/(2*bin+1);
			}//for Teilchen Typ 2
			
			//schreibe in Datei
			for(int bin=0; bin<rhoring_bins; bin++){
				fprintf(out, "%g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \n", zeit*obs_dt, bin*rhoring_dr, rho1ring[job][zeit][bin], rho2ring[job][zeit][bin], rho1zen1[job][zeit][bin], rho1zen2[job][zeit][bin], rho2zen1[job][zeit][bin], rho2zen2[job][zeit][bin]);
			}//for bin
			fprintf(out, "\n");
			
		}//for zeit
	}//for job
	
}//void berechne_dichteprofile

//berechnet für jeden Run und jede Zeit die Paarkorrelationsfunktionen g11, g12, g22. Schreibt sie in gij[][][]
void berechne_korrelationsfunktionen(){
	using namespace Observablen;
        using std::cout; using std::endl;
	
	g11 = new vector<double>*[jobs];
	g12 = new vector<double>*[jobs];
	g22 = new vector<double>*[jobs];
	
	for(int job=0; job<jobs; job++){
		g11[job] = new vector<double>[obs_anzahl];
		g12[job] = new vector<double>[obs_anzahl];
		g22[job] = new vector<double>[obs_anzahl];
		
		string dateiname = "korrfunk_run.txt";
		dateiname.insert(12, int_to_string(job+1));
		FILE* out = fopen(dateiname.c_str(), "w");
		fprintf(out, "#Format: t r g11 g12 g22 \n\n");
		
		for(int zeit=0; zeit<obs_anzahl; zeit++){
			g11[job][zeit].assign(korr_bins, 0.0);
			g12[job][zeit].assign(korr_bins, 0.0);
			g22[job][zeit].assign(korr_bins, 0.0);
                        
                        //if (498 == zeit) cout << "z365 zeit=max" << endl;
			
			//g11 zuerst
			const double rho1 = N1/L/L;
			const double vorfaktor11 = 2.0/((N1-1)*M_PI*rho1*korr_dr*korr_dr);
			const double rho2 = N2/L/L;
                        //const double vorfaktor12 = 1.0/((N1-1)*M_PI*rho2*korr_dr*korr_dr);
                        const double vorfaktor12 = V/(N2*N1*M_PI*korr_dr);
			const double vorfaktor22 = 2.0/((N2-1)*M_PI*rho2*korr_dr*korr_dr);
                        
                        
			
			for(int i=0; i<N1; i++){
				for(int j=0; j<i; j++){
					double a2 = abstand2(r1_abs[job][zeit][i], r1_abs[job][zeit][j]);
					
					int bin = (int)(sqrt(a2)/korr_dr);
					if(0 > bin || bin > korr_bins) bin=korr_bins-1;
					
					g11[job][zeit][bin] += vorfaktor11/(2.0*bin+1);
				}//for j, Typ 1
				//if (499 == zeit) cout << "z386 zeit=max" << endl;
				for(int j=0; j<N2; j++){
					double a2 = abstand2(r1_abs[job][zeit][i], r2_abs[job][zeit][j]);
					int bin = (int)(sqrt(a2)/korr_dr);
					if(0 > bin || bin > korr_bins) bin=korr_bins-1;
					g12[job][zeit][bin] += vorfaktor12/(bin*korr_dr);
				}//for j, Typ 2
			}//for i
			
			//if (498 == zeit) cout << "z394 zeit=max" << endl;
			
			for(int i=0; i<N2; i++){
				for(int j=0; j<i; j++){
					double a2 = abstand2(r2_abs[job][zeit][i], r2_abs[job][zeit][j]);
					int bin = (int)(sqrt(a2)/korr_dr);
					if(0 > bin || bin > korr_bins) bin=korr_bins-1;
					g22[job][zeit][bin] += vorfaktor22/(2.0*bin+1);
				}//for j
			}//for i
			
			//std::cout << "z400, zeit=" << zeit << std::endl;
			
			for(int bin=0; bin<korr_bins; bin++)
				fprintf(out, "%g \t %g \t %g \t %g \t %g \n", zeit*obs_dt, bin*korr_dr, g11[job][zeit][bin], g12[job][zeit][bin], g22[job][zeit][bin]);
			
			fprintf(out, "\n");
		}//for zeit
		//std::cout << "z407" << std::endl;
	fclose(out);
	}//for job
	
	//std::cout << "berechne_korrelationsfunktionen() fertig" << std::endl;
	
}//void berechne_korrelationsfunktionen



//Dichteprofile auswerten, in Dateien schreiben
void auswerten_dichteprofile(){
	using namespace Observablen;
	
	//lokale Felder: Mittelwerte und Fehler von Typ1- und Typ2-Dichteprofilen
	vector<double>* mittel1 = new vector<double>[obs_anzahl];
	vector<double>* fehler1 = new vector<double>[obs_anzahl];
	vector<double>* mittel2 = new vector<double>[obs_anzahl];
	vector<double>* fehler2 = new vector<double>[obs_anzahl];
	//vector<double>* mittel3 = new vector<double>[obs_anzahl];
	//vector<double>* fehler3 = new vector<double>[obs_anzahl];
	
	
	
	//Mittelwerte und Fehler von rho1ring und rho2ring, dh Dichten, vom Mobischwerpunkt gemessen
	FILE* out = fopen("dichteprofil_mobischwerpunkt.txt", "w");
	fprintf(out, "#Format: t TAB r TAB rho1 TAB fehler TAB rho2 TAB fehler \n #Abstände werden vom Mobilitätsschwerpunkt gemessen \n\n");
	
	//berechne Mittelwerte und Fehler von rho1
	statistik_1(rho1ring, mittel1, fehler1, obs_anzahl);
	statistik_1(rho2ring, mittel2, fehler2, obs_anzahl);
	
	
	for(int zeit=0; zeit<obs_anzahl; zeit++){
		double t = zeit*obs_dt;
		
		for(int r_bin=0; r_bin<rhoring_bins; r_bin++){
			double r = (r_bin+0.5)*rhoring_dr;
	
			double r1m = mittel1[zeit][r_bin]; //rho1, Mittelwert
			double r1f = fehler1[zeit][r_bin]; //rho1, Fehler
			double r2m = mittel2[zeit][r_bin];
			double r2f = fehler2[zeit][r_bin];
			fprintf(out, "%g \t %g \t %g \t %g \t %g \t %g \n", t, r, r1m, r1f, r2m, r2f);
		}//for zeit
		fprintf(out, "\n");
	}//for zeit
	
	fclose(out);
	
	
	//dasselbe für rho1zen1 und rho2zen1
	out = fopen("rhozen1.txt", "w");
	fprintf(out, "#Format: t TAB r TAB rho1 TAB fehler TAB rho2 TAB fehler \n #Abstände werden vom Mittelpunkt der Typ1-Teilchen gemessen \n\n");
	
	statistik_1(rho1zen1, mittel1, fehler1, obs_anzahl);
	statistik_1(rho2zen1, mittel2, fehler2, obs_anzahl);
	
	for(int zeit=0; zeit<obs_anzahl; zeit++){
		double t = zeit*obs_dt;
		
		for(int r_bin=0; r_bin<rhoring_bins; r_bin++){
			double r = (r_bin+0.5)*rhoring_dr;
	
			double r1m = mittel1[zeit][r_bin]; //rho1, Mittelwert
			double r1f = fehler1[zeit][r_bin]; //rho1, Fehler
			double r2m = mittel2[zeit][r_bin];
			double r2f = fehler2[zeit][r_bin];
			fprintf(out, "%g \t %g \t %g \t %g \t %g \t %g \n", t, r, r1m, r1f, r2m, r2f);
		}//for zeit
		fprintf(out, "\n");
	}//for zeit
	
	fclose(out);
	
	
	//dasselbe für rho1zen2 und rho2zen2
	out = fopen("rhozen2.txt", "w");
	fprintf(out, "#Format: t TAB r TAB rho1 TAB fehler TAB rho2 TAB fehler \n #Abstände werden vom Mittelpunkt der Typ2-Teilchen gemessen \n\n");
	
	statistik_1(rho1zen2, mittel1, fehler1, obs_anzahl);
	statistik_1(rho2zen2, mittel2, fehler2, obs_anzahl);
	
	for(int zeit=0; zeit<obs_anzahl; zeit++){
		double t = zeit*obs_dt;
		
		for(int r_bin=0; r_bin<rhoring_bins; r_bin++){
			double r = (r_bin+0.5)*rhoring_dr;
	
			double r1m = mittel1[zeit][r_bin]; //rho1, Mittelwert
			double r1f = fehler1[zeit][r_bin]; //rho1, Fehler
			double r2m = mittel2[zeit][r_bin];
			double r2f = fehler2[zeit][r_bin];
			fprintf(out, "%g \t %g \t %g \t %g \t %g \t %g \n", t, r, r1m, r1f, r2m, r2f);
		}//for zeit
		fprintf(out, "\n");
	}//for zeit
	
	fclose(out);
	
	
}//void auswerten_dichteprofile


//Korrelationsfunktionen auswerten, in Dateien schreiben
void auswerten_korrelationsfunktionen(){
	using namespace Observablen;
	//Indices von g11: erster Index Job, zweiter Index t/obs_dt, dritter Index r/korr_dr
	
	//lokale Felder
	vector<double>* mittel11 = new vector<double>[obs_anzahl];
	vector<double>* mittel12 = new vector<double>[obs_anzahl];
	vector<double>* mittel22 = new vector<double>[obs_anzahl];
	vector<double>* fehler11 = new vector<double>[obs_anzahl];
	vector<double>* fehler12 = new vector<double>[obs_anzahl];
	vector<double>* fehler22 = new vector<double>[obs_anzahl];
	
	statistik_1(g11, mittel11, fehler11, obs_anzahl);
	statistik_1(g12, mittel12, fehler12, obs_anzahl);
	statistik_1(g22, mittel22, fehler22, obs_anzahl);
	
	FILE* out = fopen("korr.txt", "w");
	fprintf(out, "#Format: t r g11 fehler g12 fehler g22 fehler \n\n");
	
	for(int zeit=0; zeit<obs_anzahl; zeit++){
		double t = zeit*obs_dt;
		
		for(int bin=0; bin<korr_bins; bin++){
			double r = bin*korr_dr;
			
			double g11m = mittel11[zeit][bin];
			double g11f = fehler11[zeit][bin];
			double g12m = mittel12[zeit][bin];
			double g12f = fehler12[zeit][bin];
			double g22m = mittel22[zeit][bin];
			double g22f = fehler22[zeit][bin];
			
			fprintf(out, "%g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \n", t, r, g11m, g11f, g12m, g12f, g22m, g22f);
		}//for bin
		fprintf(out, "\n");
		
	}//for zeit
	
	fclose(out);
}//void auswerten_korrelationsfunktionen


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
			for(ar=0; ar<jobs; ar++)
				mittelwerte[i][j] += input[ar][i][j];
			mittelwerte[i][j] /= jobs; //Summe durch Anzahl
			
			//Fehler bestimmen
			varianzen[i].push_back(0.0);
			for(ar=0; ar<jobs; ar++)
				varianzen[i][j] += (input[ar][i][j]-mittelwerte[i][j])*(input[ar][i][j]-mittelwerte[i][j]);
			varianzen[i][j] /= jobs; //Summe der Fehlerquadrate durch Anzahl. Ist Varianz^2.
			varianzen[i][j] = sqrt(varianzen[i][j]); //Varianz
			varianzen[i][j] /= sqrt((double)jobs); //Standardfehler
			
		}//for j, dritter Index im input bzw zweiter in mittelwerte/varianzen
		
	}//for i, zweiter Index im Input bzw erster in mittelwerte
	
}//void statistik_1



//berechnet das Quadrat des Abstands zwischen Teilchen i und Teilchen j. Berücksichtigt periodische Randbedingungen.
double abstand2(double* ri, double* rj){
	
	//Absolutpositionen
	double xi_abs = ri[0];
	double xj_abs = rj[0];
	double yi_abs = ri[1];
	double yj_abs = ri[1];
	
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

