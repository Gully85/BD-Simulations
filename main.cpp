#include "parameter.h"
#include "signaturen2.h"


#include <cstdlib>
#define _USE_MATH_DEFINES
#include <math.h>
#include <fftw3.h>
#include <algorithm>
#include <iostream>
#include <omp.h>
#include <ctime>



using std::cout; using std::endl; using std::flush; using std::vector; using std::min; using std::string;


// importiere Variablen aus parameter.h, siehe dort, was die Variable tut.
extern const int N1;
extern const double L;

extern const int densGrid_Schema, densGrid_Zellen;
extern const double densGrid_Breite;

extern const int nachList_Zellen;
extern const double nachList_Breite;

extern const double lambda_kapillar;
extern const double kapillar_vorfaktor;
extern const double T;

extern const double dq_rhoFFTW;
extern const int rhoFFTW_bins;

extern const int obs_anzahl;
extern const double obs_dt;

extern const int runs;
extern const int maxThreads;

extern const bool auswerten_korrfunk;
extern const bool auswerten_korrfunk_mixed;
extern const bool auswerten_rhovonk;
extern const bool auswerten_rhoviaFFTW;
extern const bool auswerten_rhoFT_normjerun;
extern const bool auswerten_animation;
extern const bool debugmode;
extern const bool krafttestmode;

extern const double max_reisedistanz;

/// Globale Felder. Mit Null initialisiert.
/*
// Teilchenpositionen Typ 1, Gittervektor. Erster Index TeilchenNr, zweiter Index 0 fuer x und 1 fuer y-Richtung
int** r1_git = NULL;
// dasgleiche, Vektor innerhalb der Zelle
double** r1_rel = NULL;

//dasgleiche für Typ 2
int** r2_git = NULL;
double** r2_rel = NULL;
*/

/////////// TIMING INFO START /////////////////////

// times are measured in seconds, except when stated otherwise
time_t program_started_time;
extern const int max_write_interval;
/////////// TIMING INFO ENDE //////////////////////


int main(){

    program_started_time = time(0);
	init_rng(); //Random Seed
	//fftw_init_threads();

	///Felder, in die die Ergebnisse der Runs geschrieben werden
	//Paarkorrelationsfunktion. Erster Index Durchlauf, zweiter Index Zeit in Einheiten obs_dt, dritter Index r in Einheiten korr_dr
	vector<double>** g11 = NULL;
	vector<double>** g12 = NULL;
	vector<double>** g22 = NULL;
	
	//fouriertransformierte Dichte, Realteil und Imaginärteil. Erster Index Durchlauf, zweiter Index Zeit/obs_dt, dritter Index q (nachschlagen in order[] und qabs2[])
	vector<double>** ftrho1_re = NULL;
	vector<double>** ftrho1_im = NULL;
	vector<double>** ftrho2_re = NULL;
	vector<double>** ftrho2_im = NULL;
	
	//Fouriertransformierte Dichte, berechnet via FFTW. Erster Index Durchlauf, zweiter Index obs/dt, dritter Index q/dq_rhoFFTW. Kein Index-Wrapping, Verschiebung oder Vorzeichen. Binning, q/dq_rhoFFTW läuft von 0 bis rhoFFTW_bins-1
	vector<double>** rho1FFTW_re = NULL;
	vector<double>** rho1FFTW_im = NULL;
	vector<double>** rho2FFTW_re = NULL;
	vector<double>** rho2FFTW_im = NULL;
	
	{ //reserviere Speicher für diese Felder
	g11 = new vector<double>*[runs];
	g12 = new vector<double>*[runs];
	g22 = new vector<double>*[runs];
	ftrho1_re = new vector<double>*[runs];
	ftrho1_im = new vector<double>*[runs];
	ftrho2_re = new vector<double>*[runs];
	ftrho2_im = new vector<double>*[runs];
	rho1FFTW_re = new vector<double>*[runs];
	rho1FFTW_im = new vector<double>*[runs];
	rho2FFTW_re = new vector<double>*[runs];
	rho2FFTW_im = new vector<double>*[runs];
	
	for(int run=0; run<runs; run++){
		g11[run] = new vector<double>[obs_anzahl];
		g12[run] = new vector<double>[obs_anzahl];
		g22[run] = new vector<double>[obs_anzahl];
		ftrho1_re[run] = new vector<double>[obs_anzahl];
		ftrho1_im[run] = new vector<double>[obs_anzahl];
		ftrho2_re[run] = new vector<double>[obs_anzahl];
		ftrho2_im[run] = new vector<double>[obs_anzahl];
		rho1FFTW_re[run] = new vector<double>[obs_anzahl];
		rho1FFTW_im[run] = new vector<double>[obs_anzahl];
		rho2FFTW_re[run] = new vector<double>[obs_anzahl];
		rho2FFTW_im[run] = new vector<double>[obs_anzahl];
	}//for runs
	}

//if(maxThreads != 0)
//    omp_set_num_threads(maxThreads);

//#pragma omp parallel for
for(int run=0; run<runs; run++){
    
    RunZustand theRun;
    #pragma omp critical
    {
	
	
	theRun.ausgabe_jeansgroessen();
	//main_init();
	theRun.init(run);
	theRun.obs_point(0); //initial positions
	
	cout << "Run nr " << run << ": init complete after " << time(NULL)-program_started_time <<" seconds." << endl;
    }
        
	
	/*
	cout << "record obs-Punkt 1 von " << obs_anzahl << endl;
	if (auswerten_korrfunk) record_korrelationsfunktion11(run, 0);
	if (auswerten_korrfunk) record_korrelationsfunktion22(run, 0);
	if (auswerten_korrfunk_mixed) record_korrelationsfunktion12(run,0);
	if (auswerten_rhovonk) record_ftrho_unkorrigiert(run, 0);
	if (auswerten_rhoviaFFTW) record_rhoFFTW(run, 0);
	if (auswerten_animation && run==0) pos_schreiben(0.0, pos1, pos2);
	*/
        
        if(debugmode){
            theRun.obs_point(0);
            RunZustand::RunDynamik::TimestepInfo info;
            
            FILE* debugfile1 = fopen("mrd_data1.txt", "w");
            FILE* debugfile2 = fopen("mrd_data2.txt", "w");
            fprintf(debugfile1, "# Mittelwerte über 50 Zeitschritte. Früh = erste 50 Zeitschritte. \n");
            fprintf(debugfile2, "# Mittelwerte über 50 Zeitschritte. Spät = nach 1000 Zeitschritten \n");
            
            fprintf(debugfile1, "# Format: mrd TAB <F_WCA_max> TAB |F_WCA_highest| TAB dt \n\n");
            fprintf(debugfile2, "# Format: mrd TAB <F_WCA_max> TAB |F_WCA_highest| TAB dt \n\n");
            
            double mrd = 1.0;
            for(; mrd > 0.01; mrd *= 0.95){
                
                cout << "mrd=" << mrd << endl;
                
                double F_mittel = 0.0;
                double F_highest = 0.0;
                double dt_mittel = 0.0;
                theRun.init(0);
                
                // 50 Zeitschritte früh
                for(int count=0; count<50; count++){
                    theRun.dyn.zeitschritt_debug(mrd, info, true);
                    dt_mittel += info.dt/50.0;
                    double F = sqrt( info.maxWCAx*info.maxWCAx + info.maxWCAy*info.maxWCAy );
                    F_mittel += F/50.0;
                    if (F > F_highest){
                        F_highest = F;
                    }
                }//for 50
                
                fprintf(debugfile1, "%g \t %g \t %g \t %g \n", mrd, F_mittel, F_highest, dt_mittel);
                
                // einige Zeitschritte überspringen
                int num_skip = (int) (1.5/mrd);
                if (num_skip < 50) num_skip=50;
                
                for(int count=0; count<num_skip; count++)
                    theRun.dyn.zeitschritt(0.1);
                
                //50 Zeitschritte spät
                F_mittel=0.0;
                F_highest=0.0;
                dt_mittel=0.0;
                for(int count=0; count<50; count++){
                    theRun.dyn.zeitschritt_debug(mrd, info, true);
                    dt_mittel += info.dt/50.0;
                    double F = sqrt( info.maxWCAx*info.maxWCAx + info.maxWCAy*info.maxWCAy );
                    F_mittel += F/50.0;
                    if (F > F_highest){
                        F_highest = F;
                    }
                }//for 50
                fprintf(debugfile2, "%g \t %g \t %g \t %g \n", mrd, F_mittel, F_highest, dt_mittel);
                
            }//for mrd
            fclose(debugfile1);
            fclose(debugfile2);
            return 0;
            
        }
        
        if(krafttestmode){
            
            // nur Typ 1 benutzt: r2_git und r2_rel nichtmal initialisiert.
            // Teilchen mit Index 0 sitzt auf Position (0.5*densGrid_Breite, 0.5*densGrid_Breite), 
            // also genau in der Mitte einer Zelle. Teilchen mit Index 1 wird in der Nähe platziert 
            // (y-Koordinate größer, Abstand ist Schleifenvariable t), Kraft berechnen, Kraft lesen 
            // und in Datei schreiben.
            
            FILE* outfile = fopen("kapkraft.txt", "w");
            fprintf(outfile, "#Format: r Fkap Ftotal\n\n");
            
            for(double t=0.001; t<5.0; t+=0.01){
                const double origin = theRun.r1_git[0][0]*nachList_Breite + theRun.r1_rel[0][0];
                
                double y = origin + t;
                theRun.r1_git[1][1] = (int) (y/nachList_Breite);
                theRun.r1_rel[1][1] = y - theRun.r1_git[1][1]*nachList_Breite;
                
                theRun.dyn.berechne_kapkraefte();
                theRun.dyn.berechne_WCA11();
                
                double cap = theRun.dyn.getFKap();
                double tot = theRun.dyn.getF();
                //cout << "Kapkraft bei Abstand "<<t<<": " << theRun.dyn.getFKap() << endl;
                //cout << "Gesamtkraft bei Abstand "<<t<<": " << theRun.dyn.getF() << endl << endl;
                fprintf(outfile, "%g \t %g \t %g \n", t, cap, tot);
            }//for t
            return 0;
        }//if krafttest
        
	cout << "Run " << run << ": initial positions written, starting timesteps." << endl;
        time_t last_write = time(0);
	
	for(int obs_nr=1; obs_nr<obs_anzahl; obs_nr++){
		
		/*
		double t=0.0;
		int schritte_seit_obs = 0;
		//Zeitschritte bis t=obs_dt
		while(t < obs_dt){
			t += zeitschritt(obs_dt-t);
			schritte_im_run++;
			schritte_seit_obs++;
// 			cout << "t="<<t<<endl;
		}//while obs-Punkt noch nicht erreicht
		*/
		theRun.zeitschritte_bis_obs(last_write);
                if(obs_nr > 1){
                    int numSteps = theRun.schritte_seit_obs;
                    time_t now = time(0);
                    double steps_per_second = numSteps / (now-last_write);
                    
                    cout << "\nRun " << run << ": recording obs-Point " << theRun.obs_nr+1 << " / " << obs_anzahl  << ".\n";
                    cout << "Run " << run << ": " << numSteps << " steps since CP. Avg dt: " << obs_dt/numSteps << ". Performance: " << steps_per_second <<" steps/s.\n";
                    if (0 == run) cout << "Time since program launch: " << time(0) - program_started_time << " seconds." << endl;
                }
		/*
		if(auswerten_korrfunk) record_korrelationsfunktion11(run, obs_nr);
		if(auswerten_korrfunk) record_korrelationsfunktion22(run, obs_nr);
		if(auswerten_korrfunk_mixed) record_korrelationsfunktion12(run, obs_nr);
		if(auswerten_rhovonk) record_ftrho_unkorrigiert(run, obs_nr);
		if(auswerten_rhoviaFFTW) record_rhoFFTW(run, obs_nr);
		if(auswerten_animation && run==0) pos_schreiben(obs_nr*obs_dt, pos1, pos2);
		*/
		theRun.obs_point(obs_nr);
                last_write = time(0);
	}//for obs_nr
	/*
	if(auswerten_animation && run==0){
		fclose(pos1);
		fclose(pos2);
	}//if erster Run
	*/
	
	//theRun.obs.normalize(); //TODO: rhoFFTW mit Anzahl Beiträgen normieren
	theRun.schreibe_obs();
	//TODO zurückschreiben (Adresse übergeben)

}//for run

/*
//Statistik, und Ergebnis in Datei schreiben
if(auswerten_korrfunk) auswerten_korrelationsfunktion();
if(auswerten_korrfunk_mixed) auswerten_korrelationsfunktion_mixed(); //HIER WEITER! 
if(auswerten_rhovonk) auswerten_ftrho();
if(auswerten_rhoviaFFTW) auswerten_rhoFFTW();
if(auswerten_rhoFT_normjerun) auswerten_rhoFFTW_normjerun();

pos_schreiben_einedatei(); //schreibe Positionen am Ende des letzten Runs in Datei pos.txt
*/

//mittele alle Observablen, schreibe Mittelwerte und Fehler in Dateien
//alle_auswerten(); //TODO


return 0;

}//int main




