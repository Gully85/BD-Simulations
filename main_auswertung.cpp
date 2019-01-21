// #include "signaturen_auswertung.h"
#include <vector>
#include "parameter.h"
#include <iostream>
#include <cstdio>
#include <algorithm>

using std::cout; using std::endl; using std::flush; using std::vector; using std::min;

extern const double obs_dt;
extern const int obs_anzahl;
extern const int N1, N2;
extern const int jobs;
//Teilchenpositionen. Erster Index Teilchen (bis N1 bzw N2), zweiter Index 0/1 für x/y.
double** r1_abs = NULL;
double** r2_abs = NULL;

double abstand2(double*, double*);

int main(){
	
	//reserviere Speicher
	r1_abs = new double*[N1];
	if(0 != N2) r2_abs = new double*[N2];
        
        for(int teilchen=0; teilchen<N1; teilchen++)
            r1_abs[teilchen] = new double[2];
        for(int teilchen=0; teilchen<N2; teilchen++)
            r2_abs[teilchen] = new double[2];
	
        FILE* in1 = fopen("pos1_1.txt", "r");
        FILE* in2 = fopen("pos2_1.txt", "r");
        
        
        // skip first line in both files
        char buf[100];
        if(NULL == fgets(buf, 100, in1)) cout << "error skipping empty line" << endl;
        if(NULL == fgets(buf, 100, in2)) cout << "error skipping empty line" << endl;
        
        float tmpt, tmpx, tmpy;
        
        // all correlation functions are written into this file
        FILE* out = fopen("korrfunk_run1.txt", "w");
        fprintf(out, "# Format: t r g11 g12 g22 \n\n");
        
        // vector<double> to hold the histogram.
        vector<double> g11;
        vector<double> g12;
        vector<double> g22;
        
        
        
	for(int t=1; t<obs_anzahl; t++){
            //cout << "reading t=" << t*obs_dt << endl;
            
            
            for(int teilchen=0; teilchen<N1; teilchen++){
                if(3 != fscanf(in1, "%g \t %g \t %g \n", &tmpt, &tmpx, &tmpy)){
                    cout << "Problem in pos1_1.txt: t=" <<t*obs_dt<< ", teilchen="<<teilchen << endl;
                    return 1;
                }
                if(0.9*obs_dt*t > tmpt || 1.1*obs_dt*t < tmpt){
                    cout << "Problem in pos1_1.txt: t=" <<t*obs_dt<< ", teilchen="<<teilchen<< ", time is "<<tmpt<< " instead of expected "<<t*obs_dt<<endl;
                    return 1;
                }
                
                r1_abs[teilchen][0] = (double) tmpx;
                r1_abs[teilchen][1] = (double) tmpy;
            }//for teilchen bis N1
            
            for(int teilchen=0; teilchen<N2; teilchen++){
                if(3 != fscanf(in2, "%g \t %g \t %g \n", &tmpt, &tmpx, &tmpy)){
                    cout << "Problem in pos2_1.txt: t=" <<t*obs_dt<< ", teilchen="<<teilchen << endl;
                    return 1;
                }
                if(0.9*obs_dt*t > tmpt || 1.1*obs_dt*t < tmpt){
                    cout << "Problem in pos2_1.txt: t=" <<t*obs_dt<< ", teilchen="<<teilchen<< ", time is "<<tmpt<< " instead of expected "<<t*obs_dt<<endl;
                    return 1;
                }
                
                r2_abs[teilchen][0] = (double) tmpx;
                r2_abs[teilchen][1] = (double) tmpy;
            }//for teilchen bis N1
            
            // g11
            g11.assign(korr_bins, 0.0);
            g12.assign(korr_bins, 0.0);
            g22.assign(korr_bins, 0.0);
            
            const double vorfaktor11 = L*L/(N1*(N1-1)*M_PI*korr_dr*korr_dr);
            const double vorfaktor22 = L*L/(N2*(N2-1)*M_PI*korr_dr*korr_dr);
            const double vorfaktor12 = L*L/(N1*N2  *2*M_PI*korr_dr*korr_dr);
            
            for(int i=0; i<N1; i++){
                for(int j=0; j<i; j++){
                    double a2 = abstand2(r1_abs[i], r1_abs[j]);
                    
                    int bin = (int)(sqrt(a2)/korr_dr);
                    if(0>bin || bin > korr_bins) bin=korr_bins-1;
                    g11[bin] += vorfaktor11/bin;
                }//for j, Typ 1
                
                for(int j=0; j<N2; j++){
                    double a2 = abstand2(r1_abs[i], r2_abs[j]);
                    
                    int bin = (int)(sqrt(a2)/korr_dr);
                    if(0>bin || bin > korr_bins) bin=korr_bins-1;
                    g12[bin] += vorfaktor12/bin;
                }//for j, Typ 2
            }//for i bis N1
            
            for(int i=0; i<N2; i++){
                for(int j=0; j<i; j++){
                    double a2 = abstand2(r2_abs[i], r2_abs[j]);
                    
                    int bin = (int)(sqrt(a2)/korr_dr);
                    if(0>bin || bin>korr_bins) bin = korr_bins-1;
                    g22[bin] += vorfaktor22/bin;
                }//for j
            }//for i
            
            //output to file
            for(int bin=0; bin<korr_bins; bin++)
                fprintf(out, "%g \t %g \t %g \t %g \t %g \n", t*obs_dt, bin*korr_dr, g11[bin], g12[bin], g22[bin]);
            fprintf(out, "\n");
        }//for t bis obs_anzahl
        
        fclose(out);
        fclose(in1);
        fclose(in2);
        
        return 0;
        
}//int main

//berechnet das Quadrat des Abstands zwischen Teilchen i und Teilchen j. Berücksichtigt periodische Randbedingungen.
double abstand2(double* ri, double* rj){
	
	//Absolutpositionen
	double xi_abs = ri[0];
	double xj_abs = rj[0];
	double yi_abs = ri[1];
	double yj_abs = rj[1];
	
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