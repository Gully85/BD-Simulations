// #include "signaturen_auswertung.h"
#include <vector>
#include "parameter.h"
#include <iostream>
#include <cstdio>
#include <sstream>
#include <algorithm>

using std::cout; using std::endl; using std::flush; using std::vector; using std::min;

extern const double obs_dt;
extern const int obs_anzahl, korr_bins;
extern const int N1, N2;
extern const int jobs;
//Teilchenpositionen. Erster Index Teilchen (bis N1 bzw N2), zweiter Index 0/1 für x/y.
double** r1_abs = NULL;
double** r2_abs = NULL;

extern const int snapgrid_num;
extern const double snapgrid_width;

extern const int rhohist_bins;
extern const double rhohist_drho;

extern const double rho_thresh1, rho_thresh2, rho_thresh3, drho_thresh1, drho_thresh2, drho_thresh3, drho_thresh4;

// a single particle contributes this much to the density
const double drho = 1./(snapgrid_width*snapgrid_width);

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
        
        // all correlation functions are written to this file
        FILE* korrfile = fopen("korrfunk_run1.txt", "w");
        fprintf(korrfile, "# Format: t r g11 g12 g22 \n\n");
        
        // vector<double> to hold the g(r)
        vector<double> g11;
        vector<double> g12;
        vector<double> g22;
        
        //histogram over local density and local deltarho is written to these files
        FILE* rhoHistFile = fopen("rhoHist.txt", "w");
        fprintf(rhoHistFile, "# Format: t TAB dens TAB fraction of space with <rho1 TAB rho2 TAB rho1+rho2> = dens \n\n");
        FILE* drhoHistFile= fopen("drhoHist.txt","w");
        fprintf(drhoHistFile,"# Format: t TAB rho1-rho2 TAB fraction of space with this rho1-rho2\n\n");
        
        // time evolution of (#cells with rho>thresh) is written here
        FILE* rhoThreshFile = fopen("rhoThresh.txt", "w");
        fprintf(rhoThreshFile, "# Format: t TAB #cells over thresh1 TAB #cells over thresh2 TAB #cells over thresh3\n");
        fprintf(rhoThreshFile, "# where thresh1=%g, thresh2=%g, thresh3=%g\n\n", rho_thresh1, rho_thresh2, rho_thresh3);
        // time evolution of (#cells with |rho1-rho2|/(rho1+rho2)>thresh) is written here
        FILE* drhoThreshFile= fopen("drhoThresh.txt", "w");
        fprintf(drhoThreshFile,"# Format: t TAB #cells with |rho1-rho2|/(rho1+rho2) over thresh1 TAB same thresh2 TAB same thresh3 TAB same thresh4\n");
        fprintf(drhoThreshFile,"# where thresh1=%g, thresh2=%g, thresh3=%g, thresh4=%g\n\n", drho_thresh1, drho_thresh2, drho_thresh3, drho_thresh4);
        
        
        // vector<double> to hold the rho-histograms and drho-histogram
        int rhoHist[rhohist_bins];
        int rho1Hist[rhohist_bins];
        int rho2Hist[rhohist_bins];
        int drhoHist[rhohist_bins];
        // bin assignment in rhohist:  x = rhohist_drho*(bin+0.5)
        // bin assignment in drhohist: x = rhohist_drho*(bin+0.5-rhohist_bins/2)
        // solved for bin: in rhohist bin=(int)(x/rhohist_drho - 0.5)
        // and in drhohist:           bin=(int)(x/rhohist_drho - 0.5 + rhohist_bin/2)
        
        // array of doubles to hold the gridded density
        double** snapgrid1 = new double*[snapgrid_num];
        double** snapgrid2 = new double*[snapgrid_num];
        for(int i=0; i<snapgrid_num; i++){
            snapgrid1[i] = new double[snapgrid_num];
            snapgrid2[i] = new double[snapgrid_num];
        }//for i to snapgrid_num
        FILE* snapgridfile;
        std::stringstream ss;
        
	for(int t=0; t<obs_anzahl; t++){
            //cout << "reading t=" << t*obs_dt << endl;
            
            ///// READ ALL PARTICLE POSITIONS
            {
                for(int teilchen=0; teilchen<N1; teilchen++){
                    if(3 != fscanf(in1, "%g \t %g \t %g \n", &tmpt, &tmpx, &tmpy)){
                        cout << "Problem in pos1_1.txt: t=" <<t*obs_dt<< ", teilchen="<<teilchen << endl;
                        return 0;
                    }
                    if(0.9*obs_dt*t > tmpt || 1.1*obs_dt*t < tmpt){
                        cout << "Problem in pos1_1.txt: t=" <<t*obs_dt<< ", teilchen="<<teilchen<< ", time is "<<tmpt<< " instead of expected "<<t*obs_dt<<endl;
                        return 0;
                    }
                    
                    r1_abs[teilchen][0] = (double) tmpx;
                    r1_abs[teilchen][1] = (double) tmpy;
                }//for teilchen bis N1
                
                for(int teilchen=0; teilchen<N2; teilchen++){
                    if(3 != fscanf(in2, "%g \t %g \t %g \n", &tmpt, &tmpx, &tmpy)){
                        cout << "Problem in pos2_1.txt: t=" <<t*obs_dt<< ", teilchen="<<teilchen << endl;
                        return 0;
                    }
                    if(0.9*obs_dt*t > tmpt || 1.1*obs_dt*t < tmpt){
                        cout << "Problem in pos2_1.txt: t=" <<t*obs_dt<< ", teilchen="<<teilchen<< ", time is "<<tmpt<< " instead of expected "<<t*obs_dt<<endl;
                        return 0;
                    }
                    
                    r2_abs[teilchen][0] = (double) tmpx;
                    r2_abs[teilchen][1] = (double) tmpy;
                }//for teilchen bis N2
            }
            ///// READ ALL PARTICLE POSITIONS FINISHED
            
            
            ///// CALCULATE g(r) AND WRITE TO FILES
            {
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
                        if(0>bin || bin >= korr_bins) {
                            bin=korr_bins-1;
                            //cout << bin;
                        }
                        g11[bin] += vorfaktor11/bin;
                    }//for j, Typ 1
                    
                    for(int j=0; j<N2; j++){
                        double a2 = abstand2(r1_abs[i], r2_abs[j]);
                        
                        int bin = (int)(sqrt(a2)/korr_dr);
                        if(0>bin || bin >= korr_bins) bin=korr_bins-1;
                        g12[bin] += vorfaktor12/bin;
                    }//for j, Typ 2
                }//for i bis N1
                
                for(int i=0; i<N2; i++){
                    for(int j=0; j<i; j++){
                        double a2 = abstand2(r2_abs[i], r2_abs[j]);
                        
                        int bin = (int)(sqrt(a2)/korr_dr);
                        if(0>bin || bin>=korr_bins) bin = korr_bins-1;
                        g22[bin] += vorfaktor22/bin;
                    }//for j
                }//for i
                
                //output to file
                for(int bin=0; bin<korr_bins; bin++)
                    fprintf(korrfile, "%g \t %g \t %g \t %g \t %g \n", t*obs_dt, bin*korr_dr, g11[bin], g12[bin], g22[bin]);
                fprintf(korrfile, "\n");
            }
            ///// CALCULATE g(r) AND WRITE TO FILES FINISHED
            
            ///// CALCULATE GRIDDED DENSITIES AND WRITE TO FILES
            {
                int i; //particle index
                int k,l; //lattice cell indices
                int kl,kr,ll,lr; // neighboring lattice cell indices. "kl" is left to k
                double x,y; //particle coords
                double d; // distance particle to center-of-cell
                double xcenter,xleft,xright; // fraction of particle in center,left,right cell
                double ycenter,yleft,yright; // same for y-direction
                double tmp; //distance particle to center of neighboring cell
                
                /*
                 * //Zahl zu String, keine führenden Nullen
                    string int_to_string(int zahl){
                            std::stringstream ss;
                            ss << zahl;
                            //ss << setw(n) << setfill('0') << zahl;
                            return ss.str();
                    }//string int_to_string
                */
                
                ss.str("");
                ss << "localDens/rho" << t << ".txt";
                string filename = ss.str();
                
                snapgridfile = fopen(filename.c_str(), "w");
                fprintf(snapgridfile, "# Format: t x y rho1 rho2 \n");
                fprintf(snapgridfile, "# where x and y are the center of the cell. \n\n");
                
                // start with lattice of all zeroes
                for(k=0; k<snapgrid_num; k++)
                    for(l=0; l<snapgrid_num; l++){
                        snapgrid1[k][l] = 0.0;
                        snapgrid2[k][l] = 0.0;
                    }//for k,l
                
                // loop particles type 1
                for(i=0; i<N1; i++){
                    x = r1_abs[i][0];
                    y = r1_abs[i][1];
                    
                    // indices of cells containing a fraction of the particle. 
                    // Modulo snapgrid_num to handle periodic boundaries
                    k = (int) (x/snapgrid_width);
                    k = k % snapgrid_num;
                    kl= (k-1+snapgrid_num) % snapgrid_num;
                    kr= (k+1) % snapgrid_num;
                    l = (int) (y/snapgrid_width);
                    l = l % snapgrid_num;
                    ll= (l-1+snapgrid_num) % snapgrid_num;
                    lr= (l+1) % snapgrid_num;
                    
                    // x-direction, fraction of particle in the three cells
                    d = x-(k+0.5)*snapgrid_width;
                    xcenter = 0.75 - d*d/(snapgrid_width*snapgrid_width);
                    tmp = fabs((d+snapgrid_width)/snapgrid_width);
                    xleft  = 0.5*(1.5-tmp)*(1.5-tmp);
                    tmp = fabs((d-snapgrid_width)/snapgrid_width);
                    xright = 0.5*(1.5-tmp)*(1.5-tmp);
                    
                    // same for y-direction
                    d = y-(l+0.5)*snapgrid_width;
                    ycenter = 0.75 - d*d/(snapgrid_width*snapgrid_width);
                    tmp = fabs((d+snapgrid_width)/snapgrid_width);
                    yleft  = 0.5*(1.5-tmp)*(1.5-tmp);
                    tmp = fabs((d-snapgrid_width)/snapgrid_width);
                    yright = 0.5*(1.5-tmp)*(1.5-tmp);
                    
                    // collected all info. Increase the nine cells
                    snapgrid1[kl][ll] += xleft  *yleft   * drho;
                    snapgrid1[kl][l ] += xleft  *ycenter * drho;
                    snapgrid1[kl][lr] += xleft  *yright  * drho;
                    
                    snapgrid1[k ][ll] += xcenter*yleft   * drho;
                    snapgrid1[k ][l ] += xcenter*ycenter * drho;
                    snapgrid1[k ][lr] += xcenter*yright  * drho;
                    
                    snapgrid1[kr][ll] += xright *yleft   * drho;
                    snapgrid1[kr][l ] += xright *ycenter * drho;
                    snapgrid1[kr][lr] += xright *yright  * drho;
                    
                }//for i to N1
                
                //loop particles type 2
                for(i=0; i<N2; i++){
                    x = r2_abs[i][0];
                    y = r2_abs[i][1];
                    
                    // indices of cells containing a fraction of the particle. 
                    // Modulo snapgrid_num to handle periodic boundaries
                    k = (int) (x/snapgrid_width);
                    k = k % snapgrid_num;
                    kl= (k-1+snapgrid_num) % snapgrid_num;
                    kr= (k+1) % snapgrid_num;
                    l = (int) (y/snapgrid_width);
                    l = l % snapgrid_num;
                    ll= (l-1+snapgrid_num) % snapgrid_num;
                    lr= (l+1) % snapgrid_num;
                    
                    // x-direction, fraction of particle in the three cells
                    d = x-(k+0.5)*snapgrid_width;
                    xcenter = 0.75 - d*d/(snapgrid_width*snapgrid_width);
                    tmp = fabs((d+snapgrid_width)/snapgrid_width);
                    xleft  = 0.5*(1.5-tmp)*(1.5-tmp);
                    tmp = fabs((d-snapgrid_width)/snapgrid_width);
                    xright = 0.5*(1.5-tmp)*(1.5-tmp);
                    
                    // same for y-direction
                    d = y-(l+0.5)*snapgrid_width;
                    ycenter = 0.75 - d*d/(snapgrid_width*snapgrid_width);
                    tmp = fabs((d+snapgrid_width)/snapgrid_width);
                    yleft  = 0.5*(1.5-tmp)*(1.5-tmp);
                    tmp = fabs((d-snapgrid_width)/snapgrid_width);
                    yright = 0.5*(1.5-tmp)*(1.5-tmp);
                    
                    // collected all info. Increase the nine cells
                    snapgrid2[kl][ll] += xleft  *yleft   * drho;
                    snapgrid2[kl][l ] += xleft  *ycenter * drho;
                    snapgrid2[kl][lr] += xleft  *yright  * drho;
                    
                    snapgrid2[k ][ll] += xcenter*yleft   * drho;
                    snapgrid2[k ][l ] += xcenter*ycenter * drho;
                    snapgrid2[k ][lr] += xcenter*yright  * drho;
                    
                    snapgrid2[kr][ll] += xright *yleft   * drho;
                    snapgrid2[kr][l ] += xright *ycenter * drho;
                    snapgrid2[kr][lr] += xright *yright  * drho;
                    
                }//for i to N2
                
                //write to file
                for(int k=0; k<snapgrid_num; k++)
                    for(int l=0; l<snapgrid_num; l++){
                        fprintf(snapgridfile, "%g \t %g \t %g \t %g \t %g \n", t*obs_dt, (k+0.5)*snapgrid_width, (l+0.5)*snapgrid_width, snapgrid1[k][l], snapgrid2[k][l]);
                    }//for k,l
                fprintf(snapgridfile, "\n");
                //fflush(snapgridfile);
                fclose(snapgridfile);
            }
            ///// CALCULATE GRIDDED DENSITIES AND WRITE TO FILES FINISHED
            
            ///// CALCULATE HISTOGRAM OF LOCAL DENSITIES AND WRITE TO FILE
            ///// and #cells over rho-threshold and over drho-threshold
            {
                int i,k,l; //loop indices
                const double prefactor = 1./(snapgrid_num*snapgrid_num);
                double rho1,rho2, rho, drho;
                int bin;
                for(i=0; i<rhohist_bins; i++){
                    rhoHist [i]=0;
                    rho1Hist[i]=0;
                    rho2Hist[i]=0;
                    drhoHist[i]=0;
                }
                int rho_overthresh1 = 0;
                int rho_overthresh2 = 0;
                int rho_overthresh3 = 0;
                int drho_overthresh1= 0;
                int drho_overthresh2= 0;
                int drho_overthresh3= 0;
                int drho_overthresh4= 0;
                
                for(k=0; k<snapgrid_num; k++) for(l=0; l<snapgrid_num; l++){
                    rho1 = snapgrid1[k][l];
                    rho2 = snapgrid2[k][l];
                    rho = rho1+rho2;
                    drho= rho1-rho2;
                    
                    if(rho > rho_thresh1) rho_overthresh1++;
                    if(rho > rho_thresh2) rho_overthresh2++;
                    if(rho > rho_thresh3) rho_overthresh3++;
                    if(fabs(drho)/rho > drho_thresh1) drho_overthresh1++;
                    if(fabs(drho)/rho > drho_thresh2) drho_overthresh2++;
                    if(fabs(drho)/rho > drho_thresh3) drho_overthresh3++;
                    if(fabs(drho)/rho > drho_thresh4) drho_overthresh4++;
                    
                    bin = (int)(rho/rhohist_drho - 0.5);
                    if (0>bin) bin=0;
                    if (bin >= rhohist_bins) bin=rhohist_bins-1;
                    rhoHist[bin]++;
                    
                    bin = (int) (rho1/rhohist_drho- 0.5);
                    if (0>bin) bin=0;
                    if (bin >= rhohist_bins) bin=rhohist_bins-1;
                    rho1Hist[bin]++;
                    
                    bin = (int) (rho2/rhohist_drho- 0.5);
                    if (0>bin) bin=0;
                    if (bin >= rhohist_bins) bin=rhohist_bins-1;
                    rho2Hist[bin]++;
                    
                    bin = (int)(drho/rhohist_drho - 0.5 + rhohist_bins/2);
                    if (0>bin) bin=0;
                    if (bin >= rhohist_bins) bin=rhohist_bins-1;
                    drhoHist[bin]++;
                }//for k,l
                
                fprintf(rhoThreshFile, "%g \t %g \t %g \t %g \n",       t*obs_dt, rho_overthresh1*prefactor, rho_overthresh2*prefactor, rho_overthresh3*prefactor);
                fprintf(drhoThreshFile,"%g \t %g \t %g \t %g \t %g \n", t*obs_dt,drho_overthresh1*prefactor,drho_overthresh2*prefactor,drho_overthresh3*prefactor,drho_overthresh4*prefactor);
                
                for(bin=0; bin<rhohist_bins; bin++){
                    fprintf(rhoHistFile, "%g \t %g \t %g \t %g \t %g \n", t*obs_dt, rhohist_drho*(bin+0.5), prefactor*rho1Hist[bin], prefactor*rho2Hist[bin], prefactor*rhoHist[bin]);
                    fprintf(drhoHistFile, "%g \t %g \t %g \n", t*obs_dt, rhohist_drho*(bin+0.5-rhohist_bins/2), prefactor*drhoHist[bin]);
                }//for bin
                fprintf(rhoHistFile, "\n");
                fprintf(drhoHistFile, "\n");
            }
            ///// CALCULATE HISTOGRAM OF LOCAL DENSITIES AND WRITE TO FILE FINISHED
            
        }//for t bis obs_anzahl
        
        fclose(korrfile);
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