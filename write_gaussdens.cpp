#include "parameter.h"

#include <cstdlib>
#include <cstdio>
#include <iostream>

// importiere Variablen aus parameter.h, siehe dort, was die Variable tut.
extern const int N1, N2;
extern const double L, obs_dt;

extern const int snapgrid_num;
extern const double snapgrid_width;

extern const double gaussdens_width_x, gaussdens_width_y;

void zufall_gaussverteilt(double&, double&, double sigma=1.0);

using std::cout; using std::endl;

int main(){
    
    FILE* out1 = fopen("pos1_1.txt", "w");
    FILE* out2 = fopen("pos2_1.txt", "w");
    
    fprintf(out1, "# Format: t TAB x TAB y\n\n");
    fprintf(out2, "# Format: t TAB x TAB y\n\n");
    
    // Species 1, to output out1
    for(int i=0; i<N1; i++){
        double x=0.0; double y=0.0;
        
        // randomly draw two coords. Move them into the box.
        zufall_gaussverteilt(x,y, gaussdens_width_x);
        x+=L/2.;
        y*=gaussdens_width_y/gaussdens_width_x;
        y+=L/2.;
        while(x<0) x += L;
        while(x>L) x -= L;
        while(y<0) y += L;
        while(y>L) y -= L;
        
        fprintf(out1, "%g \t %g \t %g \n", obs_dt, x, y);
    }//for i bis N1
    fprintf(out1, "\n");
    for(int i=0; i<N1; i++){
        double x=0.0; double y=0.0;
        
        // randomly draw two coords. Move them into the box.
        zufall_gaussverteilt(x,y, gaussdens_width_x);
        x+=L/2.;
        y*=gaussdens_width_y/gaussdens_width_x;
        y+=L/2.;
        while(x<0) x += L;
        while(x>L) x -= L;
        while(y<0) y += L;
        while(y>L) y -= L;
        
        fprintf(out1, "%g \t %g \t %g \n", 2*obs_dt, x, y);
    }//for i bis N1
    fclose(out1);
    
    // Species 2, to output out2
    for(int i=0; i<N2; i++){
        double x=0.0; double y=0.0;
        
        // randomly draw two coords. Move them into the box.
        zufall_gaussverteilt(x,y, gaussdens_width_x);
        x+=L/2.;
        y*=gaussdens_width_y/gaussdens_width_x;
        y+=L/2.;
        while(x<0) x += L;
        while(x>L) x -= L;
        while(y<0) y += L;
        while(y>L) y -= L;
        
        fprintf(out2, "%g \t %g \t %g \n", obs_dt, x, y);
    }//for i bis N2
    fprintf(out2, "\n");
    for(int i=0; i<N2; i++){
        double x=0.0; double y=0.0;
        
        // randomly draw two coords. Move them into the box.
        zufall_gaussverteilt(x,y, gaussdens_width_x);
        x+=L/2.;
        y*=gaussdens_width_y/gaussdens_width_x;
        y+=L/2.;
        while(x<0) x += L;
        while(x>L) x -= L;
        while(y<0) y += L;
        while(y>L) y -= L;
        
        fprintf(out2, "%g \t %g \t %g \n", 2*obs_dt, x, y);
    }//for i bis N2
    fclose(out2);
    
    
    // write a minimal dummy out.txt
    FILE* out = fopen("out.txt", "w");
    fprintf(out, "Run 0: recording obs-Point 2 / 5000.\n");
    fclose(out);
    
    
    // write noiseless gridded local density to file localDens/test_noiseless.txt
    out = fopen("localDens/test_noiseless.txt", "w");
    fprintf(out, "# Format: 0 x y rho1 rho2 \n# \n\n");
    
    cout << "gauss width: " << gaussdens_width_x << "," << gaussdens_width_y << endl;
    // variables double snapgrid_width and int snapgrid_num from parameter.h control the output grid
    for(int k=0; k<snapgrid_num; k++){
        double x = (k+0.5)*snapgrid_width;
        for(int l=0; l<snapgrid_num; l++){
            double y = (l+0.5)*snapgrid_width;
            // possible error source: Is the norm sigma_x*sigma*y or sigma_x**2+sigma_y**2 ?
            //double gaussnorm = (N1+N2)/(snapgrid_width*snapgrid_width) / sqrt(2. * M_PI * gaussdens_width_x*gaussdens_width_y);
            double gaussnorm = 1.0; //normalized to constant peak height. 
            double gaussvalue = gaussnorm * exp(-0.5*(
                                    (x-L/2.)*(x-L/2.)/(gaussdens_width_x*gaussdens_width_x) +
                                    (y-L/2.)*(y-L/2.)/(gaussdens_width_y*gaussdens_width_y)
                                                      ));
            
            fprintf(out, "0 \t %g \t %g \t %g \t %g \n", x, y, gaussvalue, gaussvalue);
        }
        }//for k,l
    fclose(out);
    return 0;
}//int main