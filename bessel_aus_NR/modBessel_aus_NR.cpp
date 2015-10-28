#include "nr.h"
#include "bessk0.cpp"
#include "bessk1.cpp"
#include "bessi0.cpp"
#include "bessi1.cpp"

#include <cstdio>
#include <math.h>


int main(){

double dx = 0.01;

FILE* out0 = fopen("K0.txt", "w");
FILE* out1 = fopen("K1.txt", "w");

double xmax = 10.0;


double x, tmp;

for(x=0.0; x<xmax; x += dx){

	tmp = NR::bessk0(x);
	fprintf(out0, "%g \t %g \n", x, tmp);

	tmp = NR::bessk1(x);
	fprintf(out1, "%g \t %g \n", x, tmp);
}//for x

fclose(out0);
fclose(out1);


return 0;}
