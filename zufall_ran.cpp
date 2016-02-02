// stellt Methoden für die Erzeugung von Zufallszahlen bereit

#pragma once

#include <cstdlib>
#include <math.h>
#include <time.h>

void init_rng(){
	srand(time(NULL));
}//void init_rng

//gleichverteilte Zufallszahl aus [0.0 , 1.0]
inline double zufall_gleichverteilt(void){
	const int max = RAND_MAX;

	return (double)rand()/max;	
}

//gleichverteilte Zufallszahl aus [min , max]
inline double zufall_gleichverteilt_vonbis(double min, double max){
	double ret = zufall_gleichverteilt(); 	//0.0 bis 1.0
	
	ret *= (max - min);			// 0.0 bis (max-min)

	return min + ret;
}//zufall_gleichverteilt_vonbis


//zwei (unabhängige) gaußverteilte Zufallszahlen mit Mittelwert Null und vorgegebener Varianz (Standardwert: Eins)
void zufall_gaussverteilt(double &z1, double &z2, double sigma=1.0){
	
	double x1 = zufall_gleichverteilt();
	double x2 = zufall_gleichverteilt();

	z1 = sqrt(-2.0 * log(x1)) * cos(2.0*M_PI*x2)*sigma;
	z2 = sqrt(-2.0 * log(x1)) * sin(2.0*M_PI*x2)*sigma;
	
	
	return;
	
}//void zufall_gaussverteilt




//zwei gaußverteile Zufallszahlen z1, z2 mit Varianzen sigma1, sigma2 und Korrelation c12
void zufall_gaussverteilt_korreliert(double &z1, double &z2, double sigma1, double sigma2, double c12){

	//erzeuge zwei unabhängige Zufallszahlen x1, x2, mit Varianzen 1.0 und Korrelation 0.0
	double x1, x2;
	
	zufall_gaussverteilt(x1, x2);
	
	z1 = sigma1 * x1;
	z2 = sigma2*(c12*x1 + sqrt(1-c12*c12)*x2);

	return;	
	
}//void zufall_gaussverteilt_korreliert
