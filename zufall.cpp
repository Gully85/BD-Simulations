// stellt Methoden für die Erzeugung von Zufallszahlen bereit

// Gaußverteilte Zufallszahlen sind nicht mehr drin. Können zB in ~/Diplomarbeit/weichcos/Phi0_2/zufall.cpp nachgelesen werden


#include <cstdlib>
#include <math.h>
#include <random>
// #include <chrono>
#include <functional>

//typedef std::chrono::high_resolution_clock myclock;
//myclock::time_point beginning = myclock::now();

using std::mt19937_64;  //generator
using std::uniform_real_distribution; 

mt19937_64 generator(time(NULL));
// mt19937_64 generator;

uniform_real_distribution<double> dist(0.0,1.0);

//gleichverteilte Zufallszahl aus [0.0 , 1.0]
auto zufall_gleichverteilt = std::bind(dist, generator);
//Typ ist: Methode (void) -> double


void init_rng(){
	
// 	generator.seed( time(NULL) );
}



double zufall_gleichverteilt_vonbis(double min, double max){
	
    return min + (max-min)*zufall_gleichverteilt();
	
}//double zufall_gleichverteilt_vonbis


//Zufallszahl aus [0,N-1] ganzzahlig
int zufall_Int(int N){
    return ((int) (N*zufall_gleichverteilt()))%N;
}//int zufall_Int

//zwei (unabhängige) gaußverteilte Zufallszahlen mit Mittelwert Null und vorgegebener Varianz (Standardwert: Eins)
void zufall_gaussverteilt(double &z1, double &z2, double sigma=1.0){
    
    double x1 = zufall_gleichverteilt();
    double x2 = zufall_gleichverteilt();
    
    z1 = sqrt(-2.0 * log(x1)) * cos(2.0*M_PI*x2)*sigma;
    z2 = sqrt(-2.0 * log(x1)) * sin(2.0*M_PI*x2)*sigma;
    
    return;
	
}//void zufall_gaussverteilt
