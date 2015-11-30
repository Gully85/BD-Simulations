// stellt Methoden für die Erzeugung von Zufallszahlen bereit

// Gaußverteilte Zufallszahlen sind nicht mehr drin. Können zB in ~/Diplomarbeit/weichcos/Phi0_2/zufall.cpp nachgelesen werden

#pragma once

#include <cstdlib>
#include <math.h>
#include <random>
#include <chrono>
#include <functional>

//typedef std::chrono::high_resolution_clock myclock;
//myclock::time_point beginning = myclock::now();

using std::mt19937_64;  //generator
using std::uniform_real_distribution; 

mt19937_64 generator;
uniform_real_distribution<double> dist(0.0,1.0);

void init_random(){
	generator.seed( time(NULL) );
}

//gleichverteilte Zufallszahl aus [0.0 , 1.0]
auto zufall_gleichverteilt = std::bind(dist, generator);
//Typ ist allerdings: Methode (void) -> double


double zufall_gleichverteilt_vonbis(double min, double max){
	
	return min + (max-min)*zufall_gleichverteilt();
	
}//double zufall_gleichverteilt_vonbis


//Zufallszahl aus [0,N-1] ganzzahlig
int zufall_Int(int N){
	return ((int) (N*zufall_gleichverteilt()))%N;
}//int zufall_Int





