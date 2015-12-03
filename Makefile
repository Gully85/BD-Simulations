# Variablen

FLAGS_IMMER = -lfftw3 -lm -std=c++11
FLAGS_DEBUG = -O0 -g
FLAGS_NORMAL = -O3

FILES = main.cpp gridRoutinen.cpp dynamik_methoden.cpp zufall.cpp auswertung.cpp signaturen.h parameter.h

BD :  $(FILES)
	g++ main.cpp $(FLAGS_IMMER) $(FLAGS_NORMAL)

debug : $(FILES)
	g++ main.cpp $(FLAGS_IMMER) $(FLAGS_DEBUG)