# Variablen

FLAGS_IMMER = -I$$FFTW_INC_DIR -L$$FFTW_LIB_DIR -lfftw3 -lm -std=c++0x
FLAGS_DEBUG = -O0 -g
FLAGS_NORMAL = -O3

FILES = main.cpp gridRoutinen.cpp dynamik_methoden.cpp zufall.cpp auswertung.cpp signaturen.h parameter.h positionen_speichernladen.cpp

BD :  $(FILES)
	g++ main.cpp $(FLAGS_IMMER) $(FLAGS_NORMAL)

debug : $(FILES)
	g++ main.cpp $(FLAGS_IMMER) $(FLAGS_DEBUG)