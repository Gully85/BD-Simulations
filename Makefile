# Variablen

FLAGS_IMMER = -I$$FFTW_INC_DIR -L$$FFTW_LIB_DIR -lfftw3 -lm -std=c++0x -fopenmp 
FLAGS_DEBUG = -O0 -g
FLAGS_NORMAL = -O3

FILES = main.cpp auswertung.cpp dynamik_methoden.cpp gridRoutinen.cpp allgemeine_methoden.cpp zufall.cpp 
FILES_AUSWERTUNG = main_auswertung.cpp parameter.h zufall.cpp signaturen_auswertung.h positionen_speichernladen.cpp observablen_auswertung.cpp

BD :  $(FILES)
	g++ $(FILES) $(FLAGS_IMMER) $(FLAGS_NORMAL) -o BDkap

debug : $(FILES)
	g++ $(FILES) $(FLAGS_IMMER) $(FLAGS_DEBUG) -o debug

auswertung: $(FILES_AUSWERTUNG)
	g++ $(FILES_AUSWERTUNG) $(FLAGS_IMMER) $(FLAGS_NORMAL)


auswertung_debug: $(FILES_AUSWERTUNG)
	g++ $(FILES_AUSWERTUNG) $(FLAGS_IMMER) $(FLAGS_DEBUG)


animation_graphen:
	gnuplot plot_dichteprofile.txt
