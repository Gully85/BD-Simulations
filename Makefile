# Variablen

FLAGS_IMMER = -I$$FFTW_INC_DIR -L$$FFTW_LIB_DIR -lfftw3 -fopenmp -lm -std=c++0x
FLAGS_DEBUG = -O0 -g
FLAGS_NORMAL = -O3

FILES = main.cpp auswertung.cpp dynamik_methoden.cpp gridRoutinen.cpp allgemeine_methoden.cpp zufall.cpp
FILES_AUSWERTUNG = main_auswertung.cpp parameter.h zufall.cpp signaturen_auswertung.h positionen_speichernladen.cpp observablen_auswertung.cpp

BD :  $(FILES)
	g++ $(FILES) $(FLAGS_IMMER) $(FLAGS_NORMAL)

debug : $(FILES)
	g++ $(FILES) $(FLAGS_IMMER) $(FLAGS_DEBUG)

auswertung_debug: $(FILES_AUSWERTUNG)
	g++ main_auswertung.cpp positionen_laden_auswertung.cpp observablen_auswertung.cpp -O0 -g

auswertung: $(FILES_AUSWERTUNG)
	g++ main_auswertung.cpp positionen_laden_auswertung.cpp observablen_auswertung.cpp -O3

animation_graphen:
	gnuplot plot_dichteprofile.txt