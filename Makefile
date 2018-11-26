# Variablen

FLAGS_IMMER = -I$$FFTW_INC_DIR -L$$FFTW_LIB_DIR -lfftw3 -lm -std=c++0x -fopenmp 
FLAGS_DEBUG = -O0 -g
FLAGS_NORMAL = -O3

FILES = main.cpp auswertung.cpp dynamik_methoden.cpp gridRoutinen.cpp allgemeine_methoden.cpp zufall.cpp 
FILES_AUSWERTUNG = main_auswertung.cpp zufall.cpp observablen_auswertung.cpp positionen_laden_auswertung.cpp
FILES_PROFILE = main_profile.cpp auswertung.cpp dynamik_methoden.cpp gridRoutinen.cpp allgemeine_methoden.cpp zufall.cpp 

BD :  $(FILES)
	g++ $(FILES) $(FLAGS_IMMER) $(FLAGS_NORMAL) -o BDkap

debug : $(FILES)
	g++ $(FILES) $(FLAGS_IMMER) $(FLAGS_DEBUG) -o debug

auswertung: $(FILES_AUSWERTUNG)
	g++ $(FILES_AUSWERTUNG) $(FLAGS_IMMER) $(FLAGS_NORMAL) -g -o auswertung

profile: $(FILES_PROFILE)
	g++ $(FILES_PROFILE) $(FLAGS_IMMER) $(FLAGS_NORMAL) -pg -o BDprofile

auswertung_debug: $(FILES_AUSWERTUNG)
	g++ $(FILES_AUSWERTUNG) $(FLAGS_IMMER) $(FLAGS_DEBUG) -o auswertung


animation_graphen:
	gnuplot plot_dichteprofile.txt
