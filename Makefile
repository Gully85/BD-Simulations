# Variablen

FLAGS_IMMER = -I$$FFTW_INC_DIR -L$$FFTW_LIB_DIR -lfftw3 -lm -std=c++0x -fopenmp 
FLAGS_DEBUG = -O0 -g
FLAGS_NORMAL = -O3

FILES = main.cpp auswertung.cpp dynamik_methoden.cpp gridRoutinen.cpp allgemeine_methoden.cpp zufall.cpp 
# FILES_AUSWERTUNG = main_auswertung.cpp zufall.cpp observablen_auswertung.cpp positionen_laden_auswertung.cpp
FILES_PROFILE = main_profile.cpp auswertung.cpp dynamik_methoden.cpp gridRoutinen.cpp allgemeine_methoden.cpp zufall.cpp 

fullsim: $(FILES) main_auswertung.cpp parameter.h ParProgress2tmp.py plot_animation.plt
	make sim
	echo "Sim complete. Calculating g(r)..."
	make corrfunc
	echo "g(r) written to korrfunk_run1.txt. Rendering snapshots..."
	make snapshots

sim :  $(FILES)
	g++ $(FILES) $(FLAGS_IMMER) $(FLAGS_NORMAL) -o BDkap
	./BDkap

simcorr: $(FILES) main_auswertung.cpp parameter.h
	make sim
	echo "Sim complete. Calculating g(r)..."
	make corrfunc

debug : $(FILES)
	g++ $(FILES) $(FLAGS_IMMER) $(FLAGS_DEBUG) -o debug

corrfunc : main_auswertung.cpp parameter.h 
	g++ main_auswertung.cpp -O3 -o auswertung_binary
	./auswertung_binary

snapshots: ParProgress2tmp.py plot_animation.plt parameter.h 
	python3 ParProgress2tmp.py
	gnuplot plot_animation.plt 

auswertung: main_auswertung.cpp parameter.h ParProgress2tmp.py plot_animation.plt
	make corrfunc
	echo "Calculation of g(r) complete. Rendering snapshots..."
	make snapshots

profile: $(FILES_PROFILE)
	g++ $(FILES_PROFILE) $(FLAGS_IMMER) $(FLAGS_NORMAL) -pg -o BDprofile

auswertung_debug: $(FILES_AUSWERTUNG)
	g++ $(FILES_AUSWERTUNG) $(FLAGS_IMMER) $(FLAGS_DEBUG) -o auswertung


animation_graphen:
	gnuplot plot_dichteprofile.txt
