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
	make obs
	echo "g(r) written to korrfunk_run1.txt. Rendering snapshots..."
	make snapshots

sim :  $(FILES)
	g++ $(FILES) $(FLAGS_IMMER) $(FLAGS_NORMAL) -o BDkap
	time ./BDkap > out.txt


obs : main_auswertung.cpp parameter.h tmp
	g++ main_auswertung.cpp -O3 -o auswertung_binary
	mkdir -p localDens
	time ./auswertung_binary
	mkdir -p rhok_iso
	time python3 calc_rhoIso.py

snapshots: ParProgress2tmp.py plot_animation.plt parameter.h tmp
	gnuplot plot_animation.plt 

	
auswertung: main_auswertung.cpp parameter.h ParProgress2tmp.py plot_animation.plt
	make obs
	echo "Calculation of g(r) complete. Rendering snapshots..."
	make snapshots


simobs: $(FILES) main_auswertung.cpp parameter.h
	make sim
	echo "Sim complete. Calculating g(r)..."
	make obs

tmp: ParProgress2tmp.py parameter.h
	python3 ParProgress2tmp.py


### These are old and probably not working anymore.
profile: $(FILES_PROFILE)
	g++ $(FILES_PROFILE) $(FLAGS_IMMER) $(FLAGS_NORMAL) -pg -o BDprofile

auswertung_debug: $(FILES_AUSWERTUNG)
	g++ $(FILES_AUSWERTUNG) $(FLAGS_IMMER) $(FLAGS_DEBUG) -o auswertung


animation_graphen:
	gnuplot plot_dichteprofile.txt

debug : $(FILES)
	g++ $(FILES) $(FLAGS_IMMER) $(FLAGS_DEBUG) -o debug