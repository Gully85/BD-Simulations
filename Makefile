# Variablen
FLAGS_IMMER = -I$$FFTW_INC_DIR -L$$FFTW_LIB_DIR -lfftw3 -lm -std=c++0x -fopenmp -Wno-div-by-zero
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
	./auswertung_binary
	#mkdir -p rhok_iso
	#python3 calc_rhoIso.py

snapshots: plot_animation.plt tmp
	gnuplot plot_animation.plt 

	
auswertung: main_auswertung.cpp parameter.h ParProgress2tmp.py plot_animation.plt rhokt_multifit.py
	make obs
	echo "Calculation of g(r) complete. Rendering snapshots..."
	make snapshots
	echo "Rendering of snapshots complete. Fitting exponential to all rho(k)"
	mkdir -p rhok_iso
	python3 calc_rhot_multik.py


simobs: $(FILES) main_auswertung.cpp parameter.h
	make sim
	echo "Sim complete. Calculating g(r)..."
	make obs

tmp: ParProgress2tmp.py parameter.h
	python3 ParProgress2tmp.py

# Test of the procedure that calculates rho(k), without time dependence. Tests are a perfect gaussian function, and
# particle coordinates that follow a gaussian distribution. In both cases, x-width and y-width of the gaussian can be different, 
# controlled by gaussdens_width_x and _y in parameter.h
# write_gaussdens.cpp will write 
# - a perfect gaussian density profile (with two diferent widths in x- and y-direction) to file localDens/test_noiseless.txt
# - particle coordinates, same probability distribution, centered at (L/2,L/2) to files pos1_1.txt and pos2_1.txt
# make obs will then apply gridding to the positions, writing results to file localDens/rho1.txt
# Result of noiseless test in rhok_testresult_noiseless.png,
# result of test with particle positions in rhok_testresult.png
# For this, the variable obs_anzahl in parameter.h should be 2
rhoktest: parameter.h write_gaussdens.cpp zufall.cpp tmp plot_FTcheck.plt FTcheck_noiseless.plt
	g++ write_gaussdens.cpp zufall.cpp -O3 -o write_gaussdens -std=c++0x
	./write_gaussdens
	make obs
	gnuplot plot_FTcheck.plt
	gnuplot FTcheck_noiseless.plt
	
# Test of the procedure that calculates the time evolution of rho(k) for all k.
# Set obs_anzahl in parameter.h to around 10..15 
rhokttest: parameter.h $(FILES) main_auswertung.cpp calc_rhot_multik.py
	make sim
	make obs
	mkdir -p rhok_timeEv
	python3 rhokt_fit_test.py
	


### These are old and probably not working anymore.

auswertung_debug: $(FILES_AUSWERTUNG)
	g++ $(FILES_AUSWERTUNG) $(FLAGS_IMMER) $(FLAGS_DEBUG) -o auswertung


animation_graphen:
	gnuplot plot_dichteprofile.txt

debug : $(FILES)
	g++ $(FILES) $(FLAGS_IMMER) $(FLAGS_DEBUG) -o debug

