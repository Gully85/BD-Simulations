#include <cstdio>
#include <iostream>
#include <sstream>

#include "parameter.h"
#include "signaturen_auswertung.h"

extern const int jobs;

extern const int N1, N2;

extern double**** r1_abs;
extern double**** r2_abs;

using std::string; using std::cout; using std::endl;

//lese Teilchenpositionen (mit zugehöriger Zeit) aus Dateien. Rückgabewert true bei Erfolg, false bei fehlender Datei oder falschem Dateiinhalt
bool rvont_lesen(){
	
	FILE* in1; //Datei, in der die Positionen der Typ1-Teilchen stehen
	FILE* in2; //Datei mit Positionen Typ2
	
	for(int job=0; job<jobs; job++){
		
		string datei1 = "pos1_.txt";
		string datei2 = "pos2_.txt";
		
		string job_str = int_to_string(job+1);
		
		datei1.insert(5, job_str);
		datei2.insert(5, job_str);
		
		if(!(in1 = fopen(datei1.c_str(), "r"))){
			cout << "Fehler: Datei " << datei1 << " nicht gefunden!" << endl;
			return false;
		}
		if(!(in2 = fopen(datei2.c_str(), "r"))){
			cout << "Fehler: Datei " << datei2 << " nicht gefunden!" << endl;
			return false;
		}
		
		// skip first two lines. One is comment, one is empty
		char buf[100];
                if(NULL == fgets(buf,100,in1)) cout << "error skipping empty line"<<endl;
                //if(NULL == fgets(buf,100,in1)) cout << "error skipping empty line"<<endl;
                //if(NULL == fgets(buf,100,in2)) cout << "error skipping empty line"<<endl;
                if(NULL == fgets(buf,100,in2)) cout << "error skipping empty line"<<endl;
		
		//zeit ist in Einheiten obs_dt, geht von Null bis obs_anzahl. In den Dateien stehen obs_anzahl viele Blöcke, und jeder Block hat N (Teilchenzahl des Typs) Einträge, und jeder Eintrag ist eine Zeile "Zeit TAB x TAB y \n"
		for(int zeit=1; zeit<obs_anzahl; zeit++){
			
			float tmpt, tmpx, tmpy;
		//Typ 1
			for(int teilchen=0; teilchen<N1; teilchen++){
				if(3 != fscanf(in1, "%g \t %g \t %g \n", &tmpt, &tmpx, &tmpy)){
					cout << "Fehler in Datei " << datei1 << ": t="<<zeit*obs_dt<<", teilchen="<<teilchen << endl;
					return false;
				}
				if(0.9*obs_dt*zeit > tmpt || 1.1*obs_dt*zeit < tmpt){
					cout << "Fehler in Datei " << datei1 << ": t="<<zeit*obs_dt<<", teilchen="<<teilchen <<", Zeit nicht wie erwartet "<<zeit*obs_dt<<", sondern " << tmpt << endl;
					return false;
				}
				
				//in r1 schreiben
				double x = (double) tmpx;
				double y = (double) tmpy;
				
				r1_abs[job][zeit][teilchen][0] = x;
				r1_abs[job][zeit][teilchen][1] = y;
			}//for teilchen bis N1
			//Leerzeile in Inputdatei!
			
		//Typ 2 analog
			for(int teilchen=0; teilchen<N2; teilchen++){
				if(3 != fscanf(in2, "%g \t %g \t %g \n", &tmpt, &tmpx, &tmpy)){
					cout << "Fehler in Datei " <<datei2<<": t="<<zeit*obs_dt<< ", teilchen="<<teilchen << endl;
					return false;
				}
				if(0.9*zeit*obs_dt > tmpt || 1.1*zeit*obs_dt < tmpt){
					cout << "Fehler in Datei "<<datei2<<": t="<<zeit*obs_dt<<", teilchen="<<teilchen<<", Zeit nicht wie erwartet "<<zeit*obs_dt<<", sondern "<<tmpt << endl;
					return false;
				}
				
				double x = (double) tmpx;
				double y = (double) tmpy;
				
				r2_abs[job][zeit][teilchen][0] = x;
				r2_abs[job][zeit][teilchen][1] = y;
			}//for teilchen bis N2
			//Leerzeile in Inputdatei!
		}//for zeit bis obs_anzahl
		
		fclose(in1);
		fclose(in2);
	}//for job
	
}//void rvont_lesen






//TODO fixe Länge mit führenden Nullen, wie in ~/Diplomarbeit/weichcos/ausgaben_methoden.cpp
//Zahl zu String, keine führenden Nullen
string int_to_string(int zahl){
	std::stringstream ss;
	ss << zahl;
	//ss << setw(n) << setfill('0') << zahl;
	return ss.str();
}//string int_to_string