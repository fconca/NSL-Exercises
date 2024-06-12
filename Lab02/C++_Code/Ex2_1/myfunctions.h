#ifndef MY_FUNCTIONS_H
#define MY_FUNCTIONS_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h" 

using namespace std;

struct pos {
    double x = 0.; 
    double y = 0.; 
    double z = 0.;
};

// Inizializza generatore numeri casuali
void InitializeRandomGenerator(const char* primesFile, const char* seedFile, Random& rnd);

// Media a blocchi dati vettori medie e medie quadratiche
vector<vector<double>> BlockAverage (const vector<double>& ave, const vector<double>& ave2);

// Realizza spostamento di una particella di lunghezza 'step' nella direzione e verso indicati
pos NewPositionLattice (const pos& part, int dir, int ver, int step);

// Realizza spostamento di una particella di lunghezza 'step' lungo la direzione indicata dagli angoli
pos NewPositionSpherical (const pos& part, double theta, double phi, double step);

#endif