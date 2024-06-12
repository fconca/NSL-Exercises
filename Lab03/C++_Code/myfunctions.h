#ifndef MY_FUNCTIONS_H
#define MY_FUNCTIONS_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h" 

using namespace std;

// Inizializza generatore numeri casuali
void InitializeRandomGenerator(const char* primesFile, const char* seedFile, Random& rnd);

// Media a blocchi dati vettori medie e medie quadratiche
vector<vector<double>> BlockAverage (const vector<double>& ave, const vector<double>& ave2);

// Calcola massimo tra due numeri
double max (const double& a, const double& b);

#endif