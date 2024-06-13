#ifndef MY_FUNCTIONS_H
#define MY_FUNCTIONS_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
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

// Legge parametri di input (tipo simulazione e generazione random, step Metropolis, punto partenza, numero e dimensione blocchi)
void ReadParameters(string& simulation_type, string& random_type, double& dim_step, int& num_step_acc, double& xcoord, double &ycoord, double& zcoord, int& num_block, int& dim_block, bool& print);

// Calcola distanza di un punto dall'origine
double dist (const pos& point);

// Calcola probabilità di una particella di trovarsi in un certo punto per groundstate e primo stato eccitato
double ProbGroundState (const pos& point);
double ProbFirstExcState (const pos& point);

// Determina minimo fra due punti
double min (double a, double b);

// Genera un nuovo punto a partire da un punto assegnato secondo distrib uniforme o gaussiana centrata in quest'ultimo
pos GetNewPointUnif (Random& rnd, const pos& point, double dim_step);
pos GetNewPointNorm (Random& rnd, const pos& point, double dim_step);

// Calcola parametro di accettazione per groundstate e primo stato eccitato
double GetAlphaGroundState (const pos& old_point, const pos& new_point);
double GetAlphaFirstExcState (const pos& old_point, const pos& new_point);

// Realizza singolo step dell'algoritmo di Metropolis, si possono scrivere su file le posizioni ricoperte
void Metropolis (Random& rnd, pos& point, double dim_step, const string& random_type, const string& simulation_type, bool write_config);

// Calcola l'accettazione media su un ciclo di <num_step> 
double Acceptance (Random& rnd, pos& point, double dim_step, int num_step, const string& random_type, const string& simulation_type);

#endif