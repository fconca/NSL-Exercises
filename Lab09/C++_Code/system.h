#ifndef __System__
#define __System__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

class System{
private:
    Random _rnd; // Generatore numeri casuali
    int _num_cities; // Numero di città
    int _square_circum; // Generazione su circonferenza o nel quadrato
    int _num_paths; // Numero di percorsi (popolazione)
    int _metric; // Metrica L1 o L2
    int _num_step; // Numero di step dell'algoritmo genetico
    double _expo; // Esponente algoritmo selezione
    double _prob_perm; // Probabilità permutazione singola
    double _prob_shift_group; // Probabilità shift gruppo
    double _prob_perm_group; // Probabilità permutazione gruppo
    double _prob_inv; // Probabilità inversione gruppo
    double _prob_cross; // Probabilità di crossover

public:
    void InitializeRandomGenerator (); // Inizializza generatore numeri casuali
    void Initialize (); // Inizializza il sistema leggendo parametri di input
    double GetRandom (const double& inf, const double& sup); // Generatore uniforme di numeri casuali
    int GetNumCities (); // Ottieni umero di città
    int GetSquareCircum (); // 0 per generare città nel quadrato, 1 sulla circonferenza
    int GetNumPaths (); // Ottieni numero di percorsi (popolazione)
    int GetMetric (); // 0 per usare metrica L1, 1 per usare metrica L2
    int GetNumStep (); // Ottieni numero step algoritmo genetico
    double GetExponent (); // Ottieni esponente per algoritmo selezione
    double GetProbPerm (); // Ottieni probabilità permutazione singola
    double GetProbShiftGroup (); // Ottieni probabilità shift gruppo
    double GetProbPermGroup (); // Ottieni probabilità permutazione gruppo
    double GetProbInversion (); // Ottieni probabilità inversione gruppo
    double GetProbCross (); // Ottieni probabilità di crossover
};

#endif // __System__