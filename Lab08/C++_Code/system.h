#ifndef __System__
#define __System__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

class System {

private:
    //PARAMETRI LETTI DA FILE INPUT
    double _x;              // Coordinata x funzione d'onda
    double _x_step;         // Dimensione step Metropolis per generare nuova coordinata
    int _num_x;             // Numero di coordinate da generare per valutare <H>
    double _mu;             // Parametro media funzione d'onda
    double _mu_step;        // Dimensione step Metropolis per generare nuova media
    double _sigma;          // Parametro dev. stnd. funzione d'onda
    double _sigma_step;     // Dimensione step Metropolis per generare nuova dev. stnd.
    int _num_block;         // Numero di blocchi per media a blocchi 
    int _dim_block;         // Dimensione blocco per media a blocchi 
    double _beta;           // Temperatura simulated annealing
    double _beta_step;      // Step temperature da simulare 
    int _num_beta;          // Numero temperature da simulare
    int _num_sim;           // Numero di simulazioni per temperatura fissata

    //PARAMETRI INTERNI 
    Random _rnd;                // Generatore numeri casuali
    double _acc_x = 0.0;        // Accettazione Metropolis coordinata x
    int _actual_step_x = 0;     // Numero step corrente Metropolis coordinata x 
    double _acc_boltz = 0.0;    // Accettazione Metropolis mu e sigma
    double _L;                  // Funzione costo: valore aspettazione hamiltoniana
    double _L_err;              // Errore funzione costo

public:
    void InitializeRandomGenerator (); // Inizializza generatore numeri casuali
    void Initialize (); // Lettura parametri di input e inizializzaazione <H>
    void InitializeCostFunction (); // Calcolo <H> assegnati mu e sigma di input
    pair<double, double> BlockAverage (const double& sum_ave, const double& sum_ave2, const int& block); // Calcolo media a blocchi con errore
    double ProbFunction (const double& x); // Funzione d'onda di test calcolata in x
    double DensityFunction (const double& x); // Densità di probabilità calcolata in x
    double Hamiltonian (const double& x); // Valore dell'hamiltoniana in x
    double min (const double& a, const double& b); // Minimo tra due valori
    void Metro (); // Step Metropolis per campionamento densità di probabilità
    pair<double, double> ExpectedHamiltonian (); // Valore di <H> con errore con media a blocchi
    void Boltzmann (); // Step Metropolis-Boltzmann per SA
    void SimulatedAnnealing (); // Simulated Annealing al variare di beta
    void PrintCostFunction (double mu_min, double mu_max, double sigma_min, double sigma_max, int num_point); // Stampa in output <H> calcolato al variare di mu e sigma
    void FindBestHamiltonian (int num_block, int dim_block, int num_x, double mu, double sigma); // Valore di <H> con errore con media a blocchi per assegnati valori di mu e sigma
    void PrintDensityFunction (double mu, double sigma, int num); // Stampa in output valori di |Psi(x)|^2
};

#endif // __System__