#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <iomanip>
#include "random.h"
#include "system.h"

using namespace std;

int main (){

    System SYS;

    SYS.InitializeRandomGenerator();                            // Inizializza generatore numeri casuali
    SYS.Initialize();                                           // Inizializza parametri simulazione

    SYS.SimulatedAnnealing();                                   // Esegue algoritmo simulated annealing
    //SYS.FindBestHamiltonian(100, 10000, 1000, 0.804, 0.617);  // Valore di <H> con media a blocchi fissati mu, sigma

    //SYS.PrintCostFunction(0.3, 1.3, 0.2, 1.0, 50);            // Funzione costo nello spazio di mu, sigma
    //SYS.PrintDensityFunction(0.804, 0.617, 1000000);          // Funzione densità di probabilità con Metropoli fissati mu, sigma

    return 0;
}