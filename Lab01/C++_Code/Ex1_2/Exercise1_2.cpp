#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"
#include "myfunctions.h"

using namespace std;
 
int main (){

    Random rnd;
    InitializeRandomGenerator("./INPUT/Primes", "./INPUT/seed.in", rnd); // Inizilizza generatore numeri casuali

    int M = 10000;                      // Numero di somme S per ciascun valore di N
    int N[4] = {1,2,10,100};            // Numeri usati per calcolo media
    ofstream out;

    

    //-----------------------------------------------------------
    //ESERCIZIO 1.2.1   Distribuzione uniforme
    //-----------------------------------------------------------

    out.open("./OUTPUT/Results_ex_1_2_1.out");

    for (int i = 0; i < 4; i++){                    // 4 cicli per valori di N = 1, 2, 10, 100
        for (int j = 0; j < M; j++){                // Ciclo su numero di somme S    
            double sum = 0.0;
            for (int k = 0; k < N[i]; k++) {        // Media calcolata su N numeri
                double r = rnd.Rannyu();            // Numero casuale da distrib. uniforme
                sum += r;
            }
            out << sum / N[i] << endl;
        }
    }

    out.close();



    //-----------------------------------------------------------
    //ESERCIZIO 1.2.2   Distribuzione esponenziale
    //-----------------------------------------------------------

    out.open("./OUTPUT/Results_ex_1_2_2.out");

    for (int i = 0; i < 4; i++){                    // 4 cicli per valori di N = 1, 2, 10, 100
        for (int j = 0; j < M; j++){                // Ciclo su numero di somme S  
            double sum = 0.;
            for (int k=0; k<N[i]; k++) {            // Media calcolata su N numeri
                double r = rnd.Expo(1.);            // Numero casuale da distrib. esponenziale
                sum += r;
            }
            out << sum / N[i] << endl;
        }
    }

    out.close();



    //-----------------------------------------------------------
    //ESERCIZIO 1.2.3   Distribuzione Cauchy-Lorentz
    //-----------------------------------------------------------

    out.open("./OUTPUT/Results_ex_1_2_3.out");

    for (int i = 0; i < 4; i++){                    // 4 cicli per valori di N = 1, 2, 10, 100
        for (int j = 0; j < M; j++){                // Ciclo su numero di somme S 
            double sum = 0.;
            for (int k=0; k<N[i]; k++) {            // Media calcolata su N numeri
                double r = rnd.Cauchy(0., 1.);      // Numero casuale da distrib. Cauchy-Lorentz
                sum += r;
            }
            out << sum / N[i] << endl;
        }
    }

    out.close();



    return 0;
}