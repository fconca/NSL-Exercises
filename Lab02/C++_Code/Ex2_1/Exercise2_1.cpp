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

    //------------------------------------------------------------
    //ESERCIZIO 2.1.1   Metodo della media (semplice)
    //------------------------------------------------------------

    int M = 100000;             // Numeri random generati in tutto 
    int N = 100;                // Numero blocchi
    int L = M/N;                // Numero elementi per blocco
    vector<double> ave(N);      // Valori medi
    vector<double> ave2(N);     // Valori quadratici medi

    ofstream out;

    for (int i = 0; i < N; i++){                                        // Ciclo su numero blocchi
        double sum = 0.0;
        for (int j = 0; j < L; j++){                                    // Ciclo su numero elementi per blocco
            double rand = rnd.Rannyu();                                 // Numero random da distrib. uniforme
            double f_value = M_PI/2 * cos(M_PI/2.0 * double(rand));     // Valore di f(x)
            sum += f_value;
        }
        ave[i] = sum / L;                                               // Media
        ave2[i] = pow(ave[i], 2);                                       // Media quadratica
    }

    vector<vector<double>> ValErr = BlockAverage(ave, ave2);            // Media a blocchi e incertezza

    out.open("./OUTPUT/Results_ex_2_1_1.out"); 
    out << setw(7) << "#BLOCK:" << setw(15) << "AVERAGE:" << setw(15) << "ERROR:" << endl;

    for (int i = 0; i < N; i++) {
        out << setw(7) << i+1 
            << setw(15) << fixed << setprecision(7) << ValErr[i][0] - 1.0
            << setw(15) << fixed << setprecision(7) << ValErr[i][1] << endl;
    }
    out.close();                                          



    //------------------------------------------------------------
    //ESERCIZIO 2.1.2   Importance Sampling
    //------------------------------------------------------------

    // Approssimazione p(x) = -2*x + 2
    // Funzione cumulativa F(x) = -x**2 +2*x
    // Funzione inversa della cumulativa G(x) = 1-sqrt(1-x)

    for (int i = 0; i < N; i++){                                                    // Ciclo su numero blocchi
        double sum = 0;
        for (int j = 0; j < L; j++){                                                // Ciclo su numero elementi per blocco
            double rand = 1.0 - sqrt(1.0 - rnd.Rannyu());                           // Numero random da distrib. lineare
            double g_value = M_PI/2.0 * cos(M_PI/2.0 * rand) / (-2.0 * rand + 2.0); // Valore di g(x) = f(x) / p(x)
            sum += g_value;
        }
        ave[i] = sum / L;                                                           // Media
        ave2[i] = pow(ave[i], 2);                                                   // Media quadratica
    }
 
    ValErr = BlockAverage(ave, ave2);                                               // Media a blocchi e incertezza

    out.open("./OUTPUT/Results_ex_2_1_2.out");
    out << setw(7) << "#BLOCK:" << setw(15) << "AVERAGE:" << setw(15) << "ERROR:" << endl;

    for (int i = 0; i < N; i++) {
        out << setw(7) << i+1 
            << setw(15) << fixed << setprecision(7) << ValErr[i][0] - 1.0
            << setw(15) << fixed << setprecision(7) << ValErr[i][1] << endl;
    }
    out.close();  



    return 0;
}


