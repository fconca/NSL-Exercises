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
    InitializeRandomGenerator("./INPUT/Primes", "./INPUT/seed.in", rnd);

    //------------------------------------------------------------
    //ESERCIZIO 3.1
    //------------------------------------------------------------

    double S0 = 100.;           // Asset price S(0) al tempo t=0
    double T = 1.;              // Delivery time
    double K = 100.;            // Strike price
    double r = 0.1;             // Interest rate
    double sigma = 0.25;        // Volatility

    int M = 100000;             // Numero totale di simulazioni
    int N = 100;                // Numero blocchi
    int L = M/N;                // Numero elementi per blocco

    vector<double> ave_C(N);    // Vettore medie opzione call
    vector<double> ave2_C(N);   // Vettore medie quadratiche opzione call
    vector<double> ave_P(N);    // Vettore medie opzione put
    vector<double> ave2_P(N);   // Vettore medie quadratiche opzione put

    ofstream out;

    //-----------------------------------------------------------
    //ESERCIZIO 3.1.1   Calcolo diretto di S(T)
    //-----------------------------------------------------------

    for (int i = 0; i < N; i++) {                                                       // Ciclo su numero blocchi
        double sum_C = 0.0;
        double sum_P = 0.0;
        for (int j = 0; j < L; j++) {                                                   // Ciclo su numero elementi per blocco
            double Z = rnd.Gauss(0, 1);                                                 // Numero random con distrib. normale
            double S_t = S0 * exp((r - 0.5 * sigma * sigma) * T + sigma * Z * sqrt(T)); // Prezzo al tempo t=T
            double C = exp(-r * T) * max(0, S_t - K);                                   // Prezzo opzione call
            double P = exp(-r * T) * max(0, K - S_t);                                   // Prezzo opzione put
            sum_C += C;
            sum_P += P;
        }
        ave_C[i] = sum_C / L;                                                           // Media call
        ave2_C[i] = pow(ave_C[i], 2);                                                   // Media quadratica call
        ave_P[i] = sum_P / L;                                                           // Media put
        ave2_P[i] = pow(ave_P[i], 2);                                                   // Media quadratica put
    }

    vector<vector<double>> ValErr = BlockAverage(ave_C, ave2_C);                        // Media a blocchi con incertezza call

    out.open("./OUTPUT/Results_ex_3_1_1_call.out");
    out << setw(7) << "#BLOCK:" << setw(15) << "AVERAGE:" << setw(15) << "ERROR:" << endl;

    for (int i = 0; i < N; i++) {
        out << setw(7) << i+1
            << setw(15) << fixed << setprecision(7) << ValErr[i][0] - 14.975790778311286 
            << setw(15) << fixed << setprecision(7) << ValErr[i][1] << endl;
    }

    out.close();

    ValErr = BlockAverage(ave_P, ave2_P);                                               // Media a blocchi con incertezza put

    out.open("./OUTPUT/Results_ex_3_1_1_put.out");
    out << setw(7) << "#BLOCK:" << setw(15) << "AVERAGE:" << setw(15) << "ERROR:" << endl;

    for (int i = 0; i < N; i++) {
        out << setw(7) << i+1
            << setw(15) << fixed << setprecision(7) << ValErr[i][0] - 5.4595325819072364 
            << setw(15) << fixed << setprecision(7) << ValErr[i][1] << endl;
    }

    out.close();



    //-----------------------------------------------------------
    //ESERCIZIO 3.1.2   Calcolo per step temporali discreti
    //-----------------------------------------------------------

    for (int i = 0; i < N; i++) {                                                       // Ciclo su numero blocchi
        double sum_C = 0.0;
        double sum_P = 0.0;
        for (int j = 0; j < L; j++) {                                                   // Ciclo su numero elementi per blocco
            double S_t = S0;    
            for (int k = 0; k < 100; k++) {                                             // Ciclo su numero step temporali
                double Z = rnd.Gauss(0, 1);                                             // Numero random con distrib. normale
                S_t = S_t * exp((r - 0.5 * sigma * sigma) * (T / 100.0) + sigma * Z * sqrt(T / 100.0)); // Prezzo al tempo t
            }
            double C = exp(-r * T) * max(0, S_t - K);                                   // Prezzo opzione call
            double P = exp(-r * T) * max(0, K - S_t);                                   // Prezzo opzione put
            sum_C += C;
            sum_P += P;
        }
        ave_C[i] = sum_C / L;                                                           // Media call
        ave2_C[i] = pow(ave_C[i], 2);                                                   // Media quadratica call
        ave_P[i] = sum_P / L;                                                           // Media put
        ave2_P[i] = pow(ave_P[i], 2);                                                   // Media quadratica put
    }

    ValErr = BlockAverage(ave_C, ave2_C);                                               // Media a blocchi con incertezza call

    out.open("./OUTPUT/Results_ex_3_1_2_call.out");
    out << setw(7) << "#BLOCK:" << setw(15) << "AVERAGE:" << setw(15) << "ERROR:" << endl;

    for (int i = 0; i < N; i++) {
        out << setw(7) << i+1
            << setw(15) << fixed << setprecision(7) << ValErr[i][0] - 14.975790778311286 
            << setw(15) << fixed << setprecision(7) << ValErr[i][1] << endl;
    }

    out.close();

    ValErr = BlockAverage(ave_P, ave2_P);                                               // Media a blocchi con incertezza put

    out.open("./OUTPUT/Results_ex_3_1_2_put.out");
    out << setw(7) << "#BLOCK:" << setw(15) << "AVERAGE:" << setw(15) << "ERROR:" << endl;

    for (int i = 0; i < N; i++) {
        out << setw(7) << i+1
            << setw(15) << fixed << setprecision(7) << ValErr[i][0] - 5.4595325819072364 
            << setw(15) << fixed << setprecision(7) << ValErr[i][1] << endl;
    }

    out.close();



    return 0;

}