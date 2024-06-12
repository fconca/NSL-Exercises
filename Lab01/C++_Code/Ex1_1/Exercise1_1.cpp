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

    //-----------------------------------------------------------
    //ESERCIZIO 1.1.1   Valor medio distribuzione uniforme
    //-----------------------------------------------------------

    int M = 100000;                             // Numero totale elementi
    int N = 100;                                // Numero blocchi
    int L = M/N;                                // Numero elementi per blocco
    vector<double> ave(N);                      // Medie per blocco
    vector<double> ave2(N);                     // Medie quadratiche per blocco

    ofstream out;

    for (int i = 0; i < N; i++){                // Ciclo su numero blocchi  
        double sum = 0;
        for (int j = 0; j < L; j++){            // Ciclo su numero elementi per blocco  
            double r = rnd.Rannyu();            // Numero pseudo-casuale da distrib. uniforme
            sum += r;
        }
        ave[i] = sum/L;                         // Media del blocco i-esimo
        ave2[i] = pow(ave[i],2);                // Media quadratica del blocco i-esimo
    }

    vector<vector<double>> ValErr = BlockAverage(ave, ave2); // Media a blocchi con incertezza

    out.open("./OUTPUT/Results_ex_1_1_1.out");  
    out << setw(7) << "#BLOCK:" << setw(15) << "AVERAGE:" << setw(15) << "ERROR:" << endl;
    for (int i = 0; i < N; i++) {
        out << setw(7) << i+1 
            << setw(15) << fixed << setprecision(7) << ValErr[i][0] - 0.5 
            << setw(15) << fixed << setprecision(7) << ValErr[i][1] << endl;
    }
    out.close();



    //-----------------------------------------------------------
    //ESERCIZIO 1.1.2   Valor medio della varianza
    //-----------------------------------------------------------

    for (int i = 0; i < N; i++){                // Ciclo su numero blocchi        
        double sum = 0;
        for (int j = 0; j < L; j++){            // Ciclo su numero elementi per blocco     
            double r = rnd.Rannyu();            // Numero pseudo-casuale da distrib. uniforme
            sum += pow((r-0.5),2);
        }
        ave[i] = sum/L;                         // Media del blocco i-esimo
        ave2[i] = pow(ave[i],2);                // Media quadratica del blocco i-esimo
    }

    ValErr = BlockAverage(ave, ave2);           // Media a blocchi con incertezza

    out.open("./OUTPUT/Results_ex_1_1_2.out");  // Stampa output
    out << setw(7) << "#BLOCK:" << setw(15) << "AVERAGE:" << setw(15) << "ERROR:" << endl;
    for (int i = 0; i < N; i++) {
        out << setw(7) << i+1 
            << setw(15) << fixed << setprecision(7) << ValErr[i][0] - 1.0/12 
            << setw(15) << fixed << setprecision(7) << ValErr[i][1] << endl;
    }
    out.close();



    //-----------------------------------------------------------
    //ESERCIZIO 1.1.3   Test di Pearson chi-quadro
    //-----------------------------------------------------------

    int num_box = 100;                              // Numero sottointervalli in [0,1]
    int num_rand = 10000;                           // Numeri random generati per ciclo
    int num_iter = 100;                             // Numero di ripetizioni

    out.open("./OUTPUT/Results_ex_1_1_3.out");
    out << setw(6) << "#ITER:" << setw(15) << "CHI_SQUARE:" << endl;

    for (int i = 0; i < num_iter; i++){                 // Ciclo su numero di iterazioni
        vector<int> num_occ(num_box, 0);                // Salva numero occorrenze in ogni sottointervallo
        for (int j = 0; j < num_rand; j++){             // Ciclo su numeri random da generare per singola iterazione
            double r = rnd.Rannyu(0, num_box);          // Numero random da distrib. uniforme tra 0 e 100
            int bin = int(r);                                       
            num_occ[bin] += 1;                          // Incrementa bin corrispondente
        }
        double chi2 = 0.0;
        for (int k = 0; k < num_box; k++){
            double chi2_box = pow(double(num_occ[k]) - double(num_rand) / num_box, 2) / (double(num_rand) / num_box);
            chi2 += chi2_box;
        }
        out << setw(6) << i+1 << setw(15) << fixed << setprecision(3) << chi2 << endl;
    }

    out.close();



    return 0;
}