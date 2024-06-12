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
    //ESERCIZIO 2.2.1   Random walk su reticolo
    //------------------------------------------------------------

    int M = 10000;                          // Numero RW simulati per numero step
    int N = 100;                            // Numero blocchi per numero step
    int L = M/N;                            // Numero RW per blocco

    int a = 1;                              // Passo reticolo
    pos part;                               // Particella in (0, 0, 0)
    int num_step_max = 100;                 // Numero step massimo

    ofstream out;
    out.open("./OUTPUT/Results_ex_2_2_1.out");
    out << setw(6) << "#STEP:" << setw(15) << "AVERAGE:" << setw(15) << "ERROR:" << endl;

    for (int i = 1; i <= num_step_max; i++){                                            // Ciclo su numero step                   
        vector<double> ave(N,0.);                                                       // Vettore medie
        vector<double> ave2(N,0.);                                                      // Vettore medie quadratiche
        vector<vector<double>> ValErr(N, vector<double>(2, 0.));                        // Per media a blocchi
        for (int j = 0; j < N; j++){                                                    // Ciclo su numero blocchi
            double sum = 0.0;
            for (int k = 0; k < L; k++){                                                // Ciclo su numero RW per blocco
                part.x = 0; part.y = 0; part.z = 0;                                     // Particella parte da origine
                for (int l = 0; l < i; l++){                                            // Numero di step da eseguire
                    int dir = int(rnd.Rannyu(0, 3));                                    // Direzione: 0 -> x   1 -> y   2 -> z
                    int ver = int(rnd.Rannyu(0, 2));                                    // Verso:     0 -> +   1 -> -
                    part = NewPositionLattice(part, dir, ver, a);                       // Nuova posizione particella
                }
                double dist2 = pow(part.x, 2) + pow(part.y, 2) + pow(part.z, 2);        // Distanza quadratica da origine
                sum += sqrt(dist2);                                             
            }
            ave[j] = sum / L;                                                           // Media
            ave2[j] = pow(sum/L,2);                                                     // Media quadratica
        }
        ValErr = BlockAverage(ave, ave2);                                               // Media a blocchi con incertezza
        out << setw(6) << i 
            << setw(15) << fixed << setprecision(7) << ValErr[N-1][0] 
            << setw(15) << fixed << setprecision(7) << ValErr[N-1][1] << endl;          // Stampa solo ultima media cumulativa 
    }

    out.close();



    //------------------------------------------------------------
    //ESERCIZIO 2.2.2    Random walk nel continuo
    //------------------------------------------------------------

    out.open("./OUTPUT/Results_ex_2_2_2.out");
    out << setw(6) << "#STEP:" << setw(15) << "AVERAGE:" << setw(15) << "ERROR:" << endl;

    for (int i = 1; i <= num_step_max; i++){                                            // Ciclo su numero step                   
        vector<double> ave(N,0.);                                                       // Vettore medie
        vector<double> ave2(N,0.);                                                      // Vettore medie quadratiche
        vector<vector<double>> ValErr(N, vector<double>(2, 0.));                        // Per media a blocchi
        for (int j = 0; j < N; j++){                                                    // Ciclo su numero blocchi
            double sum = 0.0;
            for (int k = 0; k < L; k++){                                                // Ciclo su numero RW per blocco
                part.x = 0; part.y = 0; part.z = 0;                                     // Particella parte da origine
                for (int l = 0; l < i; l++){                                            // Numero di step da eseguire
                    double theta = acos(1. - 2.*rnd.Rannyu());                          // Angolo theta
                    double phi = rnd.Rannyu(0, 2.*M_PI);                                // Angolo phi
                    part = NewPositionSpherical(part, theta, phi, double(a));           // Nuova posizione particella
                }
                double dist2 = pow(part.x, 2) + pow(part.y, 2) + pow(part.z, 2);        // Distanza quadratica da origine
                sum += sqrt(dist2);                                             
            }
            ave[j] = sum / L;                                                           // Media
            ave2[j] = pow(sum/L,2);                                                     // Media quadratica
        }
        ValErr = BlockAverage(ave, ave2);                                               // Media a blocchi con incertezza
        out << setw(6) << i 
            << setw(15) << fixed << setprecision(7) << ValErr[N-1][0] 
            << setw(15) << fixed << setprecision(7) << ValErr[N-1][1] << endl;          // Stampa solo ultima media cumulativa 
    }

    out.close();



    return 0;
}