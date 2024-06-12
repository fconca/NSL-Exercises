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

    //---------------------------------------------------------
    //ESERCIZIO 1.3   Buffon Experiment
    //---------------------------------------------------------

    // Lunghezza ago L e distanza fili d (posta = 1), posizione del filo in x = 0
    // Genero un numero random unif -0.5 < x 0.5 corrispondente a posizione centro dell'ago
    // Angolo inclinazione ago da accept/reject in un cerchio centrato nel centro dell'ago 

    int N_thr = 1000000;                            // Numero di lanci
    int N = 100;                                    // Numero blocchi
    int M = N_thr/N;                                // Numero di lanci per blocco
    vector<double> pi_val(N);                       // Valori medi di pi
    vector<double> pi_val2(N);                      // Valori quadratici medi di pi
    double d = 1.0;                                 // Distanza tra fili
    double L = 0.75;                                // Lunghezza ago 

    for (int i = 0; i < N; i++) {                                                   // Ciclo su numero blocchi
        int N_hit = 0;                                                              // Numero di intersezioni
        int N_ang = 0;                                                              // Numero di angoli generati
        while (N_ang != M) {
            double x_cm = rnd.Rannyu(-d/2, d/2);                                    // Posizione centro ago
            double x_ang = rnd.Rannyu(x_cm - L/2, x_cm + L/2);                      // Coord x per angolo
            double y_ang = rnd.Rannyu(-L/2, L/2);                                   // Coord y per angolo
            double dist = sqrt(pow(x_cm - x_ang, 2) + pow(y_ang, 2));               // Distanza punto generato da centro ago
            if (dist < L/2){                                                        // Condizione accept/reject su angolo
                N_ang ++;
                double cos_ang = abs((x_cm - x_ang)/dist);                          // Coseno angolo
                double x_ext_sx = x_cm - L/2 * cos_ang;                             // Coordinate x estremi ago
                double x_ext_dx = x_cm + L/2 * cos_ang;
                if ((x_cm > 0 && x_ext_sx < 0) || (x_cm < 0 && x_ext_dx > 0)) {     // Condizione di hit
                    N_hit ++;
                }
            }
        }
        pi_val[i] = 2 * L * N_ang / (N_hit * d);                                    // Valor medio di pi 
        pi_val2[i] = pow(pi_val[i], 2);                                             // Valor quadratico medio di pi
    }

    vector<vector<double>> ValErr = BlockAverage(pi_val, pi_val2);                  // Media a blocchi con incertezza

    ofstream out("./OUTPUT/Results_ex_1_3.out");
    out << setw(7) << "#BLOCK:" << setw(15) << "PI_VALUE:" << setw(15) << "ERROR:" << endl;
    for (int i = 0; i < N; i++) {
        out << setw(7) << i+1 
            << setw(15) << fixed << setprecision(7) << ValErr[i][0] - M_PI 
            << setw(15) << fixed << setprecision(7) << ValErr[i][1] << endl;
    }
    out.close();
    


    return 0;
}