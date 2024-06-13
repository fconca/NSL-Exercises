#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <iomanip>
#include "random.h"
#include "myfunctions.h"

using namespace std;
 
int main (){

    Random rnd; 
    InitializeRandomGenerator("./INPUT/Primes", "./INPUT/seed.in", rnd); // Inizializza generatore numeri random

    //------------------------------------------------------------
    //ESERCIZIO 5.1   Metropolis per atomo di idrogeno
    //------------------------------------------------------------

    string simulation_type;         // Tipo simulazione: groundstate o firstexcited
    string random_type;             // Tipo generatore random: uniform o normal
    double dim_step;                // Dimensione step Metropolis
    int num_step_acc;               // Numero step per calcolo accettazione media
    double xcoord;                  // Coordinate punto di partenza
    double ycoord;
    double zcoord;
    int num_block;                  // Numero blocchi
    int dim_block;                  // Dimensione blocco
    bool write_config;              // Per stampare configurazione

    ReadParameters(simulation_type, random_type, dim_step, num_step_acc, xcoord, ycoord, zcoord, num_block, dim_block, write_config);
    cout << "Tipo di simulazione: " << simulation_type << " " << random_type << endl;
    cout << "Step Metropolis: " << dim_step << endl;
    cout << "Punto di partenza: (" << xcoord << ", " << ycoord << ", " << zcoord << ")" << endl;  

    pos point;                                                  // Punto di partenza simulazione
    point.x = xcoord; point.y = ycoord; point.z = zcoord; 
     
    double acc = Acceptance (rnd, point, dim_step, num_step_acc, random_type, simulation_type);     // Accettazione media
    cout << "Probabilità di accettazione:" << acc << endl;



    vector<double> ave(num_block);      // Medie per blocco
    vector<double> ave2(num_block);     // Medie quadratiche per blocco

    for (int i = 0; i < num_block; i++){                                                    // Ciclo su numero blocchi
        pos point;
        point.x = xcoord; point.y = ycoord; point.z = zcoord;                               // Punto partenza
        double r_dist = 0.0;

        for (int j = 0; j < dim_block; j++){                                                // Ciclo su dimensione blocco
            Metropolis (rnd, point, dim_step, random_type, simulation_type, write_config);  // Step Metropolis
            r_dist += dist(point);                                                          // Distanza da origine
        }
        ave[i] = r_dist / dim_block;                                                        // Media
        ave2[i] = pow(ave[i],2);                                                            // Media quadratica
    }

    ofstream out;
    out.open("./OUTPUT/Results_ex_5.out");
    out << setw(6) << "#BLOCK" << setw(15) << "AVERAGE" << setw(15) << "ERROR" << endl;

    vector<vector<double>> ValErr = BlockAverage(ave, ave2);                                // Media a blocchi con incertezza

    for (int i=0; i<num_block; i++) {
        out << setw(6) << i+1 
        << setw(15) << fixed << setprecision(7) << ValErr[i][0]
        << setw(15) << fixed << setprecision(7) << ValErr[i][1] << endl;
    }
    out.close();

    return 0;

}