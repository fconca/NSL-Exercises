#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

// Funzione per calcolare l'autocorrelazione
double autocorrelation(int, const vector<int>&, const vector<double>&);

int main() {

    // Lettura dei valori di tempo e energia potenziale da un file
    ifstream in;
    in.open("potential_energy_sol.dat");

    int t;
    double ener, ave, err;
    string row;
    getline(in, row);  // Legge l'intestazione del file e la ignora

    vector<int> time;
    vector<double> pot_ener;

    // Legge i dati dal file e li memorizza nei vettori time e pot_ener
    while(in >> t >> ener >> ave >> err){
        time.push_back(t-1);  // Memorizza il tempo decrementato di 1
        pot_ener.push_back(ener);  // Memorizza l'energia potenziale
    }

    in.close();

    // Calcolo dell'autocorrelazione e scrittura dei risultati su file di output
    ofstream out;
    out.open("autocorrelation_sol.dat");
    out << setw(6) << "#TIME:" << setw(15) << "AUTOCORR:" << endl;

    // Calcola e scrive l'autocorrelazione per i primi 10000 intervalli di tempo
    for(int i = 0; i < 10000; i++) {
        out << setw(6) << i + 1 << setw(15) << fixed << setprecision(6) << autocorrelation(i, time, pot_ener) << endl;
    }

    out.close();

    return 0;
}

// Funzione che calcola l'autocorrelazione per un dato intervallo di tempo t
double autocorrelation(int t, const vector<int>& time, const vector<double>& ener) {
    int t_max = time.size();
    double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0;

    double diff = 1.0 / double(t_max - t);
    double inv = 1.0 / double(t_max);

    // Calcola i termini necessari per l'autocorrelazione
    for(int i = 0; i < t_max - t; i++) {
        sum1 += ener[i] * ener[i + t];  
        sum2 += ener[i];  
        sum3 += ener[i + t];  
    }

    for(int i = 0; i < t_max; i++) {
        sum4 += ener[i] * ener[i];  
        sum5 += ener[i];  
    }

    // Ritorna il valore dell'autocorrelazione normalizzata
    return (diff * sum1 - diff * sum2 * diff * sum3) / (inv * sum4 - pow(inv * sum5, 2));
}
