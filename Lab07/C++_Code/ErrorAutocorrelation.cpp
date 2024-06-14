#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

// Funzione per calcolare la media a blocchi e le incertezze
vector<vector<double>> BlockAverage(vector<double>, vector<double>);

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
    while (in >> t >> ener >> ave >> err) {
        time.push_back(t - 1);  // Memorizza il tempo decrementato di 1
        pot_ener.push_back(ener);  // Memorizza l'energia potenziale
    }

    in.close();

    // Calcolo della media a blocchi e delle incertezze
    int M = time.size();  // Numero di misure
    vector<int> L = {10, 16, 20, 25, 40, 50, 80, 100, 125, 200, 250, 400, 500, 625, 800, 1000, 1250, 1600, 2000, 2500, 3125, 4000, 5000};
    int num = L.size();  // Numero di dimensioni dei blocchi
    vector<int> N(num);  // Numero di blocchi per ogni dimensione
    for (int i = 0; i < num; i++) N[i] = M / L[i];

    ofstream out;
    out.open("error_autocorrelation_sol.dat");
    out << setw(10) << "DIM_BLOCK:" << setw(15) << "ERROR:" << setw(15) << "REL_ERR:" << endl;

    vector<double> ave1;  // Vettore delle medie dei blocchi
    vector<double> ave2;  // Vettore dei quadrati delle medie dei blocchi
    vector<vector<double>> ValErr;  // Matrice per memorizzare le medie e le incertezze

    for (int i = 0; i < num; i++) {
        int iter = 0;  // Iteratore per scorrere i dati
        ValErr.resize(N[i], vector<double>(2));
        for (int n = 0; n < N[i]; n++) {
            double sum = 0.0;
            ave1.resize(N[i]);
            ave2.resize(N[i]);
            for (int l = 0; l < L[i]; l++) {
                sum += pot_ener[iter];
                iter++;
            }
            ave1[n] = sum / L[i];  // Calcola la media del blocco
            ave2[n] = pow(ave1[n], 2);  // Calcola il quadrato della media del blocco
        }
        ValErr = BlockAverage(ave1, ave2);  // Calcola la media e l'errore per i blocchi
        double value = ValErr[ValErr.size() - 1][0];
        double err = ValErr[ValErr.size() - 1][1];
        out << setw(10) << L[i]
            << setw(15) << fixed << scientific << err
            << setw(15) << fixed << scientific << fabs(err / value) << endl;  // Scrive la dimensione del blocco, l'errore e l'errore relativo
    }

    out.close();

    return 0;
}

// Funzione che calcola la media e l'errore dei blocchi
vector<vector<double>> BlockAverage(vector<double> ave, vector<double> ave2) {
    if (ave.size() != ave2.size()) {
        cerr << "ERROR: dimensions of vector 'ave' and 'ave2' are different" << endl;
        return {};
    }
    int N = ave.size();
    vector<vector<double>> matrix(N, vector<double>(2));
    for (int i = 0; i < N; i++) {
        double sum_prog = 0;
        double sum2_prog = 0;
        for (int j = 0; j < i + 1; j++) {
            sum_prog += ave[j];
            sum2_prog += ave2[j];
        }
        sum_prog /= (i + 1);  // Media del blocco
        sum2_prog /= (i + 1);  // Media dei quadrati del blocco
        double variance = sum2_prog - pow(sum_prog, 2);  // Varianza del blocco
        double err = (i == 0) ? 0 : sqrt(variance / i);  // Deviazione standard del blocco
        matrix[i][0] = sum_prog;
        matrix[i][1] = err;
    }
    return matrix;
}
