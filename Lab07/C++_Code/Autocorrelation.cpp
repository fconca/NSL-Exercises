#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

double autocorrelation(int, const vector<int>&, const vector<double>&);

int main() {

    // Lettura da file valori di tempo e energia potenziale
    ifstream in;
    in.open("potential_energy_gas.dat");

    int t;
    double ener, ave, err;
    string row;
    getline(in, row);

    vector<int> time;
    vector<double> pot_ener;

    while(in >> t >> ener >> ave >> err){
        time.push_back(t-1);
        pot_ener.push_back(ener);
    }

    in.close();



    // Calcolo autocorrelazione e stampa file output
    ofstream out;
    out.open("autocorrelation_gas.dat");
    out << setw(6) << "#TIME:" << setw(10) << "AUTOCORR:" << endl;

    for(int i=0; i < 10000; i++)
        out << setw(6) << i+1 << setw(10) << fixed << setprecision(6) << autocorrelation(i, time, pot_ener) << endl;

    out.close();



    return 0;
}

double autocorrelation(int t, const vector<int>& time, const vector<double>& ener){
    int t_max = time.size();
    double sum1, sum2, sum3, sum4, sum5 = 0;

    double diff = 1.0/double(t_max - t);
    double inv = 1.0/double(t_max);

    for(int i=0; i < t_max-t; i++){
        sum1 += ener[i] * ener[i+t];
        sum2 += ener[i];
        sum3 += ener[i+t];
    }

    for(int i=0; i < t_max; i++){
        sum4 += ener[i] * ener[i];
        sum5 += ener[i];
    }

    return (diff * sum1 - diff * sum2 * diff * sum3) / (inv * sum4 - pow(inv * sum5, 2));
}