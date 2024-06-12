#include "myfunctions.h"

using namespace std;

void InitializeRandomGenerator(const char* primesFile, const char* seedFile, Random& rnd) {
    int seed[4];
    int p1, p2;

    ifstream Primes(primesFile);        // Legge i primi due numeri primi dal file "primesFile"
    if (Primes.is_open()) {
        Primes >> p1 >> p2;
    } 
    else {
        cerr << "PROBLEM: Unable to open file " << primesFile << endl;
        exit(-1);
    }
    Primes.close();

    ifstream input(seedFile);           // Legge il seme dal file "seedFile"
    string property;
    if (input.is_open()){
        while ( !input.eof() ){
            input >> property;
            if( property == "RANDOMSEED" ){
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed,p1,p2);
            }
        }
        input.close();
    } 
    else {
        cerr << "PROBLEM: Unable to open file " << seedFile << endl;
        exit(-1);
    }
}

vector<vector<double>> BlockAverage (const vector<double>& ave, const vector<double>& ave2) {
    if (ave.size() != ave2.size()) {
        cerr << "ERROR: dimensions of vector 'ave' and 'ave2' are different" << endl;
        exit(-1);
    }

    int N = ave.size();
    vector<vector<double>> matrix (N, vector<double>(2));
    for (int i = 0; i < N; i++){                  
        double sum_prog = 0.0;
        double sum2_prog = 0.0;
        for (int j = 0; j < i+1; j++){                         
            sum_prog += ave[j];
            sum2_prog += ave2[j];
        }
        sum_prog /= (i+1);                                      // Media del blocco
        sum2_prog /= (i+1);                                     // Media dei quadrati del blocco
        double variance = sum2_prog - pow(sum_prog, 2);         // Varianza del blocco
        double err = (i == 0) ? 0 : sqrt(variance / i);         // Deviazione standard del blocco
        matrix[i][0] = sum_prog;
        matrix[i][1] = err;
    }
    return matrix;
}

pos NewPositionLattice (const pos& part, int dir, int ver, int step) {
    pos new_pos = part;
    if (dir == 0 && ver == 0)
        new_pos.x += step;
    else if (dir == 0 && ver == 1)
        new_pos.x -= step;
    else if (dir == 1 && ver == 0)
        new_pos.y += step;
    else if (dir == 1 && ver == 1)
        new_pos.y -= step;
    else if (dir == 2 && ver == 0)
        new_pos.z += step;
    else if (dir == 2 && ver == 1)
        new_pos.z -= step;
    return new_pos;
}

pos NewPositionSpherical (const pos& part, double theta, double phi, double step) {
    pos new_pos = part;
    new_pos.x += step * sin(theta) * cos(phi);
    new_pos.y += step * sin(theta) * sin(phi);
    new_pos.z += step * cos(theta);
    return new_pos;
}