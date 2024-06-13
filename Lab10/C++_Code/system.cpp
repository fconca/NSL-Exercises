#include "system.h"

using namespace std;

void System :: InitializeRandomGenerator (int line) { // Inizializza generatore numeri casuali
    if (line < 0 || line > 383) {
        cerr << "PROBLEM in function System :: InitializeRandomGenerator: line must be between 0 and 383" << endl;
        return;
    }

    int seed[4];
    int p1, p2;

    ifstream Primes("./INPUT/Primes"); // Lettura dei primi due numeri primi
    if (Primes.is_open()) {
        string line_number;
        for (int i = 0; i <= line; ++i) {
            if (!getline(Primes, line_number)) { // Scorri fino alla riga desiderata
                cerr << "PROBLEM in function System :: InitializeRandomGenerator: unable to read line " << line << " from file ./INPUT/Primes" << endl;
                exit(-1);
            }
        } 
        if (!(Primes >> p1 >> p2)) { // Leggi i numeri dalla riga specificata
            cerr << "PROBLEM in function System :: InitializeRandomGenerator: unable to read prime numbers from line " << line << " in file ./INPUT/Primes" << endl;
            exit(-1);
        }
    }
    else {
        cerr << "PROBLEM: Unable to open file ./INPUT/Primes" << endl;
        exit(-1);
    } 
    Primes.close();

    ifstream input("./INPUT/seed.in");           // Lettura del seme 
    string property;
    if (input.is_open()){
        while ( !input.eof() ){
            input >> property;
            if( property == "RANDOMSEED" ){
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                _rnd.SetRandom(seed,p1,p2);
            }
        }
        input.close();
    } 
    else {
        cerr << "PROBLEM: Unable to open file ./INPUT/seed.in" << endl;
        exit(-1);
    }
}

void System :: Initialize() { // Inizializza il sistema leggendo parametri di input
    ifstream inp;
    inp.open("./INPUT/input.in"); // File di input
    if (inp.is_open()){
        string parameter_name;
        while (inp >> parameter_name) {
            if (parameter_name == "#ENDINPUT") 
                break;
            if (parameter_name == "NUM_CITIES") {               
                inp >> _num_cities;
                if (_num_cities <= 0) {
                    cerr << "PROBLEM while reading ./INPUT.input.in: NUM_CITIES must be a positive integer" << endl;
                    exit(-1);
                }
            }
            else if (parameter_name == "SQUARE_CIRCUM") {          
                inp >> _square_circum;
                if (_square_circum != 0 && _square_circum != 1) {
                    cerr << "PROBLEM while reading ./INPUT.input.in: SQUARE_CIRCUM must be 0 or 1" << endl;
                    exit(-1);
                }
            }
            else if (parameter_name == "NUM_PATHS") {        
                inp >> _num_paths;
                if (_num_paths <= 0 || _num_paths % 2 != 0) {
                    cerr << "PROBLEM while reading ./INPUT.input.in: NUM_PATHS must be a positive even integer" << endl;
                    exit(-1);
                }
            }
            else if (parameter_name == "METRIC") {
                inp >> _metric;
                if (_metric != 0 && _metric != 1) {
                    cerr << "PROBLEM while reading ./INPUT.input.in: METRIC must be 0 or 1" << endl;
                    exit(-1);
                }
            }
            else if (parameter_name == "NUM_STEP") {
                inp >> _num_step;
                if (_num_step <= 0) {
                    cerr << "PROBLEM while reading ./INPUT.input.in: NUM_STEP must be a positive integer" << endl;
                    exit(-1);
                }
            }
            else if (parameter_name == "NUM_EXCHANGE") {
                inp >> _num_exchange;
                if (_num_exchange <= 0) {
                    cerr << "PROBLEM while reading ./INPUT.input.in: NUM_EXCHANGE must be a positive integer" << endl;
                    exit(-1);
                }
            }
            else if (parameter_name == "EXPONENT")
                inp >> _expo;
            else if (parameter_name == "PROB_PERM") {
                inp >> _prob_perm;
                if (_prob_perm < 0 || _prob_perm > 1) {
                    cerr << "PROBLEM while reading ./INPUT.input.in: PROB_PERM must be between 0 and 1" << endl;
                    exit(-1);
                }
            }
            else if (parameter_name == "PROB_SHIFT_GROUP") {
                inp >> _prob_shift_group;
                if (_prob_shift_group < 0 || _prob_shift_group > 1) {
                    cerr << "PROBLEM while reading ./INPUT.input.in: PROB_SHIFT_GROUP must be between 0 and 1" << endl;
                    exit(-1);
                }
            }
            else if (parameter_name == "PROB_PERM_GROUP") {
                inp >> _prob_perm_group;
                if (_prob_perm_group < 0 || _prob_perm_group > 1) {
                    cerr << "PROBLEM while reading ./INPUT.input.in: PROB_PERM_GROUP must be between 0 and 1" << endl;
                    exit(-1);
                }
            }
            else if (parameter_name == "PROB_INV") {
                inp >> _prob_inv;
                if (_prob_inv < 0 || _prob_inv > 1) {
                    cerr << "PROBLEM while reading ./INPUT.input.in: PROB_INV must be between 0 and 1" << endl;
                    exit(-1);
                }
            }
            else if (parameter_name == "PROB_CROSS") {
                inp >> _prob_cross;
                if (_prob_cross < 0 || _prob_cross > 1) {
                    cerr << "PROBLEM while reading ./INPUT.input.in: PROB_CROSS must be between 0 and 1" << endl;
                    exit(-1);
                }
            }
        }
    }
    else {
        cerr << "PROBLEM: Unable to open file ./INPUT/input.in" << endl;
        exit(-1);
    }
    inp.close();

    ofstream out;
    out.open("./OUTPUT/parameters.out");
    if (out.is_open()) {
        out << setw(17) << "NUM_CITIES:" << setw(10) << _num_cities << endl
            << setw(17) << "SQUARE_CIRCUM:" << setw(10) << _square_circum << endl
            << setw(17) << "NUM_PATHS:" << setw(10) << _num_paths << endl
            << setw(17) << "METRIC:" << setw(10) << _metric << endl
            << setw(17) << "NUM_STEP:" << setw(10) << _num_step << endl
            << setw(17) << "NUM_EXCHANGE:" << setw(10) << _num_exchange << endl
            << setw(17) << "EXPONENT:" << setw(10) << _expo << endl
            << setw(17) << "PROB_PERM:" << setw(10) << _prob_perm << endl
            << setw(17) << "PROB_SHIFT_GROUP:" << setw(10) << _prob_shift_group << endl
            << setw(17) << "PROB_PERM_GROUP:" << setw(10) << _prob_perm_group << endl
            << setw(17) << "PROB_INV:" << setw(10) << _prob_inv << endl
            << setw(17) << "PROB_CROSS:" << setw(10) << _prob_cross << endl;
    }
    else {
        cerr << "PROBLEM: Unable to open file ./OUTPUT/parameters.out" << endl;
        exit(-1);
    }
    out.close();
}

double System :: GetRandom (const double& inf, const double& sup) { // Generatore uniforme di numeri casuali
    return _rnd.Rannyu(inf, sup);
}

int System :: GetNumCities () { // Ottieni umero di città
    return _num_cities;
}

int System :: GetSquareCircum () { // 0 per generare città nel quadrato, 1 sulla circonferenza
    return _square_circum;
}

int System :: GetNumPaths () { // Ottieni numero di percorsi (popolazione)
    return _num_paths;
}

int System :: GetMetric () { // 0 per usare metrica L1, 1 per usare metrica L2
    return _metric;
}

int System :: GetNumStep () { // Ottieni numero step algoritmo genetico
    return _num_step;
}

int System :: GetNumExchange () { // Ottieni intervallo scambi tra continenti
    return _num_exchange;
}

double System :: GetExponent () { // Ottieni esponente per algoritmo selezione
    return _expo;
}

double System :: GetProbPerm () { // Ottieni probabilità permutazione singola
    return _prob_perm;
}

double System :: GetProbShiftGroup () { // Ottieni probabilità shift gruppo
    return _prob_shift_group;
}

double System :: GetProbPermGroup () { // Ottieni probabilità permutazione gruppo
    return _prob_perm_group;
}

double System :: GetProbInversion () { // Ottieni probabilità inversione gruppo
    return _prob_inv;
}

double System :: GetProbCross () { // Ottieni probabilità di crossover
    return _prob_cross;
}