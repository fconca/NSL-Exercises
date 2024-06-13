#include "system.h"

using namespace std;

void System :: InitializeRandomGenerator () { // Inizializza generatore numeri casuali
    int seed[4];
    int p1, p2;

    ifstream Primes("./INPUT/Primes");          // Lettura dei primi due numeri primi
    if (Primes.is_open()) {
        Primes >> p1 >> p2;
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

void System :: Initialize() { // Lettura parametri di input e inizializzaazione <H>
    ifstream inp;
    inp.open("./INPUT/input.in");                           // File di input
    if (inp.is_open()){
        string parameter_name;
        while (inp >> parameter_name) {
            if (parameter_name == "#ENDINPUT") 
                break;
            if (parameter_name == "X_START")                // Punto di partenza coordinata x funzione d'onda
                inp >> _x;
            else if (parameter_name == "X_STEP")            // Dimensione step Metropolis per generare nuova coordinata
                inp >> _x_step;
            else if (parameter_name == "DIM_X")             // Numero coordinate da generare per valutare <H>
                inp >> _num_x;
            else if (parameter_name == "MU_START")          // Media di partenza funzione d'onda
                inp >> _mu;
            else if (parameter_name == "MU_STEP")           // Dimensione step Metropolis per generare nuova media
                inp >> _mu_step;
            else if (parameter_name == "SIGMA_START")       // Dev. stnd. di partenza funzione d'onda
                inp >> _sigma;
            else if (parameter_name == "SIGMA_STEP")        // Dimensione step Metropolis per generare nuova dev. stnd.
                inp >> _sigma_step;
            else if (parameter_name == "NUM_BLOCK")         // Numero di blocchi per media a blocchi 
                inp >> _num_block;
            else if (parameter_name == "DIM_BLOCK")         // Dimensione blocco per media a blocchi 
                inp >> _dim_block;
            else if (parameter_name == "BETA_START")        // Temperatura di partenza simulated annealing
                inp >> _beta;
            else if (parameter_name == "BETA_STEP")         // Step temperature da simulare 
                inp >> _beta_step;
            else if (parameter_name == "NUM_BETA")          // Numero temperature da simulare
                inp >> _num_beta;
            else if (parameter_name == "NUM_SIM")           // Numero di simulazioni per temperatura fissata
                inp >> _num_sim; 
        }
    }
    else {
        cerr << "PROBLEM: Unable to open file ./INPUT/input.in" << endl;
        exit(-1);
    }
    inp.close();

    ofstream out;
    out.open("./OUTPUT/parameters.out");                    // Stampa in output dei parametri letti in input
    if (out.is_open()){
        out << "--------------------------------------------------" << endl
            << "PARAMETRI SIMULAZIONE" << endl 
            << "--------------------------------------------------" << endl;
        out << setw(11) << "X_START" << setw(20) << _x << endl
            << setw(11) << "X_STEP" << setw(20) << _x_step << endl
            << setw(11) << "DIM_X" << setw(20) << _num_x << endl
            << setw(11) << "MU_START" << setw(20) << _mu << endl
            << setw(11) << "MU_STEP" << setw(20) << _mu_step << endl
            << setw(11) << "SIGMA_START" << setw(20) << _sigma << endl
            << setw(11) << "SIGMA_STEP" << setw(20) << _sigma_step << endl
            << setw(11) << "NUM_BLOCK" << setw(20) << _num_block << endl
            << setw(11) << "DIM_BLOCK" << setw(20) << _dim_block << endl
            << setw(11) << "BETA_START" << setw(20) << _beta << endl
            << setw(11) << "BETA_STEP" << setw(20) << _beta_step << endl
            << setw(11) << "NUM_BETA" << setw(20) << _num_beta << endl
            << setw(11) << "NUM_SIM" << setw(20) << _num_sim << endl << endl;
    }
    else {
        cerr << "PROBLEM: Unable to open file ./OUTPUT/parameters.out" << endl;
        exit(-1);
    }
    out.close();

    this -> InitializeCostFunction();                       // Inizializzazione di <H> dati parametri di input
}

void System :: InitializeCostFunction () { // Calcolo <H> assegnati mu e sigma di input
    pair<double, double> expect = this -> ExpectedHamiltonian();
    _L = expect.first;                                              // Funzione costo
    _L_err = expect.second;                                         // Errore funzione costo
}

pair<double, double> System::BlockAverage(const double& sum_ave, const double& sum_ave2, const int& block) { // Calcolo media a blocchi con errore
    double mean = sum_ave / block;                              // Media
    double mean2 = sum_ave2 / block;                            // Media quadratica
    double variance = (mean2 - mean * mean) / (block - 1);      // Deviazione standard
    double error = (block == 1) ? 0 : sqrt(variance);           // Incertezza statistica
    return make_pair(mean, error);
}

double System :: ProbFunction (const double& x) { // Funzione d'onda di test calcolata in x
    return exp(- 0.5 * pow((x - _mu) / _sigma, 2)) + exp(- 0.5 * pow((x + _mu) / _sigma, 2));
}

double System :: DensityFunction (const double& x) { // Densità di probabilità calcolata in x
    return exp(- pow((x - _mu) / _sigma, 2)) + exp(- pow((x + _mu) / _sigma, 2)) + 2.0 * exp(- (pow(x, 2) + pow(_mu, 2)) / pow(_sigma, 2));
}

double System :: Hamiltonian (const double& x) { // Valore dell'hamiltoniana in x
    double kinetic = - 0.5 * exp(- 0.5 * pow((x - _mu) / _sigma, 2)) / pow(_sigma, 2) * (pow((x - _mu) / _sigma, 2) - 1) - 0.5 * exp(- 0.5 * pow((x + _mu) / _sigma, 2)) / pow(_sigma, 2) * (pow((x + _mu) / _sigma, 2) - 1);
    double potential = pow(x, 4) - 2.5 * pow(x, 2);
    return kinetic / ProbFunction(x) + potential;
}

double System :: min (const double& a, const double& b) { // Minimo tra due valori
    if (a <= b) return a;
    else return b; 
}

void System :: Metro () { // Step Metropolis per campionamento densità di probabilità
    double x_new = _rnd.Rannyu(_x - _x_step, _x + _x_step);                     // Proposta nuova x
    double alpha = min(1, DensityFunction(x_new) / DensityFunction(_x));        // Probabilità accettazione
    if (_rnd.Rannyu() <= alpha) {                                               // Metropolis
        _x = x_new;
    }
    _acc_x += alpha;                                                            // Accumulatore accettazione su x
    _actual_step_x ++;                                                          // Accumulatore step Metropolis 
}

pair<double, double> System :: ExpectedHamiltonian () { // Valore di <H> con errore con media a blocchi
    _acc_x = 0.0;               // Reset accumulatore accettazione 
    _actual_step_x = 0;         // Reset numero step Metropolis
    double sum_ave = 0.0;       // Somma dei valori medi dei blocchi
    double sum_ave2 = 0.0;      // Somma dei valori quadratici medi dei blocchi

    for (int i=0; i < _num_block; i++){                 // Loop sul numero di blocchi
        double block_sum = 0.0;
        for (int j=0; j < _dim_block; j++){             // Loop su dimensione blocco
            double expected_H = 0.0;
            for (int k=0; k < _num_x; k++){
                this -> Metro();                        // Generazione coordinata x
                expected_H += this -> Hamiltonian(_x);  // Calcolo valore di H(x)
            }
            block_sum += expected_H / _num_x;           // Valore di <H>
        }
        block_sum /= _dim_block;                        // Valor medio di <H> su un blocco
        sum_ave += block_sum;                           // Accumulatore medie dei blocchi
        sum_ave2 += pow(block_sum, 2);                  // Accumulatore medie quadratiche dei blocchi
    }

    return this -> BlockAverage(sum_ave, sum_ave2, _num_block);     // Media a blocchi su tutti i blocchi
}

void System :: Boltzmann () { // Step Metropolis-Boltzmann per SA
    double mu_old = _mu;                                                    // Vecchio parametro mu
    double sigma_old = _sigma;                                              // Vecchio parametro sigma
    _mu = _rnd.Rannyu(_mu - _mu_step, _mu + _mu_step);                      // Proposta nuovo parametro mu
    _sigma = _rnd.Rannyu(_sigma - _sigma_step, _sigma + _sigma_step);       // Proposta nuovo parametro sigma

    pair<double, double> expect_new = this -> ExpectedHamiltonian();        // Calcolo di <H> con media a blocchi
    double L_new = expect_new.first;
    double L_err_new = expect_new.second;

    double boltz;                                                           // Probabilità accettazione
    if (L_new > _L) boltz = exp(- _beta * (L_new - _L));                    // Metropolis-Boltzmann
    else boltz = 1;
    if (_rnd.Rannyu() < boltz){
        _L = L_new;
        _L_err = L_err_new;
    }
    else{
        _mu = mu_old;
        _sigma = sigma_old;
    }

    _acc_boltz = boltz;                                                     // Accettazione Boltzmann singolo step
}

void System :: SimulatedAnnealing () { // Simulated Annealing al variare di beta
    ofstream out;
    out.open("./OUTPUT/simulated_annealing.out");           // Stampa in output parametri mu, sigma, <H> con errore e accettazioni
    if (out.is_open()){
        out << "--------------------------------------------------" << endl
            << "SIMULATED ANNEALING" << endl
            << "--------------------------------------------------" << endl;
        out << setw(6) << "BETA:" << setw(15) << "MU:" << setw(15) << "SIGMA:" 
            << setw(20) << "<H>:" << setw(20) << "ERR:" 
            << setw(15) << "ACC_X:" << setw(15) << "ACC_BOLTZ:" << endl;
    }
    else {
        cerr << "PROBLEM: Unable to open file ./OUTPUT/simulated_annealing.out" << endl;
        exit(-1);
    }

    for (int i = 0; i < _num_beta; i ++){                   // Loop su numero temperature da simulare
        for (int j=0; j < _num_sim; j++){                   // Loop su numero di simulazioni per fissata temperatura
            this -> Boltzmann();                            // Metropolis-Boltzmann
            out << setw(6) << setprecision(1) << fixed << _beta 
                << setw(15) << setprecision(5) << fixed << _mu
                << setw(15) << setprecision(5) << fixed << _sigma
                << setw(20) << setprecision(10) << fixed << _L
                << setw(20) << setprecision(10) << fixed << _L_err
                << setw(15) << setprecision(5) << fixed << _acc_x / _actual_step_x
                << setw(15) << setprecision(5) << fixed << _acc_boltz << endl;
            }
        _beta += _beta_step;
    }

    out.close();
}

void System :: PrintCostFunction (double mu_min, double mu_max, double sigma_min, double sigma_max, int num_point) { // Stampa in output <H> al variare di mu e sigma
    ofstream out;
    out.open("./OUTPUT/costfunction.out");
    if (out.is_open())
        out << setw(15) << "MU:" << setw(15) << "SIGMA:" << setw(20) << "<H>:" << setw(20) << "ERROR:" << endl;
    else {
        cerr << "PROBLEM: Unable to open file ./OUTPUT/costfunction.out" << endl;
        exit(-1);
    }

    double mu_step = fabs(mu_max - mu_min) / num_point;                     // Numero di mu campionati
    double sigma_step = fabs(sigma_max - sigma_min) / num_point;            // Numero di sigma campionati

    for (int i=0; i <= num_point; i++){
        for (int j=0; j <= num_point; j++){
            _mu = mu_min + i * mu_step;                                     // Valore di mu
            _sigma = sigma_min + j * sigma_step;                            // Valore di sigma
            pair<double, double> expected = this->ExpectedHamiltonian();    // Calcolo di <H>
            out << setw(15) << _mu
                << setw(15) << _sigma
                << setw(20) << fixed << setprecision(7) << expected.first
                << setw(20) << fixed << setprecision(7) << expected.second << endl; 
        }
    }

    out.close();
}

void System :: FindBestHamiltonian (int num_block, int dim_block, int num_x, double mu, double sigma) { // Valore di <H> con errore con media a blocchi per assegnati valori di mu e sigma
    _num_x = num_x;                                 // Numero di valori x per calcolo hamiltoniana     
    _mu = mu;                                       // Media funzione d'onda
    _sigma = sigma;                                 // Dev. stnd. funzione d'onda
    _num_block = num_block;                         // Numero blocchi
    _dim_block = dim_block;                         // Dimensione blocchi

    _acc_x = 0.0;                                   // Reset accumulatore accettazione 
    _actual_step_x = 0;                             // Reset numero step Metropolis
    double sum_ave = 0.0;                           // Somma dei valori medi dei blocchi
    double sum_ave2 = 0.0;                          // Somma dei valori quadratici medi dei blocchi

    ofstream out;
    out.open("./OUTPUT/hamiltonian.out");
    if (out.is_open()){
        out << "--------------------------------------------------" << endl
            << "EXPECTED VALUES HAMILTONIAN" << endl
            << "--------------------------------------------------" << endl;
        out << setw(7) << "#BLOCK:" << setw(15) << "<H>_CUMUL:" << setw(15) << "<H>_ERR:" << endl;
    }
    else {
        cerr << "PROBLEM: Unable to open file ./OUTPUT/hamiltonian.out" << endl;
        exit(-1);
    }

    for (int i=0; i < _num_block; i++){                                 // Loop sul numero di blocchi
        double block_sum = 0.0;
        for (int j=0; j < _dim_block; j++){                             // Loop su dimensione blocco
            double expected_H = 0.0;
            for (int k=0; k < _num_x; k++){
                this -> Metro();                                        // Generazione coordinata x
                expected_H += this -> Hamiltonian(_x);                  // Calcolo valore di H(x)
            }
            block_sum += expected_H / _num_x;                           // Valore di <H>
        }
        block_sum /= _dim_block;                                                        // Valor medio di <H> su un blocco
        sum_ave += block_sum;                                                           // Accumulatore medie dei blocchi
        sum_ave2 += pow(block_sum, 2);                                                  // Accumulatore medie quadratiche dei blocchi
        pair<double, double> mean_err = this -> BlockAverage(sum_ave, sum_ave2, i+1);   // Media a blocchi su tutti i blocchi
        out << setw(7) << i+1 
            << setw(15) << setprecision(7) << fixed << mean_err.first
            << setw(15) << setprecision(7) << fixed << mean_err.second << endl;
    }
    
    out.close(); 
    
}

void System :: PrintDensityFunction (double mu, double sigma, int num) { // Stampa in output valori di |Psi(x)|^2
    ofstream out;
    out.open("./OUTPUT/density_function.out");
    if (out.is_open())
        out << setw(8) << "X_COORD:" << setw(10) << "ACC_X:" << endl;
    else {
        cerr << "PROBLEM: Unable to open file ./OUTPUT/density_function.out" << endl;
        exit(-1);
    }

    _x = 0.0;
    _acc_x = 0.0;
    _mu = mu;
    _sigma = sigma;
    for (int i=0; i < num; i++) {
        this -> Metro();
        out << setw(8) << setprecision(5) << fixed << _x << setw(10) << setprecision(5) << fixed << _acc_x / (i+1) << endl;
    }

   out.close();
}