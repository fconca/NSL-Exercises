#include "myfunctions.h"

using namespace std;

// Inizializza generatore numeri casuali
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

// Media a blocchi dati vettori medie e medie quadratiche
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

// Legge parametri di input (tipo simulazione e generazione random, step Metropolis, punto partenza, numero e dimensione blocchi)
void ReadParameters(string& simulation_type, string& random_type, double& dim_step, int& num_step_acc, double& xcoord, double& ycoord, double& zcoord, int& num_block, int& dim_block, bool& print) {
    ifstream inp;
    inp.open("./INPUT/input.in");
    if (inp.is_open()){
        string parameter_name;
        while (inp >> parameter_name) {
            if (parameter_name == "#ENDINPUT") 
                break;
            if (parameter_name == "SIMULATION_TYPE")        // "groundstate" o "firstexcited" per stato fondamentale e primo stato eccitato
                inp >> simulation_type;
            else if (parameter_name == "RANDOM_TYPE")       // "uniform" o "normal" per generare nuovo punto con distrib uniforme o normale
                inp >> random_type;
            else if (parameter_name == "DIM_STEP")          // step Metropolis per generare nuovo punto
                inp >> dim_step;
            else if (parameter_name == "NUM_STEP_ACC")      // numero step per calcolo accettazione media
                inp >> num_step_acc;
            else if (parameter_name == "XCOORD")            // 3 coordinate punto di partenza
                inp >> xcoord;
            else if (parameter_name == "YCOORD")
                inp >> ycoord;
            else if (parameter_name == "ZCOORD")
                inp >> zcoord;
            else if (parameter_name == "NUM_BLOCK")         // numero blocchi
                inp >> num_block;
            else if (parameter_name == "DIM_BLOCK")         // dimensione blocco
                inp >> dim_block;
            else if (parameter_name == "PRINT_CONFIG")      // per stampare configurazione
                inp >> print;
        }
    }
    else {
        cerr << "PROBLEM: Unable to open file ./INPUT/input.in" << endl;
        exit(-1);
    }
    inp.close();
}

// Calcola distanza di un punto dall'origine
double dist (const pos& point) {
    return sqrt(pow(point.x, 2) + pow(point.y, 2) + pow(point.z, 2));
}

// Calcola probabilità di una particella di trovarsi in un certo punto per groundstate
double ProbGroundState (const pos& point) {
    return 1.0 / M_PI * exp(-2.0 * dist(point));
}

// Calcola probabilità di una particella di trovarsi in un certo punto per primo stato eccitato
double ProbFirstExcState (const pos& point) {
    double cos_theta = point.z / dist(point);
    return 1.0 / (32.0 * M_PI) * pow(dist(point), 2) * exp(-dist(point)) * pow(cos_theta, 2);
}

// Determina minimo fra due punti
double min (double a, double b) {
    if (a <= b)
        return a;
    else
        return b; 
}

// Genera un nuovo punto a partire da un punto assegnato secondo distrib uniforme centrata in quest'ultimo
pos GetNewPointUnif (Random& rnd, const pos& point, double dim_step) {
    pos new_point;
    new_point.x = rnd.Rannyu(point.x - dim_step, point.x + dim_step);
    new_point.y = rnd.Rannyu(point.y - dim_step, point.y + dim_step);
    new_point.z = rnd.Rannyu(point.z - dim_step, point.z + dim_step);
    return new_point;
}

// Genera un nuovo punto a partire da un punto assegnato secondo distrib gaussiana centrata in quest'ultimo
pos GetNewPointNorm (Random& rnd, const pos& point, double dim_step) {
    pos new_point;
    new_point.x = rnd.Gauss(point.x, dim_step);
    new_point.y = rnd.Gauss(point.y, dim_step);
    new_point.z = rnd.Gauss(point.z, dim_step);
    return new_point;
}

// Calcola parametro di accettazione per groundstate
double GetAlphaGroundState (const pos& old_point, const pos& new_point) {
    double alpha = min(1.0, ProbGroundState(new_point) / ProbGroundState(old_point));
    return alpha;
}

// Calcola parametro di accettazione per primo stato eccitato
double GetAlphaFirstExcState (const pos& old_point, const pos& new_point) {
    double alpha = min(1.0, ProbFirstExcState(new_point) / ProbFirstExcState(old_point));
    return alpha;
}

// Realizza singolo step dell'algoritmo di Metropolis, si possono scrivere su file le posizioni ricoperte
void Metropolis (Random& rnd, pos& point, double dim_step, const string& random_type, const string& simulation_type, bool write_config) {

    static bool stampa_legenda = true;                                              //possibilità di scrivere configurazione su file
    ofstream outc;                     

    if (write_config) {
        if (stampa_legenda) {
            outc.open("./OUTPUT/positions.dat");
            if (outc.is_open()) {
                outc << setw(9) << "#NUM_PART" << setw(9) << "XCOORD" << setw(9) << "YCOORD" << setw(9) << "ZCOORD" << setw(9) << "DIST" << endl;
                stampa_legenda = false;
            }
            else {
                cerr << "PROBLEM: Unable to open file ./OUTPUT/positions.dat" << endl;
                exit(-1);
            }
        }
        outc.open("./OUTPUT/positions.dat", ios::app);
    }

    pos new_point;                                                                  //proposta punto successivo
    if (random_type == "uniform") 
        new_point = GetNewPointUnif (rnd, point, dim_step);
    else if (random_type == "normal") 
        new_point = GetNewPointNorm (rnd, point, dim_step);
    else {
        cerr << "PROBLEM: unknown <random_type>: 'uniform' or 'normal'" << endl;
        exit(-1);
    }

    double alpha;                                                                   //accettazione Metropolis
    if (simulation_type == "groundstate")
        alpha = GetAlphaGroundState (point, new_point);
    else if (simulation_type == "firstexcited")
        alpha = GetAlphaFirstExcState (point, new_point);
    else {
        cerr << "PROBLEM: unknown <simulation_type>: 'groundstate' or 'firstexcited'" << endl;
        exit(-1);
    }

    if (rnd.Rannyu() <= alpha) {                                                    //aggiornamento punto successivo
        point = new_point;
    }

    static int num_part = 0;
    if (write_config) {                                                             //stampa configurazione su file
        outc << setw(8) << num_part 
            << setw(9) << fixed << setprecision(2) << point.x
            << setw(9) << fixed << setprecision(2) << point.y
            << setw(9) << fixed << setprecision(2) << point.z
            << setw(9) << fixed << setprecision(2) << dist(point) << endl;
        num_part ++;
    }
    outc.close();

}

// Calcola l'accettazione media su un ciclo di <num_step> 
double Acceptance (Random& rnd, pos& point, double dim_step, int num_step, const string& random_type, const string& simulation_type) {

    pos new_point;
    double alpha = 0.0;
    double acc = 0.0;

    for (int i = 0; i < num_step; i++){

        if (random_type == "uniform")                                                //proposta punto successivo
            new_point = GetNewPointUnif (rnd, point, dim_step);
        else if (random_type == "normal") 
            new_point = GetNewPointNorm (rnd, point, dim_step);
        else {
            cerr << "PROBLEM: unknown <random_type>: 'uniform' or 'normal'" << endl;
            exit(-1);
        }

        if (simulation_type == "groundstate")                                        //accettazione Metropolis
            alpha = GetAlphaGroundState (point, new_point);
        else if (simulation_type == "firstexcited")
            alpha = GetAlphaFirstExcState (point, new_point);
        else {
            cerr << "PROBLEM: unknown <simulation_type>: 'groundstate' or 'firstexcited'" << endl;
            exit(-1);
        }

        if (rnd.Rannyu() <= alpha) {                                                 //aggiornamento punto successivo
            point = new_point;
        }

        acc += alpha;                                                                //accettazione media

    }

    return acc / num_step;

}