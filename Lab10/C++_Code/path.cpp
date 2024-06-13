#include "path.h"

using namespace std;

vector<int> Path :: GetSequence () { // Restituisce la sequenza rappresentativa di un percorso
    return _sequence;
}

void Path :: SetSequence (const vector<int>& sequence) { // Imposta la sequenza rappresentativa di un percorso
    _sequence = sequence;
}

double Path :: GetCostFunction () const { // Restituisce valore funzione costo per una sequenza
    return _L;
}

void Path :: CostFunction (System* system, const vector<City>& cities) { // Calcolo funzione costo con metrica L1/L2
    int dim_sequence = _sequence.size();
    int num_cities = cities.size();
    if (dim_sequence != num_cities) { // Controllo che dimensione sequenza percorso = numero città
        cerr << "PROBLEM in function Path :: CostFunction: sequence and cities have different sizes" << endl;
        exit(-1);
    }

    double dist = 0.0;
    for (int i = 0; i < dim_sequence; i++) {
        double x1 = cities[_sequence[i]].GetX();                            // Coordinata x prima città
        double x2 = cities[_sequence[(i+1) % dim_sequence]].GetX();         // Coordinata x seconda città inclusa 0
        double y1 = cities[_sequence[i]].GetY();                            // Coordinata y prima città
        double y2 = cities[_sequence[(i+1) % dim_sequence]].GetY();         // Coordinata y seconda città inclusa 0

        if (system -> GetMetric() == 0)                                     // Usa metrica L1
            dist += sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));                
        else if (system -> GetMetric() == 1)                                // Usa metrica L2
            dist += pow(x1 - x2, 2) + pow(y1 - y2, 2);                      
    }
    
    _L = dist;
}

void Path :: Check (System* system) { // Controllo che sequenza soddisfi i vincoli
    int size = _sequence.size();

    if (size != system -> GetNumCities()) { // Se dimensione vettore è diverso dal numero di città
        cerr << "PROBLEM in function Path :: Check: the size of the sequence has changed" << endl;
        exit(-1);
    }

    if (_sequence[0] != 0) { // Se primo elemento della sequenza è diverso da 0
        cerr << "PROBLEM in function Path :: Check: the first element of the sequence must be 0" << endl;
        exit(-1);
    }

    vector<int> sorted_sequence = _sequence;
    sort(sorted_sequence.begin(), sorted_sequence.end()); // Copia ordinata della sequenza
    for (int i = 0; i < size; ++i) 
        if (sorted_sequence[i] != i) { // Se elementi da 0 a size non tutti presenti nella sequenza
            cerr << "PROBLEM: all items in each individual must be different" << endl;
            exit(-1);
        }
}

vector<Path> Path :: StartingPopulation (System* system, const vector<City>& cities) { // Inizializza percorsi
    int num_cities = cities.size();                                 // Numero di città
    int num_path = system -> GetNumPaths();                         // Numero percorsi 
    vector<Path> paths(num_path);                                   // Vettore di percorsi

    vector<int> sequence(num_cities);
    iota(sequence.begin(), sequence.end(), 0);                      // Sequenza 0, 1, 2, ... 
    paths[0].SetSequence(sequence);                                 // Inizializza primo percorso

    for (int i = 1; i < num_path; i++) {
        random_shuffle(sequence.begin() + 1, sequence.end());       // Sequenza ordine random con 0 fissato
        paths[i].SetSequence(sequence);                             // Imposta percorso
        paths[i].Check(system);                                     // Controllo vincoli
    }

    for (int i = 0; i < num_path; i++)
        paths[i].CostFunction(system, cities);                      // Calcolo funzione costo   

    Order(paths);                                                   // Ordina vettore percorsi secondo la funzione costo

    ofstream out;
    out.open("./OUTPUT/costfunction.out"); // Per scrivere valore funzione costo in output
    if (out.is_open())
        out << setw(14) << "COST_FUNCTION:" << setw(25) << "MEAN_COST_FUNCTION:" << endl;
    else {
        cerr << "PROBLEM: Unable to open file ./OUTPUT/costfunction.out" << endl;
        exit(-1);
    }   
    out.close();

    out.open("./OUTPUT/sequences.out"); // Per scrivere sequenza ottimale in output
    if (!out.is_open()) {
        cerr << "PROBLEM: Unable to open file ./OUTPUT/sequences.out" << endl;
        exit(-1);
    }
    out.close();

    return paths;
}

void Path :: Order (vector<Path>& paths) { // Ordina vettore secondo funzione costo in ordine crescente
    sort(paths.begin(), paths.end(), [](const Path& a, const Path& b) {
        return a.GetCostFunction() < b.GetCostFunction();
    });
}

int Path :: Selection (System* system, const vector<Path>& paths) { // Algoritmo selezione di un percorso, restituisce indice
    double exp = system -> GetExponent();                   // Esponente usato da legge estrazione indice
    int num_rows = paths.size();
    double rand = system -> GetRandom(0, 1);
    int index = int((num_rows - 1) * pow(rand, exp)) + 1;   // Estrai indice random secondo la legge specificata
    return index;
}

void Path :: Permutation (System* system) { // Permuta due elementi di una sequenza
    int vec_size = _sequence.size();
    if (vec_size < 3) { // Controllo dimensione maggiore di 3
        cerr << "PROBLEM in function Path :: Permutation: vector size must be at least 3" << endl;
        exit(-1);
    }
        
    for (int i = 1; i < vec_size; i++) {                            // Scorri su tutti elementi sequenza
        if (system -> GetRandom(0,1) <= system -> GetProbPerm()) {  // Esegui con data probabilità
            int index = i;
            while (index == i)                                      // Controllo indici diversi
                index = int(system -> GetRandom(1, vec_size));      // Indice random

            swap(_sequence[i], _sequence[index]);                   // Permutazione dei due elementi
        }
    }

    this -> Check(system); // Controllo vincoli su nuova sequenza
}

void Path :: GroupShift (System* system) { // Esegue shift di un gruppo di elementi
    int vec_size = _sequence.size();

    if (vec_size < 3) { // Controllo dimensione maggiore di 3
        cerr << "PROBLEM in function Path :: GroupShift: vector size must be at least 3" << endl;
        exit(-1);
    }

    for (int i = 1; i < vec_size; i++) {                                                        // Scorri su tutti elementi sequenza
        if (system -> GetRandom(0, 1) <= system -> GetProbShiftGroup()) {                       // Esegui con data probabilità
            int dim_group = int(system -> GetRandom(1, vec_size - 1));                          // Dimensione gruppo da shiftare
            int shift = int(system -> GetRandom(1, vec_size - dim_group));                      // Dimensione shift

            vector<int> duplicate = _sequence;
            duplicate.insert(duplicate.end(), _sequence.begin() + 1, _sequence.end());          // Seguenza raddoppiata escluso 0

            vector<int> group(dim_group + shift);              
            for (int j = 0; j < dim_group + shift; j++) 
                group[j] = duplicate[i + j];                                    // Copia del gruppo da shiftare
            rotate(group.begin(), group.begin() + dim_group, group.end());      // Esegui shift sulla copia del gruppo

            for (int j = 0; j < dim_group + shift; j++) {                       // Shift del gruppo riportato nella sequenza originaria        
                if (i + j < vec_size)
                    _sequence[i + j] = group[j];
                else if (i + j >= vec_size)
                    _sequence[(i + j) % vec_size + 1] = group[j];               
            }

            this -> Check(system); // Controllo vincoli su nuova sequenza
        }
    }
}

void Path :: GroupPermutation(System* system) { // Permuta due gruppi di elementi
    int vec_size = _sequence.size();

    if (vec_size < 3) { // Controllo dimensione maggiore di 3
        cerr << "PROBLEM in function Path :: GroupPermutation: vector size must be at least 3" << endl;
        exit(-1);
    }

    for (int i = 1; i < vec_size; i++) {                                                        // Scorri su tutti elementi sequenza
        if (system -> GetRandom(0, 1) <= system -> GetProbPermGroup()) {                        // Esegui con data probabilità
            int dim_group = int(system -> GetRandom(1, (vec_size + 1) / 2));                    // Dimensione del gruppo da permutare

            vector<int> duplicate = _sequence;
            duplicate.insert(duplicate.end(), _sequence.begin() + 1, _sequence.end());          // Seguenza raddoppiata escluso 0

            int index = int(system -> GetRandom(i + dim_group, i + vec_size - dim_group));      // Indice secondo gruppo da permutare

            vector<int> first_group(dim_group);         
            vector<int> second_group(dim_group);        
            for (int j = 0; j < dim_group; j++) {
                first_group[j] = duplicate[i + j];          // Copia del primo gruppo da permutare
                second_group[j] = duplicate[index + j];     // Copia del secondo gruppo da permutare
            }

            for (int j = 0; j < dim_group; j++) {           // Permutazione dei gruppi riportata nella sequenza originaria
                if (i + j < vec_size)
                    _sequence[i + j] = second_group[j];
                else if (i + j >= vec_size)
                    _sequence[(i + j) % vec_size + 1] = second_group[j];
                if (index + j < vec_size)
                    _sequence[index + j] = first_group[j];  
                else if (index + j >= vec_size)
                _sequence[(index + j) % vec_size + 1] = first_group[j];
            }

            this -> Check(system); // Controllo vincoli su nuova sequenza
        }
    }
}

void Path :: Inversion (System* system) { // Inverte ordine di un gruppo di elementi
    int vec_size = _sequence.size();

    if (vec_size < 3) { // Controllo dimensione maggiore di 3
        cerr << "PROBLEM in function Path :: Inversion: vector size must be at least 3" << endl;
        exit(-1);
    }

    for (int i = 1; i < vec_size; i++) {                                                        // Scorri su tutti elementi sequenza
        if (system -> GetRandom(0, 1) <= system -> GetProbInversion()) {                        // Esegui con data probabilità
            int dim_group = int(system -> GetRandom(1, vec_size));                              // Dimensione del gruppo da invertire

            vector<int> duplicate = _sequence;
            duplicate.insert(duplicate.end(), _sequence.begin() + 1, _sequence.end());          // Seguenza raddoppiata escluso 0

            vector<int> group(dim_group); 
            for (int j = 0; j < dim_group; j++) 
                group[j] = duplicate[i + j];            // Copia del gruppo da invertire
            reverse(group.begin(), group.end());        // Inversione del gruppo

            for (int j = 0; j < dim_group; j++) {       // Inversione del gruppo riportata nella sequenza originaria
                if (i + j < vec_size)
                    _sequence[i + j] = group[j];
                else if (i + j >= vec_size)
                    _sequence[(i + j) % vec_size + 1] = group[j];
            }

            this -> Check(system); // Controllo vincoli su nuova sequenza
        }
    }
}

void Path :: Crossover(System* system, Path& path1, Path& path2) { // Crossover di due sequenze
    vector<int> vec1 = path1.GetSequence();     // Prima sequenza
    vector<int> vec2 = path2.GetSequence();     // Seconda sequenza
    int dim = vec1.size(); 

    if (vec1.size() != vec2.size()) { // Controlla che sequenze abbiano stessa dimensione
        cerr << "PROBLEM in function Path :: CrossOver: the two vectors must have the same size" << endl;
        exit(-1);
    }

    if (system -> GetRandom(0, 1) <= system -> GetProbCross()) {                // Esegui con data probabilità
        int index = system -> GetRandom(2, dim - 1);                            // Indice di taglio

        vector<int> first_part_vec1(vec1.begin(), vec1.begin() + index);        // Sequenza 1 prima del taglio
        vector<int> first_part_vec2(vec2.begin(), vec2.begin() + index);        // Sequenza 2 prima del taglio    

        unordered_set<int> present_vec1(first_part_vec1.begin(), first_part_vec1.end());      // Elementi già presenti in first_part_vec1
        unordered_set<int> present_vec2(first_part_vec2.begin(), first_part_vec2.end());      // Elementi già presenti in first_part_vec2

        for (int i = 0; i < dim; ++i) {
            if (present_vec1.find(vec2[i]) == present_vec1.end()) {     // Se elemento non presente nella nuova sequenza
                first_part_vec1.push_back(vec2[i]);                     // Aggiungi elemento mancante in fondo alla sequenza 
                present_vec1.insert(vec2[i]);                           // Aggiungi elemento all'insieme degli elementi già presenti
            }
        }

        for (int i = 0; i < dim; ++i) {
            if (present_vec2.find(vec1[i]) == present_vec2.end()) {     // Se elemento non presente nella nuova sequenza
                first_part_vec2.push_back(vec1[i]);                     // Aggiungi elemento mancante in fondo alla sequenza 
                present_vec2.insert(vec1[i]);                           // Aggiungi elemento all'insieme degli elementi già presenti
            }
        }

        path1.SetSequence(first_part_vec1);                             // Setta nuova sequenza 1
        path2.SetSequence(first_part_vec2);                             // Setta nuova sequenza 2
        path1.Check(system);                                            // Controllo su nuovi vettori
        path2.Check(system);
    }
}

void Path :: PrintOut (vector<Path>& paths) { // Stampa in output sequenza ottimale e funzione costo
    double acc_cost_func = 0.0;                                     // Accumulatore funzione costo

    for (size_t i = 0; i < paths.size() / 2; i++) 
        acc_cost_func += paths[i].GetCostFunction();
    double mean_cost_func = acc_cost_func / (paths.size() / 2);     // Funzione costo media calcolata su metà della popolazione

    ofstream out;
    out.open("./OUTPUT/costfunction.out", ios::app); // File output per funzione costo
    if (out.is_open())
        out << setw(13) << fixed << setprecision(7) << paths[0].GetCostFunction()
            << setw(25) << fixed << setprecision(7) << mean_cost_func << endl;
    else {
        cerr << "PROBLEM: Unable to open file ./OUTPUT/costfunction.out" << endl;
        exit(-1);
    }   
    out.close();

    out.open("./OUTPUT/sequences.out", ios::app); // File output per sequenza ottimale
    if (out.is_open()) {
        for (size_t i = 0; i < paths[0].GetSequence().size(); i++)
            out << setw(4) << paths[0].GetSequence()[i];
        out << endl;
    }
    else {
        cerr << "PROBLEM: Unable to open file ./OUTPUT/sequences.out" << endl;
        exit(-1);
    }  
}