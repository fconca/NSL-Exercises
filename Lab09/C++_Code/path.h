#ifndef __Path__
#define __Path__

#include "city.h"
#include "system.h"
#include <algorithm>
#include <unordered_set>
#include <numeric>

using namespace std;

class Path {
private:
    vector<int> _sequence; // sequenza di interi rappresentativa di un percorso
    double _L; // valore funzione costo (distanza tra città)

public:
    Path () {}
    Path (System* system) {}
    vector<int> GetSequence (); // Restituisce la sequenza rappresentativa di un percorso
    void SetSequence (const vector<int>& sequence); // Imposta la sequenza rappresentativa di un percorso
    double GetCostFunction () const; // Restituisce valore funzione costo per una sequenza
    void CostFunction (System* system, const vector<City>& cities); // Calcolo funzione costo con metrica L1/L2
    void Check(System* system); // Controllo che sequenza soddisfi i vincoli
    static vector<Path> StartingPopulation (System* system, const vector<City>& cities); // Inizializza percorsi
    static void Order (vector<Path>& paths); // Ordina vettore secondo funzione costo in ordine crescente
    static int Selection (System* system, const vector<Path>& paths); // Algoritmo selezione di un percorso, restituisce indice
    void Permutation (System* system); // Permuta due elementi di una sequenza
    void GroupShift (System* system); // Esegue shift di un gruppo di elementi
    void GroupPermutation (System* system); // Permuta due gruppi di elementi
    void Inversion (System* system); // Inverte ordine di un gruppo di elementi
    static void Crossover (System* system, Path& path1, Path& path2); // Crossover di due sequenze
    static void PrintOut (vector<Path>& paths); // Stampa in output sequenza ottimale e funzione costo
};

#endif // __Path__