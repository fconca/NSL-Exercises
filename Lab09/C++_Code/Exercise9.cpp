#include "city.h"
#include "system.h"
#include "path.h"

using namespace std;

// --------------------------------------------------------------------------------------------------------------------
// Come usare il programma:
// - Creare un oggetto di tipo System per gestire l'inizializzazione del generatore di numeri random e per la lettura 
//   di tutti i parametri di input della simulazione;
// - Creare un vettore di oggetti di tipo City che racchiude l'insieme delle città;
// - Inizializzare la popolazione iniziale, ovvero le sequenze di città progressivamente visitate e rappresentate dalla 
//   classe Path. Su queste performare l'algoritmo genetico utilizzando i metodi della classe Path.
// --------------------------------------------------------------------------------------------------------------------

int main() {

    System SYS;
    SYS.InitializeRandomGenerator();    // Inizializza generatore di numeri random
    SYS.Initialize();                   // Inizializzazione dei parametri della simulazione

    vector<City> cities = City :: AllCities(&SYS);  // Insieme delle città da visitare
    City :: PrintOut(cities);                       // Stampa in output le città e le rispettive coordinate

    vector<Path> paths = Path :: StartingPopulation(&SYS, cities);  // Inizializza popolazione iniziale

    for (int i = 0; i < SYS.GetNumStep(); i++) {    // Quante volte realizzare lo step genetico
        vector<Path> new_paths(paths.size());       // Vettore per salvare nuova popolazione dopo lo step genetivo

        for (size_t i = 0; i < paths.size() / 2; i++) {

            int index1 = Path :: Selection(&SYS, paths);                        // Selezione indice 1 per crossover
            int index2 = index1;
            while (index1 == index2) index2 = Path :: Selection(&SYS, paths);   // Selezione indice 2 per crossover

            Path new_path1 = paths[index1];                 // Copia della sequenza 1 selezionata
            Path new_path2 = paths[index2];                 // Copia della sequenza 2 selezionata
            Path :: Crossover(&SYS, new_path1, new_path2);  // Crossover (con probabilità esecuzione inclusa)

            new_path1.Permutation(&SYS);                    // Mutazioni genetiche su sequenza 1 (con probabilità)
            new_path1.GroupShift(&SYS);
            new_path1.GroupPermutation(&SYS);
            new_path1.Inversion(&SYS);
            new_path2.Permutation(&SYS);                    // Mutazioni genetiche su sequenza 2 (con probabilità)
            new_path2.GroupShift(&SYS);
            new_path2.GroupPermutation(&SYS);
            new_path2.Inversion(&SYS);

            new_path1.CostFunction(&SYS, cities);           // Calcola funzione costo su nuova sequenza 1
            new_path2.CostFunction(&SYS, cities);           // Calcola funzione costo su nuova sequenza 2

            new_paths[2*i] = new_path1;                     // Inserisci sequenza 1 nella nuova popolazione
            new_paths[2*i + 1] = new_path2;                 // Inserisci sequenza 2 nella nuova popolazione
        }

        paths = move(new_paths);        // Copia nuova popolazione nel vettore originario
        Path :: Order(paths);           // Orinda vettore popolazione in base a funzione costo
        Path :: PrintOut(paths);        // Stampa in output funzione costo e sequenza ottimali
    }
    
    return 0;
}

