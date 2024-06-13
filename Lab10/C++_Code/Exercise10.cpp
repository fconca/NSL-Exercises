#include "city.h"
#include "system.h"
#include "path.h"
#include "mpi.h"

using namespace std;

// --------------------------------------------------------------------------------------------------------------------
// Come usare il programma:
// - Creare un oggetto di tipo System per gestire l'inizializzazione del generatore di numeri random e per la lettura 
//   di tutti i parametri di input della simulazione;
// - Creare un vettore di oggetti di tipo City che racchiude l'insieme delle città;
// - Inizializzare la popolazione iniziale, ovvero le sequenze di città progressivamente visitate e rappresentate dalla 
//   classe Path. Su queste performare l'algoritmo genetico utilizzando i metodi della classe Path.
// --------------------------------------------------------------------------------------------------------------------

int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    System SYS;
    SYS.InitializeRandomGenerator(rank); // Inizializza generatore di numeri random (diverso per ogni nodo)
    SYS.Initialize(); // Inizializzazione dei parametri della simulazione

    vector<City> cities = City :: ReadCitiesFromFile("./INPUT/cap_prov_ita.dat"); // Leggi città da file
    vector<Path> paths = Path :: StartingPopulation(&SYS, cities); // Inizializza popolazione iniziale

    for (int i = 0; i < SYS.GetNumStep(); i++) { // Quante volte realizzare lo step genetico
        vector<Path> new_paths(paths.size()); // Vettore per salvare nuova popolazione dopo lo step genetico

        for (size_t i = 0; i < paths.size() / 2; i++) {

            int index1 = Path :: Selection(&SYS, paths); // Selezione indice 1 per crossover
            int index2 = index1;
            while (index1 == index2) index2 = Path :: Selection(&SYS, paths); // Selezione indice 2 per crossover

            Path new_path1 = paths[index1]; // Copia della sequenza 1 selezionata
            Path new_path2 = paths[index2]; // Copia della sequenza 2 selezionata
            Path :: Crossover(&SYS, new_path1, new_path2); // Crossover (con probabilità esecuzione inclusa)

            new_path1.Permutation(&SYS); // Mutazioni genetiche su sequenza 1 (con probabilità)
            new_path1.GroupShift(&SYS);
            new_path1.GroupPermutation(&SYS);
            new_path1.Inversion(&SYS);
            new_path2.Permutation(&SYS); // Mutazioni genetiche su sequenza 2 (con probabilità)
            new_path2.GroupShift(&SYS);
            new_path2.GroupPermutation(&SYS);
            new_path2.Inversion(&SYS);

            new_path1.CostFunction(&SYS, cities); // Calcola funzione costo su nuova sequenza 1
            new_path2.CostFunction(&SYS, cities); // Calcola funzione costo su nuova sequenza 2

            new_paths[2*i] = new_path1; // Inserisci sequenza 1 nella nuova popolazione
            new_paths[2*i + 1] = new_path2; // Inserisci sequenza 2 nella nuova popolazione
        }

        paths = move(new_paths); // Copia nuova popolazione nel vettore originario
        Path :: Order(paths); // Orinda vettore popolazione in base a funzione costo



	    if (i % SYS.GetNumExchange() == 0) { // Scambio dei migliori individui tra continenti
            int num_exchange = int(paths.size() * 0.2); // Quanti individui scambiare

            vector<int> order(size); // Ordine casuale dei nodi
            if (rank == 0) { // Solo nodo 0 genera sequenza
                for (int k = 0; k < size; ++k) {
                    order[k] = k; // Sequenza 0, 1, 2, ...
                }
                random_shuffle(order.begin(), order.end()); // Sequenza casuale
            }

            MPI_Bcast(&order[0], size, MPI_INT, 0, MPI_COMM_WORLD); // Condivisione dell'ordine con tutti i nodi

            for (int j = 0; j < num_exchange; j++) {
                vector<int> sequence = paths[j].GetSequence(); // Ottieni sequenza del percorso
                MPI_Request send_request, recv_request;

                int current_index = -1;
                for (int k = 0; k < size; k++) {
                    if (order[k] == rank) {
                        current_index = k; // Determinazione dell'indice della sequenza casuale
                        break;
                    }
                }

                int source_index = (current_index + size - 1) % size; // Indice del nodo da cui riceve
                int destination_index = (current_index + 1) % size; // Indice del nodo a cui invia

                int source_rank = order[source_index]; // Nodo mittente
                int destination_rank = order[destination_index]; // Nodo destinatario

                MPI_Isend(&sequence[0], sequence.size(), MPI_INT, destination_rank, j, MPI_COMM_WORLD, &send_request); // Spedisci sequenza a nodo destinatario

                vector<int> received_sequence(cities.size()); // Sequenza ricevuta
                MPI_Irecv(&received_sequence[0], received_sequence.size(), MPI_INT, source_rank, j, MPI_COMM_WORLD, &recv_request); // Ricevi sequenza da nodo mittente

                MPI_Wait(&send_request, MPI_STATUS_IGNORE); // Attendi che processo invio completo
                MPI_Wait(&recv_request, MPI_STATUS_IGNORE); // Attendi che processo ricezione completo

                Path received_path; // Percordo ricevuto
                received_path.SetSequence(received_sequence); // Assegna sequenza ricevuta al percorso
                received_path.CostFunction(&SYS, cities); // Calcola funzione costo per il nuovo percorso
                paths.push_back(received_path); // Colloca nuovo percorso in fondo
            }

            Path::Order(paths); // Ordina vettore paths in base a funzione costo
            paths.erase(paths.end() - num_exchange, paths.end()); // Rimuovi percorsi in eccesso

            for (int r = 0; r < size; r++) {
                if (rank == r)
                    Path :: PrintOut(paths); // Ogni processo stampa su file
                MPI_Barrier(MPI_COMM_WORLD); // Stampa ordinata dai nodi
            }
        }
    }
    
    MPI_Finalize();
    return 0;
}   