#ifndef __City__
#define __City__

#include "system.h"
#include <vector>

using namespace std;

class City {
private:
    double _x; // Coordinata x
    double _y; // Coordinata y

public:
    City() {}
    City(System* system) {}
    void CityInSquare (System* system); // Crea una città in un quadrato di lato 2
    void CityOnCircum (System* system); // Crea una città su una circonferenza di raggio 1
    static vector<City> AllCities (System* system); // Vettore di città
    double GetX () const; // Restituisce coordinata x
    double GetY () const; // Restituisce coordinata y
    static void PrintOut (vector<City>& cities); // Stampa in output numero e posizione delle città
};

#endif // __City__