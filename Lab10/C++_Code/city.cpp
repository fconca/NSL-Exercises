#include "city.h"

void City :: CityInSquare (System* system) { // Crea una città in un quadrato di lato 2
    _x = system -> GetRandom(-1, 1);
    _y = system -> GetRandom(-1, 1);
}

void City :: CityOnCircum (System* system) { // Crea una città su una circonferenza di raggio 1
    double phi = system -> GetRandom(0, 2.0 * M_PI);
    _x = cos(phi);
    _y = sin(phi);
}

vector<City> City :: AllCities (System* system) { // Vettore di città
    int num_cities = system -> GetNumCities();                  // Numero di città
    int square_circum = system -> GetSquareCircum();            // 0: nel quadrato, 1: sulla circonferenza

    vector<City> cities;
    for (int i = 0; i < num_cities; i++) {
        City city(system);
        if (square_circum == 0) city.CityInSquare(system);      // Crea città nel quadrato
        if (square_circum == 1) city.CityOnCircum(system);      // crea città sulla circonferenza
        cities.push_back(city);
    }

    ofstream out;
    out.open("./OUTPUT/cities.out");
    if (out.is_open())
        out << setw(10) << "#NUM_CITY:" << setw(10) << "X_COORD:" << setw(10) << "Y_COORD:" << endl;
    else {
        cerr << "PROBLEM: Unable to open file ./OUTPUT/cities.out" << endl;
        exit(-1);
    }
    out.close();

    return cities;
}

vector<City> City :: ReadCitiesFromFile (string filename) {
    vector<City> cities;

    ifstream in;
    in.open(filename);
    if (in.is_open()) {
        double x, y;
        while (in >> x >> y) {
            City city;
            city.SetX(x);
            city.SetY(y);
            cities.push_back(city);
        }
    }
    else {
        cerr << "PROBLEM: Unable to open file " << filename << endl;
        exit(-1);
    }
    in.close();

    ofstream out;
    out.open("./OUTPUT/cities.out");
    if (out.is_open())
        out << setw(10) << "#NUM_CITY:" << setw(10) << "X_COORD:" << setw(10) << "Y_COORD:" << endl;
    else {
        cerr << "PROBLEM: Unable to open file ./OUTPUT/cities.out" << endl;
        exit(-1);
    }
    out.close();

    return cities;
}

double City :: GetX () const { // Restituisce coordinata x
    return _x;
}

void City :: SetX (double const& x) {
    _x = x;
}

double City :: GetY () const { // Restituisce coordinata y
    return _y;
}

void City :: SetY (double const& y) {
    _y = y;
}

void City :: PrintOut(vector<City>& cities) { // Stampa in output numero e posizione delle città
    ofstream out;
    out.open("./OUTPUT/cities.out", ios::app);
    if (out.is_open()) {
        for (size_t i = 0; i < cities.size(); i++)
            out << setw(10) << i 
                << setw(10) << fixed << setprecision(3) << cities[i].GetX() 
                << setw(10) << fixed << setprecision(3) << cities[i].GetY() << endl;
    }
    else {
        cerr << "PROBLEM: Unable to open file ./OUTPUT/cities.out" << endl;
        exit(-1);
    }
    out.close();
}