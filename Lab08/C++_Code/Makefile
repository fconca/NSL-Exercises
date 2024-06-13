# Specifica il compilatore
CXX = g++

# Opzioni del compilatore
CXXFLAGS = -std=c++11 -Wall

# Nome dell'eseguibile
EXECUTABLE = Exercise8.exe

# File sorgente (usiamo una wildcard per includere tutti i file .cpp)
SOURCES = $(wildcard *.cpp)

# File oggetto (sostituiamo l'estensione .cpp con .o)
OBJECTS = $(SOURCES:.cpp=.o)

# File header (usiamo una wildcard per includere tutti i file .h)
HEADERS = $(wildcard *.h)

# Regola principale: compilazione dell'eseguibile
$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Regola generica per la compilazione dei file oggetto
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Pulisce i file oggetto e l'eseguibile
clean:
	rm -f $(EXECUTABLE) $(OBJECTS)