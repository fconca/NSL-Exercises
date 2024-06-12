# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -Wall -O3 --std=c++11

# Source files
SRCS = Exercise1_1.cpp random.cpp myfunctions.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Header files
INCS = $(wildcard *.h)

# Executable name
EXE = Exercise1_1.exe

# Rule to compile .cpp files to .o files
%.o: %.cpp $(INCS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to build the executable
$(EXE): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(EXE)

# Clean rule
clean:
	rm -f $(OBJS) $(EXE)
