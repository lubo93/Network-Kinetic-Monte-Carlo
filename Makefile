CXX = g++
CXXFLAGS = -std=gnu++11 -O3 -Wall -D_NDEBUG

all: dynamicNet

clean:
	rm *.o dynamicNet

dynamicNet: main.o kmc.o network.o
	$(CXX) $(CXXFLAGS) -o $@ main.o kmc.o network.o

%.o: %.cpp *.hpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<