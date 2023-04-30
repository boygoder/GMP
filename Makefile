CXX=g++
CXXFLAGS=-g -std=c++17
GMPLINKS=-lgmpxx -lgmp

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(GMPLINKS)

.PHONY: clean
clean:
	rm -f *.o

