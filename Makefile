CXX=g++
CXXFLAGS=-g -std=c++17
GMPLINKS=-lgmpxx -lgmp

#%: %.cpp
#	$(CXX) $(CXXFLAGS) -o $@ $^ $(GMPLINKS)
test_gmp: test_gmp.cpp *.h
	$(CXX) $(CXXFLAGS) test_gmp.cpp -o test_gmp $(GMPLINKS)


test_Gauss: test_Gauss.cpp *.h
	$(CXX) $(CXXFLAGS) test_Gauss.cpp -o test_Gauss $(GMPLINKS)

test_Newton: test_Newton.cpp *.h
	$(CXX) $(CXXFLAGS) test_Newton.cpp -o test_Newton $(GMPLINKS)

test_polynomial: test_polynomial.cpp *.h
	$(CXX) $(CXXFLAGS) test_polynomial.cpp -o test_polynomial $(GMPLINKS)

.PHONY: clean
clean:
	rm -f *.o

