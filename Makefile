mycprog: mycprog.c
	gcc -g mycprog.c -o mycprog -lgmp

mycxxprog: mycxxprog.cpp
	g++ -g -std=c++17 mycxxprog.cpp -o mycxxprog -lgmpxx -lgmp

.Phony: clean
clean:
	rm mycprog mycxxprog
