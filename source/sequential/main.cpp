//gmp library
#include <gmpxx.h>

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include <ctime>

#include "batchGCD.hpp"

using namespace std;

//read input and return a vector of mpz_class 
vector<mpz_class> readInput (const char *filename) {
    vector<mpz_class> input;
    ifstream file(filename);
    int index;
    while (file >> index) {
        mpz_class x;
        file >> x;
        input.push_back(x);
    }
    return input;
}

int main (int argc, char* argv[]) {
    
    //argv[1] is the number of line to be read
    if (argc < 3) {
        cout << "too few arguments to program\n";
        return 1;
    }

    //number of lines
    long dataSize = atol(argv[1]);

    //input vector
    vector<mpz_class> input;

    //read input
    input = readInput(argv[2]);

    clock_t begin = clock();
    //Pointer to object solution
    BatchGCD* solution = new BatchGCD(input);
    
    //vector of factors. Each Factor has the fields (Index, P, Q)
    vector<Factor> factors = solution->getFactorization();
    clock_t end = clock();
    double exec_secs = double(end - begin) / CLOCKS_PER_SEC;

    //print N, P and Q
    for (int i=0; i<factors.size(); i++) {
        long index = factors[i].getIndex();
        mpz_class p = factors[i].getP();
        mpz_class q = factors[i].getQ();
        mpz_class N = p * q;
        //if (p!= 1 && q!=1)
        cout << "N = " << N << " P = " << p << " Q = " << q << "\n";
    }

    cout << "TIME EXECUTION " << exec_secs << "\n";
    delete solution;
    return 0;
}



