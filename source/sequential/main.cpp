#include <gmpxx.h>
#include <iostream>
#include <vector>
#include <stdlib.h>

#include "batchGCD.hpp"

using namespace std;

int main (int argc, char* argv[]) {
    
    if (argc < 2) {
        cout << "too few arguments to program\n";
        return 1;
    }

    long dataSize = atol(argv[1]);

    vector<mpz_class> input;

    for (long i = 0; i < dataSize; i++) {
        int index;
        cin >> index;
        mpz_class x;
        cin >> x;
        input.push_back(x);
    }

    BatchGCD* solution = new BatchGCD(input);

    vector<Factor> factors = solution->getFactorization();

    //solution->printKeys();
    for (int i=0; i<factors.size(); i++) {
        long index = factors[i].getIndex();
        mpz_class p = factors[i].getP();
        mpz_class q = factors[i].getQ();
        //if (p!= 1 && q!=1)
            cout << index << " " << p <<  " || " << q << "\n";
    }

    delete solution;
    return 0;
}



