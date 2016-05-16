#include <gmpxx.h>
#include <iostream>
#include "batchGCD.hpp"
#include <vector>

using namespace std;

int main () {
    vector< vector<mpz_class> > tree;
    vector<mpz_class> v;
    mpz_class x = 1;
    mpz_class i = 1;
    for (i=0; i<6; i++) {
        mpz_class x;
        cin >> x;
        v.push_back(x);
    }

    BatchGCD* solution = new BatchGCD(v);

//    for (int i = 0; i < tree.size(); i++) {
//        for (int j = 0; j<tree[i].size(); j++)
//            cout << tree[i][j] << " ";
//        cout << "\n";
//    }

    vector<Factor*> factors = solution->getFactorization();

    solution->printKeys();
    for (int i=0; i<factors.size(); i++) {
        long index = factors[i]->getIndex();
        mpz_class p = factors[i]->getP();
        mpz_class q = factors[i]->getQ();
        cout << index << " " << p <<  " " << q << "\n";
    }
    return 0;
}



