
#include<gmpxx.h>

#include <iostream>
#include<vector>
#include <math.h>

#include "batchGCD.hpp"
#include "factors.hpp"

using namespace std;


vector< vector<mpz_class> > BatchGCD::productTree () {
    vector< vector<mpz_class> > tree;
    
    vector<mpz_class> v (this->keys);
    tree.push_back(v);
    long N = keys.size();
    long h = 0; //height tree

    while (N > 0) {
        vector<mpz_class> x;
        x.resize(N % 2 ? N/2 + 1 : N/2);
        long k = 0;
        for (long i = 0; i < tree[h].size(); i += 2) {
            if (i < tree[h].size() - 1) x[k++] = tree[h][i] * tree[h][i+1];
            else x[k++] = tree[h][i];
        }
        tree.push_back(x);
        h++; N /= 2;
    }
    return tree;
}

vector<mpz_class> BatchGCD::getRemainders () {
    vector < vector<mpz_class> > tree = productTree ();
    for (long i = tree.size() - 2; i >= 0; i--) {
        for (long j = 0; j < tree[i].size(); j++) {
                tree[i][j] = tree[i+1][j/2] % (tree[i][j] * tree[i][j]);
        }
    }
    return tree[0];
    //debug
    //for (int i=0; i<remainders.size(); i++)
    //    cout << remainders[i] << " ";
    //cout << endl;
}

//constructor
BatchGCD::BatchGCD (vector<mpz_class>& keys) {
    this->keys = keys;
}

// factorization of the public keys
vector<Factor*> BatchGCD::getFactorization () {
    vector <Factor*> factors;
    factors.resize(this->keys.size());
    vector<mpz_class> remainders = getRemainders();
    for (long i = 0; i < remainders.size(); i++) {
        mpz_class N = this->keys[i];
        mpz_class p = gcd(N, remainders[i]/N);
        mpz_class q = N/p;
        factors[i] = new Factor(i, p, q);
    }
    return factors;
}

void BatchGCD::printKeys () {
    for (long i = 0; i < this->keys.size(); i++)
        cout << this->keys[i] << " ";
    cout << "\n";
}
