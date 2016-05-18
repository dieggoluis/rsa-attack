
#ifndef BATCH_GCDH
#define BATCH_GCDH

//gmp library
#include<gmpxx.h>

//STL 
#include<vector>

#include "factors.hpp"

#define debug(x) { std::cout << #x << ": "  << x << std::endl; } 

using namespace std;

class BatchGCD {
    private:
        vector<mpz_class> keys; //vector of keys
    
        vector< vector<mpz_class> > productTree ();
        vector< mpz_class > getRemainders ();
        mpz_class gcdCPP (mpz_class p1, mpz_class p2);
    public:
        //constructor
        BatchGCD (vector<mpz_class>& keys);
        // factorization of the public keys
        vector<Factor> getFactorization ();
        // print keys
        void printKeys ();
};

#endif
