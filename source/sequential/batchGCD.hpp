
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
        //vector of keys
        vector<mpz_class> keys;
        
        //product tree, initiated in the constructor
        vector< vector<mpz_class> > tree;
        
    public:
        //fill product tree
        void initTree ();

        //change root value (useful for the distributed algorithm)
        void changeValueRoot(mpz_class x);

        //get root value (useful for the distributed algorithm)
        mpz_class getValueRoot();

        //get remainders to find factors in the batchGCD algorithm
        vector< mpz_class > getRemainders ();

        //c++ wrapper to gcd
        mpz_class gcdCPP (mpz_class p1, mpz_class p2);

        //constructor
        BatchGCD (vector<mpz_class>& keys);

        // factorization of the public keys
        vector<Factor> getFactorization ();

        // print keys
        void printKeys ();
};

#endif
