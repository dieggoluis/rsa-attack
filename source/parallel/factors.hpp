
#ifndef FACTORS
#define FACTORS

//gmp library
#include<gmpxx.h>

class Factor {
    private:
        long i;
        mpz_class p, q;
    public:
        //constructor
        Factor (long i, mpz_class p, mpz_class q);
        long getIndex ();
        mpz_class getP();
        mpz_class getQ();
};

#endif
