
//gmp library
#include <gmpxx.h>

//factors hpp
#include "factors.hpp"

//constructor
Factor::Factor (long i, mpz_class p, mpz_class q) {
    this->i = i;
    this->p = p;
    this->q = q;
}

//get index
long Factor::getIndex () {
    return this->i;
}

//get factor p
mpz_class Factor::getP() {
    return this->p;
}

//get factor q
mpz_class Factor::getQ() {
    return this->q;
}
