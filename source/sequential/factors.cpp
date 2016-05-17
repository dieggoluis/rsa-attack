
#include <gmpxx.h>
#include "factors.hpp"

Factor::Factor (long i, mpz_class p, mpz_class q) {
    this->i = i;
    this->p = p;
    this->q = q;
}

long Factor::getIndex () {
    return this->i;
}

mpz_class Factor::getP() {
    return this->p;
}

mpz_class Factor::getQ() {
    return this->q;
}
