#ifndef TRANSFAC_H
#define TRANSFAC_H

#include "factory/factory.h" // maybe less?
#include "polys/monomials/ring.h"
#include "polys/monomials/p_polys.h"

#include <list>

class transFac
{
  private:
    ring          f_r;         // polynomial ring where coefficient polynomials live
    CanonicalForm numerator;   // numerator of fraction
    CanonicalForm denominator; // factorized denominator of fraction
    int           complexity;

  public:
    transFac (const ring r);
    transFac (long n, const ring r);

    transFac (CanonicalForm n, const ring r);
    transFac (CanonicalForm n, CanonicalForm, int c, const ring r);

    transFac (poly n, const ring r);
    transFac (poly n, poly d, int c, const ring r);

    transFac (const transFac *l);

    ring getRing() const;
    CanonicalForm const& getNum()   const;
    CanonicalForm const& getDenom() const;
    int getComp() const;

    poly getNumPoly() const;
    poly getDenomPoly() const;

    bool numeratorIsZero() const;
    bool numeratorIsOne() const;
    bool denominatorIsOne() const;

    void negateInplace();

    void normalize(bool);
};

typedef transFac* pTransFac;

BOOLEAN n_transFacInitChar(coeffs, void*);

#endif
