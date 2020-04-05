#ifndef TRANSFAC_H
#define TRANSFAC_H

#include "factory/factory.h" // maybe less?
#include "polys/monomials/ring.h"
#include "polys/monomials/p_polys.h"

class transFac
{
  private:
    CanonicalForm numerator;   // numerator of fraction
    CanonicalForm denominator; // factorized denominator of fraction

  public:
    transFac();
    transFac (long n);

    transFac (CanonicalForm n);
    transFac (CanonicalForm n, CanonicalForm d);

    transFac (poly n, const ring r);
    transFac (poly n, poly d, const ring r);

    transFac (const transFac *l);

    CanonicalForm const& getNum()   const;
    CanonicalForm const& getDenom() const;

    poly getNumPoly (const ring r) const;
    poly getDenomPoly (const ring r) const;

    bool numeratorIsZero() const;
    bool numeratorIsOne() const;
    bool denominatorIsOne() const;

    void negateInplace();

    void normalize();
};

typedef transFac* pTransFac;

// needed at some points, as n_Init does not overload
number nfInit(poly, const coeffs);
int nftIsParam (number, const coeffs);
BOOLEAN n_transFacInitChar (coeffs, void*);

#endif
