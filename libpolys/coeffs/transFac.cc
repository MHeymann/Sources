#include <coeffs/transFac.h>

#include "coeffs/coeffs.h"
#include "coeffs/numbers.h"
#include "omalloc/omalloc.h"
#include "misc/auxiliary.h"

#define TRANSEXT_PRIVATES
#include "polys/ext_fields/transext.h"
#include "coeffs/longrat.h"

#include "Singular/links/ssiLink.h"

#include "polys/clapsing.h"
#include "polys/clapconv.h"
#include "polys/prCopy.h"

// member functions
transFac::transFac (const ring r) :
  f_r (r),
  numerator (CanonicalForm (0)),
  denominator (CanonicalForm (1))
{}

transFac::transFac (long n, const ring r) :
  f_r (r),
  numerator (CanonicalForm (n)),
  denominator (CanonicalForm (1))
{}

transFac::transFac (CanonicalForm n, const ring r) :
  f_r (r),
  numerator (n),
  denominator (CanonicalForm (1))
{}

transFac::transFac (CanonicalForm n, CanonicalForm d, const ring r) :
  f_r (r),
  numerator (n),
  denominator (d)
{}

transFac::transFac (poly n, const ring r) :
  f_r (r),
  denominator (CanonicalForm (1))
{
  numerator = convSingPFactoryP (n, f_r);
  p_Delete (&n, f_r);
}

transFac::transFac (poly n, poly d, const ring r) :
  f_r (r)
{
  numerator = convSingPFactoryP (n, f_r);
  p_Delete (&n, f_r);
  denominator = convSingPFactoryP (d, f_r);
  p_Delete (&d, f_r);
}

transFac::transFac (const transFac *l) :
  f_r (l->getRing()),
  numerator (l->getNum()),
  denominator (l->getDenom())
{}

ring transFac::getRing() const
{
  return f_r;
}

CanonicalForm const& transFac::getNum() const
{
  return numerator;
}

CanonicalForm const& transFac::getDenom() const
{
  return denominator;
}

poly transFac::getNumPoly() const
{
  return convFactoryPSingP (numerator, f_r);
}

poly transFac::getDenomPoly() const
{
  return convFactoryPSingP (denominator, f_r);
}

bool transFac::numeratorIsZero() const
{
  return numerator.isZero();
}

bool transFac::numeratorIsOne() const
{
  return numerator.isOne();
}

bool transFac::denominatorIsOne() const
{
  return denominator.isOne();
}

void transFac::negateInplace()
{
  // TODO probably easier
  numerator = -numerator;
}

void transFac::normalize ()
{
  if (numeratorIsZero())
  {
    denominator = CanonicalForm (1);
    return;
  }

  Off (SW_RATIONAL);
  CanonicalForm G = gcd (numerator, denominator);
  if (getCharacteristic() == 0)
  {
    On (SW_RATIONAL);
  }

  if (!G.isOne())
  {
    numerator   /= G;
    denominator /= G;
  }

  if (getCharacteristic() == 0)
  {
    CanonicalForm denN = bCommonDen (numerator);
    CanonicalForm denD = bCommonDen (denominator);
    if (!denN.isOne())
      numerator   *= denN;
    if (!denD.isOne())
      denominator *= denD;

    Off (SW_RATIONAL);
    CanonicalForm gcdDen = gcd (denN, denD);
    if (!gcdDen.isOne())
    {
      denN /= gcdDen;
      denD /= gcdDen;
    }
    On (SW_RATIONAL);
    if (!denN.isOne())
      denominator *= denN;
    if (!denD.isOne())
      numerator   *= denD;

    Off (SW_RATIONAL);
    CanonicalForm gcon = gcd (icontent (numerator), icontent (denominator));
    if (!gcon.isOne())
    {
      numerator   /= gcon;
      denominator /= gcon;
    }
  }

  // FIXME: canonicalform.lc() is in general not the lc w.r.t. ring
  // order of parameters (thus obtain different term when printing)
  if (denominator.lc() < 0)
  {
    numerator   *= -1;
    denominator *= -1;
  }
}

// coef field functions

// proper copying is done by constructor
static number nfCopy (number f, const coeffs cf)
{
  n_Test(f, cf);
  pTransFac ff = new transFac ((pTransFac) f);
  n_Test ((number) ff, cf);
  return (number) ff;
}

static void nfDelete (number *f, coeffs cf)
{
  n_Test(*f, cf);
  delete (pTransFac) (*f);
  *f = NULL;
}

static void nfNormalize (number &f, const coeffs cf)
{
  n_Test(f, cf);
  ((pTransFac) f)->normalize();
}

static BOOLEAN nfIsZero (number f, const coeffs cf)
{
  n_Test(f, cf);
  return ((pTransFac) f)->numeratorIsZero();
}

static BOOLEAN nfIsOne (number f, const coeffs cf)
{
  n_Test(f, cf);
  pTransFac ff = (pTransFac) f;
  ff->normalize ();
  return (ff->numeratorIsOne() && ff->denominatorIsOne());
}

// is -1?
static BOOLEAN nfIsMOne (number f, const coeffs cf)
{
  n_Test(f, cf);
  pTransFac ff = (pTransFac) f;
  ff->normalize ();

  if (ff->numeratorIsZero() || !ff->denominatorIsOne())
  {
    return FALSE;
  }

  const CanonicalForm& num = ff->getNum();
  if (!num.inCoeffDomain())
  {
    return FALSE;
  }

  return (-num).isOne();
}

static BOOLEAN nfGreaterZero (number f, const coeffs cf)
{
  n_Test(f, cf);
  pTransFac ff = (pTransFac) f;

  if (ff->numeratorIsZero())
  {
    return FALSE;
  }

  ff->normalize();
  const CanonicalForm& num = ff->getNum();

  return num > 0;
}

static BOOLEAN nfGreater (number f, number g, const coeffs cf)
{
  n_Test(f, cf);
  n_Test(g, cf);

  pTransFac ff = (pTransFac) f;
  pTransFac gg = (pTransFac) g;

  if (ff->numeratorIsZero())
  {
    // have to exclude case g==0
    if (gg->numeratorIsZero())
    {
      return FALSE;
    }
    return !nfGreaterZero(g, cf);
  }

  if (gg->numeratorIsZero())
  {
    return nfGreaterZero(f, cf);
  }

  // now have f!=0, g!=0
  const CanonicalForm& fnum = ff->getNum();
  const CanonicalForm& fden = ff->getDenom();
  CanonicalForm fNumCoeff   = fnum.lc();
  int           fNumDeg     = totaldegree (fnum);

  CanonicalForm fDenCoeff   = fden.lc();
  int           fDenDeg     = totaldegree (fden);

  const CanonicalForm& gnum = gg->getNum();
  const CanonicalForm& gden = gg->getDenom();
  CanonicalForm gNumCoeff   = gnum.lc();
  int           gNumDeg     = totaldegree (gnum);

  CanonicalForm gDenCoeff   = gden.lc();
  int           gDenDeg     = totaldegree (gden);

  // in case of degree inequality lead coefficients are unnecessary
  if (fNumDeg - fDenDeg != gNumDeg - gDenDeg)
  {
    return (fNumDeg - fDenDeg > gNumDeg - gDenDeg);
  }

  // degrees are equal, have to compare coefficients
  CanonicalForm fg = fNumCoeff * gDenCoeff,
                gf = gNumCoeff * fDenCoeff;

  return fg > gf;
}

static BOOLEAN nfEqual (number f, number g, const coeffs cf)
{
  n_Test(f, cf);
  n_Test(g, cf);

  // simple tests
  if (f == g)
  {
    return TRUE;
  }

  pTransFac ff = (pTransFac) f;
  pTransFac gg = (pTransFac) g;

  // check whether one is zero but not the other
  bool fz = ff->numeratorIsZero();
  bool gz = gg->numeratorIsZero();

  // if at least one number is zero, they are equal iff both are zero
  if (fz || gz)
  {
    return fz && gz;
  }

  if (ff->denominatorIsOne() && gg->denominatorIsOne())
  {
    return ff->getNum() == gg->getNum();
  }

  // else do product test: fn/fd == gn/gd <=> fn*gd == gn*fd
  // TODO Optimize
  return ff->getNum() * gg->getDenom() == gg->getNum() * ff->getDenom();
}

number nfInit (long i, const coeffs cf)
{
  pTransFac f = new transFac (i, cf->extRing);
  n_Test ((number) f, cf);
  return (number) f;
}

// takes p, does not copy
number nfInit (poly p, const coeffs cf)
{
  if (p != NULL)
  {
    p_Test (p, cf->extRing);
    pTransFac f = new transFac (p, cf->extRing);
    n_Test ((number) f, cf);
    return (number) f;
  }
  return (number) new transFac (cf->extRing);
}

static number nfGetNumerator (number &f, const coeffs cf)
{
  n_Test (f, cf);
  pTransFac ff = (pTransFac) f;
  if (ff->numeratorIsZero())
  {
    return nfInit(0l, cf);
  }

  ff->normalize ();
  pTransFac num = new transFac(ff->getNum(), cf->extRing);
  n_Test((number) num, cf);
  return (number) num;
}

static number nfGetDenom (number &f, const coeffs cf)
{
  n_Test (f, cf);
  pTransFac ff = (pTransFac) f;
  if (ff->numeratorIsZero())
  {
    return nfInit(1l, cf);
  }

  ff->normalize ();
  pTransFac num = new transFac(ff->getDenom(), cf->extRing);
  n_Test((number) num, cf);
  return (number) num;
}

static long nfInt(number &f, const coeffs cf)
{
  n_Test (f, cf);

  pTransFac ff = (pTransFac) f;
  if (ff->numeratorIsZero())
  {
    return 0;
  }

  ff->normalize();
  if (!ff->denominatorIsOne())
  {
    return 0;
  }

  const CanonicalForm& num = ff->getNum();

  if (!num.inCoeffDomain())
  {
    return 0;
  }

  return num.intval();
}

static number nfNeg (number f, const coeffs cf)
{
  n_Test(f, cf);
  pTransFac ff = (pTransFac) f;
  if (!ff->numeratorIsZero())
  {
    ff->negateInplace();
  }
  n_Test (f, cf);
  return f;
}

static number nfAdd (number f, number g, const coeffs cf)
{
  n_Test(f, cf);
  n_Test(g, cf);

  pTransFac ff = (pTransFac) f,
            gg = (pTransFac) g;

  // check if one summand is zero
  if (ff->numeratorIsZero())
  {
    return nfCopy (g, cf);
  }
  if (gg->numeratorIsZero())
  {
    return nfCopy(f, cf);
  }

  const CanonicalForm& df = ff->getDenom();
  const CanonicalForm& dg = gg->getDenom();
  const CanonicalForm& nf = ff->getNum();
  const CanonicalForm& ng = gg->getNum();

  CanonicalForm newNum, newDen;

  if (df == dg)
  {
    newNum = nf + ng;
    newDen = df;
  }
  else if (df.isOne())
  {
    newNum = nf * dg + ng;
    newDen = dg;
  }
  else if (dg.isOne())
  {
    newNum = nf + ng*df;
    newDen = df;
  }
  else
  {
    CanonicalForm gd = gcd (df, dg);
    if (gd.isOne())
    {
      newNum = nf * dg + ng * df;
      newDen = df * dg;
    }
    else
    {
      CanonicalForm qf = df / gd,
                    qg = dg / gd;
      newNum = qg * nf + qf * ng;
      newDen = qg * df;
    }
  }

  pTransFac sum = new transFac (newNum, newDen, cf->extRing);
  sum->normalize();
  n_Test ((number) sum, cf);
  return (number) sum;
}

static number nfSub (number f, number g, const coeffs cf)
{
  n_Test(f, cf);
  n_Test(g, cf);

  pTransFac ff = (pTransFac) f,
            gg = (pTransFac) g;

  // check if one summand is zero
  if (ff->numeratorIsZero())
  {
    return nfNeg (nfCopy (g, cf), cf);
  }
  if (gg->numeratorIsZero())
  {
    return nfCopy (f, cf);
  }

  const CanonicalForm& df = ff->getDenom();
  const CanonicalForm& dg = gg->getDenom();
  const CanonicalForm& nf = ff->getNum();
  const CanonicalForm& ng = gg->getNum();

  CanonicalForm newNum, newDen;

  if (df == dg)
  {
    newNum = nf - ng;
    newDen = df;
  }
  else if (df.isOne())
  {
    newNum = nf * dg - ng;
    newDen = dg;
  }
  else if (dg.isOne())
  {
    newNum = nf - ng*df;
    newDen = df;
  }
  else
  {
    CanonicalForm gd = gcd (df, dg);
    if (gd.isOne())
    {
      newNum = nf * dg - ng * df;
      newDen = df * dg;
    }
    else
    {
      CanonicalForm qf = df / gd,
                    qg = dg / gd;
      newNum = qg * nf - qf * ng;
      newDen = qg * df;
    }
  }

  pTransFac diff = new transFac (newNum, newDen, cf->extRing);
  diff->normalize();
  n_Test ((number) diff, cf);
  return (number) diff;
}

static number nfMult (number f, number g, const coeffs cf)
{
  n_Test(f, cf);
  n_Test(g, cf);

  pTransFac ff = (pTransFac) f,
            gg = (pTransFac) g;

  if (ff->numeratorIsZero() || gg->numeratorIsZero())
  {
    return nfInit (0l, cf);
  }

  CanonicalForm N, D;

  CanonicalForm const& fn = ff->getNum();
  CanonicalForm const& fd = ff->getDenom();
  CanonicalForm const& gn = gg->getNum();
  CanonicalForm const& gd = gg->getDenom();

  if (fd.isOne())
  {
    CanonicalForm G = gcd (fn, gd);
    if (G.isOne())
    {
      N = fn * gn;
      D = gd;
    }
    else
    {
      N = (fn / G) * gn;
      D = gd / G;
    }
  }
  else if (gd.isOne())
  {
    CanonicalForm G = gcd (gn, fd);
    if (G.isOne())
    {
      N = fn * gn;
      D = fd;
    }
    else
    {
      N = fn * (gn / G);
      D = fd / G;
    }
  }
  else
  {
    CanonicalForm G1 = gcd (fn, gd);
    CanonicalForm G2 = gcd (fd, gn);
    if (G1.isOne())
    {
      N = fn;
      D = gd;
    }
    else
    {
      N = fn / G1;
      D = gd / G1;
    }
    if (G2.isOne())
    {
      N *= gn;
      D *= fd;
    }
    else
    {
      N *= gn / G2;
      D *= fd / G2;
    }
  }

  pTransFac prod = new transFac (N, D, cf->extRing);
  prod->normalize();
  n_Test ((number) prod, cf);
  return (number) prod;
}

static number nfDiv (number f, number g, const coeffs cf)
{
  n_Test(f, cf);
  n_Test(g, cf);

  pTransFac ff = (pTransFac) f,
            gg = (pTransFac) g;

  if (gg->numeratorIsZero())
  {
    WerrorS (nDivBy0);
  }
  if (ff->numeratorIsZero())
  {
    return nfInit (0l, cf);
  }

  CanonicalForm N, D;

  CanonicalForm const& fn = ff->getNum();
  CanonicalForm const& fd = ff->getDenom();
  // flip gg here
  CanonicalForm const& gd = gg->getNum();
  CanonicalForm const& gn = gg->getDenom();


  // from here identical to nfMult
  if (fd.isOne())
  {
    CanonicalForm G = gcd (fn, gd);
    if (G.isOne())
    {
      N = fn * gn;
      D = gd;
    }
    else
    {
      N = (fn / G) * gn;
      D = gd / G;
    }
  }
  else if (gd.isOne())
  {
    CanonicalForm G = gcd (gn, fd);
    if (G.isOne())
    {
      N = fn * gn;
      D = fd;
    }
    else
    {
      N = fn * (gn / G);
      D = fd / G;
    }
  }
  else
  {
    CanonicalForm G1 = gcd (fn, gd);
    CanonicalForm G2 = gcd (fd, gn);
    if (G1.isOne())
    {
      N = fn;
      D = gd;
    }
    else
    {
      N = fn / G1;
      D = gd / G1;
    }
    if (G2.isOne())
    {
      N *= gn;
      D *= fd;
    }
    else
    {
      N *= gn / G2;
      D *= fd / G2;
    }
  }

  pTransFac div = new transFac (N, D, cf->extRing);
  div->normalize();
  n_Test ((number) div, cf);
  return (number) div;
}

static number nfInvers (number f, const coeffs cf)
{
  pTransFac ff = (pTransFac) f;

  if (ff->numeratorIsZero())
  {
    WerrorS (nDivBy0);
  }

  pTransFac inv = new transFac (ff->getDenom(), ff->getNum(), cf->extRing);
  inv->normalize();
  n_Test ((number) inv, cf);
  return (number) inv;
}

static void nfPower (number f, int exp, number *p, const coeffs cf)
{
  n_Test(f, cf);
  pTransFac ff = (pTransFac) f;

  if (ff->numeratorIsZero())
  {
    // get here that 0^0 == 1, but eh ...
    if (exp >= 0)
    {
      *p = nfInit (0l, cf);
    }
    else
    {
      WerrorS (nDivBy0);
    }
  }

  // might reduce amount of work to do
  ff->normalize ();

  // destroy F later on
  number F;
  if (exp < 0)
  {
    exp = -exp;
    F = nfInvers (f, cf);
  }
  else
  {
    F = nfCopy (f, cf);
  }

  pTransFac FF = (pTransFac) F;

  CanonicalForm npow = power (FF->getNum(), exp);
  CanonicalForm dpow = power (FF->getDenom(), exp);

  nfDelete(&F, cf);
  pTransFac pow = new transFac (npow, dpow, cf->extRing);
  *p = (number) pow;
  n_Test(*p, cf);
}

// TODO
number nfFarey (number f, number n, const coeffs cf)
{
  // n is bigint
  extern coeffs coeffs_BIGINT;
  pTransFac ff = (pTransFac) f;
  // does this work?
  CanonicalForm N = n_convSingNFactoryN (n, TRUE, coeffs_BIGINT);
  CanonicalForm num = Farey (ff->getNum(), N);

  CanonicalForm denom = Farey (ff->getDenom(), N);

  return (number) new transFac(num, denom, cf->extRing);
}

static int nfSize (number f, const coeffs cf)
{
  n_Test(f, cf);
  pTransFac ff = (pTransFac) f;

  if (ff->numeratorIsZero())
  {
    return 0;
  }

  unsigned long size_n = 0UL, size_d = 0UL;

  poly num = ff->getNumPoly();
  for (poly it = num; it != NULL; pIter (it))
  {
    number co = p_GetCoeff (it, cf->extRing);
    size_n += (unsigned long) cf->extRing->cf->cfSize (co, cf->extRing->cf) + 1UL;
  }
  p_Delete (&num, cf->extRing);

  poly den = ff->getDenomPoly();
  for (poly it = den; it != NULL; pIter (it))
  {
    number co = p_GetCoeff (it, cf->extRing);
    size_d += (unsigned long) cf->extRing->cf->cfSize (co, cf->extRing->cf) + 1UL;
  }
  p_Delete (&den, cf->extRing);

  unsigned long t = size_n + size_d + 1;

  if (t > INT_MAX)
  {
    return INT_MAX;
  }
  else
  {
    return (int) t;
  }
}

#if 1
static number nfMap00 (number a, const coeffs src, const coeffs dst)
{
  n_Test(a, src);

  if (n_IsZero (a, src))
  {
    return nfInit ((poly) NULL, dst);
  }
  assume (src->rep == dst->extRing->cf->rep);

  On (SW_RATIONAL);
  CanonicalForm N = n_convSingNFactoryN (a, TRUE, src);
  Off (SW_RATIONAL);

  // normalization will take care of nested fractions later on
  number res = (number) new transFac (N, dst->extRing);

  n_Test(res, dst);
  return res;
}

static number nfMapPP(number a, const coeffs src, const coeffs dst)
{
  n_Test(a, src) ;
  if (n_IsZero(a, src))
  {
    return nfInit ((poly) NULL, dst);
  }
  assume(src == dst->extRing->cf);

  number res = (number) new transFac (n_convSingNFactoryN (a, TRUE, src), dst->extRing);
  n_Test(res, dst);
  return res;
}

// copy between same coeffs
// TODO probably very inefficient
static number nfCopyMap (number f, const coeffs cf, const coeffs dst)
{
  n_Test (f, cf);
  pTransFac ff = (pTransFac) f;

  const ring rSrc = cf->extRing;
  const ring rDst = dst->extRing;

  if (ff->numeratorIsZero())
  {
    pTransFac cp = new transFac (rDst);
    n_Test ((number) cp, dst);
    return (number) cp;
  }

  if (rSrc == rDst)
  {
    return nfCopy (f, dst);
  }

  poly num  = ff->getNumPoly();
  poly nNew = prCopyR (num, rSrc, rDst);
  p_Delete (&num, rSrc);

  if (ff->denominatorIsOne())
  {
    return nfInit(nNew, dst);
  }

  poly d = ff->getDenomPoly();
  poly dNew = prCopyR (d, rSrc, rDst);
  p_Delete (&d, rSrc);

  pTransFac fNew = new transFac (nNew, dNew, rDst);
  n_Test ((number) fNew, dst);

  return (number) fNew;
}

number nfSubMap (number f, const coeffs src, const coeffs dst)
{
  n_Test (f, src);
  pTransFac ff = (pTransFac) f;

  int save = getCharacteristic();
  //const ring rSrc = src->extRing;
  const ring rDst = dst->extRing;

  setCharacteristic (dst->ch);
  pTransFac result = new transFac
    ( ff->getNum().mapinto()
    , ff->getDenom().mapinto()
    , rDst
    ); // does this work?

  setCharacteristic (save);
  n_Test ((number) result, dst);
  return (number) result;
}

// TODO
// mapping from transExt into transFac
number nfTransMap (number f, const coeffs src, const coeffs dst)
{
  n_Test (f, src);
  fraction ff = (fraction) f;

  const ring rSrc = src->extRing;
  const ring rDst = dst->extRing;


  const nMapFunc nMap = n_SetMap (rSrc->cf, rDst->cf);
  poly num  = NUM (ff);
  poly g = prMapR (num, nMap, rSrc, rDst);
  p_Delete (&num, rSrc);

  /* g may contain summands with coeff 0 */
  poly hh = g;
  poly prev = NULL;
  while(hh != NULL)
  {
    if (n_IsZero(pGetCoeff(hh), rDst->cf))
    {
      if (prev == NULL)
      {
        g = p_LmFreeAndNext(g, rDst);
        hh = g;
      }
      else
      {
        prev->next = p_LmFreeAndNext(prev->next, rDst);
        hh = prev->next;
      }
    }
    else
    {
      prev = hh;
      pIter(hh);
    }
  }

  if (g == NULL)
  {
    return nfInit((poly) NULL, dst);
  }

  poly denNew;

  if (!DENIS1(ff))
  {
    poly d = DEN(ff);

    poly h = prMapR(d, nMap, rSrc, rDst);
    p_Delete(&d, rSrc);

    /* h may contain summands with coeff 0 */
    hh = h;
    prev = NULL;
    while(hh != NULL)
    {
      if (n_IsZero(pGetCoeff(hh), rDst->cf))
      {
        if (prev == NULL)
        {
          h = p_LmFreeAndNext(h, rDst);
          hh = h;
        }
        else
        {
          prev->next = p_LmFreeAndNext(prev->next, rDst);
          hh = prev->next;
        }
      }
      else
      {
        prev = hh;
        pIter(hh);
      }
    }
    if (h == NULL)
    {
      WerrorS("mapping to */0");
    }
    denNew = h;
  }
  else
  {
    denNew = p_ISet (1, rDst);
  }

  pTransFac result = new transFac (g, denNew, rDst);
  // denominator must be factorized
  result->normalize ();
  n_Test ((number) result, dst);
  return (number) result;
}

nMapFunc nfSetMap (const coeffs src, const coeffs dst)
{
  // easiest case:
  if (src == dst)
  {
    return ndCopyMap;
  }

  if (src->extRing == NULL)
  {
    if ((src->rep==n_rep_gap_rat) && nCoeff_is_Q(dst->extRing->cf))
    {
      return nfMap00;
    }
    if (nCoeff_is_Zp(src) && nCoeff_is_Zp(dst->extRing->cf))
    {
      if (src->ch == dst->ch)
      {
        return nfMapPP;
      }
      return NULL;
    }
  }

  // case where we also map from facDem coeffs
  if (getCoeffType(src) == getCoeffType(dst))
  {
    // if src has more parameters than dst, then we can't do anything:
    // for any substitution there are numbers which will be sent to poly/0
    if (rVar(src->extRing) > rVar(dst->extRing))
    {
      return NULL;
    }

    // pars should match by name
    for (int i = 0; i < rVar(src->extRing); i++)
    {
      if (strcmp(rRingVar(i, src->extRing), rRingVar(i, dst->extRing)) != 0)
      {
        return NULL;
      }
    }

    if (src->extRing->cf==dst->extRing->cf)
    {
      return nfCopyMap;
    }
    else
    {
      return nfSubMap;
    }
  }
  if (nCoeff_is_transExt(src))
  {
    if (rVar(src->extRing) > rVar(src->extRing))
    {
      return NULL;
    }

    // pars should match by name
    for (int i = 0; i < rVar(src->extRing); i++)
    {
      if (strcmp(rRingVar(i, src->extRing), rRingVar(i, dst->extRing)) != 0)
      {
        return NULL;
      }
    }
    return nfTransMap;
  }

  return NULL;
}
#endif // map stuff

// read as poly in extRing
static const char* nfRead (const char *s, number *f, const coeffs cf)
{
  poly p;
  setCharacteristic (cf->ch); // important
  const char * result = p_Read(s, p, cf->extRing);
  if (p == NULL)
  {
    *f = (number) new transFac(cf->extRing);
  }
  else
  {
    *f = nfInit(p, cf);
  }
  n_Test(*f, cf);
  return result;
}

static void nfWriteLong (number f, const coeffs cf)
{
  n_Test(f, cf);
  pTransFac ff = (pTransFac) f;
  ff->normalize();

  if (ff->numeratorIsZero())
  {
    StringAppendS("0");
  }
  else
  {
    poly fn = ff->getNumPoly();
    BOOLEAN omitParenth = p_IsConstant (fn, cf->extRing);

    if (!omitParenth) { StringAppendS("("); }
    p_String0Long(fn, cf->extRing, cf->extRing);
    if (!omitParenth) { StringAppendS(")"); }
    p_Delete(&fn, cf->extRing);

    if (!ff->denominatorIsOne())
    {
      StringAppendS("/");
      poly fd = ff->getDenomPoly();
      omitParenth = p_IsConstant (fd, cf->extRing);

      if (!omitParenth) { StringAppendS("("); }
      p_String0Long(fd, cf->extRing, cf->extRing);
      if (!omitParenth) { StringAppendS(")"); }
      p_Delete(&fd, cf->extRing);

      p_Delete (&fd, cf->extRing);
    }
  }
}

static void nfWriteShort(number f, const coeffs cf)
{
  n_Test(f, cf);
  pTransFac ff = (pTransFac) f;
  ff->normalize();

  if (ff->numeratorIsZero())
  {
    StringAppendS("0");
  }
  else
  {
    poly fn = ff->getNumPoly();
    BOOLEAN omitParenth = p_IsConstant (fn, cf->extRing);

    if (!omitParenth) { StringAppendS("("); }
    p_String0Short(fn, cf->extRing, cf->extRing);
    if (!omitParenth) { StringAppendS(")"); }
    p_Delete(&fn, cf->extRing);

    if (!ff->denominatorIsOne())
    {
      StringAppendS("/");
      poly fd = ff->getDenomPoly();
      omitParenth = p_IsConstant (fd, cf->extRing);

      if (!omitParenth) { StringAppendS("("); }
      p_String0Short(fd, cf->extRing, cf->extRing);
      if (!omitParenth) { StringAppendS(")"); }
      p_Delete(&fd, cf->extRing);

      p_Delete (&fd, cf->extRing);
    }
  }
}

// correct?
static number nfConvFactoryNSingN(CanonicalForm n, const coeffs cf)
{
  if (n.isZero())
  {
    return nfInit(0l, cf);
  }
  pTransFac p = new transFac(n, cf->extRing);
  n_Test((number) p, cf);
  return (number) p;
}

static CanonicalForm nfConvSingNFactoryN(number n, BOOLEAN setChar, const coeffs cf)
{
  n_Test(n, cf);
  pTransFac nn = (pTransFac) n;
  if (setChar)
  {
    setCharacteristic (cf->ch);
  }

  CanonicalForm nnum = nn->getNum();
  return nnum;
}

static void nfKillChar(const coeffs cf)
{
  if ((--cf->extRing->ref) == 0)
  {
    rDelete(cf->extRing);
  }
}

// adapted from algExt.cc
// prints transFac(cf->ch,"par1",...,"parN")
char* nfCoeffString(const coeffs cf)
{
  const char* const* p=n_ParameterNames(cf);
  int l = 0;
  int i;
  for(i = 0; i < n_NumberOfParameters(cf); i++)
  {
    l += (strlen(p[i]) + 1);
  }
  char *s = (char *)omAlloc(10+l+10+1);
  s[0] = '\0';
  snprintf(s, 10+10+1, "%s%d", "transFac(", cf->ch); /* Fp(a) or Q(a) */
  char tt[2];
  tt[0] = ',';
  tt[1] = '\0';
  for(i = 0; i < n_NumberOfParameters(cf); i++)
  {
    strcat(s, tt);
    strcat(s, p[i]);
  }
  strcat(s, ")");
  return s;
}

char* nfCoeffName(const coeffs cf)
{
  const char* const* p=n_ParameterNames(cf);
  int i;
  static char s[200];
  s[0] = '\0';
  snprintf(s, 10+10+1, "%s%d", "transFac(", cf->ch); /* Fp(a) or Q(a) */
  char tt[2];
  tt[0] = ',';
  tt[1] = '\0';
  for(i = 0; i < n_NumberOfParameters(cf); i++)
  {
    strcat(s, tt);
    strcat(s, p[i]);
  }
  strcat(s, ")");
  return s;
}

static void nfCoeffWrite(const coeffs cf, BOOLEAN details)
{
  const ring R = cf->extRing;
  n_CoeffWrite(R->cf, details);

  const int N = rVar(R);
  assume( N > 0 );

  PrintS("(");

  for (int nop=0; nop < N; nop ++)
  {
    Print("%s", rRingVar(nop, R));
    if (nop != N-1)
    {
      PrintS(",");
    }
  }

  PrintS(") in factory representation");
}

static void nfWriteFd(number f, const ssiInfo* d, const coeffs cf)
{
  // format: <numerator> <denominator>
  n_Test(f, cf);
  pTransFac ff = (pTransFac) f;
  // force normalization
  ff->normalize ();

  poly num = ff->getNumPoly();
  ssiWritePoly_R (d, POLY_CMD, num, cf->extRing);
  p_Delete (&num, cf->extRing);

  poly den = ff->getDenomPoly();

  ssiWritePoly_R (d, POLY_CMD, den, cf->extRing);
  p_Delete (&den, cf->extRing);
}

number nfReadFd(const ssiInfo* d, const coeffs cf)
{
  poly num = ssiReadPoly_R(d, cf->extRing);
  poly denom = ssiReadPoly_R(d, cf->extRing);

  pTransFac F = new transFac(num, denom, cf->extRing);
  n_Test((number) F, cf);
  return (number) F;
}

static number nfParameter(const int iParameter, const coeffs cf)
{
  const ring R = cf->extRing;
  assume( R != NULL );
  assume( 0 < iParameter && iParameter <= rVar(R) );

  pTransFac p = new transFac(CanonicalForm(Variable(iParameter)), cf->extRing);

  n_Test((number) p, cf);
  return (number) p;
}

static int nfParDeg (number f, const coeffs cf)
{
  n_Test (f, cf);
  pTransFac ff = (pTransFac) f;
  return totaldegree (ff->getNum());
}

static number nfGcd(number f, number g, const coeffs cf)
{
  pTransFac ff = (pTransFac) f,
            gg = (pTransFac) g;

  pTransFac p = new transFac(gcd(ff->getNum(),gg->getNum()),cf->extRing);
  return (number)p;
}

static number nfNormalizeHelper(number f, number g, const coeffs cf)
{
  /* NUM(f)*DEN(g)/gcd(NUM(f),DEN(g)) */
  pTransFac ff = (pTransFac) f,
            gg = (pTransFac) g;
  if (gg->denominatorIsOne()) return nfCopy(f,cf);
  CanonicalForm p = ff->getNum();
  CanonicalForm const& d = gg->getDenom();
  CanonicalForm GCD = gcd (p, d);

  pTransFac res=new transFac ((p/GCD)*d, cf->extRing);
  return (number) res;
}

static BOOLEAN nfCoeffIsEqual(const coeffs cf, n_coeffType n, void *param)
{
  if (getCoeffType(cf) != n)
  {
    return FALSE;
  }
  ring r = (ring) param;
  if (cf->extRing == r)
    return TRUE;

  // NOTE: Q(a)[x] && Q(a)[y] should better share the _same_ Q(a)...
  if( rEqual(cf->extRing, r, TRUE) )
  {
    rDelete(r);
    return TRUE;
  }

  return FALSE;
}

// auxillary stuff
int nftIsParam (number f, const coeffs cf)
{
  pTransFac ff = (pTransFac) f;
  if (!ff->denominatorIsOne())
  {
    return 0;
  }
  poly fn = ff->getNumPoly();
  int v = p_Var (fn, cf->extRing);
  p_Delete (&fn, cf->extRing);
  return v;
}

// initialization function, assume we are given the polynomial ring in which base
// polys shall lie
BOOLEAN n_transFacInitChar(coeffs cf, void* parInfo)
{
  assume(parInfo != NULL);
  ring extRing = (ring) parInfo;

  assume( extRing         != NULL );
  assume( extRing->cf     != NULL );
  assume( extRing->qideal == NULL );

  assume( cf != NULL );

  extRing->ref++;

  cf->extRing          = extRing;
  cf->ch               = extRing->cf->ch;
  cf->factoryVarOffset = extRing->cf->factoryVarOffset + rVar(extRing);
  cf->has_simple_Alloc = FALSE;

  cf->is_field  = TRUE;
  cf->is_domain = TRUE;
  cf->rep       = n_rep_transFac;

  cf->cfCopy         = nfCopy;
  cf->cfDelete       = nfDelete;
  cf->cfNormalize    = nfNormalize;

  cf->cfRePart       = nfCopy;

  cf->cfIsZero       = nfIsZero;
  cf->cfIsOne        = nfIsOne;
  cf->cfIsMOne       = nfIsMOne;
  cf->cfGreater      = nfGreater;
  cf->cfGreaterZero  = nfGreaterZero;
  cf->cfEqual        = nfEqual;

  cf->cfInit         = nfInit;
  cf->cfInt          = nfInt;

  cf->cfGetNumerator = nfGetNumerator;
  cf->cfGetDenom     = nfGetDenom;

  cf->cfInpNeg       = nfNeg;
  cf->cfAdd          = nfAdd;
  cf->cfSub          = nfSub;
  cf->cfMult         = nfMult;
  cf->cfDiv          = nfDiv;
  cf->cfExactDiv     = nfDiv;
  cf->cfInvers       = nfInvers;
  cf->cfPower        = nfPower;

  cf->cfFarey        = nfFarey;
  cf->cfSize         = nfSize;

  cf->cfSetMap       = nfSetMap;
  cf->cfRead         = nfRead;
  cf->cfWriteLong    = nfWriteLong;
  cf->cfWriteShort   = rCanShortOut( cf->extRing ) ? nfWriteShort : nfWriteLong;

  cf->convFactoryNSingN = nfConvFactoryNSingN;
  cf->convSingNFactoryN = nfConvSingNFactoryN;

  cf->cfKillChar     = nfKillChar;

  cf->cfCoeffString  = nfCoeffString;
  cf->cfCoeffName    = nfCoeffName;
  cf->cfCoeffWrite   = nfCoeffWrite;

  cf->cfWriteFd      = nfWriteFd;
  cf->cfReadFd       = nfReadFd;

  cf->iNumberOfParameters = rVar(extRing);
  cf->pParameterNames     = (const char**) extRing->names;
  cf->cfParameter         = nfParameter;
  cf->cfParDeg            = nfParDeg;
  cf->has_simple_Inverse  = FALSE;

  cf->cfSubringGcd        = nfGcd;
  cf->cfNormalizeHelper   = nfNormalizeHelper;

  cf->nCoeffIsEqual  = nfCoeffIsEqual;
  // copied from transext.cc: write functions step by step
  /*
  cf->cfChineseRemainder = ntChineseRemainder;
#ifdef LDEBUG
  cf->cfDBTest       = ntDBTest;
#endif



  if( nCoeff_is_Q(R->cf) )
    cf->cfClearContent = ntClearContent;

  cf->cfClearDenominators = ntClearDenominators;
  */

  return FALSE;
}
// vim: spelllang=en
