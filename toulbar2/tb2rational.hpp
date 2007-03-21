/** \file tb2rational.hpp
 *  \brief Rational numbers with basic operations.
 * 
 */

#ifndef TB2RATIONAL_HPP_
#define TB2RATIONAL_HPP_

#include <gmp.h>

struct Rational {
  mpq_t rational;     // the number

  Rational() {
    mpq_init(rational);
  }

  Rational(Long p_, Long q_) {
    assert(q_ != 0);
    mpq_init(rational);
    mpq_set_si(rational, p_, q_);
    mpq_canonicalize(rational);
  }

  Rational(Long p_) {
    mpq_init(rational);
    mpq_set_si(rational, p_, 1);
    mpq_canonicalize(rational);
  }

  Rational(const Rational &r) {
    mpq_init(rational);
    mpq_set(rational, r.rational);
  }

  ~Rational () {
    mpq_clear(rational);
  }

  double to_double() const {
    return mpq_get_d(rational);
  }

  Long getNumerator() const {
    return mpz_get_si(mpq_numref(rational));
  }

  Long getDenominator() const {
    return mpz_get_si(mpq_denref(rational));
  }

  Rational &operator=(const Rational &r) {
    mpq_set(rational, r.rational);
    return *this;
  }
  Rational &operator+=(const Rational &r) {
    mpq_add(rational, rational, r.rational);
    return *this;
  }
  Rational &operator-=(const Rational &r) {
    mpq_sub(rational, rational, r.rational);
    return *this;
  }

  const Rational operator-() const {
    Rational r;
    mpq_neg(r.rational, rational);
    return r;
  }

  friend const Rational operator+(const Rational& left, const Rational& right) {
    Rational r;
    mpq_add(r.rational, left.rational, right.rational);
    return r;
  }
  friend const Rational operator-(const Rational& left, const Rational& right) {
    Rational r;
    mpq_sub(r.rational, left.rational, right.rational);
    return r;
  }
  friend const Rational operator*(const Rational& left, const Rational& right) {
    Rational r;
    mpq_mul(r.rational, left.rational, right.rational);
    return r;
  }
  friend const Rational operator/(const Rational& left, const Rational& right) {
    Rational r;
    assert(right != 0);
    mpq_div(r.rational, left.rational, right.rational);
    return r;
  }
  friend bool operator==(const Rational& left, const Rational& right) {
    return (mpq_equal(left.rational, right.rational));
  }
  friend bool operator!=(const Rational& left, const Rational& right) {
    return (!mpq_equal(left.rational, right.rational));
  }
  friend bool operator<=(const Rational& left, const Rational& right) {
    return (mpq_cmp(left.rational, right.rational) <= 0);
  }
  friend bool operator>=(const Rational& left, const Rational& right) {
    return (mpq_cmp(left.rational, right.rational) >= 0);
  }
  friend bool operator<(const Rational& left, const Rational& right) {
    return (mpq_cmp(left.rational, right.rational) < 0);
  }
  friend bool operator>(const Rational& left, const Rational& right) {
    return (mpq_cmp(left.rational, right.rational) > 0);
  }
  void print(ostream& os) const {
    assert(mpq_denref(rational) != 0);
    if (mpz_get_si(mpq_denref(rational)) == 1) {
      os << mpz_get_si(mpq_numref(rational));
      return;
    }
    if (mpq_denref(rational) == 0) {
      os << "NaN";
      return;
    }
    os << mpz_get_si(mpq_numref(rational)) << '/';
    os << mpz_get_si(mpq_denref(rational));
  }
  friend ostream& operator<<(ostream& os, const Rational &r) {
    r.print(os);
    return os;
  }
  friend istream& operator>>(istream& is, Rational& r) {
    // READ ONLY NON-NEGATIVE INTEGER, NOT RATIONAL !!!!!!!! 
    Long p;
    is >> p;
    mpq_set_si(r.rational, p, 1);
    mpq_canonicalize(r.rational);
    return is;
  }
};

const Rational RATIONAL_MAX = Rational((Long) (sqrt((double) LONG_LONG_MAX)/2.), (Long) 1);

inline double to_double(const int cost) {return (double) cost;}
inline double to_double(const Long cost) {return (double) cost;}
inline double to_double(const Rational r) {return r.to_double();}

inline int ceil(const int e) {return e;}
inline int floor(const int e) {return e;}
inline Long ceil(const Long e) {return e;}
inline Long floor(const Long e) {return e;}
inline Long ceil(const Rational r) {return (Long) ceil(r.to_double());}
inline Long floor(const Rational r) {return (Long) floor(r.to_double());}

#endif /*TB2RATIONAL_HPP_*/
