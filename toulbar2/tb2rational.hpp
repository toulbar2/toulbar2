/** \file tb2rational.hpp
 *  \brief Rational numbers with basic operations.
 * 
 */

#ifndef TB2RATIONAL_HPP_
#define TB2RATIONAL_HPP_

struct Rational {
	Long p; // numerator
	Long q; // denominator
	
	Rational() : p(0), q(1) {}
	Rational(Long p_, Long q_) : p(p_), q(q_) {}
	Rational(Long p_) : p(p_), q(1) {}

    double to_double() const {return (double) p / (double) q;}
     
    Rational(const Rational &r) : p(r.p), q(r.q) {}
    	
	Rational &operator=(const Rational &r) {
        p = r.p;
        q = r.q;
        return *this;
    }
    Rational &operator+=(const Rational &r) {
        p = (p * r.q + q * r.p);
        q = (q * r.q);
        return *this;
    }
    Rational &operator-=(const Rational &r) {
        p = (p * r.q - q * r.p);
        q = (q * r.q);
        return *this;
    }
    const Rational operator-() const {return Rational(-p,q);}
	friend const Rational operator+(const Rational& left, const Rational& right) {
		return Rational(left.p * right.q + left.q * right.p, left.q * right.q);
	}
	friend const Rational operator-(const Rational& left, const Rational& right) {
		return Rational(left.p * right.q - left.q * right.p, left.q * right.q);
	}
	friend bool operator==(const Rational& left, const Rational& right) {
		return left.p * right.q == left.q * right.p;
	}
	friend bool operator!=(const Rational& left, const Rational& right) {
		return left.p * right.q != left.q * right.p;
	}
	friend bool operator<=(const Rational& left, const Rational& right) {
		return left.p * right.q <= left.q * right.p;
	}
	friend bool operator>=(const Rational& left, const Rational& right) {
		return left.p * right.q >= left.q * right.p;
	}
	friend bool operator<(const Rational& left, const Rational& right) {
		return left.p * right.q < left.q * right.p;
	}
	friend bool operator>(const Rational& left, const Rational& right) {
		return left.p * right.q > left.q * right.p;
	}
	void print(ostream& os) const { os << p << '/' << q; }
	friend ostream& operator<<(ostream& os, const Rational &r) {
		os << r.p << '/' << r.q;
		return os;
	}
	friend istream& operator>>(istream& is, Rational& r) {
		is >> r.p;
		r.q = 1;		// READ ONLY INTEGER, NOT RATIONAL !!!!!!!! 
		return is;
	}
};

const Rational RATIONAL_MAX = Rational((Long) (sqrt((double) LONG_LONG_MAX)/2.), (Long) 1);

inline double to_double(const int cost) {return (double) cost;}
inline double to_double(const Long cost) {return (double) cost;}
inline double to_double(const Rational r) {return r.to_double();}

inline double ceil(const Rational r) {return ceil(r.to_double());}
inline double floor(const Rational r) {return floor(r.to_double());}

#endif /*TB2RATIONAL_HPP_*/
