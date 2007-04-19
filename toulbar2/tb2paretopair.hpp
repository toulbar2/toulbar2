/** \file tb2paretopair.hpp
 *  \brief ParetoPair numbers with basic operations.
 * 
 */

#ifndef TB2PARETOPAIR_HPP_
#define TB2PARETOPAIR_HPP_

struct ParetoPair {
	int p; // 
	int q; // 
	
    ParetoPair() : p(0), q(0) {}

    ParetoPair(int p_, int q_) : p(p_), q(q_) {}

    ParetoPair(int p_) : p(p_), q(p_) {}
     
    double to_double() const {cerr << "to_double not implemented on Paretopair"; exit(EXIT_FAILURE);}

    ParetoPair(const ParetoPair &r) : p(r.p), q(r.q) {}
    	
	ParetoPair &operator=(const ParetoPair &r) {
        p = r.p;
        q = r.q;
        return *this;
    }
    ParetoPair &operator+=(const ParetoPair &r) {
        p = (p +  r.p);
        q = (q + r.q);
        return *this;
    }
    ParetoPair &operator-=(const ParetoPair &r) {
        p = (p  - r.p);
        q = (q  - r.q); 
        return *this;
    }
    const ParetoPair operator-() const {return ParetoPair(-p,-q);}

  friend const ParetoPair operator+(const ParetoPair& left, const ParetoPair& right) {
		return ParetoPair(left.p + right.p, left.q + right.q);
	}
     
  friend const ParetoPair operator-(const ParetoPair& left, const ParetoPair& right) {
		return ParetoPair(left.p - right.p, left.q - right.q);
	}


  friend const ParetoPair operator*(const ParetoPair& left, const ParetoPair& right) {
        return ParetoPair(left.p * right.p, left.q * right.q);
    }
    
  friend const ParetoPair operator/(const ParetoPair& left, const ParetoPair& right) {
        return ParetoPair(left.p / right.p, left.q / right.q);
    }
	
  friend bool operator==(const ParetoPair& left, const ParetoPair& right) {
		return (left.p == right.p) & (left.q == right.q);
	}
	
  friend bool operator!=(const ParetoPair& left, const ParetoPair& right) {
		return (left.p != right.p)  | (left.q != right.q);
	}

  friend bool operator<=(const ParetoPair& left, const ParetoPair& right) {
		return (left.p <= right.p) & (left.q <= right.q);
	}

  friend bool operator>=(const ParetoPair& left, const ParetoPair& right) {
		return (left.p >= right.p) & ( left.q >= right.q);
	}

  friend bool operator<(const ParetoPair& left, const ParetoPair& right) {
		return ((left.p <= right.p) & ( left.q < right.q)) |
		  ((left.p < right.p) & ( left.q <= right.q));
	}

  friend bool operator>(const ParetoPair& left, const ParetoPair& right) {
		return ((left.p >= right.p) & ( left.q > right.q)) |
		  ((left.p > right.p) & ( left.q >= right.q));
	}

  void print(ostream& os) const { os << '(' << p << ',' << q << ')'; }
	
  friend ostream& operator<<(ostream& os, const ParetoPair &r) {
    os << '(' << r.p << ',' << r.q << ')';
		return os;
	}

// friend istream& operator>>(istream& is, ParetoPair& r) {
//   int c1,c2;
//   is >> c1;
//   is >> c2;
//      r.p = c1;
//      r.q = c2;       
//      return is;
//  }

   friend istream& operator>>(istream& is, ParetoPair& r) {
  		is >> r.p;
  		r.q = r.p;		// READ ONLY INTEGER, NOT PARETOPAIR !!!!!!!! 
  		return is;
  	}
};

const ParetoPair PARETOPAIR_MAX = ParetoPair(INT_MAX/2, INT_MAX/2);

const ParetoPair PARETOPAIR_MIN = ParetoPair(0,0);

inline double to_double(const ParetoPair r) {cerr << "to_double not implemented on Paretopair"; exit(EXIT_FAILURE);}
inline ParetoPair ceil(const ParetoPair r) {return r;}
inline ParetoPair floor(const ParetoPair r){return r;}
inline ParetoPair randomCost(ParetoPair min, ParetoPair max) {return ParetoPair(min.p + (myrand() % (max.p - min.p + 1)), min.q + (myrand() % (max.q - min.q + 1)));}
inline ParetoPair string2Cost(char *ptr) {return atoi(ptr);}

inline int cost2log2(int x)
{
        if (x==0) return -1;
        register int l2 = 0;
        x>>=1;
        for (; x != 0; x >>=1)
        {
                ++ l2;
        }
        return (l2);
}
inline int cost2log2(const ParetoPair &r)
{  return  min(cost2log2(r.p), cost2log2(r.p));
}

#endif /*TB2PARETOPAIR_HPP_*/
