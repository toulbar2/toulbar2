/*
 solution.h

 Author: 
  Sebastien Verel, 
  Univ. du Littoral Côte d'Opale, France.
 
  David Simoncini
  Univ. Toulouse 1 Capitole, France.

*/

#ifndef _solution_h
#define _solution_h

#include <iostream>
#include <vector>
#include <cfloat>

class Solution : public std::vector<int> {
public:
  Solution() : std::vector<int>() {
    valid = false;
  }

  Solution(unsigned int _n) : std::vector<int>(_n) {
    valid = false;
  }

  Solution(const Solution & _s) : std::vector<int>(_s) {
    _fitness = _s.fitness();
    valid = _s.isValid();
  }

  Solution& operator=(const Solution & _s) {
    this->resize(_s.size());
    for(unsigned int i = 0; i < _s.size(); i++)
      this->operator[](i) = _s[i];

    _fitness = _s.fitness();
    valid = _s.isValid();

    return *this;
  }

  /**
   * set the fitness
   */
  void fitness(Cost _fit) {
    valid    = true;
    _fitness = _fit;
  }

  /**
   * get the fitness
   */
  Cost fitness() const {
    return _fitness; 
  }

  /**
   * print the solution
   */
  void print() {
      printOn(std::cout);
  }

  virtual void printOn(std::ostream& _os) const {
    if (valid) {
      Cost res = this->fitness();
      long long resi = (long long)res;
      if (resi < (res - DBL_EPSILON) || resi > (res + DBL_EPSILON)) {
    	_os << res;
      } else {
        _os << resi;
      }
    } else {
      _os << "I" ;
    }

    // print size
    _os << " " << this->size();

    for(unsigned int i = 0; i < this->size(); i++)
      _os << " " << this->operator[](i) ;
  }

  virtual void readFrom(std::istream& _is) {
    std::string fit_str;
    int pos = _is.tellg();
    _is >> fit_str;

    if (fit_str == "I") 
      invalidate();
    else {
      valid = true;
      _is.seekg(pos); 
      _is >> _fitness;
    }

    // size
    unsigned s;
    _is >> s;

    this->resize(s);

    for(unsigned i = 0; i < s; i++)
      _is >> this->operator[](i);

  }

  void invalidate() {
    valid = false;
  }

  bool isValid() const {
    return valid;
  }

private:
  // quality of the solution
  Cost _fitness;

  // valid fitness value if true
  bool valid;
};

std::ostream & operator<<(std::ostream& _os, const Solution& _s) {
    _s.printOn(_os);
    return _os;
}

std::istream & operator >> (std::istream& _is, Solution& _s) {
  _s.readFrom(_is);
  return _is;
}

#endif
