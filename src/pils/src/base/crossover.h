
#ifndef _crossover_h
#define _crossover_h


class Crossover {
public:

	/*
		return true when the child solution is different to both parent solutions (improvement)
	*/
	virtual bool operator()(Solution & x1, Solution & x2, Solution & child) = 0;

};

#endif