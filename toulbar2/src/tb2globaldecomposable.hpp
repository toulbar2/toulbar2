/** \file tb2globaldecomposable.hpp
 *  \brief Decomposable global cost functions : WeightedRegular, WeightedAmong
 */
#ifndef TB2GLOBALDECOMPOSABLE_HPP_
#define TB2GLOBALDECOMPOSABLE_HPP_

#include <iostream>
#include <fstream>
#include <set>
#include <list>
#include "tb2wcsp.hpp"
#include "tb2types.hpp"
#include "tb2enumvar.hpp"

using namespace std;

class DecomposableGlobalCostFunction {
	protected:
			int 		arity;		
			int*	 	scope;
	public:
		DecomposableGlobalCostFunction();
		DecomposableGlobalCostFunction(unsigned int _arity, int* _scope);
		//~DecomposableGlobalCostFunction(); 
		static DecomposableGlobalCostFunction* FactoryDGCF(string type, unsigned int _arity, int* _scope, ifstream &file);
		
		virtual void addToCostFunctionNetwork(WCSP* wcsp)=0;
		virtual void display()=0;
};

class WeightedAmong : public DecomposableGlobalCostFunction {
	private:
		set<int> 	values;
		string 		semantics;
		Cost 		baseCost;
		unsigned 	int lb;
		unsigned 	int ub;
	public:
		WeightedAmong();
		WeightedAmong(unsigned int _arity, int* _scope);
		WeightedAmong(unsigned int _arity, int* _scope, ifstream &file);
		//~WeightedAmong();
		
		inline void addValue(int _value)							{ values.insert(_value); }
		inline void setSemantics(string _semantics) 				{ semantics=_semantics; }
		inline void setBaseCost(Cost _baseCost) 					{ baseCost=_baseCost; }
		inline void setBounds(unsigned int _lb, unsigned int _ub) 	{ lb=_lb; ub=_ub; }
		
		void addToCostFunctionNetwork(WCSP* wcsp);
		void display();
};

struct WTransition {
	unsigned int start;
	unsigned int end;
	unsigned int symbol;
	Cost weight;
	
	WTransition(unsigned int _start, unsigned int _end, unsigned int _symbol, Cost  _weight) {
		start = _start;
		end = _end;
		symbol = _symbol;
		weight = _weight;
	}
	
	void display() {
		cout << start << " x " << symbol << " --("<<weight <<")--> " << end << endl;
	}
};
struct WFA{
	unsigned int nbStates;
	list<pair<int,Cost> > initialStates;
	list<pair<int,Cost> > acceptingStates;
	list<WTransition*> transitions;
	
	WFA(unsigned int _nbStates) {
		nbStates = _nbStates;
	}
	/*~WFA() {
		initialStates.clear();
		acceptingStates.clear();
		transitions.clear();
	}*/
	void display() {
		cout << "Number of states = " << nbStates << endl;
		cout << "Initial States : " << endl;
		for (list<pair<int,Cost> >::iterator it = initialStates.begin() ; it != initialStates.end() ; it++) {
			pair<int,int> initial = *it;
			cout << initial.first << "(" << initial.second << ")" << endl;
		}
		cout << "Accepting States : " << endl;
		for (list<pair<int,Cost> >::iterator it = acceptingStates.begin() ; it != acceptingStates.end() ; it++) {
			pair<int,int> accepting = *it;
			cout << accepting.first << "(" << accepting.second << ")" << endl;
		}
		cout << "Transition : " << endl;
		for (list<WTransition*>::iterator it = transitions.begin() ; it != transitions.end() ; it++) {
			WTransition* transition = *it;
			transition->display();
		}
	}
};

class WeightedRegular : public DecomposableGlobalCostFunction{
	private:
		WFA* automaton;
	public:
		WeightedRegular();
		WeightedRegular(unsigned int _arity, int* _scope);
		WeightedRegular(unsigned int _arity, int* _scope, ifstream &file);
		//~WeightedRegular();
		
		inline void setWFA(WFA* _automaton) {automaton=_automaton;}
		
		void addToCostFunctionNetwork(WCSP* wcsp);
		void display(); 
};

class WeightedSum : public DecomposableGlobalCostFunction{
	private:
		string 	comparator;
		int 	rightRes;
		string 	semantics;
		Cost 	baseCost;
	public:
		WeightedSum();
		WeightedSum(unsigned int _arity, int* _scope);
		WeightedSum(unsigned int _arity, int* _scope, ifstream &file);
		//~WeightedSum();
		
		inline void setSemantics(string _semantics) 				{ semantics=_semantics; }
		inline void setBaseCost(Cost _baseCost) 					{ baseCost=_baseCost; }
		inline void setComparator(string _comparator) 				{ comparator=_comparator;}
		inline void setRightRes(int _rightRes)						{ rightRes=_rightRes;}
		
		void addToCostFunctionNetwork(WCSP* wcsp);
		void display(); 
};

#endif
