#include "tb2globaldecomposable.hpp"

/// DECOMPOSABLE COST FUNCTION /////////////////////////////////////////

DecomposableGlobalCostFunction::DecomposableGlobalCostFunction() : arity(0), label("empty") {
	ToulBar2::Berge_Dec=1;
}

DecomposableGlobalCostFunction::DecomposableGlobalCostFunction(unsigned int _arity, int* _scope) : arity(_arity), label("empty"){
	scope = new int[arity];
	for (unsigned int variable = 0 ; variable < _arity ; ++variable) {
		scope[variable] = _scope[variable];
	}
	ToulBar2::Berge_Dec=1;
}

DecomposableGlobalCostFunction::~DecomposableGlobalCostFunction() {
	delete [] scope;
}

DecomposableGlobalCostFunction* 
DecomposableGlobalCostFunction::FactoryDGCF(string type, unsigned int _arity, int* _scope, istream &file) {
	//cout << "Creating a " << type << " global cost function " << endl;
	if (type == "wamong") 				return new WeightedAmong(_arity,_scope,file);
	if (type == "wregular")				return new WeightedRegular(_arity,_scope,file);
	if (type == "wsum")					return new WeightedSum(_arity,_scope,file);
	if (type == "woverlap")				return new WeightedOverlap(_arity,_scope,file);
	cout << type << " unknown decomposable global cost function" << endl;
	return 0;
}


void 
DecomposableGlobalCostFunction::color(int i) {
	switch (i) {
		case 1 : cout << "\033[41m"; break;
		case 2 : cout << "\033[42m"; break;
		case 3 : cout << "\033[43m"; break;
		case 4 : cout << "\033[44m"; break;
		case 5 : cout << "\033[45m"; break;
		case 6 : cout << "\033[46m"; break;
		case 7 : cout << "\033[47m"; break;
		case 8 : cout << "\033[40m\033[37m"; break;
		default : cout << "\033[0m"; break;
	};
}
/// WEIGHTED AMONG /////////////////////////////////////////////////////

WeightedAmong::WeightedAmong() : DecomposableGlobalCostFunction() {}

WeightedAmong::WeightedAmong(unsigned int _arity, int* _scope) : DecomposableGlobalCostFunction(_arity,_scope) {
}

WeightedAmong::WeightedAmong(unsigned int _arity, int* _scope, istream &file) : DecomposableGlobalCostFunction(_arity,_scope) {
	file >> semantics >> baseCost;
	unsigned int nbValue;
	file >> nbValue;
	for (unsigned int value = 0 ; value < nbValue ; ++value) {
		int valueRead;
		file >> valueRead;
		values.insert(valueRead);
	}
	file >> lb >> ub;
}

WeightedAmong::~WeightedAmong() {
	values.clear();
}

void
WeightedAmong::addToCostFunctionNetwork(WCSP* wcsp) {
	bool VERBOSE = false;
	int nbVariableCFN = wcsp->numberOfVariables();
	//cout << nbVariableCFN << endl;
	
	// -- new variables : counters -- //
	int addVariablesIndex[arity+1];
	for (int newVariable = 0 ; newVariable <= arity ; newVariable++) {
		string varname="WAmong" + to_string(nbVariableCFN)+"_"+ to_string(newVariable);
		addVariablesIndex[newVariable] = wcsp->makeEnumeratedVariable(varname,0,newVariable);
		if (VERBOSE) {color(5); cout << "new variable " << addVariablesIndex[newVariable] << "("<< ((EnumeratedVariable *) wcsp->getVar(addVariablesIndex[newVariable]))->getDomainInitSize()<< ")"; color(0); cout << endl;}
		wcsp->getListSuccessors()->push_back(vector<int>());
	}
	
	Cost top = wcsp->getUb();
	// -- ternary constraints : partial sum -- //
	for (int variable = 0 ; variable < arity ; ++variable) {
		int indexCi = addVariablesIndex[variable];
		int indexCj = addVariablesIndex[variable+1];
		int indexXi  = scope[variable];
		if (VERBOSE) {color(5); cout << indexCi << "--" << indexXi << "--" << indexCj; color(0); cout << endl;}
		EnumeratedVariable* varCi = (EnumeratedVariable *) wcsp->getVar(indexCi); 
		EnumeratedVariable* varCj = (EnumeratedVariable *) wcsp->getVar(indexCj); 
		EnumeratedVariable* varXi = (EnumeratedVariable *) wcsp->getVar(indexXi);
		wcsp->getListSuccessors()->at(indexCi).push_back(indexXi);
		wcsp->getListSuccessors()->at(indexXi).push_back(indexCj); 


		
		unsigned long tableSize = long (varCi->getDomainInitSize() * varCj->getDomainInitSize() * varXi->getDomainInitSize());
		vector<Cost> ternaryCosts(tableSize,top);
		
		for (unsigned long value = 0; value < varXi->getDomainInitSize() ; value++) {
				for (unsigned long counter = 0; counter < varCi->getDomainInitSize() ; counter++) {
					int nextCounter = counter;
					if (values.find(value) != values.end()) {
						nextCounter++;						
					}
					unsigned long position =  (counter) 		* (varXi->getDomainInitSize()*varCj->getDomainInitSize()) 
											+ (value) 			* (varCj->getDomainInitSize()) 
											+ (nextCounter);
					ternaryCosts[position] = 0;
			}
		}
		wcsp->postTernaryConstraint(indexCi, indexXi, indexCj, ternaryCosts);			
	}
	
	// -- unary constraints : final variable -- // 
	if (VERBOSE) {color(5); cout << "post unary constraint on " << addVariablesIndex[arity]; color(0); cout << endl;}
	vector<Cost> unaryCosts(arity+1,top);
	for (int i = 0 ; i <= arity ; i++) {
		int gap = max( 0 , max( int(lb - i), int(i - ub) ) );
		if (semantics == "hard") {
			if (((unsigned int) i) >= lb && ((unsigned int) i)  <= ub) unaryCosts[i] = 0;							
			else unaryCosts[i] = min(top,baseCost);
		}
		if (semantics == "lin")  unaryCosts[i] = min (top, baseCost * gap);							
		if (semantics == "quad") unaryCosts[i] = min (top,baseCost * gap * gap);
	}	
	wcsp->postUnary(addVariablesIndex[arity],unaryCosts);
}

Cost 
WeightedAmong::evaluate(int* tuple) {
	int occurency = 0;
	for (int var = 0 ; var < arity ; var++) {
		if(values.find(tuple[var]) != values.end()) occurency++;
	}
	int gap = max( 0 , max( int(lb - occurency), int(occurency - ub) ) );
	if (gap) {
		if (semantics == "hard") return baseCost;
		if (semantics == "lin")  return baseCost * gap;	
		if (semantics == "quad") return baseCost * gap * gap;	
	}
	return 0;
}

void 
WeightedAmong::display() {
	cout << "WAmong (" << arity << ") : ";
	for (int variable = 0 ; variable < arity ; ++variable) {
		cout << scope[variable] << " ";
	}
	cout << endl;
	cout << "sem : " << semantics << " " << baseCost << endl;
	cout << "val : ";
	for (set<int>::iterator value = values.begin() ; value != values.end() ; value++) {
		cout << *value << " ";
	}
	cout << endl;
	cout << "bounds [" << lb << ":" << ub << "]" << endl;
}


/// WEIGHTED REGULAR ///////////////////////////////////////////////////

WeightedRegular::WeightedRegular() : DecomposableGlobalCostFunction(), automaton(0) {}

WeightedRegular::WeightedRegular(unsigned int _arity, int* _scope) : DecomposableGlobalCostFunction(_arity,_scope),automaton(0) {}

WeightedRegular::WeightedRegular(unsigned int _arity, int* _scope, istream &file) : DecomposableGlobalCostFunction(_arity,_scope) {
	automaton = new WFA(file);
}

WeightedRegular::~WeightedRegular() {
	delete automaton;
}

void
WeightedRegular::addToCostFunctionNetwork(WCSP* wcsp) {
	ToulBar2::Berge_Dec=1; 
	//automaton->display();
	Cost top = wcsp->getUb();
	int unsigned current_var_number = wcsp->numberOfVariables();
	int unsigned q0 = current_var_number;
	if ( ToulBar2::verbose > 1 ) 	{
		cout << "DEBUG>> wregular found : inital number of variables before creation = " <<  wcsp->numberOfVariables() <<endl;
		cout << "DEBUG>> wregular Automatum Total number of states: " <<  automaton->getNbStates() <<endl;
		cout << "DEBUG>> wregular Initial states number: " << automaton->getInitialStates().size()  <<endl;
		cout << "DEBUG>> wregular add new variable from " << q0 <<" to "<<  current_var_number+arity+1 <<endl;
	}
	if(  current_var_number > 0 ) { 
		int unsigned domsize  = automaton->getNbStates()-1;
		string varname = "WR" + to_string(current_var_number);
		if ( ToulBar2::verbose > 1 ) cout << "DEBUG>> wregular q0 index "<< q0 << " domain = " << domsize+1 << endl; 
		int theindex = -1;
		theindex=wcsp->makeEnumeratedVariable(varname,0,domsize);			// add q0 variable
		wcsp->getListSuccessors()->push_back(vector<int>()); 				// add the new variable in the topological order;
		if ( ToulBar2::verbose > 1 ) cout << "wregular add varname =" << varname <<"=> var index "<<  wcsp->numberOfVariables() << " domain size = " << domsize+1 << endl;
		} else { exit(EXIT_FAILURE); 
	} 		
	//################################################initial state ##############################################
	if (automaton->getInitialStates().size() > 0 ) {
		vector<Cost> initial_states_costs(automaton->getNbStates(),top);
		list< pair<int,Cost> > initialStates = automaton->getInitialStates();
		for (list<pair<int,Cost> >::iterator it = initialStates.begin(); it !=  initialStates.end() ; ++it){
			pair<int,Cost> initial = *it;
			//cout << initial.first << " " << initial.second << endl;
			initial_states_costs[initial.first]=initial.second;
		}
		wcsp->postUnary(q0,initial_states_costs);
		if ( ToulBar2::verbose > 1 ) {
			cout << "DEBUG>> wregular initial state (q0) vector size ( nbre value) = "<<  initial_states_costs.size() << endl;
			cout << "DEBUG>> wregular var q0= "<< q0 <<" number of constraint wregular initialisation ==>" << wcsp->numberOfConstraints()  << endl;
		}
	}
	//################################################accepting state ##############################################
	for( int v = 1 ; v < arity+1 ; v++ )  {
	 	int unsigned domsize = automaton->getNbStates()-1;
		string varname = to_string(v+q0);
					
		int theindex = -1;
		theindex=wcsp->makeEnumeratedVariable(varname,0,domsize);			// add qi variable
		wcsp->getListSuccessors()->push_back(vector<int>()); 				// add new variable in the topological order;
		assert(theindex=v+current_var_number);
		if ( ToulBar2::verbose > 1 ) cout << "DEBUG>> wregular add varname =" << varname <<"=> rank "<<  wcsp->numberOfVariables() << " domain = " << domsize+1 << endl;
	}	
	int unsigned q_last = wcsp->numberOfVariables() -1 ;
	if ( ToulBar2::verbose > 1 ) cout << "DEBUG>> wregular Final number of variables : " << wcsp->numberOfVariables() << endl;
	vector<Cost>final_states_costs(automaton->getNbStates(),top);
	
	list< pair<int,Cost> > acceptingStates = automaton->getAcceptingStates();
	if (acceptingStates.size()>=0) {
					
		for (list<pair<int,Cost> >::iterator it = acceptingStates.begin(); it !=  acceptingStates.end() ; ++it ) {
			pair<int,Cost> accept = *it;
			int unsigned t_index = accept.first; 
			Cost ucost = accept.second;
			
		 	EnumeratedVariable* Qv  = (EnumeratedVariable *) wcsp->getVar(q_last); // get domaine size of last qi var
			unsigned long DomVar=Qv->getDomainInitSize();
	
			if(t_index < DomVar){
				final_states_costs[t_index]=ucost;
			} else { 
				cout <<"wregular tuple error " << t_index << "out of domain" << DomVar << endl; 
				exit(EXIT_FAILURE); 
			} 		  
		}			  
		wcsp->postUnary(q_last,final_states_costs);
					
		if ( ToulBar2::verbose > 1 )  {
			cout << "DEBUG>> wregular last q varname = "<<  q_last <<endl;
		}
	}
				/*
				//################################################### lecture des transition ????
				if ( ToulBar2::verbose > 1 ) 
				cout << "DEBUG>>wregular final number of Unary constraint Post after q0 and qi post : " << numberOfConstraints()  << endl;
				//==================
				// transition stat reading 
				//==================
				int nb_transition;
				vector<unsigned int> VQi;
				vector<unsigned int> VQj;
				vector<unsigned int> VXi;
				vector<Cost> transition_costs;
				file >> nb_transition;
				if ( ToulBar2::verbose > 1 ) cout << "DEBUG>> wregular transitions number :  " <<  nb_transition <<endl;
				for( int s = 0 ; s < nb_transition ; s++)
				{
					int qi;
					int xi;
					int qj;
					Cost transition_COST;
					file >> qi;
					file >> xi;
					file >> qj;
					file >> transition_COST;
					VQi.push_back(qi);
					VXi.push_back(xi);
					VQj.push_back(qj);
					transition_costs.push_back(transition_COST);

					if ( ToulBar2::verbose > 1 ) {
					cout << "DEBUG>> wregular read transition table qi =" << qi << " xi =" << xi << " qj =" << qj << " cost=" << transition_COST << endl;
					cout << "DEBUG>> wregular scope xi " << scopeIndex[xi] <<endl;
					}
				}
*/
//##################################################ajout ternaire#############################
	for ( int q = 0 ; q < arity ; q++) {
	    int qi = q0+q;
		int	xi = scope[q] ;
		int qj = qi+1;
		if ( ToulBar2::verbose > 1 ) cout << "DEBUG>>wregular  post ternary on  qi =" << qi << " xi =" << xi << " qj =" << qj << endl;
		// poiner on qi , xi , qj varibale 	
		EnumeratedVariable* Qi = (EnumeratedVariable *) wcsp->getVar(qi); //current qi variable;
		EnumeratedVariable* Xi = (EnumeratedVariable *) wcsp->getVar(xi); //current Xi variable;
		EnumeratedVariable* Qj = (EnumeratedVariable *) wcsp->getVar(qj); //current qj variable;
		// domain definition 
		unsigned long DomQi=Qi->getDomainInitSize();
		unsigned long DomQj=Qi->getDomainInitSize();
		unsigned long DomXi=Xi->getDomainInitSize();
		unsigned long Domsize = long (Qi->getDomainInitSize() * Xi->getDomainInitSize() * Qj->getDomainInitSize());
		
		vector<Cost>tmp_ternary_costs(Domsize,top);
		list<WTransition*> transitions = automaton->getTransitions();
		for(list<WTransition*>::iterator it = transitions.begin() ; it != transitions.end() ; ++it)  {
			WTransition* transition = *it;
			int start = transition->start;
			int end = transition->end;
			int symbol = transition->symbol;
			int weight = transition->weight;
			
			if(symbol < DomXi) {
				unsigned long cindex= start*DomXi*DomQj + symbol*DomQj + end;
				tmp_ternary_costs[cindex]=weight;
			}
			if(ToulBar2::verbose > 1) {
					cout << "DEBUG>> wregular init cost vector for ternary rank= " << q << endl;
					cout << "DEBUG>> wregular add ternary table qi =" << qi << " xi =" << xi << " qj =" << qj << endl; 
					cout <<"DEBUG>> wregular Ternary const DOMAIN SIZE = " << Domsize <<" Dom qi ="<< DomQi << " Dom Xi=" << DomXi << " Dom Qj =" << DomQj <<" -------"<<endl;
					cout << "DEBUG>> wregular initial COST SIZE = " << tmp_ternary_costs.size() << endl;
					cout << "DEBUG>> wregular number of constraint before ternary cost function : " << wcsp->numberOfConstraints()  << endl;
				}
		}
		wcsp->postTernaryConstraint(qi, xi, qj, tmp_ternary_costs);
		wcsp->getListSuccessors()->at(qi).push_back(xi);
		wcsp->getListSuccessors()->at(xi).push_back(qj);
	}
	if ( ToulBar2::verbose >=1 ) cout << "DEBUG>> wregular Total number of constraints after wregular posting " << wcsp->numberOfConstraints()  << endl;	
}

/*void WeightedRegular::addToCostFunctionNetwork(WCSP* wcsp) {
	//display();
	
	int nbVariableCFN = wcsp->numberOfVariables();
	
	// -- new variables : counters -- //
	int addVariablesIndex[arity+1];
	for (int newVariable = 0 ; newVariable <= arity ; newVariable++) {
		string varname="WRegular" + to_string(nbVariableCFN)+"_"+ to_string(newVariable);
		addVariablesIndex[newVariable] = wcsp->makeEnumeratedVariable(varname,0,automaton->nbStates-1);
		//cout << "\033[45m" << "new variable " << addVariablesIndex[newVariable] << "("<< ((EnumeratedVariable *) wcsp->getVar(addVariablesIndex[newVariable]))->getDomainInitSize()<< ")" << "\033[0m" << endl;
		wcsp->getListSuccessors()->push_back(vector<int>());
	}
	
	Cost top = wcsp->getUb();
	// -- ternary constraints : transition table -- //
	for (int variable = 0 ; variable < arity ; ++variable) {
		int indexCi = addVariablesIndex[variable];
		int indexCj = addVariablesIndex[variable+1];
		int indexXi  = scope[variable];
		EnumeratedVariable* varCi = (EnumeratedVariable *) wcsp->getVar(indexCi); 
		EnumeratedVariable* varCj = (EnumeratedVariable *) wcsp->getVar(indexCj); 
		EnumeratedVariable* varXi = (EnumeratedVariable *) wcsp->getVar(indexXi);
		wcsp->getListSuccessors()->at(indexCi).push_back(indexXi);
		wcsp->getListSuccessors()->at(indexXi).push_back(indexCj); 
		cout << "\033[45m" << "post ternary constraint on " << indexCi <<"("<<varCi->getDomainInitSize()<<")" << ","  << indexXi <<"("<<varXi->getDomainInitSize()<<")" << ","  << indexCj <<"("<<varCj->getDomainInitSize()<<")" << "\033[0m" << endl;		
		
		unsigned long tableSize = long (varCi->getDomainInitSize() * varCj->getDomainInitSize() * varXi->getDomainInitSize());
		vector<Cost> ternaryCosts(tableSize,top);
		
		for (list<WTransition*>::iterator it = automaton->transitions.begin() ; it != automaton->transitions.end() ; ++it) {
			WTransition* transition = *it;
			unsigned int start  = transition->start;
			unsigned int end    = transition->end;
			unsigned int symbol = transition->symbol;
			Cost weight         = transition->weight;
			
			unsigned long position =  long(start) 		* long(varXi->getDomainInitSize()*varCj->getDomainInitSize()) 
									+ long(symbol) 		* long(varCj->getDomainInitSize()) 
									+ long(end);
			if (start >= varCi->getDomainInitSize()) cout << "WARNING !!! WARNING !!!" << endl;
			if (end >= varCj->getDomainInitSize()) cout << "WARNING !!! WARNING !!!" << endl;
			if (symbol >= varXi->getDomainInitSize()) cout << "WARNING !!! WARNING !!!" << endl;
			ternaryCosts[position] = weight;			
		}
		wcsp->postTernaryConstraint(indexCi, indexXi, indexCj, ternaryCosts);		
	}
	
	// -- unary constraints : initial and accepting state -- //
	//cout << "\033[45m" << "post unary constraint on " << addVariablesIndex[0] << "\033[0m" << endl;
	 vector<Cost> unaryCostsInit(arity+1,top);
	 for (list<pair<int,Cost> >::iterator it = automaton->initialStates.begin() ; it != automaton->initialStates.end() ; it++) {
		 pair<int,Cost> init = *it;
		 unaryCostsInit[init.first] = init.second;
	 }
	 wcsp->postUnary(addVariablesIndex[0],unaryCostsInit);
	 
	 //cout << "\033[45m" << "post unary constraint on " << addVariablesIndex[arity] << "\033[0m" << endl;
	 vector<Cost> unaryCostsAccept(arity+1,top);
	 for (list<pair<int,Cost> >::iterator it = automaton->acceptingStates.begin() ; it != automaton->acceptingStates.end() ; it++) {
		 pair<int,Cost> accept = *it;
		 unaryCostsAccept[accept.first] = accept.second;
	 }
	 wcsp->postUnary(addVariablesIndex[arity],unaryCostsAccept);
}*/

void 
WeightedRegular::display() {
	cout << "WRegular (" << arity << ") : ";
	for (int variable = 0 ; variable < arity ; ++variable) {
		cout << scope[variable] << " ";
	}
	cout << endl;
	if (automaton) {
		automaton->display();
	}
	else {
		cout << "no automaton associated" << endl;
	}
}
			
/// WEIGHTED SUM ///////////////////////////////////////////////////////

WeightedSum::WeightedSum() : DecomposableGlobalCostFunction() {}

WeightedSum::WeightedSum(unsigned int _arity, int* _scope) : DecomposableGlobalCostFunction(_arity,_scope) {}

WeightedSum::WeightedSum(unsigned int _arity, int* _scope, istream &file) : DecomposableGlobalCostFunction(_arity,_scope) {
	file >> semantics >> baseCost;
	file >> comparator >> rightRes;
}

WeightedSum::~WeightedSum() {}

void 
WeightedSum::addToCostFunctionNetwork(WCSP* wcsp) {
	int nbVariableCFN = wcsp->numberOfVariables();
	//cout << nbVariableCFN << endl;
	
	// -- new variables : counters -- //
	int addVariablesIndex[arity+1];
	int cumul = 0;
	for (int newVariable = 0 ; newVariable <= arity ; newVariable++) {
		string varname="WSum" + to_string(nbVariableCFN)+"_"+ to_string(newVariable);
		addVariablesIndex[newVariable] = wcsp->makeEnumeratedVariable(varname,0,cumul);
		//cout << "\033[45m" << "new variable " << addVariablesIndex[newVariable] << "("<< ((EnumeratedVariable *) wcsp->getVar(addVariablesIndex[newVariable]))->getDomainInitSize()<< ")" << "\033[0m" << endl;
		wcsp->getListSuccessors()->push_back(vector<int>());
		if (newVariable < arity) cumul += ((EnumeratedVariable *) wcsp->getVar(scope[newVariable]))->getDomainInitSize()-1;
	}
	
	Cost top = wcsp->getUb();
	// -- ternary constraints : partial sum -- //
	for (int variable = 0 ; variable < arity ; ++variable) {
		int indexCi = addVariablesIndex[variable];
		int indexCj = addVariablesIndex[variable+1];
		int indexXi  = scope[variable];
		EnumeratedVariable* varCi = (EnumeratedVariable *) wcsp->getVar(indexCi); 
		EnumeratedVariable* varCj = (EnumeratedVariable *) wcsp->getVar(indexCj); 
		EnumeratedVariable* varXi = (EnumeratedVariable *) wcsp->getVar(indexXi);
		wcsp->getListSuccessors()->at(indexCi).push_back(indexXi);
		wcsp->getListSuccessors()->at(indexXi).push_back(indexCj); 

		//cout << "\033[45m" << "post ternary constraint on " << indexCi <<"("<<varCi->getDomainInitSize()<<")" << ","  << indexXi <<"("<<varXi->getDomainInitSize()<<")" << ","  << indexCj <<"("<<varCj->getDomainInitSize()<<")" << "\033[0m" << endl;		
		
		unsigned long tableSize = long (varCi->getDomainInitSize() * varCj->getDomainInitSize() * varXi->getDomainInitSize());
		vector<Cost> ternaryCosts(tableSize,top);
		
		for (unsigned long value = 0; value < varXi->getDomainInitSize() ; value++) {
				for (unsigned long counter = 0; counter < varCi->getDomainInitSize() ; counter++) {
					int nextCounter = counter + value;
					unsigned long position =  (counter) 		* (varXi->getDomainInitSize()*varCj->getDomainInitSize()) 
											+ (value) 			* (varCj->getDomainInitSize()) 
											+ (nextCounter);
					ternaryCosts[position] = 0;
					//cout << counter << " + " << value << " => " << nextCounter << endl;
			}
		}
		wcsp->postTernaryConstraint(indexCi, indexXi, indexCj, ternaryCosts);			
	}
	
	// -- unary constraints : final variable -- // 
	//cout << "\033[45m" << "post unary constraint on " << addVariablesIndex[arity] << "\033[0m" << endl;
	vector<Cost> unaryCosts(cumul+1 ,top);
	if (comparator == "==") {
		for (int i = 0 ; i <= cumul ; i++) {
			if (i == rightRes) unaryCosts[i] = 0;
			else {
				int gap = (i < rightRes)?  rightRes - i: i - rightRes;
				if (semantics == "hard") unaryCosts[i] = baseCost;
				if (semantics == "lin")  unaryCosts[i] = gap*baseCost;
				if (semantics == "quad")  unaryCosts[i] = gap*gap*baseCost;
			}
		}
	} else
	if (comparator == "!=") {
		for (int i = 0 ; i <= cumul ; i++) {
			if (i != rightRes) unaryCosts[i] = 0;
			else {
				unaryCosts[i] = baseCost;
			}
		}
	} else 
	if (comparator == "<" || comparator == "<=") {
		int newRightRes = rightRes;
		if (comparator == "<") newRightRes--;
		
		for (int i = 0 ; i <= cumul ; i++) {
			if (i <= newRightRes) unaryCosts[i] = 0;
			else {
				int gap = max(0,i - newRightRes);
				if (semantics == "hard") unaryCosts[i] = baseCost;
				if (semantics == "lin")  unaryCosts[i] = gap*baseCost;
				if (semantics == "quad")  unaryCosts[i] = gap*gap*baseCost;
			}
		}
	} else
		if (comparator == ">" || comparator == ">=") {
		int newRightRes = rightRes;
		if (comparator == ">") newRightRes++;
		
		for (int i = 0 ; i <= cumul ; i++) {
			if (i >= newRightRes) unaryCosts[i] = 0;
			else {
				int gap = max(0,newRightRes - i);
				if (semantics == "hard") unaryCosts[i] = baseCost;
				if (semantics == "lin")  unaryCosts[i] = gap*baseCost;
				if (semantics == "quad")  unaryCosts[i] = gap*gap*baseCost;
			}
		}
	}
	else cout << "comparator " << comparator << " not handle yet" << endl;
	wcsp->postUnary(addVariablesIndex[arity],unaryCosts);
}

Cost 
WeightedSum::evaluate(int* tuple) {
	int sum = 0;
	for (int var = 0 ; var < arity ; var++) {
		sum += tuple[var];
	}
	if (comparator == "==") {
		int gap = (sum < rightRes)?  rightRes - sum: sum - rightRes;
		if (semantics == "hard") return min(gap*baseCost,baseCost);
		if (semantics == "lin")  return gap*baseCost;
		if (semantics == "quad") return gap*gap*baseCost;
	}
	if (comparator == "!=") {
		if (sum != rightRes) return baseCost;
	}
	if (comparator == "<" || comparator == "<=") {
		int newRightRes = rightRes; if (comparator == "<") newRightRes--;
		int gap = max(0,sum - newRightRes);
		if (semantics == "hard") return min(gap*baseCost,baseCost);
		if (semantics == "lin")  return gap*baseCost;
		if (semantics == "quad") return gap*gap*baseCost;
	}
		if (comparator == ">" || comparator == ">=") {
		int newRightRes = rightRes; if (comparator == ">") newRightRes++;
		int gap = max(0,newRightRes - sum);
		if (semantics == "hard") return min(gap*baseCost,baseCost);
		if (semantics == "lin")  return gap*baseCost;
		if (semantics == "quad") return gap*gap*baseCost;
	}
	return 0;
}

void 
WeightedSum::display() {
	cout << "WSum (" << arity << ") : ";
	for (int variable = 0 ; variable < arity ; ++variable) {
		cout << scope[variable] << " ";
	}
	cout << endl;
	cout << comparator << " " << rightRes << endl;
	cout << semantics << " " << baseCost << endl;
}

/// WEIGHTED SUM ///////////////////////////////////////////////////////

WeightedOverlap::WeightedOverlap() : DecomposableGlobalCostFunction() {}

WeightedOverlap::WeightedOverlap(unsigned int _arity, int* _scope) : DecomposableGlobalCostFunction(_arity,_scope) {}

WeightedOverlap::WeightedOverlap(unsigned int _arity, int* _scope, istream &file) : DecomposableGlobalCostFunction(_arity,_scope) {
	file >> semantics >> baseCost;
	file >> comparator >> rightRes;
	//display();
}

WeightedOverlap::~WeightedOverlap() {}

void 
WeightedOverlap::addToCostFunctionNetwork(WCSP* wcsp) {
	int nbVariableCFN = wcsp->numberOfVariables();
	
	// -- new variables : counters OVERLAP -- //
	int addVariablesOverlap[arity/2];
	for (int newVariable = 0 ; newVariable < arity/2 ; newVariable++) {
		string varname="WOVERL_OVER_" + to_string(nbVariableCFN)+"_"+ to_string(newVariable);
		addVariablesOverlap[newVariable] = wcsp->makeEnumeratedVariable(varname,0,1);
		wcsp->getListSuccessors()->push_back(vector<int>());
		//cout << "add overlap " << newVariable << endl;
	}
	
	Cost top = wcsp->getUb();
	// -- ternary -- //
	for (int newVariable = 0 ; newVariable < arity/2 ; newVariable++) {
		int X = scope[newVariable];
		int Y = scope[newVariable+arity/2];
		int O = addVariablesOverlap[newVariable];
		EnumeratedVariable* varX = (EnumeratedVariable *) wcsp->getVar(X); 
		EnumeratedVariable* varY = (EnumeratedVariable *) wcsp->getVar(Y); 
		EnumeratedVariable* varO = (EnumeratedVariable *) wcsp->getVar(O);
		wcsp->getListSuccessors()->at(X).push_back(O);
		wcsp->getListSuccessors()->at(Y).push_back(O); 
		
		unsigned long tableSize = long (varX->getDomainInitSize() * varY->getDomainInitSize() * varO->getDomainInitSize());
		vector<Cost> ternaryCosts(tableSize,top);
		
		//cout << X << "--" << Y  << "--" << O << endl;
		for (unsigned int vX = 0 ; vX < varX->getDomainInitSize() ; vX++) {
			for (unsigned int vY = 0 ; vY < varY->getDomainInitSize() ; vY++) {
				unsigned int vO = 0;
				if (vX == vY && vY >= 1)  vO = 1;
				//cout << "X=" << vX << " Y=" << vY << " O=" << vO << endl;
				unsigned long position =  (vX) 		* (varY->getDomainInitSize()*varO->getDomainInitSize()) 
										+ (vY) 		* (varO->getDomainInitSize()) 
										+ (vO);
				ternaryCosts[position] = 0;
			}
		}
		wcsp->postTernaryConstraint(X,Y,O,ternaryCosts);
	}
	
	// -- new variables : counters Among -- //
	int addVariablesAmong[arity/2+1];
	for (int newVariable = 0 ; newVariable < arity/2+1 ; newVariable++) {
		string varname="WOVERL_AMONG_" + to_string(nbVariableCFN)+"_"+ to_string(newVariable);
		addVariablesAmong[newVariable] = wcsp->makeEnumeratedVariable(varname,0,newVariable);
		wcsp->getListSuccessors()->push_back(vector<int>());
		//cout << "add among " << newVariable << endl;
	}
	
	// -- ternary -- //
	for (int newVariable = 0 ; newVariable < arity/2 ; newVariable++) {
		int indexCi = addVariablesAmong[newVariable];
		int indexCj = addVariablesAmong[newVariable+1];
		int indexXi  = addVariablesOverlap[newVariable];
		
		//cout << indexCi << " -- " << indexXi << " -- "  << indexCj << endl;
		EnumeratedVariable* varCi = (EnumeratedVariable *) wcsp->getVar(indexCi); 
		EnumeratedVariable* varCj = (EnumeratedVariable *) wcsp->getVar(indexCj); 
		EnumeratedVariable* varXi = (EnumeratedVariable *) wcsp->getVar(indexXi);
		wcsp->getListSuccessors()->at(indexCi).push_back(indexXi);
		wcsp->getListSuccessors()->at(indexXi).push_back(indexCj); 
		
		unsigned long tableSize = long (varCi->getDomainInitSize() * varCj->getDomainInitSize() * varXi->getDomainInitSize());
		vector<Cost> ternaryCosts(tableSize,top);
		
		for (unsigned long value = 0; value < varXi->getDomainInitSize() ; value++) {
				for (unsigned long counter = 0; counter < varCi->getDomainInitSize() ; counter++) {
					int nextCounter = counter;
					if (value == 1) {
						nextCounter++;						
					}
					unsigned long position =  (counter) 		* (varXi->getDomainInitSize()*varCj->getDomainInitSize()) 
											+ (value) 			* (varCj->getDomainInitSize()) 
											+ (nextCounter);
					ternaryCosts[position] = 0;
			}
		}
		
		wcsp->postTernaryConstraint(indexCi,indexXi,indexCj,ternaryCosts);
	}
	
	// -- unary constraints : final variable -- // 
	//cout << "\033[45m" << "post unary constraint on " << addVariablesIndex[arity] << "\033[0m" << endl;
	vector<Cost> unaryCosts(arity/2+1 ,top);
	if (comparator == "==") {
		for (int i = 0 ; i <= arity/2 ; i++) {
			if (i == rightRes) unaryCosts[i] = 0;
			else {
				int gap = (i < rightRes)?  rightRes - i: i - rightRes;
				if (semantics == "hard") unaryCosts[i] = baseCost;
				if (semantics == "lin")  unaryCosts[i] = gap*baseCost;
				if (semantics == "quad")  unaryCosts[i] = gap*gap*baseCost;
			}
		}
	} else
	if (comparator == "!=") {
		for (int i = 0 ; i <= arity/2 ; i++) {
			if (i != rightRes) unaryCosts[i] = 0;
			else {
				unaryCosts[i] = baseCost;
			}
		}
	} else 
	if (comparator == "<" || comparator == "<=") {
		int newRightRes = rightRes;
		if (comparator == "<") newRightRes--;
		
		for (int i = 0 ; i <= arity/2 ; i++) {
			if (i <= newRightRes) unaryCosts[i] = 0;
			else {
				int gap = max(0,i - newRightRes);
				if (semantics == "hard") unaryCosts[i] = baseCost;
				if (semantics == "lin")  unaryCosts[i] = gap*baseCost;
				if (semantics == "quad")  unaryCosts[i] = gap*gap*baseCost;
			}
		}
	} else
		if (comparator == ">" || comparator == ">=") {
		int newRightRes = rightRes;
		if (comparator == ">") newRightRes++;
		
		for (int i = 0 ; i <= arity/2 ; i++) {
			if (i >= newRightRes) unaryCosts[i] = 0;
			else {
				int gap = max(0,newRightRes - i);
				if (semantics == "hard") unaryCosts[i] = baseCost;
				if (semantics == "lin")  unaryCosts[i] = gap*baseCost;
				if (semantics == "quad")  unaryCosts[i] = gap*gap*baseCost;
			}
		}
	}
	else cout << "comparator " << comparator << " not handle yet" << endl;
	//cout << "control " << addVariablesAmong[arity/2] << endl;
	wcsp->postUnary(addVariablesAmong[arity/2],unaryCosts);
}

Cost 
WeightedOverlap::evaluate(int* tuple) {
	int occurency = 0;
	for (int var = 0 ; var < arity/2 ; var++) {
		if(tuple[var] && tuple[var] == tuple[var+arity/2]) {
			//cout << var << " && " << var+arity/2;
			occurency++;
		}
	}
	//cout << " => " << occurency << " " << comparator << " " << rightRes << endl;
	if (comparator == "==") {
		int gap = (occurency < rightRes)?  rightRes - occurency: occurency - rightRes;
		if (semantics == "hard") return min(gap*baseCost,baseCost);
		if (semantics == "lin")  return gap*baseCost;
		if (semantics == "quad") return gap*gap*baseCost;
	}
	if (comparator == "!=") {
		if (occurency != rightRes) return baseCost;
	}
	if (comparator == "<" || comparator == "<=") {
		int newRightRes = rightRes; if (comparator == "<") newRightRes--;
		int gap = max(0,occurency - newRightRes);
		if (semantics == "hard") return min(gap*baseCost,baseCost);
		if (semantics == "lin")  return gap*baseCost;
		if (semantics == "quad") return gap*gap*baseCost;
	}
		if (comparator == ">" || comparator == ">=") {
		int newRightRes = rightRes; if (comparator == ">") newRightRes++;
		int gap = max(0,newRightRes - occurency);
		if (semantics == "hard") return min(gap*baseCost,baseCost);
		if (semantics == "lin")  return gap*baseCost;
		if (semantics == "quad") return gap*gap*baseCost;
	}
	return 0;
}

void 
WeightedOverlap::display() {
	cout << "WOverlap (" << arity << ") : ";
	for (int variable = 0 ; variable < arity ; ++variable) {
		cout << scope[variable] << " ";
	}
	cout << endl;
	cout << semantics << " " << baseCost << endl;
	int i = 0;
	cout << "{ " ; 
	for ( ; i < arity/2 ; i++) cout << scope[i]<< " ";
	cout << "}" << endl;
	cout << "{ " ; 
	for ( ; i < arity ; i++) cout << scope[i]<< " ";
	cout << "}" << endl;
	cout << comparator << " " << rightRes << endl;
}

