
#include "tb2naryconstr.hpp"
#include "tb2wcsp.hpp"


NaryConstraint::NaryConstraint(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in, Cost defval)
			:nonassigned(arity_in, &wcsp->getStore()->storeValue), AbstractNaryConstraint(wcsp, scope_in, arity_in)
{
	default_cost = defval;
	store_top = default_cost < wcsp->getUb();
	xy = new BinaryConstraint(wcsp, &wcsp->getStore()->storeCost );
}


// for adding a tuple in f
// scope_in contains the order of the values in string tin 
void NaryConstraint::insertTuple( string tin, Cost c, EnumeratedVariable** scope_in = NULL)
{
	string t(tin);

	if(scope_in) {
		for(int i = 0; i < arity_; i++) {
			int pos = getIndex(scope_in[i]);
			t[pos] = tin[i];
		}  
	}
	f[t] = c;
}  


Cost NaryConstraint::eval( string s ) {
	TUPLES::iterator  it = f.find(s);
	if(it != f.end()) return it->second;
	else return default_cost;  
}


// USED ONLY DURING SEARCH to project the nary constraint 
// to the binary constraint xy, when all variables but 2 are assigned 
void NaryConstraint::projectNaryBinary()
{
	int indexs[2];
	EnumeratedVariable* unassigned[2];

	char* tbuf = new char [arity_ + 1]; 
	string t;

	int i,j = 0;
	for(i=0;i<arity_;i++) { 		
		if(getVar(i)->unassigned()) { 
			unassigned[j] = (EnumeratedVariable*) getVar(i); 
			indexs[j++] = i;
			tbuf[i] = CHAR_FIRST;
		} else tbuf[i] = getVar(i)->getValue() + CHAR_FIRST;
	}
	tbuf[arity_] =  '\0';
	t = tbuf;
	delete [] tbuf;
	
	EnumeratedVariable* x = unassigned[0];
	EnumeratedVariable* y = unassigned[1];
	xy->fillElimConstr(x,y);

	for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
    for (EnumeratedVariable::iterator itery = y->begin(); itery != y->end(); ++itery) {	
		Value xval = *iterx;
		Value yval = *itery;
		t[indexs[0]] =  xval + CHAR_FIRST;			
		t[indexs[1]] =  yval + CHAR_FIRST;					
		xy->setcost(xval,yval,eval(t));
    }}

	BinaryConstraint* ctr = x->getConstr(y);   			
	if(ctr) {
		ctr->reconnect();
		ctr->addCosts(xy);
		ctr->propagate();
	}
	else {
		xy->reconnect();
		xy->propagate();
	}	
}


// USED ONLY DURING SEARCH 
void NaryConstraint::assign(int varIndex) {
	assert(nonassigned > 2);	
	nonassigned = nonassigned - 1;
	if(nonassigned == 2) {
			deconnect();
			projectNaryBinary();
	}
}





void NaryConstraint::sum( AbstractNaryConstraint* nary )
{
	
}

// Projection of variable x of the nary constraint 
// complexity O(2|f|)
// this function is independent of the search 
void NaryConstraint::project( EnumeratedVariable* x )
{
	if(arity_ <= 1) return;
	if(getIndex(x) < 0) return;
	
	string t,tnext,tproj;
	Cost c,cnext;
	
	TUPLES fproj;
	
	int xindex = getIndex(x);
	
	// First part of the projection: complexity O(|f|)
	// we swap positions between the projected variable
	// and the last variable
	
	TUPLES::iterator  it = f.begin();
	while(it != f.end()) {
		t = it->first;
		c =  it->second;
		
		string tswap(t);
		char a = tswap[arity_-1];
		tswap[arity_-1] = tswap[xindex];
		tswap[xindex] = a;
				
		fproj[tswap] = c;		
		f.erase(t);		
				
		it++;
	} 

	// Second part of the projection: complexity O(|f|)
	// as the projected variable is in the last position,
	// it is sufficient to look for tuples with the same
	// arity-1 posotions. If there are less than d (domain of
	// the projected variable) tuples, we have also to perform 
	// the minimum with default_cost
	// this is only true when the tuples are LEXICOGRAPHICALY ordered	
    it = fproj.begin();
	while(it != fproj.end()) {
		t = it->first;
		c =  it->second;
		
		bool sameprefix = true;
		bool end = false;

		int ntuples = 1;

		while(!end) {		
			it++;
			end = (it == fproj.end());
			
			if(!end) { 
				tnext = it->first;
				cnext = it->second;
				sameprefix = (t.compare(0,arity_-2,tnext) == 0);
				if(!sameprefix) {
					if((ntuples < x->getDomainInitSize()) && (cnext < c)) c = cnext;
					if(c != default_cost) f[ t.substr(0,arity_-2) ] = c;
					t = tnext;
					c = cnext;
					ntuples = 1;
				} else {
					ntuples++;
					if(cnext < c) c = cnext;
				}			
			}
		}
	}
	
	// Update of internal structures: scope, indexes, links, arity_.
	varElem    x_ve = {x, xindex};
	varElem last_ve = {(EnumeratedVariable*) getVar(arity_-1), xindex};
		
	scope_inv.erase(&x_ve);
	x->deconnect(links[xindex]);
	if(!links[xindex]->removed) nconnected_links--;
	links[xindex] = links[arity_-1];

	scope_inv.insert(&last_ve);
	scope[ xindex ] = (EnumeratedVariable*) getVar(arity_-1);
	arity_--;	
}



void NaryConstraint::print(ostream& os)
{
	os << endl << "f(";
	for(int i = 0; i < arity_;i++) 
	{
		os << scope[i]->wcspIndex;
		if(i < arity_-1) os << ",";
	}
	os << ")    ";
	os << " |f| = " << f.size();
	os << "    default_cost: " << default_cost << endl;
	os << "tuples: {";
	
	TUPLES::iterator  it = f.begin();
	while(it != f.end()) {
		string t = it->first;
		Cost c =  it->second;		
		it++;
		os << "<" << t << "," << c << ">";
		if(it != f.end()) os << " "; 
	} 
	os << "} " << endl;
}
