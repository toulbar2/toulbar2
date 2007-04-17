/*
 * ****** Random WCSP generator *******
 */
 

#include "tb2randomgen.hpp"
#include "tb2constraint.hpp"
#include "tb2variable.hpp"
#include "tb2enumvar.hpp"


bool naryRandom::connected() {
	return true;
}

void naryRandom::generateTernCtr( int i, int j, int k, long nogoods, Cost costMin, Cost costMax )
{
 	int a,b,c,dice;
 	EnumeratedVariable* x = (EnumeratedVariable*) wcsp.getVar(i);
    EnumeratedVariable* y = (EnumeratedVariable*) wcsp.getVar(j);
    EnumeratedVariable* z = (EnumeratedVariable*) wcsp.getVar(k);
	int mx = x->getDomainInitSize();
	int my = y->getDomainInitSize();
	int mz = z->getDomainInitSize();
    long total_nogoods = mx*my*mz;

    vector<Cost> costs;
	for (a = 0; a < mx; a++)
	    for (b = 0; b < my; b++) 
 		    for (c = 0; c < mz; c++) 
		        costs.push_back(0);

    while(nogoods>0) {
      dice = myrand() % total_nogoods;
      for(a=0;a<mx;a++)
	  for(b=0;b<my;b++)
	  for(c=0;c<mz;c++) {
	     if(costs[my*mz*a + b*mz + c] == 0) {
		    if(dice == 0) {
		    	costs[my*mz*a + b*mz + c] = randomCost(costMin, costMax);
			    nogoods--;
			    total_nogoods--;
			    a=mx;b=my;c=mz;
			 }
			 dice--;
	      }
	  }
    }
	wcsp.postTernaryConstraint(i,j,k,costs);
}




void naryRandom::generateBinCtr( int i, int j, long nogoods, Cost costMin, Cost costMax )
{
 	int a,b,dice;
 	EnumeratedVariable* x = (EnumeratedVariable*) wcsp.getVar(i);
    EnumeratedVariable* y = (EnumeratedVariable*) wcsp.getVar(j);
	int mx = x->getDomainInitSize();
	int my = y->getDomainInitSize();
    long total_nogoods = mx*my;

    vector<Cost> costs;
	for (a = 0; a < mx; a++)
	    for (b = 0; b < my; b++) 
	        costs.push_back(0);

    while(nogoods>0) {
      dice = myrand() % total_nogoods;
      for(a=0;a<mx;a++)
	  for(b=0;b<my;b++) {
	     if(costs[my*a+b] == 0) {
		    if(dice == 0) {
		    	costs[my*a+b] = randomCost(costMin, costMax);
			    nogoods--;
			    total_nogoods--;
			    a=mx;b=my;
			 }
			 dice--;
	      }
	  }
    }
	wcsp.postBinaryConstraint(i,j,costs);
}


long naryRandom::toIndex( vector<int>& index )
{
	long result = 1;
	for(int i=0;i<(int)index.size();i++) result += (long) pow((double)n,i)*index[i];
	return result;
}


void naryRandom::ini( vector<int>& index, int arity )
{
	index.clear();
	for(int i=0;i<arity;i++) index.push_back(i);
}


bool naryRandom::inc( vector<int>& index )
{
	int res = inc(index, index.size()-1);
	if(res < 0) return false;
	else return true;
}


int naryRandom::inc( vector<int>& index, int i )
{
	if(i < 0) return i;
	assert( i < (int) index.size());
	
	index[i]++;
	if(index[i] == n - ((int) index.size() - i - 1) ) {
		if(i>=0) {
			int val = inc(index,i-1);
			if(val < 0) return -1;
			index[i] = val+1;
			if(index[i] == n) return -1;					
		}
	}
	return index[i];
}


void naryRandom::Input( int in_n, int in_m, vector<int>& p )
{
  n = in_n;	
  m = in_m;
  	
  assert(p.size() >= 2);

  int i,arity;
  vector<int>  indexs;
  vector<long> totalCtrs;
  vector<long> numCtrs;

  int maxa = p.size();
  
  for(arity=0; arity <= maxa; arity++) {
	    if(arity < 2) numCtrs.push_back(0);
	    else 	      numCtrs.push_back(p[arity-1]);
  }
    
  for(i=0;i<n;i++) {
  	string varname = to_string(i);
  	wcsp.makeEnumeratedVariable(varname,0,m-1);
  }

  for(arity=maxa;arity>1;arity--)        
  {
  	  long nogoods =  (long) (((double)p[0] / 100.) * pow((double)m, arity));
  	  long totalarraysize = (long) pow( (double)n, arity);
	  vector<bool> existCtr;
      long tCtrs = 1;
	  for(i=0; i < totalarraysize; i++) existCtr.push_back(false);
	  for(i=0; i < arity; i++)  tCtrs *= (n - i);
	  for(i=2; i <= arity; i++)  tCtrs /= i;	
		
	  if(numCtrs[arity] > tCtrs) { 
	  	cout << numCtrs[arity] << "  " << arity << "ary constraints and the maximum is " << tCtrs << endl;
	  	numCtrs[arity] = tCtrs;
	  } 
	
  	  while(numCtrs[arity])
  	  {	
		  bool oneadded = false;
	      int dice = myrand() % tCtrs;
	
		  ini(indexs,arity);
		  do {  
 
		  		 if(!existCtr[toIndex(indexs)]) {
					 if(dice == 0) {
					    existCtr[toIndex(indexs)] = true;
					    switch(arity) {
					    	case 2:  generateBinCtr(indexs[0],indexs[1],nogoods); break;
					    	case 3:  generateTernCtr(indexs[0],indexs[1],indexs[2],nogoods); break;
					    	default: /*generateNaryCtr(indexs,nogoods)*/ ;
					    }
					    tCtrs--;
					    numCtrs[arity]--;
					    oneadded = true;
					 }
					 dice--;
		  		 }
		  } while(inc(indexs) && !oneadded);
  	   }
    }
}
	  
 


