/** \file tb2rangeminquery.hpp
 *  \brief Special query for Tree global constraint
 */

#ifndef TB2RANGEMINQUERY_HPP
#define TB2RANGEMINQUERY_HPP

#include <vector>

using namespace std;

template <class T>
class RangeMinQuery {

private:
	vector<T> A;

	int n;
	vector<int> pow2array;
	vector<int> log2array;	
	vector<vector<int> > M;

public:
	RangeMinQuery(): n(0) {}
	
	~RangeMinQuery() {	
	}

	T& operator[](int i) {return A[i];}

	void push_back(const T &t){A.push_back(t);}

	int size() {return A.size();}

	void clear() {A.clear();}

	void pre_compute(){						
		if (n != (int) A.size())
		{			
			pow2array.clear();
			log2array.clear();
			M.clear();
                        
			n = A.size();			
			pow2array.resize(n+1);
			log2array.resize(n+1);
			
			for (int i=0;i<n+1;i++) log2array[i] = -1;
			pow2array[0] = 1;		
			log2array[1] = 0;		
			for (int i=1;i<n+1;i++) {
				pow2array[i] = pow2array[i-1]*2;
				if (pow2array[i] < n+1) log2array[pow2array[i]] = i;
			}

			int logVal = 0;
			for (int i=1;i<n+1;i++) {
				if (log2array[i] == -1) log2array[i] = logVal;
				else logVal = log2array[i];			
			}

			M.resize(n);
			for (int i=0;i<n;i++) {
				M[i].resize(n);			
			}
		} 		
		
		for (int i=0;i<n;i++) M[i][0] = i;
		for (int j=1;pow2array[j]<=n;j++) {
			for (int i=0;i < n - pow2array[j] + 1;i++) {				
				int minL = M[i][j-1];                               
				int minR = M[i + pow2array[j-1]][j-1];				
				if (A[minL] < A[minR]) M[i][j] = minL;
				else M[i][j] = minR;
			}
		}
	}

	int query(int start, int end) {
		int logWidth = log2array[end - start + 1];		
		int minL = M[start][logWidth];
		int minR = M[end - pow2array[logWidth] + 1][logWidth];
		return ((A[minL] < A[minR])?minL:minR);
	}

};

#endif //TB2RANGEMINQUERY_HPP
