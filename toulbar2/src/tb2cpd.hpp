#ifndef TB2CPD_HPP_
#define TB2CPD_HPP_

#include "tb2wcsp.hpp"
#include "tb2trie.hpp"

class Cpd
{
public:
	Cpd();
	~Cpd();
	void read_rotamers2aa(ifstream &file, vector<Variable *> &vars) throw (int);   ///< \brief read rotamer to amino acid correspondence
	void readPSMatrix(const char *filename);
	void readPSSMatrix(const char *filename);
	void storeSequence(const vector<Variable *> &vars, Cost _cost);
	void printSequences();
	void printSequence(const vector<Variable *> &vars, Cost _cost);
	void printSequence(TAssign &vars);
	int getTotalSequences() { return cpdtrie->getTotalSequences();}
	vector< vector<char> > &getRotamers2AA() { return rotamers2aa; }
	char getAA(int varIndex, Value value) { return rotamers2aa[varIndex][value];}
	Value getLeft(int varIndex, Value value) { return LeftAA[varIndex][value];}
	Value getRight(int varIndex, Value value) { return RightAA[varIndex][value];}
	size_t rot2aaSize(int varIndex) { return rotamers2aa[varIndex].size();}
	string nativeSequence;
	double PSMBias = 0.0;
	double PSSMBias = 0.0;
private:
	const static map<char,int> PSMIdx; // converts AA char to indices in PSMatrix
	const static map<char,int> PSSMIdx; // converts AA char to indices in PsiBlast PSSMatrix
	vector< vector<char> > rotamers2aa;
	vector< vector<Value> > LeftAA;
	vector< vector<Value> > RightAA;
	int PSM[24][24] = {{1, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
							 {0, 1, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
							 {0, 0, 1, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
							 {0, 0, 0, 1, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
							 {0, 0, 0, 0, 1, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
							 {0, 0, 0, 0, 0, 1, 0 ,0 ,0 ,0, 0, 0, 0, 0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
							 {0, 0, 0, 0, 0, 0, 1 ,0 ,0 ,0, 0, 0, 0, 0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
							 {0, 0, 0, 0, 0, 0, 0 ,1 ,0 ,0, 0, 0, 0, 0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
							 {0, 0, 0, 0, 0, 0, 0 ,0 ,1 ,0, 0, 0, 0, 0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
							 {0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,1, 0, 0, 0, 0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
							 {0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 1, 0, 0, 0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
							 {0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 1, 0, 0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
							 {0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 1, 0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
							 {0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 1 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
							 {0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0 ,1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
							 {0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0 ,0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
							 {0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0 ,0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
							 {0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0 ,0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
							 {0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0 ,0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
							 {0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0 ,0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
							 {0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0 ,0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
							 {0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0 ,0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
							 {0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0 ,0, 0, 0, 0, 0, 0, 0, 0, 1, 0},
							 {0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 1}};
	vector <vector<int>> PSSM;
	TrieCpd *cpdtrie;
};


#endif
