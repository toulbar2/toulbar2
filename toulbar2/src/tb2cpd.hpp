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
	void storeSequence(const vector<Variable *> &vars, Cost _cost);
	void printSequences();
	void printSequence(const vector<Variable *> &vars, Cost _cost);
	void printSequence(TAssign &vars);
	int getTotalSequences() { return cpdtrie->getTotalSequences();}
	vector< vector<char> > &getRotamers2AA() { return rotamers2aa; }
	char getAA(int varIndex, Value value) { return rotamers2aa[varIndex][value];}
	size_t rot2aaSize(int varIndex) { return rotamers2aa[varIndex].size();}
private:
	vector< vector<char> > rotamers2aa;
	TrieCpd *cpdtrie;
};


#endif
