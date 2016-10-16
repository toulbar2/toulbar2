#ifndef TB2TRIE_HPP_
#define TB2TRIE_HPP_

#include "tb2wcsp.hpp"

using namespace std;

class TrieCpd
{
public:
	TrieCpd();
	~TrieCpd();
	int aa2int(char aa);
	char int2aa(int i);
	bool present(char aa);
	void insert(char aa);
	void insert_sequence(string seq, Cost _cost);
	void print_tree();
	void print_tree(string acc);
	size_t getTotalSequences() {return total_sequences;}
private:
	vector<TrieCpd *> sons;
	size_t sequence_count;
	static size_t total_sequences;
	Cost maxc;
	Cost minc;
};

#endif
