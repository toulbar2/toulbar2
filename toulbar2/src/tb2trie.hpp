#ifndef TB2TRIE_HPP_
#define TB2TRIE_HPP_

#include "tb2wcsp.hpp"

using namespace std;

class TrieLeaf;
class TrieNode {
public:
    TrieNode();
    ~TrieNode();
    void insert_sequence(string seq, unsigned int pos, Cost _cost);
    vector<TrieNode*> sons;
    bool present(char aa);
    void insertNode(char aa);
    void insertLeaf(char aa);
    int aa2int(char aa);
    char int2aa(unsigned int i);
    void print_tree(string acc);
private:
    TrieLeaf* getLeaf(char aa);
    static const string i2a;
};

class TrieLeaf : public TrieNode {
public:
    TrieLeaf();
    ~TrieLeaf();
    unsigned long sequence_count;
    Cost minc;
    Cost maxc;
};

class TrieCpd {
public:
    TrieCpd() {};
    ~TrieCpd() {};
    void insert_sequence(string seq, Cost _cost);
    void print_tree();
    size_t getTotalSequences() { return total_sequences; }

private:
    static size_t total_sequences;
    TrieNode root;
};

#endif

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
