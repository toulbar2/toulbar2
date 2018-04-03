#ifndef TB2TRIE_HPP_
#define TB2TRIE_HPP_

#include "tb2wcsp.hpp"

using namespace std;

class TrieLeaf;
class TrieNode {
public:
    TrieNode();
    ~TrieNode();
    void insert_sequence(string seq, unsigned int pos, Double energy);
    vector<TrieNode*> sons;
    bool present(char aa);
    void insertNode(char aa);
    void insertLeaf(char aa);
    int aa2int(char aa);
    char int2aa(unsigned int i);
    void print_tree(string acc);
    static size_t total_sequences;

private:
    TrieLeaf* getLeaf(char aa);
    static const string i2a;
    static const size_t numAA;
};

class TrieLeaf : public TrieNode {
public:
    TrieLeaf();
    ~TrieLeaf();
    unsigned long sequence_count;
    Double mine;
    Double maxe;
};

class TrieCpd {
public:
    TrieCpd(){};
    ~TrieCpd(){};
    void init() { root.sons.clear(); };
    void insert_sequence(string seq, Double energy);
    void print_tree();
    size_t getTotalSequences() { return root.total_sequences; }

private:
    TrieNode root;
};

#endif

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
