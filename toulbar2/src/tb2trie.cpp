#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include "tb2trie.hpp"

using namespace std;

const string TrieNode::i2a = "ACDEFGHIKLMNPQRSTVWYZ";
const size_t TrieNode::numAA = i2a.length();

size_t TrieNode::total_sequences = 0;

int TrieNode::aa2int(char aa)
{
    for (std::size_t i = 0; i < i2a.size(); ++i)
        if (aa == i2a[i])
            return i;

    cout << "Error, unrecognized residue: " << aa << endl;
    exit(1);
}

char TrieNode::int2aa(unsigned int i)
{
    if (i < i2a.length()) {
        return i2a[i];
    } else {
        cout << "Error, index out of range: " << i << endl;
        exit(1);
    }
}

void TrieCpd::insert_sequence(string seq, Double energy)
{
    root.insert_sequence(seq, 0, energy);
}

void TrieNode::insert_sequence(string seq, unsigned int pos, Double energy)
{
    if (pos + 1 < seq.length()) {
        if (!present(seq[pos])) {
            insertNode(seq[pos]);
        }
        sons[aa2int(seq[pos])]->insert_sequence(seq, pos + 1, energy);
    } else {
        if (!present(seq[pos])) {
            insertLeaf(seq[pos]);
        }
        TrieLeaf* leaf = getLeaf(seq[pos]);
        leaf->sequence_count++;
        leaf->maxe = std::max(leaf->maxe, energy);
        leaf->mine = std::min(leaf->mine, energy);
    }
}

void TrieCpd::print_tree()
{
    root.print_tree("");
}

void TrieNode::print_tree(string acc)
{
    if (sons.size() == 0) {
        TrieLeaf* imaleaf = static_cast<TrieLeaf*>(this);
        cout << acc << " min cost: " << imaleaf->mine << " max cost: " << imaleaf->maxe
             << " occurred " << imaleaf->sequence_count << " times" << endl;
    } else {
        for (size_t i = 0; i < TrieNode::numAA; i++)
            if (sons[i] != NULL) {
                acc.push_back(int2aa(i));
                sons[i]->print_tree(acc);
                acc.pop_back();
            }
    }
}

TrieNode::TrieNode()
{
    sons.resize(TrieNode::numAA, NULL);
}

TrieLeaf::TrieLeaf()
{
    sons.clear();
    sequence_count = 0;
    mine = -log(0.0);
    maxe = -mine;
}

TrieNode::~TrieNode()
{
    for (size_t i = 0; i < sons.size(); i++)
        delete sons[i];
}

bool TrieNode::present(char aa)
{
    return (sons[aa2int(aa)] != NULL);
}

void TrieNode::insertNode(char aa)
{
    sons[aa2int(aa)] = new TrieNode();
}

void TrieNode::insertLeaf(char aa)
{
    sons[aa2int(aa)] = new TrieLeaf();
    total_sequences++;
}

TrieLeaf* TrieNode::getLeaf(char aa)
{
    return static_cast<TrieLeaf*>(sons[aa2int(aa)]);
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
