#include <iostream>
#include <vector>
#include <algorithm>    // std::sort
#include "tb2trienum.hpp"
using namespace std;



VarNode *VarNode::find_child(int c)
{
    for (auto it : m_Children) {
        VarNode *tmp = it;
        if (tmp->get_Content() == c) {
            return tmp;
        }
    }
    return NULL;
}

TrieNum::TrieNum()
{
    root = new VarNode();
    m_logz = MAX_COST;
}

TrieNum::~TrieNum()
{
    // Free memory
}

void TrieNum::add_Sol(vector<int> Sol, const Cost ct)
{
    VarNode *current = root; // Init to root node

    //assert(Sol.size()  != 0);

    for (auto i : Sol) {
        VarNode *child = current->find_child(i);
        if (child != NULL) { // If the child is already in the list go to the next one
            current = child;
        } else {
            VarNode *tmp = new VarNode();
            tmp->set_Content(i);
            current->append_child(tmp);
            current = tmp;
        }
        if (i == Sol.back()) {
            current->true_Marker();
            current->set_cost(ct);
        }
    }
}


bool TrieNum::search_Sol(vector<int> Sol)
{
    VarNode *current = root;

    while (current != NULL) { // Until tree is finished
        for (auto i : Sol) {
            VarNode *tmp = current->find_child(i);
            if (tmp == NULL) return false;
            current = tmp;
        }
        if (current->get_Marker()) return true;
        else return false;
    }

    return false;
}
