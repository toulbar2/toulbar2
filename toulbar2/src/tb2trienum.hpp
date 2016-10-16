#include <iostream>
#include <vector>
#include "toulbar2lib.hpp"

// Class Node : A node matches with a variable
class VarNode
{
public:
	VarNode() { m_Content = -1; m_Marker = false; } // Constructor
	~VarNode() {} //Destructor
	int get_Content() { return m_Content; } // Get the content of the NODE
	void set_Content(int rot) { m_Content = rot; } // Set the content of the NODE
	Cost get_cost() { return m_cost; } // Get the cost of the NODE
	void set_cost(Cost ct) { m_cost = ct; } // Set the cost of the NODE
	bool get_Marker() { return m_Marker; } // Get bool Marker
	void true_Marker() { m_Marker = true; } // Set the marker to true
	VarNode *find_child(int c); // find child c in Children's list
	void append_child(VarNode *child) { m_Children.push_back(child); } // Append a child to the children's list
	std::vector<VarNode *> get_Children() { return m_Children; } // Get the children's list

private:
	int m_Content; // Content of the NODE (Numero of the rotamer for Variable)
	Cost m_cost;
	bool m_Marker; // Marker: if true the "word" is over (use to find partial assigment)
	std::vector<VarNode *> m_Children; // Vector with all the child of the node
};

class TrieNum
{
public:
	TrieNum(); //Constructor
	~TrieNum(); // Destructor
	std::vector<Cost> get_sols() { return sol_costs; }// Get the solution's list
	void add_costs(Cost c) { sol_costs.push_back(c);} // add solution to solution's list
	Cost get_logz() { return m_logz; }
	void set_logz(Cost logz) {m_logz = logz;}
	void add_Sol(std::vector<int> Sol, const Cost ct);
	bool search_Sol(std::vector<int> Sol);
private:
	VarNode *root; // Root of the trie list
	Cost m_logz;
	std::vector<Cost> sol_costs;
};
