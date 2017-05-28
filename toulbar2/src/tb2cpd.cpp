#include "tb2cpd.hpp"
#include <sstream>

Cpd::Cpd()
{
	cpdtrie = new TrieCpd();
}

Cpd::~Cpd()
{
    delete cpdtrie;
}

void Cpd::read_rotamers2aa(ifstream &file, vector<Variable *> &vars) throw (int)
{
	istringstream line;
	string s;

	file.unget();
	while (file) {
		vector<char> rot2aa_var;
        
		char current_char;
		getline(file, s, '\n');
		line.str(s);
		line.clear();
		if (line.str().empty()) {
			line.str().clear();
			continue;
		}
		
		while (line >> current_char) {
			if (!isspace(current_char)) rot2aa_var.push_back(current_char);
		}
		if (rot2aa_var.size() != 0) rotamers2aa.push_back(rot2aa_var);
	}
	//~ for (int i=0;i<rotamers2aa.size();i++){
		//~ for(int j=0;j<rotamers2aa[i].size();j++){
			//~ cout<< rotamers2aa[i][j];
		//~ }
		//~ cout<<endl;
	//~ }
	if (rotamers2aa.size() != vars.size()) {
		cout << "Wrong variable number " << rotamers2aa.size() << " " << vars.size() << endl;
		throw 1;
	} else {
		for (size_t i = 0; i < rotamers2aa.size(); i++) {
            unsigned int initsize = dynamic_cast<EnumeratedVariable *>(vars[i])->getDomainInitSize();
			if (rotamers2aa[i].size() != initsize) {
				cout << "Wrong domain size " << rotamers2aa[i].size() << " " << initsize << " of variable "<<dynamic_cast<EnumeratedVariable *>(vars[i])->getName()<<endl;
				throw 2;
			}
        }
	}

    
    for (auto &rv : rotamers2aa) {
        
        vector<Value> leftidx_var;
        vector<Value> rightidx_var;

        char prev_char = '0';
        size_t pos = 0;

        for (size_t i = 0; i < rv.size(); i++) {
            if (rv[i] != prev_char) {
                prev_char = rv[i];
                pos = i;
            }
            leftidx_var.push_back(pos);
        }
        
        prev_char = '0';
        for (int i = rv.size()-1; i >= 0; i--) {
            if (rv[i] != prev_char) {
                prev_char = rv[i];
                pos = i;
            }
            rightidx_var.push_back(pos);
        }

        LeftAA.push_back(leftidx_var);
        reverse(rightidx_var.begin(), rightidx_var.end());
        RightAA.push_back(rightidx_var);
    }
}


void Cpd::storeSequence(const vector<Variable *> &vars, Cost _cost)
{
	string sequence;
	for (size_t i = 0; i < vars.size(); i++) {
		sequence.push_back(rotamers2aa[i][vars[i]->getValue()]);
	}
	cpdtrie->insert_sequence(sequence, _cost);
}

void Cpd::printSequences()
{
	cpdtrie->print_tree();
}

void Cpd::printSequence(const vector<Variable *> &vars, Cost _cost)
{
	string sequence;
	for (size_t i = 0; i < vars.size(); i++) {
		sequence.push_back(rotamers2aa[i][vars[i]->getValue()]);
	}
	cout << "New sequence: " << sequence << " Cost: " << _cost << endl;
}

void Cpd::printSequence(TAssign &vars)
{
	//  cpdtrie->print_tree();
	string sequence;
	// TAssign::iterator it = vars.begin();
	// while(it != vars.end())
	//   {
	//     sequence.push_back(rotamers2aa[it->first][it->second]);
	//     ++it;
	//   }
	for (size_t i = 0; i < vars.size(); i++) {
		sequence.push_back(rotamers2aa[i][vars[i]]);
	}

	cout << "New sequence: " << sequence << endl;
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
