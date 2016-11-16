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
			//       line >> ws;
			rot2aa_var.push_back(current_char);
			//       cout << "Push: " << current_char << endl;getchar();
		}
		rotamers2aa.push_back(rot2aa_var);
	}
	if (rotamers2aa.size() != vars.size()) {
		cout << "Wrong variable number " << rotamers2aa.size() << " " << vars.size() << endl;
		throw 1;
	} else {
		for (size_t i = 0; i < rotamers2aa.size(); i++) {
            unsigned int initsize = dynamic_cast<EnumeratedVariable *>(vars[i])->getDomainInitSize();
			if (rotamers2aa[i].size() != initsize) {
				cout << "Wrong domain size " << rotamers2aa[i].size() << " " << initsize << endl;
				throw 2;
			}
        }
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
	//  cpdtrie->print_tree();
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
