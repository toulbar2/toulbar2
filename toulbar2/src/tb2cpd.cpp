#include "tb2cpd.hpp"
#include <sstream>

const map<char,int> Cpd::PSMIdx = {{'A', 0}, {'R' , 1},  {'N', 2}, {'D', 3}, {'C', 4}, {'Q', 5},
                                   {'E', 6},  {'G', 7},  {'H', 8},  {'I', 9},  {'L', 10},  {'K', 11},
                                   {'M', 12},  {'F', 13},  {'P', 14},  {'S', 15},  {'T', 16},  {'W', 17},
                                   {'Y', 18},  {'V', 19},  {'B', 20},  {'Z', 21},  {'X', 22},  {'*', 23}};

const map<char,int> Cpd::PSSMIdx = {{'A', 0}, {'G' , 1},  {'I', 2}, {'L', 3}, {'V', 4}, {'M', 5},
                                    {'F', 6},  {'W', 7},  {'P', 8},  {'C', 9},  {'S', 10},  {'T', 11},
                                    {'Y', 12},  {'N', 13},  {'Q', 14},  {'H', 15},  {'K', 16},  {'R', 17},
                                    {'D', 18},  {'E', 19}};

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

void Cpd::readPSMatrix(const char* filename)
{
    ifstream file;
    file.open(filename);
    
    if (!file.is_open()) {
        cerr << "Could not open PSSM file, aborting." << endl;
        exit(EXIT_FAILURE);
    }

    string s;
	istringstream line;
    int minscore = std::numeric_limits<int>::max();
        
    do getline(file, s); //Skip comments and AA line
    while (s[0] == '#');
    
    for (int i = 0; i < 24; i++) {
        file >> s; // skip AA
        for (int j = 0; j < 24; j++) {
            file >> PSM[i][j];
            PSM[i][j] = -PSM[i][j];
            minscore = min(minscore,PSM[i][j]);
        }
    }

    // renormalize to have only penalties
    for (int i = 0; i < 24; i++) 
        for (int j = 0; j < 24; j++) 
            PSM[i][j] -= minscore;
}

void Cpd::fillPSMbiases(size_t varIndex, vector<Cost> &biases)
{
    for (char c : rotamers2aa[varIndex]) {
        int bias = PSMBias*PSM[PSMIdx.find(c)->second][PSMIdx.find(nativeSequence[varIndex])->second];
        biases.push_back((Cost)bias);
    }
}


void Cpd::readPSSMatrix(const char* filename)
{
    ifstream file;
    file.open(filename);

    if (!file.is_open()) {
        cerr << "Could not open PSSM file, aborting." << endl;
        exit(EXIT_FAILURE);
    }
    
    string s;
	istringstream line;
    int minscore = std::numeric_limits<int>::max();

    getline(file, s); //Skip header line
    
    while (file) {
        int pos;
        vector<int> scores;
        char cons;
        file >> pos;
        file >> cons;

        scores.clear();
        for (int j = 0; j < 20; j++) {
            int  score;
            file >> score;
            score = -score;
            minscore = min(minscore,score);
            scores.push_back(score);
        }
        PSSM.push_back(scores);
    }
    
    // renormalize to have only penalties
    for (auto &v : PSSM)
        for (auto &i : v)
            i -= minscore;
}


void Cpd::fillPSSMbiases(size_t varIndex, vector<Cost> &biases)
{
    for (char c : rotamers2aa[varIndex]) {
        int bias = PSSMBias*PSSM[varIndex][PSSMIdx.find(c)->second];
        biases.push_back((Cost)bias);
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
