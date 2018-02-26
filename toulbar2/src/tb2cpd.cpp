#include "tb2cpd.hpp"
#include <sstream>

const map<char, int> Cpd::PSMIdx = { { 'A', 0 }, { 'R', 1 }, { 'N', 2 }, { 'D', 3 }, { 'C', 4 }, { 'Q', 5 },
    { 'E', 6 }, { 'G', 7 }, { 'H', 8 }, { 'I', 9 }, { 'L', 10 }, { 'K', 11 },
    { 'M', 12 }, { 'F', 13 }, { 'P', 14 }, { 'S', 15 }, { 'T', 16 }, { 'W', 17 },
    { 'Y', 18 }, { 'V', 19 }, { 'B', 20 }, { 'Z', 21 }, { 'X', 22 }, { '*', 23 } };

const map<char, int> Cpd::PSSMIdx = { { 'A', 0 }, { 'G', 1 }, { 'I', 2 }, { 'L', 3 }, { 'V', 4 }, { 'M', 5 },
    { 'F', 6 }, { 'W', 7 }, { 'P', 8 }, { 'C', 9 }, { 'S', 10 }, { 'T', 11 },
    { 'Y', 12 }, { 'N', 13 }, { 'Q', 14 }, { 'H', 15 }, { 'K', 16 }, { 'R', 17 },
    { 'D', 18 }, { 'E', 19 } };

const map<char, int> AminoMRF::AminoMRFIdx = { { 'A', 0 }, { 'R', 1 }, { 'N', 2 }, { 'D', 3 }, { 'C', 4 }, { 'Q', 5 },
    { 'E', 6 }, { 'G', 7 }, { 'H', 8 }, { 'I', 9 }, { 'L', 10 }, { 'K', 11 },
    { 'M', 12 }, { 'F', 13 }, { 'P', 14 }, { 'S', 15 }, { 'T', 16 }, { 'W', 17 },
    { 'Y', 18 }, { 'V', 19 } };

// AminoMRF Class
// Read MRF trained from multiple alignment. The MRF (alignment) should have
// exactly the same number of variables/columns as the currently solved design problem.
AminoMRF::AminoMRF(const char* filename)
{
    constexpr int NumNatAA = 20;

    ifstream file;
    file.open(filename);

    if (!file.is_open()) {
        cerr << "Could not open alignment MRF file, aborting." << endl;
        exit(EXIT_FAILURE);
    }

    TLogProb LP;
    int nv = 0;
    // read unaries
    bool binariesReached = false;
    string s;

    while (!binariesReached) {
        getline(file, s);
        if (s[0] == '#') {
            binariesReached = true;
            break;
        }

        stringstream ss(s);
        ss.ignore();

        TLogProb minscore = std::numeric_limits<TLogProb>::max();
        for (int i = 0; i < NumNatAA; i++) {
            ss >> LP;
            unaries[nv].push_back(LP);
            minscore = min(minscore, LP);
        }
        for (int i = 0; i < NumNatAA; i++) {
            unaries[nv][i] -= minscore;
        }
        nv++;
    }

    nVar = nv;

    // read binaries
    int n, m;
    nPot = 0;
    do {
        while (s[0] != '#' && getline(file, s)) {
        }
        if (file.eof())
            break;

        stringstream ss(s);
        s.clear();
        ss.ignore();

        ss >> n;
        ss >> m;

        TLogProb minscore = std::numeric_limits<TLogProb>::max();
        pair<int, int> pv(n, m);
        for (int i = 0; i < NumNatAA; i++) {
            binaries[pv].resize(NumNatAA);
            for (int j = 0; j < NumNatAA; j++) {
                file >> LP;
                binaries[pv][i].push_back(LP);
                minscore = min(minscore, LP);
            }
        }

        for (int i = 0; i < NumNatAA; i++) {
            for (int j = 0; j < NumNatAA; j++) {
                binaries[pv][i][j] -= minscore;
            }
        }
        nPot++;
    } while (!file.eof());
    cout << "loaded evolutionary MRF with " << nVar << " residues and " << nPot << " correlated pairs\n";
}

AminoMRF::~AminoMRF()
{
}

TLogProb AminoMRF::getUnary(int var, int AAidx)
{
    return unaries[var][AAidx];
}

TLogProb AminoMRF::getBinary(int var1, int var2, int AAidx1, int AAidx2)
{
    if (var1 > var2) {
        swap(var1, var2);
        swap(AAidx1, AAidx2);
    }

    pair<int, int> pv(var1, var2);

    return binaries[pv][AAidx1][AAidx2];
}

void AminoMRF::Penalize(WeightedCSP* pb, TLogProb CMRFBias)
{
    // check residue numbers
    if (pb->numberOfVariables() < nVar) {
        cerr << "The loaded evolutionary MRF has more variables than the number of variables in the problem\nAborting\n";
        exit(1);
    }
    if (pb->numberOfVariables() > nVar) {
        cout << "WARNING: the loaded evolutionary MRF has less variables than in the problem. Extra variables won't be penalized\n";
        exit(1);
    }
    

    // process unaries
    for (size_t varIdx = 0; varIdx < pb->numberOfVariables(); varIdx++) {
        bool warn = true;
        vector<Cost> biases;

        for (char c : ToulBar2::cpd->getRotamers2AA()[varIdx]) {
            int valIdx = AminoMRFIdx.find(c)->second;
            if (unaries[varIdx][valIdx] == 0.0)
                warn = false;
            Cost bias = CMRFBias * unaries[varIdx][valIdx];

            biases.push_back(bias);
        }
        if (warn) {
            cout << "WARNING: the preferred amino acid has been excluded from design at residue " << varIdx + 1 << endl;
        }
        pb->postUnaryConstraint(varIdx, biases);
    }

    // process binaries
    for (auto const& bincf : binaries) {
        int varIdx1 = bincf.first.first;
        int varIdx2 = bincf.first.second;
        bool warn = true;
        vector<Cost> biases;

        for (char c1 : ToulBar2::cpd->getRotamers2AA()[varIdx1]) {
            for (char c2 : ToulBar2::cpd->getRotamers2AA()[varIdx2]) {
                int valIdx1 = AminoMRFIdx.find(c1)->second;
                int valIdx2 = AminoMRFIdx.find(c2)->second;
                if (bincf.second[valIdx1][valIdx2] == 0.0)
                    warn = false;
                Cost bias = CMRFBias * bincf.second[valIdx1][valIdx2];
                biases.push_back(bias);
            }
        }
        if (warn) {
            cout << "WARNING: the preferred amino acid pair has been excluded from design at residues " << varIdx1 + 1 << "-" << varIdx2 + 1 << endl;
        }
        pb->postBinaryConstraint(varIdx1, varIdx2, biases);
    }
}

// CPD Class
Cpd::Cpd()
{
}

Cpd::~Cpd()
{
}

void Cpd::init()
{
    rotamers2aa.clear();
    LeftAA.clear();
    RightAA.clear();
    cpdtrie.init();
    AminoMat = NULL;
}

void Cpd::read_rotamers2aa(ifstream& file, vector<Variable*>& vars)
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
            if (!isspace(current_char))
                rot2aa_var.push_back(current_char);
        }
        if (rot2aa_var.size() != 0)
            rotamers2aa.push_back(rot2aa_var);
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
            unsigned int initsize = dynamic_cast<EnumeratedVariable*>(vars[i])->getDomainInitSize();
            if (rotamers2aa[i].size() != initsize) {
                cout << "Wrong domain size " << rotamers2aa[i].size() << " " << initsize << " of variable " << dynamic_cast<EnumeratedVariable*>(vars[i])->getName() << endl;
                throw 2;
            }
        }
    }

    for (auto& rv : rotamers2aa) {

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
        for (int i = rv.size() - 1; i >= 0; i--) {
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
    int minscore = std::numeric_limits<int>::max();

    do
        getline(file, s); //Skip comments and AA line
    while (s[0] == '#');

    for (int i = 0; i < 24; i++) {
        file >> s; // skip AA
        for (int j = 0; j < 24; j++) {
            file >> PSM[i][j];
            PSM[i][j] = -PSM[i][j];
            minscore = min(minscore, PSM[i][j]);
        }
    }

    // renormalize to have only penalties
    for (int i = 0; i < 24; i++)
        for (int j = 0; j < 24; j++)
            PSM[i][j] -= minscore;
}

void Cpd::fillPSMbiases(size_t varIndex, vector<Cost>& biases)
{
    for (char c : rotamers2aa[varIndex]) {
        int bias = PSMBias * PSM[PSMIdx.find(c)->second][PSMIdx.find(nativeSequence[varIndex])->second];
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
            int score;
            file >> score;
            score = -score;
            minscore = min(minscore, score);
            scores.push_back(score);
        }
        PSSM.push_back(scores);
    }

    // renormalize to have only penalties
    for (auto& v : PSSM)
        for (auto& i : v)
            i -= minscore;
}

void Cpd::fillPSSMbiases(size_t varIndex, vector<Cost>& biases)
{
    for (char c : rotamers2aa[varIndex]) {
        int bias = PSSMBias * PSSM[varIndex][PSSMIdx.find(c)->second];
        biases.push_back((Cost)bias);
    }
}

void Cpd::storeSequence(const vector<Variable*>& vars, Cost _cost)
{
    string sequence;
    for (size_t i = 0; i < vars.size(); i++) {
        char aa = rotamers2aa[i][vars[i]->getValue()];
        if (aa != '*')
            sequence.push_back(aa);
    }
    cpdtrie.insert_sequence(sequence, _cost);
}

void Cpd::printSequences()
{
    cpdtrie.print_tree();
}

void Cpd::printSequence(const vector<Variable*>& vars, Cost _cost)
{
    string sequence;
    cout << "New rotamers:";
    for (size_t i = 0; i < vars.size(); i++) {
        char aa = rotamers2aa[i][vars[i]->getValue()];
        if (aa != '*') {
            sequence.push_back(aa);
            cout << " " << vars[i]->getValue();
        }
    }
    cout << "\nNew sequence: " << sequence << " Cost: " << _cost << endl;
}

void Cpd::printSequence(TAssign& vars)
{
    string sequence;
    cout << "New rotamers:";
    for (size_t i = 0; i < vars.size(); i++) {
        char aa = rotamers2aa[i][vars[i]];
        if (aa != '*') {
            sequence.push_back(aa);
            cout << " " << vars[i];
        }
    }
    cout << "\nNew sequence: " << sequence << endl;
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
