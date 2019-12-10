#include "tb2cpd.hpp"

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

const map<char, int> AminoMRF::AminoPMRFIdx = { { 'A', 0 }, { 'C', 1 }, { 'D', 2 }, { 'E', 3 }, { 'F', 4 }, { 'G', 5 },
    { 'H', 6 }, { 'I', 7 }, { 'K', 8 }, { 'L', 9 }, { 'M', 10 }, { 'N', 11 },
    { 'P', 12 }, { 'Q', 13 }, { 'R', 14 }, { 'S', 15 }, { 'T', 16 }, { 'V', 17 },
    { 'W', 18 }, { 'Y', 19 } };

const string AminoMRFs = "ARNDCQEGHILKMFPSTWYV-";
const string AminoPMRFs = "ACDEFGHIKLMNPQRSTVWY-";
constexpr int NumNatAA = 20;

// AminoMRF Class
// Read MRF trained from multiple alignment using CCMPred. The MRF (alignment) should have
// exactly the same number of variables/columns as the native protein.
AminoMRF::AminoMRF(const char* filename)
{
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

    // We don't normalize unaries as they will receive projections from binaries
    while (!binariesReached) {
        getline(file, s);
        if (s[0] == '#') {
            binariesReached = true;
            break;
        }
        stringstream ss(s);

        for (int i = 0; i < NumNatAA; i++) {
            ss >> LP;
            unaries[nv].push_back(-LP);
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

        // we just ignore the gap
        TLogProb minscore = std::numeric_limits<TLogProb>::max();
        TLogProb maxscore = std::numeric_limits<TLogProb>::min();
        pair<int, int> pv(n, m);
        for (int i = 0; i < NumNatAA + 1; i++) {
            if (i < NumNatAA)
                binaries[pv].resize(NumNatAA);

            for (int j = 0; j < NumNatAA + 1; j++) {
                file >> LP;
                if (i < NumNatAA && j < NumNatAA) {
                    binaries[pv][i].push_back(-LP);
                    minscore = min(minscore, -LP);
                    maxscore = max(maxscore, -LP);
                }
            }
        }

        // renormalize globally
        maxscore -= minscore;
        for (int i = 0; i < NumNatAA; i++) {
            for (int j = 0; j < NumNatAA; j++) {
                binaries[pv][i][j] -= minscore;
            }
        }

        if (maxscore > 1e-1)
            nPot++;
    } while (!file.eof());
    cout << "loaded evolutionary MRF with " << nVar << " residues and " << nPot << " coupled pairs (dev > 1e-1)\n";
}

// AminoMRF Class
// Read MRF trained from multiple alignment using PMRF. The MRF (alignment) should have
// exactly the same number of variables/columns as the native protein. The fmt is there
// to have different signatire and should be the expected fmt number of the file (ie. 1).
AminoMRF::AminoMRF(const char* filename, size_t fmt)
{
    static bool debug = false;

    ifstream file;
    file.open(filename);

    if (!file.is_open()) {
        cerr << "Could not open alignment MRF file, aborting." << endl;
        exit(EXIT_FAILURE);
    }

    size_t fmtnum;
    file.read((char*)&fmtnum, sizeof(size_t));
    if (fmtnum != fmt) {
        cerr << "toulbar2 only supports format 1 PMRF files\n";
        exit(EXIT_FAILURE);
    }
    if (debug)
        cout << "fmtnum = " << fmtnum << endl;

    float neff; // effective number of sequences in the alignment
    file.read((char*)&neff, sizeof(float));
    if (debug)
        cout << "neff = " << neff << endl;

    // Reading native sequence
    size_t n, n2; //n is native length, n2 is the number of edges
    file.read((char*)&n, sizeof(size_t));
    std::getline(file, ToulBar2::cpd->nativeSequence, '\0');
    // Reading MRF architecture
    file.read((char*)&n2, sizeof(size_t));
    if (debug) {
        cout << "Native: " << ToulBar2::cpd->nativeSequence << endl;
        cout << "Variables: " << n << endl;
        cout << "Edges: " << n2 << endl;
    }
    nVar = n;

    vector<pair<size_t, size_t>> edgeList;
    size_t idx1, idx2;
    for (size_t i = 0; i < n2; ++i) {
        file.read((char*)&idx1, sizeof(size_t));
        file.read((char*)&idx2, sizeof(size_t));
        if (debug)
            cout << "Edge " << idx1 << " - " << idx2 << endl;
        edgeList.push_back(make_pair(idx1, idx2));
    }
    // Skipping sequence profile for now - TODO usable for PSSM
    file.ignore(n * (NumNatAA + 1) * sizeof(float));

    // Reading unaries
    // We don't normalize unaries as they will receive projections from binaries
    float tmpUnary[NumNatAA + 1];
    for (size_t i = 0; i < n; ++i) {
        file.read((char*)tmpUnary, (NumNatAA + 1) * sizeof(float));
        unaries[i].resize(NumNatAA);
        for (size_t j = 0; j < NumNatAA + 1; ++j) {
            // We need to convert the PMRF idx to a CCMpred idx
            unaries[i][AminoMRFIdx.find(AminoPMRFs[j])->second] = -tmpUnary[j];
        }
        if (debug)
            cout << "Read unary on variable " << i << endl;
    }

    // Reading binaries
    float tmpBinary[(NumNatAA + 1) * (NumNatAA + 1)];
    float minscore = std::numeric_limits<float>::max();
    float maxscore = std::numeric_limits<float>::min();

    for (pair<size_t, size_t> e : edgeList) {
        file.read((char*)tmpBinary, (NumNatAA + 1) * (NumNatAA + 1) * sizeof(float));
        binaries[e].resize(NumNatAA);

        for (int i = 0; i < NumNatAA + 1; i++) {
            if (i < NumNatAA)
                binaries[e][i].resize(NumNatAA);

            for (int j = 0; j < NumNatAA + 1; j++) {
                if (i < NumNatAA && j < NumNatAA) {
                    float LP = -tmpBinary[(AminoPMRFIdx.find(AminoMRFs[i])->second * NumNatAA) + AminoPMRFIdx.find(AminoMRFs[j])->second];
                    binaries[e][i][j] = LP;
                    minscore = min(minscore, LP);
                    maxscore = max(maxscore, LP);
                }
            }
        }

        // renormalize globally
        maxscore -= minscore;
        for (int i = 0; i < NumNatAA; i++) {
            for (int j = 0; j < NumNatAA; j++) {
                binaries[e][i][j] -= minscore;
            }
        }

        if (maxscore > 1e-1)
            nPot++;
        if (debug)
            cout << "Read binary on " << e.first << " - " << e.second << endl;
    }

    cout << "loaded evolutionary MRF with " << nVar << " residues and " << nPot << " coupled pairs (dev > 1e-1)\n";
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

TLogProb AminoMRF::eval(const string& sequence, const vector<Variable*>& vars)
{
    static bool debug = false;

    size_t nbCFNVars = vars.size();
    TLogProb paid = 0.0;
    for (size_t varIdx1 = 0; varIdx1 < nbCFNVars; varIdx1++) {
        int pos1 = vars[varIdx1]->getPosition();
        int code1 = AminoMRFIdx.find(sequence[varIdx1])->second;
        if (debug) {
            cout << "Unary at position " << pos1 << " amino acid code " << code1 << endl;
        }
        paid += unaries[pos1][code1];
        for (size_t varIdx2 = varIdx1 + 1; varIdx2 < nbCFNVars; varIdx2++) {
            int pos2 = vars[varIdx2]->getPosition();
            int code2 = AminoMRFIdx.find(sequence[varIdx2])->second;
            if (debug)
                cout << "Binary at positions " << pos1 << " - " << pos2 << " amino acid codes " << code1 << " - " << code2 << endl;

            pair<int, int> pv(pos1, pos2);
            if (binaries.count(pv))
                paid += binaries[pv][code1][code2];
        }
    }
    return paid;
}

// Must be called after the problem is loaded, CFN format only.
void AminoMRF::Penalize(WeightedCSP* pb, TLogProb CMRFBias)
{
    const static bool debug = false;

    if (ToulBar2::cpd->nativeSequence.empty()) {
        cerr << "Error: the native sequence must be provided using --native.\n";
        exit(1);
    }

    // check residue numbers
    if (ToulBar2::cpd->nativeSequence.size() != nVar) {
        cerr << "Error: the loaded evolutionary MRF has not the size of the native protein.\n";
        exit(1);
    }

    // project and process binaries
    map<int, size_t> posList; // map of designed position to varIdx
    for (size_t varIdx = 0; varIdx < pb->numberOfVariables(); varIdx++) {
        int pos = pb->getVars()[varIdx]->getPosition();
        bool isAA = pb->getVars()[varIdx]->getName()[0] != 'Z';

        if (debug)
            cout << "Variable " << pb->getVars()[varIdx]->getName() << " has position " << pos << endl;

        if ((pos < 0 || (unsigned int)pos >= nVar) && isAA) {
            cerr << "Variable " << pb->getVars()[varIdx]->getName() << " has an out-of-range sequence position (wrt. native)" << endl;
        }
        if (isAA)
            posList[pos] = varIdx;
        if (debug)
            cout << "Position " << pb->getVars()[varIdx]->getPosition() << " on var " << varIdx << endl;
    }

    // parse all binaries with one non mutable variable and project on the corresponding unary
    for (auto const& bincf : binaries) {
        int pos1 = bincf.first.first;
        int pos2 = bincf.first.second;

        if (debug)
            cout << "MRF pair " << pos1 << "-" << pos2 << " ";

        if (posList.count(pos1) + posList.count(pos2) == 1) { // different types (designable/absent)
            if (debug)
                cout << "must be projected." << endl;
            if (posList.count(pos2)) { // project on second
                char natAA = ToulBar2::cpd->nativeSequence[pos1];
                int natIdx = AminoMRFIdx.find(natAA)->second;
                if (debug)
                    cout << "Projecting on position " << pos2 << " using AA " << natAA << " for position " << pos1 << endl;
                for (size_t i = 0; i < NumNatAA; i++) {
                    unaries[pos2][i] += bincf.second[natIdx][i];
                }
            } else { // project on first
                char natAA = ToulBar2::cpd->nativeSequence[pos2];
                int natIdx = AminoMRFIdx.find(natAA)->second;
                if (debug)
                    cout << "Projecting on position " << pos1 << " using " << natAA << " for position " << pos2 << endl;
                for (size_t i = 0; i < NumNatAA; i++) {
                    unaries[pos1][i] += bincf.second[i][natIdx];
                }
            }
        }

        if (posList.count(pos1) + posList.count(pos2) == 2) { // both designable

            bool warn = true;
            vector<Cost> biases;
            int varIdx1 = posList[pos1];
            int varIdx2 = posList[pos2];
            if (debug)
                cout << "will be normalized and posted on " << varIdx1 << "-" << varIdx2 << endl;

            for (char c1 : ToulBar2::cpd->getRotamers2AA()[varIdx1]) {
                for (char c2 : ToulBar2::cpd->getRotamers2AA()[varIdx2]) {
                    int valIdx1 = AminoMRFIdx.find(c1)->second;
                    int valIdx2 = AminoMRFIdx.find(c2)->second;
                    if (bincf.second[valIdx1][valIdx2] == 0.0)
                        warn = false;
                    Cost bias = Round(powl(10, ToulBar2::decimalPoint) * ToulBar2::costMultiplier * CMRFBias * bincf.second[valIdx1][valIdx2]);
                    biases.push_back(bias);
                }
            }
            if (warn && ToulBar2::verbose > 0) {
                cout << "WARNING: the preferred amino acid pair (";
                for (int i = 0; i < NumNatAA; i++)
                    for (int j = 0; j < NumNatAA; j++)
                        if (bincf.second[i][j] == 0.0)
                            cout << "[" << AminoMRFs[i] << AminoMRFs[j] << "]";
                cout << ") has been excluded from design at residues " << varIdx1 + 1 << "-" << varIdx2 + 1 << endl;
            }
            pb->postBinaryConstraint(varIdx1, varIdx2, biases);
        }
        if (debug && (posList.count(pos1) + posList.count(pos2) == 0))
            cout << " will be ignored" << endl;
    }

    // process unaries
    for (size_t varIdx = 0; varIdx < pb->numberOfVariables(); varIdx++) {
        bool warn = true;
        vector<Cost> biases;
        int pos = pb->getVars()[varIdx]->getPosition();
        if (pb->getVars()[varIdx]->getName()[0] == 'Z')
            break;

        if (debug)
            cout << "Processing unary MRF potential on position " << pos << ", variable " << varIdx << endl;

        // Normalize
        TLogProb minscore = std::numeric_limits<TLogProb>::max();
        for (int i = 0; i < NumNatAA; i++) {
            minscore = min(minscore, unaries[pos][i]);
        }
        for (int i = 0; i < NumNatAA; i++) {
            unaries[pos][i] -= minscore;
        }
        if (debug)
            cout << "Normalized." << endl;

        // Insert in the CFN model
        for (char c : ToulBar2::cpd->getRotamers2AA()[varIdx]) {
            int valIdx = AminoMRFIdx.find(c)->second;
            if (unaries[pos][valIdx] == 0.0)
                warn = false;
            Cost bias = Round(powl(10, ToulBar2::decimalPoint) * ToulBar2::costMultiplier * CMRFBias * unaries[pos][valIdx]);

            biases.push_back(bias);
        }
        if (warn) {
            cout << "WARNING: the preferred amino acid (";
            for (int i = 0; i < NumNatAA; i++) {
                if (unaries[pos][i] == 0.0)
                    cout << AminoMRFs[i];
            }
            cout << ") has been excluded from design at residue " << varIdx + 1 << endl;
        }
        pb->postUnaryConstraint(varIdx, biases);
        if (debug)
            cout << "Posted." << endl;
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

void Cpd::read_rotamers2aa(istream& file, vector<Variable*>& vars)
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
    computeAAbounds();
}

void Cpd::computeAAbounds()
{
    LeftAA.clear();
    RightAA.clear();

    for (auto& rv : rotamers2aa) {

        vector<Value> leftidx_var, rightidx_var;
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
    static bool debug = false;

    ifstream file;
    file.open(filename);

    if (!file.is_open()) {
        cerr << "Could not open PSM file, aborting." << endl;
        exit(EXIT_FAILURE);
    }
    int minscore = std::numeric_limits<int>::max();

    for (int i = 0; i < 20; i++) {
        for (int j = 0; j < 20; j++) {
            file >> PSM[i][j];
            if (debug)
                cout << i << " " << j << " " << PSM[i][j] << endl;
            PSM[i][j] = -PSM[i][j];
            minscore = min(minscore, PSM[i][j]);
        }
    }
    if (debug)
        cout << "Neg-shifting by " << minscore << endl;

    // renormalize to have only penalties
    // the indices 20, 21, 22 for Z/X/* are
    // left to their initial zero value for now (don't care)
    for (int i = 0; i < 20; i++)
        for (int j = 0; j < 20; j++)
            PSM[i][j] -= minscore;
}

void Cpd::fillPSMbiases(size_t varIndex, vector<Cost>& biases)
{
    const bool debug = false;

    if (debug) {
        cout << "natbias on var " << varIndex << " ";
    }

    for (char c : rotamers2aa[varIndex]) {
        Cost bias = (Cost)powl(10, ToulBar2::decimalPoint) * ToulBar2::costMultiplier * PSMBias * PSM[PSMIdx.find(c)->second][PSMIdx.find(nativeSequence[varIndex])->second];
        biases.push_back((Cost)bias);
        if (debug)
            cout << bias << " ";
    }
    if (debug)
        cout << endl;
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
    int pos;

    while (file >> pos) {
        vector<int> scores;
        char cons;
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
        for (auto& i : v) {
            i -= minscore;
        }
}

void Cpd::fillPSSMbiases(size_t varIndex, vector<Cost>& biases)
{
    for (char c : rotamers2aa[varIndex]) {
        int bias = PSSMBias * PSSM[varIndex][PSSMIdx.find(c)->second];
        biases.push_back((Cost)bias);
    }
}

void Cpd::storeSequence(const vector<Variable*>& vars, Double energy)
{
    string sequence;
    for (size_t i = 0; i < vars.size(); i++) {
        char aa = rotamers2aa[i][vars[i]->getValue()];
        if (aa != '*')
            sequence.push_back(aa);
    }
    cpdtrie.insert_sequence(sequence, energy);
}

void Cpd::printSequences()
{
    cpdtrie.print_tree();
}

void Cpd::printSequence(const vector<Variable*>& vars, Double energy)
{
    string sequence;
    size_t mutations = 0;

    cout << "New rotamers:";
    for (size_t i = 0; i < vars.size(); i++) {
        char aa = rotamers2aa[i][vars[i]->getValue()];
        if (aa != '*') {
            sequence.push_back(aa);
            cout << " " << vars[i]->getValue();
        }
        mutations += (aa != vars[i]->getNativeResidue());
    }
    cout << "\nNew sequence: " << sequence << " Mutations: " << mutations << " Energy: " << std::setprecision(ToulBar2::decimalPoint);
    if (AminoMRFBias != 0.0) {
        Double Evol = AminoMRFBias * AminoMat->eval(sequence, vars);
        cout << energy - Evol << " (evol " << Evol << ", joint " << energy << ")";
    } else {
        cout << energy;
    }

    cout << endl;
}

void Cpd::printSequence(TAssign& assig, const vector<Variable*>& vars)
{
    string sequence;
    size_t mutations = 0;

    cout << "New rotamers:";
    for (size_t i = 0; i < assig.size(); i++) {
        char aa = rotamers2aa[i][assig[i]];
        if (aa != '*') {
            sequence.push_back(aa);
            cout << " " << assig[i];
        }
        mutations += (aa != vars[i]->getNativeResidue());
    }
    cout << "\nNew sequence: " << sequence << " Mutations: " << mutations;
    if (AminoMRFBias != 0.0)
        cout << " (evol " << AminoMRFBias * AminoMat->eval(sequence, vars) << ")";
    cout << endl;
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
