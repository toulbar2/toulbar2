#ifndef TB2CPD_HPP_
#define TB2CPD_HPP_

#include "core/tb2wcsp.hpp"
#include "tb2trie.hpp"

// Storing an MRF to capture evolutionary information
// could be a WCSP but simpler: all variables have 20 values (hard encoded) and
// the number of variables should be exactly the same as in the designed protein.
// So it's not compatible with rigid positions for now.
class AminoMRF {
public:
    AminoMRF(const char* filename); // read from CCMPredfile
    AminoMRF(const char* filename, size_t fmtnum); // read from PMRF binary file. The second argument must be equal to 1
    ~AminoMRF();
    size_t nVar = 0;
    size_t nPot = 0;
    TLogProb getUnary(int var, int AAidx);
    TLogProb getBinary(int var1, int var2, int AAidx1, int AAidx2);
    TLogProb eval(const string& sequence, const vector<Variable*>& vars);
    void Penalize(WeightedCSP* pb, TLogProb biasStrength);

private:
    map<int, vector<TLogProb>> unaries;
    map<pair<int, int>, vector<vector<TLogProb>>> binaries;
    static const map<char, int> AminoMRFIdx;
    static const map<char, int> AminoPMRFIdx;
};

class Cpd {
public:
    Cpd();
    ~Cpd();
    void init();
    void read_rotamers2aa(istream& file, vector<Variable*>& vars); ///< \brief read rotamer to amino acid array from wcsp file
    void computeAAbounds(); ///< \brief compute beg/end of each AA in the rotamer array
    void newRotamerArray(const vector<char>& rots) { rotamers2aa.push_back(rots); };
    void readPSMatrix(const char* filename);
    void readPSSMatrix(const char* filename);
    void readEvolMat(const char* filename);
    void fillPSMbiases(size_t varIndex, vector<Cost>& biases);
    void fillPSSMbiases(size_t varIndex, vector<Cost>& biases);
    void storeSequence(const vector<Variable*>& vars, Double energy);
    void printSequences();
    void printSequence(const vector<Variable*>& vars, Double energy);
    void printSequence(TAssign& assig, const vector<Variable*>& vars);
    int getTotalSequences() { return cpdtrie.getTotalSequences(); }
    vector<vector<char>>& getRotamers2AA() { return rotamers2aa; }
    string getMutationDomain(size_t varIdx) { return mutationDomains[varIdx]; }
    bool isAAVariable(const Variable* var) { return (std::count(NotAnAA.begin(), NotAnAA.end(), var->getName()[0]) == 0); }
    bool isSeqVariable(const Variable* var) { return ((std::count(seqVarPrefixes.begin(), seqVarPrefixes.end(), var->getName()[0]) != 0) && isdigit(var->getName()[1])); }
    char getAA(int varIndex, Value value) { return rotamers2aa[varIndex][value]; }
    Value getLeft(int varIndex, Value value) { return LeftAA[varIndex][value]; }
    Value getRight(int varIndex, Value value) { return RightAA[varIndex][value]; }
    size_t rot2aaSize(int varIndex) { return rotamers2aa[varIndex].size(); }
    string nativeSequence;
    AminoMRF* AminoMat;
    size_t PSSMlen() { return PSSM.size(); };
    float PSMBias = 0.0;
    float PSSMBias = 0.0;
    float AminoMRFBias = 0.0;

private:
    const string seqVarPrefix1 = "Z";
    const string seqVarPrefix2 = DIVERSIFIED_VAR_TAG; // We use diversity on sequence variables.
    const string seqVarPrefixes = seqVarPrefix1 + seqVarPrefix2;
    const string NotAnAA = seqVarPrefix1 + HIDEABLE_VAR_TAGS; // A list of chars indicating a variable that does not define an AA identiy (as first char of the variable name).
    const static map<char, int> PSMIdx; // converts AA char to indices in PSMatrix
    const static map<char, int> PSSMIdx; // converts AA char to indices in PsiBlast PSSMatrix
    TrieCpd cpdtrie;
    vector<string> mutationDomains;
    vector<vector<char>> rotamers2aa;
    vector<vector<Value>> LeftAA;
    vector<vector<Value>> RightAA;
    int PSM[24][24] = {
        { 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0 },
        { 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0 },
        { 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0 },
        { 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0 },
        { 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0 },
        { 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0 },
        { 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0 },
        { 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }
    };
    vector<vector<int>> PSSM;
};

#endif

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
