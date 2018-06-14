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
    AminoMRF(const char* filename); // read from file
    ~AminoMRF();
    size_t nVar;
    size_t nPot;
    TLogProb getUnary(int var, int AAidx);
    TLogProb getBinary(int var1, int var2, int AAidx1, int AAidx2);
    void Penalize(WeightedCSP* pb, TLogProb biasStrength);

private:
    map<int, vector<TLogProb>> unaries;
    map<pair<int, int>, vector<vector<TLogProb>>> binaries;
    static const map<char, int> AminoMRFIdx;
};

class Cpd {
public:
    Cpd();
    ~Cpd();
    void init();
    void read_rotamers2aa(ifstream& file, vector<Variable*>& vars); ///< \brief read rotamer to amino acid array from wcsp file
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
    void printSequence(TAssign& vars);
    int getTotalSequences() { return cpdtrie.getTotalSequences(); }
    vector<vector<char>>& getRotamers2AA() { return rotamers2aa; }
    char getAA(int varIndex, Value value) { return rotamers2aa[varIndex][value]; }
    Value getLeft(int varIndex, Value value) { return LeftAA[varIndex][value]; }
    Value getRight(int varIndex, Value value) { return RightAA[varIndex][value]; }
    size_t rot2aaSize(int varIndex) { return rotamers2aa[varIndex].size(); }
    char* nativeSequence = NULL;
    AminoMRF* AminoMat;
    bool isPSSMlen() { return PSSM.size(); };
    int PSMBias = 0;
    int PSSMBias = 0;
    TLogProb AminoMRFBias = 0.0;

private:
    const static map<char, int> PSMIdx; // converts AA char to indices in PSMatrix
    const static map<char, int> PSSMIdx; // converts AA char to indices in PsiBlast PSSMatrix
    TrieCpd cpdtrie;
    vector<vector<char>> rotamers2aa;
    vector<vector<Value>> LeftAA;
    vector<vector<Value>> RightAA;
    int PSM[24][24] = { { 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0 } };
    vector<vector<int>> PSSM;
};

#endif

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
