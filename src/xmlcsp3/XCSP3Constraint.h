/*=============================================================================
 * parser for CSP instances represented in XCSP3 Format
 * 
 * Copyright (c) 2015 xcsp.org (contact <at> xcsp.org)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *=============================================================================
 */
#ifndef XCONSTRAINT_H
#define    XCONSTRAINT_H

#include "XCSP3Variable.h"
#include "XCSP3utils.h"
#include "XCSP3Constants.h"
#include <typeinfo>
#include <regex>
#include<map>

namespace XCSP3Core {

    using namespace std;

    class XConstraintGroup;

    class XConstraintIntension;

    /***************************************************************************
     * A constraint XValue
     **************************************************************************/
    class XConstraint {
    public :
        std::string id;
        vector<XVariable *> list;
        string classes;


        XConstraint(std::string dd, std::string c) : id(dd), classes(c) {}


        virtual ~XConstraint() {}


        virtual void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original);
    };


    /***************************************************************************
     * A constraint group
     * Only basic abstraction is implemented: All parameters occur in the same tag
     **************************************************************************/

    class XConstraintGroup : public XConstraint {
    public :
        XConstraint *constraint;                 // A Constraint group contains a constraint
        vector<vector<XVariable *> > arguments; // The list of all arguments
        ConstraintType type;                     // Use it to discover the type of constraint... and perform cast
        map<string, XVariable *> toArguments;


        XConstraintGroup(std::string idd, std::string c) : XConstraint(idd, c), constraint(NULL), type(UNKNOWN) {}


        virtual ~XConstraintGroup() { delete constraint; }


        void unfoldVector(vector<XVariable *> &toUnfold, vector<XVariable *> &args, vector<XVariable *> &initial);
        void unfoldString(string &toUnfold, vector<XVariable *> &args);
        void unfoldArgumentNumber(int, XConstraint *builtConstraint);
    };



    // All different structures related to XVariables used in constraints
    // This facilitates the unfolding of group constraint


    class XCondition {
    public :
        OrderType op;
        OperandType operandType;
        int val;
        int min, max;
        string var;
        std::vector<int> set;

        friend ostream &operator<<(ostream &f, const XCSP3Core::XCondition &ie);
    };


    class XInitialCondition {
    public :
        string condition;

        virtual void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original);
        void extractCondition(XCondition &xc);  // Create the op and the operand (which can be a value, an interval or a XVariable)
        static void extract(XCondition &xc, string &c);
    };

    class XValues {
    public :
        vector<XVariable *> values; // Only XVariable or XInteger

        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original);

    };

    class XValue {
    public :
        XVariable *value;


        XValue() : value(0) {}


        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original);
    };

    class XIndex {
    public :
        XVariable *index;


        XIndex() : index(0) {}


        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original);
    };

    class XLengths {
    public :
        vector<XVariable *> lengths;

        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original);
    };


    /***************************************************************************
     * constraint extension
     **************************************************************************/
    class XConstraintExtension : public XConstraint {

    public :
        vector<vector<int> > tuples;
        bool isSupport;
        bool containsStar;


        XConstraintExtension(std::string idd, std::string c) : XConstraint(idd, c), containsStar(false) {}


        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original) override;
    };

    /***************************************************************************
     * constraint intension
     **************************************************************************/
    class XConstraintIntension : public XConstraint {

    public :
        string function;


        XConstraintIntension(std::string idd, std::string c) : XConstraint(idd, c) {}


        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original) override;
    };






    /***************************************************************************
     ****************************************************************************
     *                  CONSTRAINTS DEFINED ON LANGAGES
     ****************************************************************************
     ***************************************************************************/

    /***************************************************************************
     * constraint regular
     **************************************************************************/
    class XTransition {
    public :
        XTransition(string f, int v, string t) : from(f), val(v), to(t) {}


        string from;
        int val;
        string to;
    };

    extern vector<XTransition> tr; // Not beautiful but remove code to fixed data in group constraint.
    extern string st;
    extern vector<string> fi;

    class XConstraintRegular : public XConstraint {
    public :
        string &start;
        vector<string> &final;
        vector<XTransition> &transitions;

        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original) override;
        XConstraintRegular(std::string idd, std::string c) : XConstraint(idd, c), start(st), final(fi), transitions(tr) {}
    };

    /***************************************************************************
     * constraint mdd
     **************************************************************************/

    class XConstraintMDD : public XConstraint {
    public :
        vector<XTransition> &transitions;


        XConstraintMDD(std::string idd, std::string c) : XConstraint(idd, c), transitions(tr) {}
    };

    /***************************************************************************
     ****************************************************************************
     *                  COMPARISON BASED CONSTRAINTS
     ****************************************************************************
     ***************************************************************************/
    extern vector<int> ___except;

    class XConstraintAllDiff : public XConstraint, public XValues {
    // Values refer to except values
    public :


        XConstraintAllDiff(std::string idd, std::string c) : XConstraint(idd, c) {}
        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original) override;

    };


    class XConstraintAllDiffMatrix : public XConstraint {
    public :
        vector<vector<XVariable *> > matrix;


        XConstraintAllDiffMatrix(std::string idd, std::string c) : XConstraint(idd, c) {}


        XConstraintAllDiffMatrix(std::string idd, std::string c, vector<vector<XVariable *> > &mat) : XConstraint(idd, c) {
            matrix.resize(mat.size());
            for(unsigned int i = 0 ; i < mat.size() ; i++)
                matrix[i].assign(mat[i].begin(), mat[i].end());
        }


        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original) override;
    };

    typedef XConstraintAllDiffMatrix XConstraintAllDiffList;

    typedef XConstraint XConstraintAllEqual;

    /***************************************************************************
     * constraint ordered and lex
     **************************************************************************/

    extern OrderType _op;

    class XConstraintOrdered : public XConstraint, public XLengths {
    public :
        OrderType &op;

        XConstraintOrdered(std::string idd, std::string c) : XConstraint(idd, c), op(_op) {}
        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original) override;
    };


    class XConstraintLex : public XConstraintOrdered {
    public :
        // list is cleard and all lists are stored in lists
        vector<vector<XVariable *> > lists;


        XConstraintLex(std::string idd, std::string c) : XConstraintOrdered(idd, c) {}


        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original) override;
    };

    class XConstraintLexMatrix : public XConstraintOrdered {
    public :
        vector<vector<XVariable *> > matrix;


        XConstraintLexMatrix(std::string idd, std::string c) : XConstraintOrdered(idd, c) {}


        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original) override;
    };

    /***************************************************************************
     ****************************************************************************
     *                  SUMMING and COUNTING CONSTRAINTS
     ****************************************************************************
     ***************************************************************************/


    /***************************************************************************
     * constraint sum
     **************************************************************************/

//static vector<int> _coeffs;

    class XConstraintSum : public XConstraint, public XInitialCondition, public XValues {
    public :
        //     vector<int> &coeffs;


        XConstraintSum(std::string idd, std::string c) : XConstraint(idd, c) {}


        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original) override;
    };


    class XConstraintNValues : public XConstraint, public XInitialCondition {
    public :
        vector<int> &except;


        XConstraintNValues(std::string idd, std::string c) : XConstraint(idd, c), except(___except) {}


        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original) override;
    };


    class XConstraintCardinality : public XConstraint, public XValues {
    public :
        vector<XVariable *> occurs;
        bool closed;


        XConstraintCardinality(std::string idd, std::string c) : XConstraint(idd, c) {}


        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original) override;
    };

    /***************************************************************************
     * constraint count
     **************************************************************************/

    class XConstraintCount : public XConstraint, public XInitialCondition, public XValues {
    public :
        XConstraintCount(std::string idd, std::string c) : XConstraint(idd, c) {}


        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original) override;
    };


    /***************************************************************************
     ****************************************************************************
     *                  CONNECTION  CONSTRAINTS
     ****************************************************************************
     ***************************************************************************/

    /***************************************************************************
     * constraint minimum and maximum
     * startIndex and rank are useful iff index!=NULL
     * in such a case condition is optionnal (=="")
     **************************************************************************/

    class XConstraintMaximum : public XConstraint, public XInitialCondition, public XIndex {
    public :
        int startIndex;
        RankType rank;


        XConstraintMaximum(std::string idd, std::string c) : XConstraint(idd, c), startIndex(0), rank(ANY) {}


        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original) override;
    };

    typedef XConstraintMaximum XConstraintMinimum;

    /***************************************************************************
     * constraint element
     **************************************************************************/

    class XConstraintElement : public XConstraint, public XIndex, public XValue, public XInitialCondition {
    public :
        int startIndex;
        RankType rank;


        XConstraintElement(std::string idd, std::string c) : XConstraint(idd, c), startIndex(0), rank(ANY) {}


        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original) override;
    };

    class XConstraintElementMatrix : public XConstraintElement {
    public :
        vector<vector<XVariable *> > matrix;
        XVariable *index2;
        int startRowIndex, startColIndex;

        XConstraintElementMatrix(std::string idd, std::string c) : XConstraintElement(idd, c) {}

        XConstraintElementMatrix(std::string idd, std::string c, vector<vector<XVariable *> > &mat) : XConstraintElement(idd, c) {
            matrix.resize(mat.size());
            for(unsigned int i = 0 ; i < mat.size() ; i++)
                matrix[i].assign(mat[i].begin(), mat[i].end());
        }


        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original) override;
    };

    /***************************************************************************
     * constraint channel
     **************************************************************************/

    class XConstraintChannel : public XConstraint, public XValue {
    public :
        vector<XVariable *> secondList;
        int startIndex1;
        int startIndex2;


        XConstraintChannel(std::string idd, std::string c) : XConstraint(idd, c), startIndex1(0), startIndex2(0) {}


        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original) override;
    };


    /***************************************************************************
     ****************************************************************************
     *                  PACKING AND SCHEDULING CONSTRAINTS
     ****************************************************************************
     ***************************************************************************/


    /***************************************************************************
     * constraint noOverlap
     **************************************************************************/

    class XConstraintNoOverlap : public XConstraint, public XLengths {
    public :
        // Be careful origins is the vector list!!!
        vector<XVariable *> &origins;
        bool zeroIgnored;


        XConstraintNoOverlap(std::string idd, std::string c) : XConstraint(idd, c), origins(list) {}


        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original) override;
    };

    /***************************************************************************
     * constraint Cumulative
     **************************************************************************/
    class XConstraintCumulative : public XConstraint, public XLengths, public XInitialCondition {
    public :
        // Be careful origin is the vector list!!!
        vector<XVariable *> &origins;
        vector<XVariable *> ends;
        vector<XVariable *> heights;


        XConstraintCumulative(std::string idd, std::string c) : XConstraint(idd, c), origins(list) {}


        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original) override;
    };

    /***************************************************************************
     * constraint BinPacking
     **************************************************************************/
    class XConstraintBinPacking : public XConstraint, public XValues, public XInitialCondition {
    public :
        vector<XVariable*> limits;
        vector<XVariable*> loads;
        std::string conditions;
        int startIndex;

        XConstraintBinPacking(std::string idd, std::string c) : XConstraint(idd, c) {}
        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original) override;
    };
    /***************************************************************************
     * constraint stretch
     **************************************************************************/

    class XConstraintStretch : public XConstraint {
    public :
        vector<int> values; // Only integers
        vector<XInterval> widths; // interval
        vector<vector<int> > patterns;

        XConstraintStretch(std::string idd, std::string c) : XConstraint(idd, c) {}

        // Group is valied with lists only.
        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original) override;
    };

    /***************************************************************************
     * constraint instantiation
     **************************************************************************/


    class XConstraintInstantiation : public XConstraint, public XValues {
    public :
        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original) override;
        XConstraintInstantiation(std::string idd, std::string c) : XConstraint(idd, c) {}
    };


    /***************************************************************************
     * constraint clause
     **************************************************************************/


    class XConstraintClause : public XConstraint {
    public :
        vector<XVariable *> positive;
        vector<XVariable *> negative;


        XConstraintClause(std::string idd, std::string c) : XConstraint(idd, c) { }
        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original) override;
    };

    /***************************************************************************
     ****************************************************************************
    *                  COMPARISON BASED CONSTRAINTS
    ****************************************************************************
    ***************************************************************************/

    class XConstraintCircuit : public XConstraint, public XValue { // value => size
    public :
        XConstraintCircuit(std::string idd, std::string c) : XConstraint(idd, c) {}
        int startIndex;

        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original) override;

    };


    class XConstraintPrecedence :  public XConstraint, public XValues {

        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original) override;

    public:
        bool covered;
        XConstraintPrecedence(std::string idd, std::string c) : XConstraint(idd, c), covered(false) {}
    };


    class XConstraintFlow :  public XConstraint, public XInitialCondition {
        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original) override;
    public:
        vector<vector<int> > arcs;
        vector<XVariable *>balance;
        vector<XVariable *> weights;

        XConstraintFlow(std::string idd, std::string c) : XConstraint(idd, c) {}
    };


    class XConstraintKnapsack :  public XConstraint, public XInitialCondition {
        void unfoldParameters(XConstraintGroup *group, vector<XVariable *> &arguments, XConstraint *original) override;
    public:
        vector<XVariable *>profits;
        vector<XVariable *> weights;
        XInitialCondition profitCondition;

        XConstraintKnapsack(std::string idd, std::string c) : XConstraint(idd, c) {}
    };



}


#endif    /* XCONSTRAINT_H */

