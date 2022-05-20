/*=============================================================================
 * parser for CSP instances represented in XCSP3 Format
 *
 * Copyright (c) 2015 xcsp.org (contact <at> xcsp.org)
 * Copyright (c) 2008 Olivier ROUSSEL (olivier.roussel <at> cril.univ-artois.fr)
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
#ifndef COSOCO_XMANAGER_H
#define COSOCO_XMANAGER_H

#include <XCSP3CoreCallbacks.h>
#include "XCSP3Constants.h"
#include "XCSP3Variable.h"
#include "XCSP3Constraint.h"
#include "XCSP3Objective.h"
#include <string>
#include <regex>
#include <map>


namespace XCSP3Core {


    class PrimitivePattern;
    class XCSP3Manager {

    public :
        XCSP3CoreCallbacks *callback;
        std::map<std::string, XEntity *> &mapping;
        std::string blockClasses;


        inline bool discardedClasses(string classes) {
            return callback->discardedClasses(blockClasses) || callback->discardedClasses(classes);
        }


    private :
        std::vector<XCSP3Core::PrimitivePattern*> patterns;
        bool recognizePrimitives(std::string id, Tree *tree);
        void createPrimitivePatterns();
        void destroyPrimitivePatterns();

        void containsTrees(vector<XVariable *>&list, vector<Tree *>&newlist);

    public :
        // XCSP3CoreCallbacks *c, std::map<std::string, XEntity *> &m, bool
        XCSP3Manager(XCSP3CoreCallbacks *c, std::map<std::string, XEntity *> &m, bool = true) : callback(c), mapping(m), blockClasses("") { }


        void beginInstance(InstanceType type) {
            callback->_arguments = nullptr;
            callback->beginInstance(type);
        }


        void endInstance() {
            callback->endInstance();
        }


        void beginVariables() {
            callback->beginVariables();
        }


        void endVariables() {
            callback->endVariables();
        }


        void buildVariable(XVariable *variable);


        void beginVariableArray(std::string id) {
            callback->beginVariableArray(id);
        }


        void endVariableArray() {
            callback->endVariableArray();
        }


        void buildVariableArray(XVariableArray *variable);


        void beginConstraints() {
            if(callback->recognizeSpecialIntensionCases)
                createPrimitivePatterns();
            callback->beginConstraints();
        }


        void endConstraints() {
            callback->endConstraints();
            if(callback->recognizeSpecialIntensionCases)
                destroyPrimitivePatterns();
        }


        void beginSlide(std::string id, bool circular) {
            callback->beginSlide(id, circular);
        }


        void endSlide() {
            callback->endSlide();
        }


        void beginObjectives() {
            callback->beginObjectives();
        }


        void endObjectives() {
            callback->endObjectives();
        }


        void beginAnnotations() {
            callback->beginAnnotations();
        }


        void endAnnotations() {
            callback->endAnnotations();
        }
        //--------------------------------------------------------------------------------------
        // Basic constraints
        //--------------------------------------------------------------------------------------

        void newConstraintExtension(XConstraintExtension *constraint);


        void newConstraintExtensionAsLastOne(XConstraintExtension *constraint);


        void newConstraintIntension(XConstraintIntension *constraint);

        //--------------------------------------------------------------------------------------
        // Languages constraints
        //--------------------------------------------------------------------------------------

        void newConstraintRegular(XConstraintRegular *constraint);


        void newConstraintMDD(XConstraintMDD *constraint);

        //--------------------------------------------------------------------------------------
        // Comparison constraints
        //--------------------------------------------------------------------------------------

        void newConstraintAllDiff(XConstraintAllDiff *constraint);


        void newConstraintAllDiffMatrix(XConstraintAllDiffMatrix *constraint);


        void newConstraintAllDiffList(XConstraintAllDiffList *constraint);


        void newConstraintAllEqual(XConstraintAllEqual *constraint);


        void newConstraintOrdered(XConstraintOrdered *constraint);


        void newConstraintLex(XConstraintLex *constraint);


        void newConstraintLexMatrix(XConstraintLexMatrix *constraint);


        //--------------------------------------------------------------------------------------
        // Summin and Counting constraints
        //--------------------------------------------------------------------------------------
    protected :
        void normalizeSum(vector < XVariable * > &list, vector<int> & coefs);

    public :
        void newConstraintSum(XConstraintSum *constraint);


        void newConstraintCount(XConstraintCount *constraint);


        void newConstraintNValues(XConstraintNValues *constraint);


        void newConstraintCardinality(XConstraintCardinality *constraint);



        //--------------------------------------------------------------------------------------
        // Connection constraints
        //--------------------------------------------------------------------------------------

        void newConstraintMinimum(XConstraintMinimum *constraint);


        void newConstraintMaximum(XConstraintMaximum *constraint);

        void newConstraintMinMaxArg(XConstraintMaximum *constraint, bool max = true);
        void newConstraintMinArg(XConstraintMaximum *constraint);
        void newConstraintMaxArg(XConstraintMaximum *constraint);

        void newConstraintElement(XConstraintElement *constraint);


        void newConstraintElementMatrix(XConstraintElementMatrix *constraint);


        void newConstraintChannel(XConstraintChannel *constraint);
        //--------------------------------------------------------------------------------------
        // packing and scheduling constraints
        //--------------------------------------------------------------------------------------

        void newConstraintStretch(XConstraintStretch *constraint);


        void newConstraintNoOverlap(XConstraintNoOverlap *constraint);

        void newConstraintNoOverlapKDim(XConstraintNoOverlap *constraint);



        void newConstraintCumulative(XConstraintCumulative *constraint);


        void newConstraintBinPacking(XConstraintBinPacking *constraint);
        //--------------------------------------------------------------------------------------
        // Instantiation  constraint
        //--------------------------------------------------------------------------------------
        void newConstraintInstantiation(XConstraintInstantiation *constraint);


        //--------------------------------------------------------------------------------------
        // Clause  constraint
        //--------------------------------------------------------------------------------------
        void newConstraintClause(XConstraintClause *constraint);


        //--------------------------------------------------------------------------------------
        // Graph  constraints
        //--------------------------------------------------------------------------------------

        void newConstraintCircuit(XConstraintCircuit *constraint);

        //--------------------------------------------------------------------------------------
        // Precedence  constraints
        //--------------------------------------------------------------------------------------

        void newConstraintPrecedence(XConstraintPrecedence *constraint);

        //--------------------------------------------------------------------------------------
        // Flow constraint
        //--------------------------------------------------------------------------------------

        void newConstraintFlow(XConstraintFlow *constraint);

        //--------------------------------------------------------------------------------------
        // Knapsack constraint
        //--------------------------------------------------------------------------------------

        void newConstraintKnapsack(XConstraintKnapsack *constraint);


        //--------------------------------------------------------------------------------------
        // block of  constraints
        //--------------------------------------------------------------------------------------
        void beginBlock(string classes) {
            blockClasses = classes;
            callback->beginBlock(classes);
        }


        void endBlock() {
            blockClasses = "";
            callback->endBlock();
        }

        //--------------------------------------------------------------------------------------
        // group constraints
        //--------------------------------------------------------------------------------------

        void beginGroup(std::string id) {
            callback->beginGroup(id);
        }


        void endGroup() {
            callback->endGroup();
        }

        template<class T> void unfoldConstraint(XConstraintGroup *group, int i, void (XCSP3Manager::*newConstraint)(T* ));


        void newConstraintGroup(XConstraintGroup *group);



        //--------------------------------------------------------------------------------------
        // Objective constraints
        //--------------------------------------------------------------------------------------

        void addObjective(XObjective *objective);

        //--------------------------------------------------------------------------------------
        // Annotation : Decision variables
        //--------------------------------------------------------------------------------------

        void buildAnnotationDecision(vector<XVariable*> &list) {
            callback->buildAnnotationDecision(list);
        }
    };

}

#endif //COSOCO_XMANAGER_H

























