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
#ifndef XCSP3CONSTANTS_H
#define	XCSP3CONSTANTS_H
#include<climits>

namespace XCSP3Core {
    typedef enum {
        CSP, COP
    } InstanceType;

    typedef enum CONSTRAINTS {
        UNKNOWN, EXTENSION, INTENSION, ALLDIFF, ALLEQUAL, SUM, ORDERED, COUNT, NVALUES, CARDINALITY, MAXIMUM, MINIMUM,
        ELEMENT, ELEMENTMATRIX, NOOVERLAP, STRETCH, LEX, CHANNEL, REGULAR, MDD, CUMULATIVE, INSTANTIATION, CIRCUIT, CLAUSE,
        PRECEDENCE, BINPACKING, FLOW, KNAPSACK, MINARG, MAXARG
    } ConstraintType;

    typedef enum order {
        LE, LT, GE, GT, IN, EQ, NE, NOTIN
    } OrderType;
    typedef enum operantype {
        INTEGER, INTERVAL, VARIABLE, SET
    } OperandType;


    typedef enum tag {
        UNKNOWNTAG, LISTTAG, FUNCTIONALTAG, VALUESTAG, VALUETAG, CONDITIONTAG, INDEXTAG, LENGTHSTAG
    } Tag;

    typedef enum ranktype {
        ANY, FIRST, LAST
    } RankType;

    typedef enum objective {
        MINIMIZE, MAXIMIZE
    } ObjectiveGoal;

    typedef enum expressionObjective {
        EXPRESSION_O, SUM_O, PRODUCT_O, MINIMUM_O, MAXIMUM_O, NVALUES_O, LEX_O
    } ExpressionObjective;

#define STAR INT_MAX
}
#endif	/* XCSP3CONSTANTS_H */

