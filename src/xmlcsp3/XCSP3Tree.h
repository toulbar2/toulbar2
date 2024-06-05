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

#ifndef TREE_H
#define    TREE_H

#include<cmath>
#include <vector>
#include <map>
#include<iostream>
#include<assert.h>
#include "XCSP3TreeNode.h"


namespace XCSP3Core {
    class Tree {
    protected:
        std::string expr;

        void createOperator(std::string currentElement, std::vector<NodeOperator *> &stack, std::vector<Node *> &params);
        void closeOperator(std::vector<NodeOperator *> &stack, std::vector<Node *> &params);
        void createBasicParameter(std::string currentElement, std::vector<NodeOperator *> &stack, std::vector<Node *> &params);
    public:
        Node *root;
        std::vector<std::string> listOfVariables;


        Tree(std::string e) : expr(e) {
            root = fromStringToTree(expr);
        }

        Tree(Node *r) : root(r) { }

        Node *fromStringToTree(std::string);

        int arity() {
            return static_cast<int>(listOfVariables.size());
        }


        int evaluate(std::map<std::string, int> &tuple) {
            return root->evaluate(tuple);
        }

        std::string toString() {
            return root->toString();
        }

        void prefixe()  {
            std::cout << root->toString();
        }

        void canonize() {
            root =  root->canonize();
        }
    };
}



#endif	/* TREE_H */

