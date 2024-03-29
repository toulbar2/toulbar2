/*=============================================================================
 * parser for CSP instances represented in XML format
 *
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
#ifndef _AST_hh_
#define _AST_hh_

#include <iostream>
#include <iomanip>
#include <sstream>
#include <exception>
#include <stdexcept>
#include <string>
#include <map>
#include <deque>
#include <cctype>
#include <cstdlib>
#include <cstring>

#include "XMLParser_constants.h"
#include "C_AST.h"

/**
 * @file AST.hh
 * @brief Defines an Abstract Syntax Tree
 */

namespace CSPXMLParser {
using namespace std;

/**
 * TODO
 *
 * - handle operator priorities (infix only)
 * - distinguish booleans from integers in the valuation of an AST
 */

/**
 * return value of a function represented in an abstract syntax tree
 * (AST).
 *
 * one day we'll switch to a union with an int/bool/...
 */
typedef int FunctionValue;

/// possible types for a variable or an expression
enum VarType { TYPE_INT,
    TYPE_BOOL };

/// user provided information about a variable
struct VarInfo {
    int id; // identifier corresponding to the variable name
    VarType type; // type of the variable
};

typedef map<string, VarInfo> VariableInfo; // list of variables (when used)

/**
 * an interface to get the value of a variable
 */
class VariableValuation {
public:
    virtual ~VariableValuation() {}

    /**
     * return the value of a variable in the current interpretation
     *
     * @param id: the identifier of the variable. Whether id is a global
     * or local identifier is dependent on what has been recorded in the
     * AST.
     */
    virtual FunctionValue getVarValue(unsigned int id) const = 0;

    // ??? could we get rid of the virtual
};

class ASTList;
class ASTDict;

/**
 * @brief An abstract node of an abstract syntax tree (AST).
 */
class AST {
public:
    virtual ~AST() {}

    /**
     * compute the value of this expression
     *
     * @param varValuation: provides the value of each variable. This
     * parameter is needed so that the code be usable from several
     * threads or contexts.
     */
    virtual FunctionValue value(const VariableValuation& varValuation) const = 0;

    void expression(ostream& s, Syntax syntax = INFIX_C)
    {
        switch (syntax) {
        case TREE:
            throw runtime_error("can't use this method to get an AST");
            break;
        case PREFIX:
            prefixExpression(s);
            break;
        case INFIX_C:
            infixExpression(s);
            break;
        case POSTFIX:
            postfixExpression(s);
            break;
        case MATHML:
            mathmlExpression(s);
            break;
        default:
            throw runtime_error("undefined syntax for an expressioin");
        }
    }

    /**
     * return the number of children of this node
     */
    virtual int nbArg() const
    {
        return 0;
    }

    int size() const { return nbArg(); }

    /**
     * return the i-th child of this node
     */
    virtual const AST& getArg(int i) const
    {
        throw runtime_error("request for an non existing child of an AST node");
    }

    /**
     * return the type of this node
     */
    virtual NodeType getType() const = 0;

    virtual void prefixExpression(ostream& s) const = 0;
    virtual void infixExpression(ostream& s) const = 0;
    virtual void postfixExpression(ostream& s) const = 0;

    void mathmlExpression(ostream& s) const
    {
        s << "<math>" << endl;
        mathmlExpressionRec(s);
        s << "</math>" << endl;
    }

    virtual void mathmlExpressionRec(ostream& s) const = 0;

    virtual void setArg(int pos, AST* subExpr) = 0;

    /**
     * return a tree of C structures (see C_AST.h) representing the AST 
     *
     * every structure is allocated with malloc. The tree is owned by
     * the caller and memory must be freed by the caller.
     */
    virtual C_AST* makeCTree() const = 0;

    bool isVar() const { return getType() == VAR; }
    bool isInteger() const { return getType() == CST_INT; }
    bool isBoolean() const { return getType() == CST_BOOL; }
    bool isSymbol() const;
    bool isOperator() const;
    bool isList() const { return getType() == LIST; }
    bool isDict() const { return getType() == DICT; }

    bool hasKey(const string& key) const;

    const string& getVarName() const;
    int getVarId() const;
    int getInteger() const;
    bool getBoolean() const;

    /**
     * will throw bad_cast if this is not a dictionary
     */
    const ASTList& list() const;

    /**
     * will throw bad_cast if this is not a dictionary
     */
    const ASTDict& dict() const;

    const AST& operator[](int i) const;
    const AST& operator[](const string& key) const;
};

/**
   * @brief a symbolic constant represented in the AST.
   */
class ASTSymb : public AST {
protected:
    const NodeType cst;

public:
    ASTSymb(NodeType c)
        : cst(c)
    {
    }

    virtual FunctionValue value(const VariableValuation& varValuation) const
    {
        throw runtime_error("no valuation defined for a symbolic constant");
    }

    virtual NodeType getType() const { return cst; }

    /**
     * return the i-th child of this node
     */
    virtual const AST& getArg(int i) const
    {
        throw runtime_error("abstract method called");
    }

    virtual void setArg(int pos, AST* subExpr)
    {
        throw runtime_error("abstract method called");
    }

    virtual void prefixExpression(ostream& s) const
    {
        infixExpression(s);
    }

    virtual void infixExpression(ostream& s) const
    {
        switch (cst) {
        case SYMB_NIL:
            s << "nil";
            break;
        case SYMB_EQ:
            s << "==";
            break;
        case SYMB_NE:
            s << "!=";
            break;
        case SYMB_GE:
            s << ">=";
            break;
        case SYMB_GT:
            s << ">";
            break;
        case SYMB_LE:
            s << "<=";
            break;
        case SYMB_LT:
            s << "<";
            break;
        default:
            throw runtime_error("unimplemented");
        }
    }

    virtual void postfixExpression(ostream& s) const
    {
        infixExpression(s);
    }

    virtual void mathmlExpressionRec(ostream& s) const
    {
        throw runtime_error("unimplemented");
    }

    virtual C_AST* makeCTree() const
    {
        C_AST* node = reinterpret_cast<C_AST*>(malloc(sizeof(C_AST)));

        node->type = cst;

        return node;
    }
};

/**
 * @brief A variable represented in the AST.
 */
class ASTVar : public AST {
protected:
    string name;
    int id;

public:
    ASTVar(const string& name, int id = -1)
    {
        this->name = name;
        this->id = id;
    }

    const string& getName() const { return name; }
    int getId() const { return id; }

    virtual FunctionValue value(const VariableValuation& varValuation) const
    {
        return varValuation.getVarValue(id);
    }

    virtual NodeType getType() const { return VAR; }

    virtual void setArg(int pos, AST* subExpr)
    {
        throw runtime_error("a variable cannot have any argument");
    }

    virtual void prefixExpression(ostream& s) const
    {
        s << name
          << "/$" << id;
    }

    virtual void infixExpression(ostream& s) const
    {
        s // << name << "/"
            //<< "$" << id
            << "parm[" << id << "]";
    }

    virtual void postfixExpression(ostream& s) const
    {
        s << name
          << "/$" << id;
    }

    virtual void mathmlExpressionRec(ostream& s) const
    {
        s << "  <ci>" << name << "</ci>" << endl;
    }

    virtual C_AST* makeCTree() const
    {
        C_AST_VarNode* node = reinterpret_cast<C_AST_VarNode*>(malloc(sizeof(C_AST_VarNode)));

        node->type = VAR;
        node->varName = strdup(name.c_str());
        node->idVar = id;

        return reinterpret_cast<C_AST*>(node);
    }
};

/**
 * @brief An integer constant represented in the AST.
 */
class ASTInteger : public AST {
protected:
    int val;

public:
    ASTInteger(int n)
    {
        val = n;
    }

    ASTInteger(const string& expr)
    {
        istringstream f(expr);
        f >> val;
    }

    int getVal() const { return val; }

    virtual FunctionValue value(const VariableValuation& varValuation) const
    {
        return val;
    }

    virtual NodeType getType() const { return CST_INT; }

    virtual void setArg(int pos, AST* subExpr)
    {
        throw runtime_error("a constant cannot have any argument");
    }

    virtual void prefixExpression(ostream& s) const
    {
        s << val;
    }

    virtual void infixExpression(ostream& s) const
    {
        s << val;
    }

    virtual void postfixExpression(ostream& s) const
    {
        s << val;
    }

    virtual void mathmlExpressionRec(ostream& s) const
    {
        s << "  <cn>" << val << "</cn>" << endl;
    }

    virtual C_AST* makeCTree() const
    {
        C_AST_CstNode* node = reinterpret_cast<C_AST_CstNode*>(malloc(sizeof(C_AST_CstNode)));

        node->type = CST_INT;
        node->val = val;

        return reinterpret_cast<C_AST*>(node);
    }
};

/**
 * @brief A boolean constant represented in the AST.
 */
class ASTBoolean : public AST {
protected:
    bool val;

public:
    ASTBoolean(bool v)
    {
        val = v;
    }

    bool getVal() const { return val; }

    virtual FunctionValue value(const VariableValuation& varValuation) const
    {
        if (val)
            return 1;
        else
            return 0;
    }

    virtual NodeType getType() const { return CST_BOOL; }

    virtual void setArg(int pos, AST* subExpr)
    {
        throw runtime_error("a constant cannot have any argument");
    }

    virtual void prefixExpression(ostream& s) const
    {
        s << boolalpha << val << noboolalpha;
    }

    virtual void infixExpression(ostream& s) const
    {
        s << boolalpha << val << noboolalpha;
    }

    virtual void postfixExpression(ostream& s) const
    {
        s << boolalpha << val << noboolalpha;
    }

    virtual void mathmlExpressionRec(ostream& s) const
    {
        if (val)
            s << "  <true/>" << endl;
        else
            s << "  <false/>" << endl;
    }

    virtual C_AST* makeCTree() const
    {
        C_AST_CstNode* node = reinterpret_cast<C_AST_CstNode*>(malloc(sizeof(C_AST_CstNode)));

        node->type = CST_BOOL;
        node->val = val;

        return reinterpret_cast<C_AST*>(node);
    }
};

/**
 * @brief a list of (heterogeneous) elements represented in the AST.
 */
class ASTList : public AST {
protected:
    deque<AST*> list;

public:
    ASTList()
    {
    }

    ~ASTList()
    {
        for (deque<AST*>::iterator it = list.begin(); it != list.end(); ++it)
            delete *it;
    }

    virtual FunctionValue value(const VariableValuation& varValuation) const
    {
        throw runtime_error("no valuation defined for a list");
    }

    virtual NodeType getType() const { return LIST; }

    virtual int nbArg() const
    {
        return list.size();
    }

    const AST& operator[](int i) const
    {
        return *list[i];
    }

    void clear()
    {
        for (deque<AST*>::size_type i = 0; i < list.size(); ++i)
            if (list[i])
                delete list[i];

        list.clear();
    }

    void push_back(AST* subExpr)
    {
        list.push_back(subExpr);
    }

    void push_front(AST* subExpr)
    {
        list.push_front(subExpr);
    }

    /**
     * return the i-th child of this node
     */
    virtual const AST& getArg(int i) const
    {
        return *list[i];
    }

    virtual void setArg(int pos, AST* subExpr)
    {
        throw runtime_error("abstract method called");
    }

    virtual void prefixExpression(ostream& s) const
    {
        s << "[";
        for (deque<AST*>::const_iterator it = list.begin(); it != list.end(); ++it) {
            s << ' ';
            (*it)->prefixExpression(s);
        }
        s << " ]";
    }

    virtual void infixExpression(ostream& s) const
    {
        s << "[";
        for (deque<AST*>::const_iterator it = list.begin(); it != list.end(); ++it) {
            s << ' ';
            (*it)->infixExpression(s);
        }
        s << " ]";
    }

    virtual void postfixExpression(ostream& s) const
    {
        s << "[";
        for (deque<AST*>::const_iterator it = list.begin(); it != list.end(); ++it) {
            s << ' ';
            (*it)->postfixExpression(s);
        }
        s << " ]";
    }

    virtual void mathmlExpressionRec(ostream& s) const
    {
        throw runtime_error("unimplemented");
    }

    virtual C_AST* makeCTree() const
    {
        C_AST_ListNode* node = reinterpret_cast<C_AST_ListNode*>(malloc(sizeof(C_AST_ListNode) + sizeof(void * [list.size()])));

        node->type = LIST;
        node->size = list.size();

        for (unsigned int i = 0; i < list.size(); ++i)
            node->items[i] = list[i]->makeCTree();

        return reinterpret_cast<C_AST*>(node);
    }
};

/**
 * @brief an associative array of (heterogeneous) elements
 * represented in the AST.
 */
class ASTDict : public AST {
protected:
    map<string, AST*> dict;

public:
    ASTDict()
    {
    }

    ~ASTDict()
    {
        for (map<string, AST*>::iterator it = dict.begin(); it != dict.end(); ++it)
            delete (*it).second;
    }

    virtual FunctionValue value(const VariableValuation& varValuation) const
    {
        throw runtime_error("no valuation defined for a dict");
    }

    virtual NodeType getType() const { return DICT; }

    virtual int nbArg() const
    {
        return dict.size();
    }

    const AST& operator[](const string& key) const
    {
        map<string, AST*>::const_iterator it = dict.find(key);

        if (it == dict.end())
            throw runtime_error("couldn't find key in dictionary");
        else
            return *(*it).second;
    }

    bool hasKey(const string& key) const
    {
        return dict.find(key) != dict.end();
    }

    /**
     * return the i-th child of this node
     */
    virtual const AST& getArg(int i) const
    {
        throw runtime_error("dictionaries can't be accessed by index");
    }

    virtual void setArg(int pos, AST* subExpr)
    {
        throw runtime_error("dictionaries can't be accessed by index");
    }

    virtual const AST& getArg(const string& key) const
    {
        map<string, AST*>::const_iterator it = dict.find(key);

        if (it == dict.end())
            throw runtime_error("couldn't find key in dictionary");
        else
            return *(*it).second;
    }

    virtual void setArg(const string& key, AST* subExpr)
    {
        map<string, AST*>::iterator it = dict.find(key);

        if (it != dict.end())
            throw runtime_error("multiple assignments to the same key in dictionary");
        else
            dict[key] = subExpr;
    }

    virtual void prefixExpression(ostream& s) const
    {
        s << "{";
        for (map<string, AST*>::const_iterator it = dict.begin(); it != dict.end(); ++it) {
            s << " /" << (*it).first << ' ';
            (*it).second->prefixExpression(s);
        }
        s << " }";
    }

    virtual void infixExpression(ostream& s) const
    {
        s << "{";
        for (map<string, AST*>::const_iterator it = dict.begin(); it != dict.end(); ++it) {
            s << " /" << (*it).first << ' ';
            (*it).second->infixExpression(s);
        }
        s << " }";
    }

    virtual void postfixExpression(ostream& s) const
    {
        s << "{";
        for (map<string, AST*>::const_iterator it = dict.begin(); it != dict.end(); ++it) {
            s << " /" << (*it).first << ' ';
            (*it).second->postfixExpression(s);
        }
        s << " }";
    }

    virtual void mathmlExpressionRec(ostream& s) const
    {
        throw runtime_error("unimplemented");
    }

    virtual C_AST* makeCTree() const
    {
        C_AST_DictNode* node = reinterpret_cast<C_AST_DictNode*>(malloc(sizeof(C_AST_DictNode) + sizeof(C_AST_DictEntry[dict.size()])));

        node->type = DICT;
        node->size = dict.size();

        int i = 0;
        for (map<string, AST*>::const_iterator it = dict.begin(); it != dict.end(); ++it) {
            node->items[i].key = strdup((*it).first.c_str());
            node->items[i].value = (*it).second->makeCTree();
            ++i;
        }

        return reinterpret_cast<C_AST*>(node);
    }
};

/**
 * @brief A generic representation of a function inside the AST
 */
class ASTAbstractFunction : public AST {
protected:
    int nbarg;
    AST** args;

    const char* prefixSymbol;
    const char* infixSymbol;
    const char* mathMLSymbol;
    NodeType type;

public:
    ASTAbstractFunction(int nbarg,
        const char* prefixSymbol,
        const char* infixSymbol,
        const char* mathMLSymbol,
        NodeType type)
    {
        if (nbarg < 0)
            throw runtime_error("number of arguments cannot be negative");

        this->nbarg = nbarg;
        this->prefixSymbol = prefixSymbol;
        this->infixSymbol = infixSymbol;
        this->mathMLSymbol = mathMLSymbol;
        this->type = type;

        if (nbarg == 0)
            args = NULL;
        else {
            args = new AST*[nbarg];
            for (int i = 0; i < nbarg; ++i)
                args[i] = NULL;
        }
    }

    virtual ~ASTAbstractFunction()
    {
        if (args != NULL) {
            for (int i = 0; i < nbarg; ++i)
                if (args[i])
                    delete args[i];

            delete[] args;
        }
    }

    virtual AST* makeNode() = 0;

    virtual void setArg(int pos, AST* subExpr)
    {
        if (pos < 0 || pos >= nbarg)
            throw runtime_error("incorrect argument number");

        args[pos] = subExpr;
    }

    virtual int nbArg() const
    {
        return nbarg;
    }

    virtual const AST& getArg(int i) const
    {
        return *args[i];
    }

    virtual NodeType getType() const { return type; }

    const char* getPrefixSymbol() { return prefixSymbol; }
    const char* getInfixSymbol() { return infixSymbol; }

    void prefixExpression(ostream& s) const
    {
        s << prefixSymbol << "(";

        if (nbarg != 0 && args[0] != NULL)
            args[0]->prefixExpression(s);

        for (int i = 1; i < nbarg; ++i) {
            s << ",";
            if (args[i] != NULL)
                args[i]->prefixExpression(s);
        }

        s << ")";
    }

    void infixExpression(ostream& s) const
    {
        if (infixSymbol != NULL) {
            if (nbarg == 2) {
                s << "("; // be careful with priorities

                if (args[0] != NULL)
                    args[0]->infixExpression(s);

                s << infixSymbol;

                if (args[1] != NULL)
                    args[1]->infixExpression(s);

                s << ")";
            } else if (nbarg == 1) {
                s << infixSymbol << "(";
                if (args[0] != NULL)
                    args[0]->infixExpression(s);
                s << ")";
            } else
                throw runtime_error("unexpected arity");
        } else {
            s << prefixSymbol << "(";

            if (nbarg != 0 && args[0] != NULL)
                args[0]->infixExpression(s);

            for (int i = 1; i < nbarg; ++i) {
                s << ",";
                if (args[i] != NULL)
                    args[i]->infixExpression(s);
            }

            s << ")";
        }
    }

    void postfixExpression(ostream& s) const
    {
        if (nbarg != 0 && args[0] != NULL)
            args[0]->postfixExpression(s);

        for (int i = 1; i < nbarg; ++i) {
            s << " ";
            if (args[i] != NULL)
                args[i]->postfixExpression(s);
        }

        s << " " << prefixSymbol;
    }

    virtual void mathmlExpressionRec(ostream& s) const
    {
        s << "  <apply>" << endl;
        s << "  <" << mathMLSymbol << "/>" << endl;
        for (int i = 0; i < nbarg; ++i) {
            if (args[i] != NULL)
                args[i]->mathmlExpressionRec(s);
        }
        s << "  </apply>" << endl;
    }

    virtual C_AST* makeCTree() const
    {
        C_AST_FxNode* node = reinterpret_cast<C_AST_FxNode*>(malloc(sizeof(C_AST_FxNode) + sizeof(void * [nbarg])));

        node->type = type;
        node->nbarg = nbarg;
        for (int i = 0; i < nbarg; ++i)
            node->args[i] = args[i]->makeCTree();

        return reinterpret_cast<C_AST*>(node);
    }

protected:
    /**
     * check that all args are defined (not NULL)
     */
    void checkargs() const
    {
        for (int i = 0; i < nbarg; ++i)
            if (args[i] == NULL)
                throw runtime_error("missing argument for a function");
    }
};

class FunctionNeg : public ASTAbstractFunction {
public:
    FunctionNeg()
        : ASTAbstractFunction(1, "neg", NULL, "minus", F_NEG)
    {
    }

    virtual AST* makeNode() { return new FunctionNeg; }

    FunctionValue value(const VariableValuation& varValuation) const
    {
        checkargs();
        return -args[0]->value(varValuation);
    }

    void infixExpression(ostream& s) const
    {
        s << "-";
        s << "("; // ??? avoid these parentheses unless we need them

        if (args[0] != NULL)
            args[0]->infixExpression(s);

        s << ")";
    }
};

class FunctionAbs : public ASTAbstractFunction {
public:
    FunctionAbs()
        : ASTAbstractFunction(1, "abs", NULL, "abs", F_ABS)
    {
    }

    virtual AST* makeNode() { return new FunctionAbs; }

    FunctionValue value(const VariableValuation& varValuation) const
    {
        checkargs();
        return abs(args[0]->value(varValuation));
    }
};

class FunctionAdd : public ASTAbstractFunction {
public:
    FunctionAdd()
        : ASTAbstractFunction(2, "add", "+", "plus", F_ADD)
    {
    }

    AST* makeNode() { return new FunctionAdd; }

    FunctionValue value(const VariableValuation& varValuation) const
    {
        checkargs();
        return args[0]->value(varValuation) + args[1]->value(varValuation);
    }
};

class FunctionSub : public ASTAbstractFunction {
public:
    FunctionSub()
        : ASTAbstractFunction(2, "sub", "-", "minus", F_SUB)
    {
    }

    AST* makeNode() { return new FunctionSub; }

    FunctionValue value(const VariableValuation& varValuation) const
    {
        checkargs();
        return args[0]->value(varValuation) - args[1]->value(varValuation);
    }
};

class FunctionMul : public ASTAbstractFunction {
public:
    FunctionMul()
        : ASTAbstractFunction(2, "mul", "*", "times", F_MUL)
    {
    }

    AST* makeNode() { return new FunctionMul; }

    FunctionValue value(const VariableValuation& varValuation) const
    {
        checkargs();
        return args[0]->value(varValuation) * args[1]->value(varValuation);
    }
};

class FunctionDiv : public ASTAbstractFunction {
public:
    FunctionDiv()
        : ASTAbstractFunction(2, "div", "/", "quotient", F_DIV)
    {
    }

    AST* makeNode() { return new FunctionDiv; }

    FunctionValue value(const VariableValuation& varValuation) const
    {
        checkargs();

        int div = args[1]->value(varValuation);

        if (div == 0)
            throw runtime_error("divide by 0");

        return args[0]->value(varValuation) / div;
    }
};

class FunctionMod : public ASTAbstractFunction {
public:
    FunctionMod()
        : ASTAbstractFunction(2, "mod", "%", "rem", F_MOD)
    {
    }

    AST* makeNode() { return new FunctionMod; }

    FunctionValue value(const VariableValuation& varValuation) const
    {
        checkargs();

        int div = args[1]->value(varValuation);

        if (div == 0)
            throw runtime_error("modulo by 0");

        // ??? implementation defined !!
        return args[0]->value(varValuation) % div;
    }
};

class FunctionPow : public ASTAbstractFunction {
public:
    FunctionPow()
        : ASTAbstractFunction(2, "pow", NULL, "power", F_POW)
    {
    }

    AST* makeNode() { return new FunctionPow; }

    FunctionValue value(const VariableValuation& varValuation) const
    {
        checkargs();

        long long prod = 1;
        long long x = args[0]->value(varValuation);
        int n = args[1]->value(varValuation);
        int test;

        if (x == 0 && n == 0)
            throw runtime_error("pow(0,0) is undefined");

        while (n != 0) {
            if ((test = prod) != prod || (test = x) != x)
                throw runtime_error("overflow in pow(x,y) computation");

            if (n & 0x01)
                prod *= x;

            x *= x;
            n >>= 1;
        }

        if ((test = prod) != prod) // ??? reliable ?
            throw runtime_error("overflow in pow(x,y) computation");

        return prod;
    }
};

class FunctionIf : public ASTAbstractFunction {
public:
    FunctionIf()
        : ASTAbstractFunction(3, "if", NULL, NULL, F_IF)
    {
    }

    AST* makeNode() { return new FunctionIf; }

    FunctionValue value(const VariableValuation& varValuation) const
    {
        checkargs();

        if (args[0]->value(varValuation))
            return args[1]->value(varValuation);
        else
            return args[2]->value(varValuation);
    }
};

class FunctionMin : public ASTAbstractFunction {
public:
    FunctionMin()
        : ASTAbstractFunction(2, "min", NULL, "min", F_MIN)
    {
    }

    AST* makeNode() { return new FunctionMin; }

    FunctionValue value(const VariableValuation& varValuation) const
    {
        checkargs();

        return std::min(args[0]->value(varValuation), args[1]->value(varValuation));
    }
};

class FunctionMax : public ASTAbstractFunction {
public:
    FunctionMax()
        : ASTAbstractFunction(2, "max", NULL, "max", F_MAX)
    {
    }

    AST* makeNode() { return new FunctionMax; }

    FunctionValue value(const VariableValuation& varValuation) const
    {
        checkargs();

        return std::max(args[0]->value(varValuation), args[1]->value(varValuation));
    }
};

class FunctionEQ : public ASTAbstractFunction {
public:
    FunctionEQ()
        : ASTAbstractFunction(2, "eq", "==", "eq", F_EQ)
    {
    }

    AST* makeNode() { return new FunctionEQ; }

    FunctionValue value(const VariableValuation& varValuation) const
    {
        checkargs();

        return args[0]->value(varValuation) == args[1]->value(varValuation);
    }
};

class FunctionNE : public ASTAbstractFunction {
public:
    FunctionNE()
        : ASTAbstractFunction(2, "ne", "!=", "neq", F_NE)
    {
    }

    AST* makeNode() { return new FunctionNE; }

    FunctionValue value(const VariableValuation& varValuation) const
    {
        checkargs();

        return args[0]->value(varValuation) != args[1]->value(varValuation);
    }
};

class FunctionGE : public ASTAbstractFunction {
public:
    FunctionGE()
        : ASTAbstractFunction(2, "ge", ">=", "geq", F_GE)
    {
    }

    AST* makeNode() { return new FunctionGE; }

    FunctionValue value(const VariableValuation& varValuation) const
    {
        checkargs();

        return args[0]->value(varValuation) >= args[1]->value(varValuation);
    }
};

class FunctionGT : public ASTAbstractFunction {
public:
    FunctionGT()
        : ASTAbstractFunction(2, "gt", ">", "gt", F_GT)
    {
    }

    AST* makeNode() { return new FunctionGT; }

    FunctionValue value(const VariableValuation& varValuation) const
    {
        checkargs();
        if (args[0]->value(varValuation) > args[1]->value(varValuation))
            return FunctionValue(1);
        else
            return FunctionValue(0);
    }
};

class FunctionLE : public ASTAbstractFunction {
public:
    FunctionLE()
        : ASTAbstractFunction(2, "le", "<=", "leq", F_LE)
    {
    }

    AST* makeNode() { return new FunctionLE; }

    FunctionValue value(const VariableValuation& varValuation) const
    {
        checkargs();

        return args[0]->value(varValuation) <= args[1]->value(varValuation);
    }
};

class FunctionLT : public ASTAbstractFunction {
public:
    FunctionLT()
        : ASTAbstractFunction(2, "lt", "<", "lt", F_LT)
    {
    }

    AST* makeNode() { return new FunctionLT; }

    FunctionValue value(const VariableValuation& varValuation) const
    {
        checkargs();
        if (args[0]->value(varValuation) < args[1]->value(varValuation))
            return FunctionValue(1);
        else
            return FunctionValue(0);
    }
};

class FunctionNot : public ASTAbstractFunction {
public:
    FunctionNot()
        : ASTAbstractFunction(1, "not", "!", "not", F_NOT)
    {
    }

    AST* makeNode() { return new FunctionNot; }

    FunctionValue value(const VariableValuation& varValuation) const
    {
        checkargs();

        return !args[0]->value(varValuation);
    }
};

class FunctionAnd : public ASTAbstractFunction {
public:
    FunctionAnd()
        : ASTAbstractFunction(2, "and", "&&", "and", F_AND)
    {
    }

    AST* makeNode() { return new FunctionAnd; }

    FunctionValue value(const VariableValuation& varValuation) const
    {
        checkargs();

        return args[0]->value(varValuation) && args[1]->value(varValuation);
    }
};

class FunctionOr : public ASTAbstractFunction {
public:
    FunctionOr()
        : ASTAbstractFunction(2, "or", "||", "or", F_OR)
    {
    }

    AST* makeNode() { return new FunctionOr; }

    FunctionValue value(const VariableValuation& varValuation) const
    {
        checkargs();

        return args[0]->value(varValuation) || args[1]->value(varValuation);
    }
};

class FunctionXor : public ASTAbstractFunction {
public:
    FunctionXor()
        : ASTAbstractFunction(2, "xor", "^", "xor", F_XOR)
    {
    }

    AST* makeNode() { return new FunctionXor; }

    FunctionValue value(const VariableValuation& varValuation) const
    {
        checkargs();

        return args[0]->value(varValuation) ^ args[1]->value(varValuation);
    }
};

class FunctionIff : public ASTAbstractFunction {
public:
    FunctionIff()
        : ASTAbstractFunction(2, "iff", "==", "equivalent", F_IFF)
    {
    }

    AST* makeNode() { return new FunctionIff; }

    FunctionValue value(const VariableValuation& varValuation) const
    {
        checkargs();

        return args[0]->value(varValuation) == args[1]->value(varValuation);
    }
};

/**
 * for internal use only
 *
 * ??? this should be moved somewhere else
 */
class ASTDictKey : public AST {
private:
    const string key;

public:
    ASTDictKey(const string& key)
        : key(key)
    {
    }

    virtual NodeType getType() const { return _DICTKEY; }

    const string& getKey() const { return key; }

    virtual void setArg(int pos, AST* subExpr)
    {
        throw runtime_error("abstract method called");
    }

    virtual FunctionValue value(const VariableValuation& varValuation) const
    {
        throw runtime_error("abstract method called");
    }

    virtual void prefixExpression(ostream& s) const
    {
        throw runtime_error("abstract method called");
    }

    virtual void infixExpression(ostream& s) const
    {
        throw runtime_error("abstract method called");
    }

    virtual void postfixExpression(ostream& s) const
    {
        throw runtime_error("abstract method called");
    }

    virtual void mathmlExpressionRec(ostream& s) const
    {
        throw runtime_error("abstract method called");
    }

    virtual C_AST* makeCTree() const
    {
        throw runtime_error("abstract method called");
    }
};

inline bool AST::isSymbol() const
{
    return dynamic_cast<const ASTSymb*>(this) == this;
}

inline bool AST::isOperator() const
{
    return dynamic_cast<const ASTAbstractFunction*>(this) == this;
}

inline const string& AST::getVarName() const
{
    return dynamic_cast<const ASTVar*>(this)->getName();
}

inline int AST::getVarId() const
{
    return dynamic_cast<const ASTVar*>(this)->getId();
}

inline int AST::getInteger() const
{
    return dynamic_cast<const ASTInteger*>(this)->getVal();
}

inline bool AST::getBoolean() const
{
    return dynamic_cast<const ASTBoolean*>(this)->getVal();
}

/**
 * will throw bad_cast if this is not a dictionary
 */
inline const ASTList& AST::list() const
{
    return dynamic_cast<const ASTList&>(*this);
}

/**
 * will throw bad_cast if this is not a dictionary
 */
inline const ASTDict& AST::dict() const
{
    return dynamic_cast<const ASTDict&>(*this);
}

inline bool AST::hasKey(const string& key) const
{
    return dynamic_cast<const ASTDict*>(this)->hasKey(key);
}

inline const AST& AST::operator[](int i) const
{
    return dynamic_cast<const ASTList&>(*this)[i];
}

inline const AST& AST::operator[](const string& key) const
{
    return dynamic_cast<const ASTDict&>(*this)[key];
}

/**
 * this class is used by the parser whenever it needs to create an
 * AST node.
 */
class DefaultASTFactory {
public:
    typedef AST ASTType;

    inline static ASTVar* mkVar(const string& name, int id = -1)
    {
        return new ASTVar(name, id);
    }

    inline static ASTInteger* mkInteger(int val)
    {
        return new ASTInteger(val);
    }

    inline static ASTInteger* mkInteger(const string& s)
    {
        return new ASTInteger(s);
    }

    inline static ASTBoolean* mkBoolean(bool val)
    {
        return new ASTBoolean(val);
    }

    inline static ASTSymb* mkSymb(NodeType c)
    {
        return new ASTSymb(c);
    }

    inline static ASTList* mkList()
    {
        return new ASTList();
    }

    inline static ASTDict* mkDict()
    {
        return new ASTDict();
    }
};

} // namespace

#endif
