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
#ifndef _ExpressionParser_hh_
#define _ExpressionParser_hh_

#include <iostream>
#include <iomanip>
#include <sstream>
#include <exception>
#include <stdexcept>
#include <string>
#include <map>
#include <deque>
#include <cctype>

#include "AST.hh"

/**
 * @file ExpressionParser.hh
 * @brief Defines a parser for prefix/infix/postfix expressions.
 */

namespace CSPXMLParser
{
  using namespace std;

  template<class ASTFactory>
  class ExpressionParser
  {
  private:
    deque<ASTAbstractFunction *> functionList;

    const VariableInfo *varInfo; // an optional map which provides some
    // information on the variables that
    // may occur in the expression

  public:
    ExpressionParser()
    {
      varInfo=NULL;

      functionList.push_back(new FunctionNeg);
      functionList.push_back(new FunctionAbs);
      functionList.push_back(new FunctionAdd);
      functionList.push_back(new FunctionSub);
      functionList.push_back(new FunctionMul);
      functionList.push_back(new FunctionDiv);
      functionList.push_back(new FunctionMod);
      functionList.push_back(new FunctionPow);
      functionList.push_back(new FunctionIf);
      functionList.push_back(new FunctionMin);
      functionList.push_back(new FunctionMax);
      functionList.push_back(new FunctionEQ);
      functionList.push_back(new FunctionNE);
      functionList.push_back(new FunctionGE);
      functionList.push_back(new FunctionGT);
      functionList.push_back(new FunctionLE);
      functionList.push_back(new FunctionLT);
      functionList.push_back(new FunctionNot);
      functionList.push_back(new FunctionAnd);
      functionList.push_back(new FunctionOr);
      functionList.push_back(new FunctionXor);
      functionList.push_back(new FunctionIff);
    }

    ~ExpressionParser()
    {
      for(deque<ASTAbstractFunction *>::iterator it=functionList.begin();
	  it!=functionList.end();++it)
	delete *it;
    }

    void setVarInfo(const VariableInfo *varInfo=NULL)
    {
      this->varInfo=varInfo;
    }

    void unsetVarInfo()
    {
      this->varInfo=NULL;
    }

    AST *prefixParser(const string &expr)
    {
      string s;

      // remove spaces
      for(unsigned int i=0;i<expr.length();++i)
	if (!isspace(expr[i]))
	  s+=expr[i];

      return recursivePrefixParser(s);
    }

    AST *infixParser(const string &expr)
    {
      return ASTFactory::mkVar("infix parser unimplemented");
    }

    AST *postfixParser(const string &expr)
    {
      return ASTFactory::mkVar("postfix parser unimplemented");
    }

  private:
    AST *recursivePrefixParser(const string &f)
    {
      int level=0;
      int argNum=0;
      int subExprStart=0;

      AST *node=NULL;

      for(unsigned int i=0;i<f.length();++i)
      {
	if (f[i]=='(')
	{
	  if (level==0)
	  {
	    node=findPrefixFunction(f.substr(0,i));
	    subExprStart=i+1;
	  }
	  ++level;
	}
	else
	{
	  if (level==1 && (f[i]==',' || f[i]==')'))
	  {
	    node->setArg(argNum,prefixParser(f.substr(subExprStart,
						      i-subExprStart)));
	    ++argNum;

	    subExprStart=i+1;
	  }

	  if (f[i]==')')
	    --level;
	}
      }

      if (level!=0)
	throw runtime_error("unbalanced parentheses");

      if (node==NULL)
      {
	// no opening parenthese found, this is a constant or a variable
	if (isalpha(f[0]))
	{
	  if (f=="true")
	    node=ASTFactory::mkBoolean(true);
	  else
	    if (f=="false")
	      node=ASTFactory::mkBoolean(false);
	    else
	    {
	      // a variable
	      if (varInfo==NULL)
		node=ASTFactory::mkVar(f);
	      else
	      {
		VariableInfo::const_iterator it=varInfo->find(f);

		if(it==varInfo->end())
		  throw runtime_error("undefined variable found in expression");

		node=ASTFactory::mkVar(f,(*it).second.id);
	      }
	    }
	}
	else
	{
	  // an int
	  node=ASTFactory::mkInteger(f);
	}
      }

      return node;
    }

    AST *findPrefixFunction(const string &name)
    {
      for(unsigned int i=0;i<functionList.size();++i)
	if (functionList[i]->getPrefixSymbol()==name)
	  return functionList[i]->makeNode();

      throw runtime_error("unknown function symbol");
    }
  };

} // namespace

#endif
