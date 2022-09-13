/*
 costFunction.h

 Author: 
  Sebastien Verel, 
  Univ. du Littoral Côte d'Opale, France.
 
  David Simoncini
  Univ. Toulouse 1 Capitole, France.

*/

#ifndef _costFunction_h
#define _costFunction_h

#include <string>
#include <iostream>
#include <streambuf>

#include <base/solution.h>
//#include "rapidjson/document.h"
//#include "rapidjson/filereadstream.h"
#include <cstdio>

using namespace std;

/*
    Basic Energy Function
*/
class CostFunction {
public:
//    CostFunction(const char * _instance_fileName) {
//        n_variables = 0;
//
//        rapidjson::Document document;
//
//        FILE* fp = fopen(_instance_fileName, "rb"); // non-Windows use "r"
//        if (!fp)
//            std::cerr << "Impossible to open " << _instance_fileName << std::endl;
//
//        // read first line which is an optional comment
//        int i, ch;
//        if (fgetc(fp) == '#') {
//        	while ((ch = fgetc(fp)) != '\n' && ch != EOF) { i++; }
//        } else {
//        	rewind(fp);
//        }
//
//        char readBuffer[65536];
//        rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));
//
//        document.ParseStream(is);
//        fclose(fp);
//
//        if (document.IsObject()) {
//            // variables
//            if (document.HasMember("variables") && document["variables"].IsObject()) {
//                for (rapidjson::Value::ConstMemberIterator itr = document["variables"].MemberBegin(); itr != document["variables"].MemberEnd(); ++itr) {
//                    std::string v = itr->name.GetString();
//                    n_variables++;
//                    // read possibles values
//                    const rapidjson::Value& a = document["variables"][v.c_str()];
//
//                    // number of possible values
//                    n_values.push_back( a.Size() );
//
//                    // save the name of the values
//                    std::vector<std::string> names;
//                    for (rapidjson::SizeType k = 0; k < a.Size(); k++) {
//                        names.push_back( a[k].GetString() );
//                    }
//                    value_name.push_back(names);
//                }
//            }
//
//            // energy functions
//            if (document.HasMember("functions") && document["functions"].IsObject()) {
//                    // linear part
//                    energy.resize(n_variables);
//                    char e[256];
//                    for(unsigned i = 0; i < n_variables; i++) {
//                        sprintf(e, "E%d", i+1);
//                        rapidjson::Value::ConstMemberIterator itr = document["functions"].FindMember(e);
//                        if (itr == document["functions"].MemberEnd()) {
//                        	sprintf(e, "F_%d", i);
//                        	itr = document["functions"].FindMember(e);
//                        }
//                        if (itr != document["functions"].MemberEnd()) {
//                    		const rapidjson::Value& ucosts = document["functions"][e]["costs"];
//                        	// two possibilities
//                        	rapidjson::Value::ConstMemberIterator itr2 = itr->value.FindMember("defaultcost");
//                        	if (itr2 == itr->value.MemberEnd()) {
//                        		// without defaultcost:                               // without defaultcost:
//                        		for (rapidjson::SizeType k = 0; k < ucosts.Size(); k++) {
//                        			energy[i].push_back(ucosts[k].GetDouble());
//                        		}
//                        	} else {
//                        		// with defaultcost
//                        		double defaultcost = itr2->value.GetDouble();
//                        		for (unsigned k = 0; k < n_values[i]; k++) {
//                        			energy[i].push_back(defaultcost);
//                        		}
//                        		for (rapidjson::SizeType k = 0; k < ucosts.Size(); k+=2) {
//                        			energy[i][ucosts[k].GetInt()] = ucosts[k+1].GetDouble();
//                        		}
//
//                        	}
//                        }
//
//                    }
//
//                    // quadratic part
//                    energy2.resize(n_variables);
//                    links.resize(n_variables);
//                    backlinks.resize(n_variables);
//                    for(unsigned i = 0; i < n_variables; i++) {
//                        energy2[i].resize(n_variables);
//                        for(unsigned j = i + 1; j < n_variables; j++) {
//                            sprintf(e, "E%d_%d", i+1, j+1);
//
//                            rapidjson::Value::ConstMemberIterator itr = document["functions"].FindMember(e);
//                            if (itr == document["functions"].MemberEnd()) {
//                            	sprintf(e, "F_%d_%d", i, j);
//                            	itr = document["functions"].FindMember(e);
//                            }
//                            rapidjson::Value::ConstMemberIterator itr2;
//                            if (itr != document["functions"].MemberEnd()) {
//                                links[i].push_back(j);
//                                backlinks[j].push_back(i);
//                                // two possibilities
//                                itr2 = itr->value.FindMember("defaultcost");
//
//                                if (itr2 == itr->value.MemberEnd()) {
//                                    // without defaultcost:
//                                    const rapidjson::Value& costs = itr->value["costs"];
//
//                                    // quadratic term : energy2[i][j][a1][a2] : position i and j, and then acide/angle a1 in position i, and a2 in position j
//                                    energy2[i][j].resize(n_values[i]);
//
//                                    for (rapidjson::SizeType k = 0; k < costs.Size(); k++) {
//                                        if (k / n_values[j] < n_values[i])
//                                            energy2[i][j][k / n_values[j]].push_back(costs[k].GetDouble());
//                                    }
//                                } else {
//                                    // with defaultcost
//                                    double defaultcost = itr2->value.GetDouble();
//                                    const rapidjson::Value& costs = itr->value["costs"];
//
//                                    // quadratic term : energy2[i][j][a1][a2] : position i and j, and then acide/angle a1 in position i, and a2 in position j
//                                    energy2[i][j].resize(n_values[i]);
//                                    for(unsigned a1 = 0; a1 < n_values[i]; a1++)
//                                        energy2[i][j][a1].resize(n_values[j], defaultcost);
//
//                                    rapidjson::SizeType k = 0;
//                                    while (k < costs.Size()) {
//                                        energy2[i][j][ costs[k].GetInt() ][ costs[k + 1].GetInt() ] = costs[k + 2].GetDouble();
//                                        k += 3;
//                                    }
//                                }
//                                // check size
//                                if (energy2[i][j].size() != n_values[i]) {
//                                    std::cout << "error instance " << i << " " << j << " " << n_values[i] << " " << energy2[i][j].size() << std::endl;
//                                }
//                                for(unsigned a = 0; a < n_values[i]; a++)
//                                    if (energy2[i][j][a].size() != n_values[j]) {
//                                        std::cout << "error instance " << i << " " << n_values[i] << " " << a << " " << energy2[i][j][a].size() << std::endl;
//                                    }
//                                }
//                            }
//                    }
//
//            	    // constant part added to first variable
//            	    rapidjson::Value::ConstMemberIterator itr = document["functions"].FindMember("F");
//            	    if (itr != document["functions"].MemberEnd() && n_variables > 0) {
//            	      const rapidjson::Value& costs = itr->value["costs"];
//            	      for (unsigned k = 0; k < n_values[0]; k++) {
//            	    	energy[0][k] += costs[0].GetDouble();
//            	      }
//            	    }
//
//                }
//        } else
//            std::cerr << "Impossible to parse the json file " << _instance_fileName << std::endl;
//    }

    CostFunction(WeightedCSP* wcsp) {
        n_variables = wcsp->numberOfUnassignedVariables();
        // energy functions
        energy.resize(n_variables);
        Cost gap = wcsp->getUb() - wcsp->getLb();

        // variables
        int nbvar = 0;
        for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
            var2index.push_back(nbvar);
            if (wcsp->unassigned(i)) {
                n_values.push_back( wcsp->getDomainSize(i));
                index2var.push_back( i );

                // save the associated values
                std::vector<Value> domvals;
                vector<pair<Value, Cost>> domcosts = wcsp->getEnumDomainAndCost(i);
                for (size_t k = 0; k < domcosts.size(); k++) {
                    domvals.push_back( domcosts[k].first );
                    // linear part: unary cost
                    energy[nbvar].push_back( min(gap, domcosts[k].second) );
                }
                values.push_back(domvals);
                nbvar++;
            }
        }

        energy0 = wcsp->getLb();

        // quadratic part
        energy2.resize(n_variables);
        links.resize(n_variables);
        backlinks.resize(n_variables);
        for(unsigned int i = 0; i < n_variables; i++) {
            energy2[i].resize(n_variables);
        }

        for (unsigned int k = 0; k < wcsp->numberOfConstraints(); k++) {
            if (((WCSP *)wcsp)->getCtr(k)->connected() &&
                    !((WCSP *)wcsp)->getCtr(k)->isSep() &&
                    ((WCSP *)wcsp)->getCtr(k)->isBinary()) {
                BinaryConstraint *ctr = (BinaryConstraint *)((WCSP *)wcsp)->getCtr(k);
                unsigned int i = var2index[ctr->getVar(0)->wcspIndex];
                assert(i < n_variables);
                unsigned int j = var2index[ctr->getVar(1)->wcspIndex];
                assert(j < n_variables);
                if (i < j) {
                    links[i].push_back(j);
                    backlinks[j].push_back(i);
                    // quadratic term
                    energy2[i][j].resize(n_values[i]);
                    vector<Value> domi = wcsp->getEnumDomain(ctr->getVar(0)->wcspIndex);
                    vector<Value> domj = wcsp->getEnumDomain(ctr->getVar(1)->wcspIndex);
                    for(unsigned a1 = 0; a1 < n_values[i]; a1++) {
                        assert(energy2[i][j][a1].size() == 0);
                        for(unsigned a2 = 0; a2 < n_values[j]; a2++) {
                            energy2[i][j][a1].push_back( min(gap, ctr->getCost(domi[a1], domj[a2])) );
                        }
                    }
                } else {
                    links[j].push_back(i);
                    backlinks[i].push_back(j);
                    // quadratic term
                    energy2[j][i].resize(n_values[j]);
                    vector<Value> domi = wcsp->getEnumDomain(ctr->getVar(0)->wcspIndex);
                    vector<Value> domj = wcsp->getEnumDomain(ctr->getVar(1)->wcspIndex);
                    for(unsigned a1 = 0; a1 < n_values[j]; a1++) {
                        assert(energy2[j][i][a1].size() == 0);
                        for(unsigned a2 = 0; a2 < n_values[i]; a2++) {
                            energy2[j][i][a1].push_back( min(gap, ctr->getCost(domi[a2], domj[a1])) );
                        }
                    }
                }
            }
        }

        for (int i = 0; i < ((WCSP *)wcsp)->getElimBinOrder(); i++) {
            BinaryConstraint* ctr = (BinaryConstraint *)((WCSP *)wcsp)->getElimBinCtr(i);
            if (ctr->connected() && !ctr->isSep()) {
                unsigned int i = var2index[ctr->getVar(0)->wcspIndex];
                assert(i < n_variables);
                unsigned int j = var2index[ctr->getVar(1)->wcspIndex];
                assert(j < n_variables);
                if (i < j) {
                    links[i].push_back(j);
                    backlinks[j].push_back(i);
                    // quadratic term
                    energy2[i][j].resize(n_values[i]);
                    vector<Value> domi = wcsp->getEnumDomain(ctr->getVar(0)->wcspIndex);
                    vector<Value> domj = wcsp->getEnumDomain(ctr->getVar(1)->wcspIndex);
                    for(unsigned a1 = 0; a1 < n_values[i]; a1++) {
                        assert(energy2[i][j][a1].size() == 0);
                        for(unsigned a2 = 0; a2 < n_values[j]; a2++) {
                            energy2[i][j][a1].push_back( min(gap, ctr->getCost(domi[a1], domj[a2])) );
                        }
                    }
                } else {
                    links[j].push_back(i);
                    backlinks[i].push_back(j);
                    // quadratic term
                    energy2[j][i].resize(n_values[j]);
                    vector<Value> domi = wcsp->getEnumDomain(ctr->getVar(0)->wcspIndex);
                    vector<Value> domj = wcsp->getEnumDomain(ctr->getVar(1)->wcspIndex);
                    for(unsigned a1 = 0; a1 < n_values[j]; a1++) {
                        assert(energy2[j][i][a1].size() == 0);
                        for(unsigned a2 = 0; a2 < n_values[i]; a2++) {
                            energy2[j][i][a1].push_back( min(gap, ctr->getCost(domi[a2], domj[a1])) );
                        }
                    }
                }
            }
        }

    }

    unsigned int size() const {
        return n_variables;
    };

    unsigned int variable_size(unsigned k) const {
        return n_values[k];
    };

    int getWCSPIndex(unsigned i) const {
        assert(i < n_variables);
        return index2var[i];
    }

    Value getValue(unsigned i, unsigned a) {
        assert(i < n_variables);
        assert(a < n_values[i]);
        return values[i][a];
    }

    int getPILSValueIndex(unsigned i, Value a) {
        assert(i < n_variables);
        unsigned pos = 0;
        for (auto iter = values[i].begin(); iter != values[i].end(); iter++) {
            if (*iter == a) {
                break;
            }
            pos++;
        }
        return (pos < n_values[i])?pos:-1;
    }

    Cost getLb() {
        return energy0;
    }

    void operator()(Solution & _x) {
        Cost fit = 0;
        unsigned j;


        for(size_t i = 0; i < _x.size(); i++) {
            fit += energy[i][ _x[i] ];

            for(size_t k = 0; k < links[i].size(); k++) {
                j = links[i][k];

                if (_x[i] < (int)energy2[i][j].size() && _x[j] < (int)energy2[i][j][_x[i]].size())
                    fit += energy2[i][j][ _x[i] ][ _x[j] ];
                else {
                    if (_x[i] >= (int)energy2[i][j].size())
                        std::cerr << "error i=" << i << " " << _x[i] << " " << energy2[i][j].size() << std::endl;
                    else
                        std::cerr << "error j=" << j << " " << _x[j] << " " << energy2[i][j][_x[i]].size() << std::endl;
                }
            }
        }

        _x.fitness(fit);
    }

    void print() {
        std::cout << n_variables << std::endl;

        for(unsigned i = 0; i < n_variables; i++) {
            std::cout << i << " " << n_values[i] ;
            for(unsigned j = 0; j < n_values[i]; j++) {
                std::cout << " " << values[i][j];
            }
            std::cout << std::endl;
        }

        for(unsigned i = 0; i < n_variables; i++) {
            std::cout << i ;
            for(size_t k = 0; k < links[i].size(); k++) {
                std::cout << " " << links[i][k];
            }
            std::cout << std::endl;
        }

        for(unsigned i = 0; i < n_variables; i++) {
            std::cout << i << " " << n_values[i] ;
            for(unsigned j = 0; j < n_values[i]; j++) {
                std::cout << " " << energy[i][j];
            }
            std::cout << std::endl;
        }

        for(unsigned i = 0; i < n_variables; i++) {
            for(unsigned j = i+1; j < n_variables; j++) {
                if (energy2[i][j].size() > 0) {
                    for(size_t a1 = 0; a1 < energy2[i][j].size(); a1++)
                        for(size_t a2 = 0; a2 < energy2[i][j][a1].size(); a2++)
                            std::cout << i << " " << j << " " << a1 << " " << a2 << " " << energy2[i][j][a1][a2] << std::endl;
                }
            }
        }
    }

  	void printOnShort(std::ostream& _os) const {
        _os << n_variables << std::endl;

        _os << n_values[0] ;
        for(unsigned i = 1; i < n_variables; i++) {
            _os << " " << n_values[i] ;
        }
        _os << std::endl;

        // l_i
        for(unsigned i = 0; i < n_variables; i++) {
        	_os << i ;
            for(size_t k = 0; k < links[i].size(); k++) {
            	if (i < links[i][k])
	                _os << " " << links[i][k] ;
            }
            _os << std::endl ;
        }
    }

    void printlinks() {
        for(size_t i; i<links.size() ; i ++)
        {
            cout << " links " << i << ": ";
                        for (unsigned j : links[i]){
                cout << j << " ";
            }
            cout << endl;
        }
        cout << endl;

    }


    void printEnergy2() {
        for(unsigned i = 0; i< n_variables; i ++)
        {
            cout << i << ": " << endl;
                for (unsigned j = i+1; j< n_variables; j++){
                cout << "      " << j << ": " << endl ;

                    for(size_t k = 0; k< energy2[i][j].size(); k++)   {
                            for (size_t l = 0; l< energy2[i][j][k].size(); l++){
                               cout << "           " << k << "    " << l << "   " << energy2[i][j][k][l] << endl ;
                            }
                    }
            cout << endl;
                }
        cout << endl;

        }
    }

//protected:

    unsigned int n_variables;

    // number of values for each variable
    std::vector<unsigned> n_values;

    // value names
    std::vector< std::vector< Value > > values;

    // list of PILS variable's indexes
    vector<int> var2index;
    // list of WCSP wcspIndex's indexes
    vector<int> index2var;

    // constant term
    Cost energy0;

    // linear term
    std::vector< std::vector<Cost> > energy;

    // interaction between variables with i < j
    std::vector< std::vector<unsigned> > links; // who am I connecting to ?
    std::vector< std::vector<unsigned> > backlinks; // who is connecting to me ?

    // quadratic term : energy2[i][j][a1][a2] : position i and j, and then acide/angle a1 in position i, and a2 in position j
    std::vector< std::vector< std::vector< std::vector<Cost> > > > energy2;
};

#endif
