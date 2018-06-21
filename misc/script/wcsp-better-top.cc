// -*- Mode: c++; c-basic-offset: 4 -*-

// This C++ script reads a CFN with enumerated variables and table cost
// functions only and computes a possibly improved starting upper
// bound (the sum of all the maximum finite costs in tables where
// finite means less than the initially provided upper bound).
// Compile with  g++
// Author: George Katsirelos (CPD version Thomas Schiex).

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <tuple>
#include <map>
#include <unordered_map>
#include <boost/format.hpp>

using namespace std;
using boost::format;

const bool debug = false;

template<typename T>
ostream& operator<<(ostream& os, vector<T> const& v)
{
    os << "[";
    for(auto& t : v)
        os << t << ' ';
    os << "]";
    return os;
}

typedef long long Cost;

struct wcsptuple {
    vector<size_t> tup;
    Cost cost;
};

struct wcspfunc {
    Cost defcost;
    vector<size_t> scope;
    vector< wcsptuple > specs;

    size_t arity() const { return scope.size(); }
};

struct wcsp {
    string name;
    Cost ub;

    vector<size_t> domains;
    size_t nvars() const { return domains.size(); }

    vector<wcspfunc> functions;
    vector<string> value_names;
};

template<typename T>
vector<T> read_vec(istream& is)
{
    vector<T> r;
    T s;
    is >> s;
    while(is) {
        r.push_back(s);
        is >> s;
    }
    return r;
}

template<typename T>
vector<T> read_vec(string const& line)
{
    istringstream iss(line);
    return read_vec<T>(iss);
}

tuple<string, size_t, size_t, size_t, Cost> read_header(string const& line)
{
    istringstream iss(line);

    string name;
    size_t nvars;
    size_t domsize;
    size_t nfun;
    Cost ub;

    iss >> name >> nvars >> domsize >> nfun >> ub;
    return make_tuple(name, nvars, domsize, nfun, ub);
}

wcspfunc read_fun(istream& is)
{
    string line;
    getline(is, line);
    vector<Cost> hd = read_vec<Cost>(line);
    size_t arity = hd[0];
    Cost defcost = hd[hd.size()-2];
    size_t nspec = hd[hd.size()-1];

    vector<wcsptuple> specs;
    for(size_t i = 0; i != nspec; ++i) {
        getline(is, line);
        vector<Cost> v = read_vec<Cost>(line);
        specs.push_back( {vector<size_t>(v.begin(), v.begin()+arity), v[v.size()-1]} );
    }

    return { defcost, vector<size_t>(hd.begin()+1, hd.begin()+1+arity),
            specs };
}

string read_value_names(istream& is)
{
    string line;
    getline(is,line);
    return line;
}

wcsp readwcsp(istream& is)
{
    wcsp w;

    size_t nvars;
    size_t domsize;
    size_t nfun;

    string line;

    getline(is, line);
    tie(w.name, nvars, domsize, nfun, w.ub) = read_header(line);

    getline(is, line);
    w.domains = read_vec<size_t>(line);

    for(size_t i = 0; i != nfun; ++i)
        w.functions.push_back(read_fun(is));

    for(size_t i = 0; i != nvars; ++i)
	w.value_names.push_back(read_value_names(is));

    return w;
}

void write_wcsp(wcsp const& w, ostream &ofs)
{
    size_t maxd = *max_element(w.domains.begin(), w.domains.end());
    ofs << w.name << ' ' << w.nvars()
        << ' ' << maxd
        << ' ' << w.functions.size() << ' ' << w.ub << "\n";

    for(auto& d : w.domains)
        ofs << d << ' ';
    ofs << "\n";

    for(auto& f: w.functions) {
        ofs << f.arity() << ' ';
        for(auto& v : f.scope)
            ofs << v << ' ';
        ofs << f.defcost << ' ' << f.specs.size() << "\n";
        for(auto& s : f.specs) {
            for(auto& v : s.tup)
                ofs << v << ' ';
            ofs << min(s.cost, w.ub) << "\n";
        }
    }

    for(auto& vn: w.value_names) {
	ofs << vn << std::endl;
    }

}

int main(int argc, char* argv[])
{
    if( argc != 3 ) {
        cout << "usage: " << argv[0] << " <wcsp-input> <wcnf-output>\n";
        return 1;
    }

    ifstream ifs(argv[1]);
    ofstream ofs(argv[2]);

    if( !ifs ) {
        cout << "could not open " << argv[1] << "\n";
        return 1;
    }

    if( !ofs ) {
        cout << "could not open " << argv[2] << "\n";
        return 1;
    }

    wcsp w = readwcsp(ifs);
    //cout << "initial top " << w.ub << "\n";
    Cost newtop = 0;
    for(auto const& f : w.functions) {
        auto me = std::max_element(f.specs.begin(), f.specs.end(),
                                   [&](wcsptuple const& m, wcsptuple const& c){
                                       return c.cost < w.ub &&
                                       (c.cost > m.cost ||
                                        m.cost >= w.ub);
                                   });
        //cout << "max of function " << f.scope << ": " << me->cost << "\n";
        newtop += me->cost;
    }
    cout << "new top " << newtop << "\n";
    w.ub = newtop;
    write_wcsp(w, ofs);

    return 0;
}
