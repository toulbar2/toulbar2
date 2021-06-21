#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <algorithm>
#include <cassert>
#include <vector>

/* normalize lb/ub local per instance but time global (i.e. 1 means TIMELIMIT) to produce anytime curves */

/* it requires two arguments:  a file with a list of problem names followed by .best extension and a solver name */
/* it assumes relevant information are available in files problem_solver.ub,  problem_solver.lb,  problem.lb0, problem.best, problem.worst */
/* change the following CPU-time limit to your setting */

#define TIMELIMIT 3600

using namespace std;

/* 3 modes of operation wrt preprocessing, set by the following two
   variables

   1. assume the lower bound is known at time 0 (false, false)
   2. assume preprocessing takes 0 time (true, false)

   3. start with lower bound 0, preprocessing reports a better one at
   a later time (false, true). Clamp normalization to [0,1], so this
   acts like option 1 most of the time but can tell when
*/
static const bool ignore_preprocessing{false};
static const bool the_right_thing{true};

using Cost = long long;

#define PV(x) "("#x")=" << (x)

/* basename("/path/to/file.ext", ".ext") returns "file" */
std::string basename(std::string fn, std::string extension)
{
    using std::string;
    size_t dotpos = fn.rfind(extension);
    size_t dirpos = fn.rfind('/');
    if( dirpos == string::npos )
        dirpos = 0;
    else
        ++dirpos;
    return fn.substr(dirpos, dotpos-dirpos);
}

std::string replace_ext(std::string fn, std::string extension, std::string newext)
{
    using std::string;
    size_t dotpos = fn.rfind(extension);
    return fn.substr(0, dotpos) + newext;
}


void read_pairs(ifstream& ifs, vector<pair<double, Cost>>& pairs)
{
    string line;
    do {
        getline(ifs, line);
        if( !ifs )
            break;
        istringstream iss(line);
        double t;
        long double b;
        iss >> t >> b;
        if( !iss )
            continue;
        pairs.push_back( make_pair(t, (Cost) b) );
    } while(true);
}

// re-implementation of std::unique, with control: we want to
// keep the earliest occurrence
template<class ForwardIt, typename Pred>
ForwardIt myunique(ForwardIt first, ForwardIt last, Pred p)
{
    if (first == last)
        return last;

    ForwardIt result = first;
    while (++first != last) {
        if (!p(*result, *first)) {
            *(++result) = std::move(*first);
        }
    }
    return ++result;
}

int main(int argc, char *argv[])
{
    if( argc != 3 ) {
        cout << "give list of .best and algname\n";
        exit(1);
    }

    assert(!ignore_preprocessing || !the_right_thing);

    map<string, Cost> lb0s;
    map<string, Cost> opts;

    map<string, vector<pair<double, double>>> alllbs, allubs;

    map<string, double> lastlb, lastub;

    const char *bestfiles = argv[1];
    string algname{argv[2]};

    ifstream bests(bestfiles);

    string bestline, lbline, ubline;
    do {
        getline(bests, bestline);
        if(!bests) break;

        lbline = replace_ext( bestline, ".best",  "_"+algname+".lb" );
        ubline = replace_ext( bestline, ".best",  "_"+algname+".ub" );

        cout << bestline << "\n";

        Cost opt;
        ifstream bfs(bestline);
        bfs >> opt;
        if( !bfs ) {
            cout << "no solution found by anyone, skipping\n";
            continue;
        } else
            cout << "best solution " << opt << "\n";

        Cost ub0;
        cout << "worst solution in " << replace_ext(bestline, ".best", ".worst") << "\n";
        ifstream wfs( replace_ext(bestline, ".best", ".worst"));
        if( !wfs ) {
            cout << "could not open worst solution file for '" << bestline << "'\n";
            exit(1);
        }
        wfs >> ub0;
        if( !wfs ) {
            cout << "something is wrong: best solution found but not worst\n";
            exit(1);
        }

        cout << "lb0 in " << replace_ext(bestline, ".best", ".lb0") << "\n";
        ifstream lb0fs( replace_ext(bestline, ".best", ".lb0"));
        if( !lb0fs ) {
            cout << "could not open lb0 file for '" << bestline << "'\n";
            exit(1);
        }
        Cost lb0;
        lb0fs >> lb0;
        if( !lb0fs ) {
            cout << "something is wrong: best solution found but not lb0\n";
            exit(1);
        }

        vector<pair<double, Cost> > instancelbs, instanceubs;

        cout << "reading lbs from " << lbline << "\n";
        ifstream ilfs(lbline);
        cout << "reading ubs from " << ubline << "\n";
        ifstream iufs(ubline);

        if( (!ilfs && !!iufs) || (!!ilfs && !iufs ) ) {
            cout << "only one of .lb, .ub exist, skipping\n";
            continue;
        }

        if( !!ilfs ) {
            read_pairs( ilfs, instancelbs );
        } else
            instancelbs.push_back(make_pair(0.0, 0));

        instanceubs.push_back( make_pair(0.0, std::numeric_limits<Cost>::max()));
        if( !!iufs )
            read_pairs( iufs, instanceubs );

        if( instancelbs.empty() ) {
            cout << "no lb0, using fake\n";
            instancelbs.push_back(make_pair(0.0,0));
            //lb0 = 0;
        } else {
            //cout << "lb0 = " << instancelbs[0].second << "\n";
            double pptime = instancelbs[0].first;
            //lb0 = instancelbs[0].second;
            if( ignore_preprocessing ) {
                for(auto & lbe : instancelbs )
                    lbe.first -= pptime;
                for(auto &ube : instanceubs )
                    ube.first -= pptime;
            } else if( the_right_thing) {
                if( instancelbs[0].first > 0 )
                    instancelbs.insert(instancelbs.begin(), make_pair(0.0, Cost(0)));
            } else {
                instancelbs[0].first = 0.0; // assume it's from the beginning
            }
        }

        cout << "lb0 = " << lb0 << "\n";

        cout << "read " << instancelbs.size() << " lbs, "
             << instanceubs.size() << " ubs" << endl;

        auto i = find_if(begin(instancelbs), end(instancelbs), [&](pair<double, Cost> x) {
                return x.second > instanceubs.back().second;
            });
        if( i != end(instancelbs) ) {
            cout << "removing " << (end(instancelbs)-i) << " bad lower bounds\n";
            instancelbs.erase(i, end(instancelbs));
        }

        cout << "lbs: ";
        for( auto l : instancelbs )
            cout << "<" << l.first << "," << l.second << "> ";
        cout << "\n";

        cout << "ubs: ";
        for( auto l : instanceubs )
            cout << "<" << l.first << "," << l.second << "> ";
        cout << endl;

        auto samebound = [&](pair<double, Cost> p1, pair<double, Cost> p2) {
            return p1.second == p2.second;
        };

        auto lbi = myunique(begin(instancelbs), end(instancelbs),
                          samebound);
        instancelbs.erase(lbi, end(instancelbs));

        auto ubi = myunique(begin(instanceubs), end(instanceubs),
                            samebound);
        instanceubs.erase(ubi, end(instanceubs));

        cout << "cleaned up: " << instancelbs.size() << " lbs, "
             << instanceubs.size() << " ubs" << endl;

        cout << "lbs: ";
        for( auto l : instancelbs )
            cout << "<" << l.first << "," << l.second << "> ";
        cout << "\n";

        cout << "ubs: ";
        for( auto l : instanceubs )
            cout << "<" << l.first << "," << l.second << "> ";
        cout << "\n";

        double totaltime;
        if( instanceubs.back().second != instancelbs.back().second )
            totaltime = TIMELIMIT;
        else
            totaltime = max(instancelbs.back().first, instanceubs.back().first);
        totaltime = TIMELIMIT;

        // normalize
        bool lbfault{false}, ubfault{false};
        vector<pair<double, double>> normlbs, normubs;
        for(auto& l : instancelbs) {
            double nlb;
            if( opt == lb0 )
                nlb = 1.0;
            else
                nlb = (l.second - lb0)/double(opt - lb0);
            cout << PV(l.second) << "\n";
            cout << PV(opt) << "\n";
            cout << PV(lb0) << "\n";
            cout << PV(l.second - lb0) << "\n";
            cout << PV(opt - lb0) << "\n";
            cout << PV( (l.second - lb0) < (opt -lb0) ) << "\n";
            if( 0.0 > nlb || nlb > 1.0)
                lbfault = true;
            auto clampednlb = min(1.0, max(nlb, 0.0));
            if( clampednlb != nlb )
                cout << "clamped from " << nlb << " to " << clampednlb << endl;
            nlb = clampednlb;
            normlbs.push_back( make_pair(l.first/totaltime, nlb));
        }

        cout << "ub0 = " << ub0 << "\n";
        for(auto& l : instanceubs) {
            double nub;
            if( ub0 == opt )
                nub = 1.0;
            else if( l.second > ub0 )
                nub = 2.0;
            else
                nub = 2 - (ub0-l.second)/double(ub0-opt);
            if( nub > 2.0 ) {
                cout << "having to normalize nub from " << nub << " to 2.0, orig ub "
                     << l.second << "\n";
                nub = 2.0;
            }
            if(1.0 > nub || nub > 2.0)
                ubfault = true;
            nub = min(2.0, max(nub, 1.0));
            if( l.first > TIMELIMIT )
                l.first = TIMELIMIT;
            normubs.push_back( make_pair(l.first/totaltime, nub));
        }

        cout << "normalized:\n";
        cout << "lbs: ";
        for( auto l : normlbs )
            cout << "<" << l.first << "," << l.second << "> ";
        cout << "\n";

        for(auto i = normlbs.begin(); i != normlbs.end(); ++i) {
            auto n = i;
            ++n;
            if( n == normlbs.end() )
                continue;
            if( i->second > n->second ) {
                cout << "at index " << distance(normlbs.begin(), i)
                     << ", <" << i->first << ", " << i->second << ">, then"
                     << ", <" << n->first << ", " << n->second << ">"
                     << endl;
            }
            assert((long double) i->second <= ((long double)n->second)*1.00001);
        }

        cout << "ubs: ";
        for( auto l : normubs )
            cout << "<" << l.first << "," << l.second << "> ";
        cout << "\n";

        if( lbfault )
            cout << "LBFAULT\n";
        if( ubfault )
            cout << "UBFAULT\n";

        // store everything
        lb0s[bestline] = lb0;
        opts[bestline] = opt;
        alllbs[bestline] = move(normlbs);
        allubs[bestline] = move(normubs);
    } while(true);

    cout << "read all points" << endl;

    int numlbpoints{1};
    for( auto & ilb : alllbs )
        numlbpoints += ilb.second.size()-1;
    cout << "will output " << numlbpoints << " lb points" << endl;
    int numubpoints{1};
    for( auto & iub : allubs )
        numubpoints += iub.second.size()-1;
    cout << "will output " << numubpoints << " ub points" << endl;

    ofstream normlbfs(algname+".LB");
    cout << "normalized average to \'" << (algname+".LB") << "\'" << endl;

    multimap<double, string> nextlbchange, nextubchange;

    double sumlb{0.0};
    int emptylbs{0};
    for(auto & ilb : alllbs ) {
        auto & v = ilb.second;
        reverse(begin(v), end(v));
        if( v.empty() ) {
            ++emptylbs;
            lastlb[ilb.first] = 0;
        } else {
            sumlb += v.back().second;
            lastlb[ilb.first] = v.back().second;
            v.pop_back();
        }
        if( !v.empty() )
            nextlbchange.insert(make_pair(v.back().first, ilb.first));
    }
    cout << PV(emptylbs) << endl;

    auto ninst = alllbs.size();
    normlbfs << "0.0 " << sumlb/ninst << "\n";
    while( !nextlbchange.empty() ) {
        auto first = *nextlbchange.begin();
        nextlbchange.erase(nextlbchange.begin());
        double nexttime = first.first; //awesome
        string& iname = first.second;
        auto & v = alllbs[iname];
        //cout << iname << " had " << lastlb[iname] << " before, now " << v.back().second << "\n";
        assert(v.back().second >= lastlb[iname]);
        sumlb += v.back().second - lastlb[iname];
        normlbfs << nexttime << " " << sumlb/ninst << "\n";
        lastlb[iname] = v.back().second;
        v.pop_back();
        if( !v.empty() )
            nextlbchange.insert(make_pair(v.back().first, iname));
    }

    normlbfs << 1.1 << " " << sumlb/ninst << "\n";

    ofstream normubfs(algname+".UB");
    cout << "normalized average to \'" << (algname+".UB") << "\'" << endl;
    double sumub{0.0};
    for(auto & iub : allubs ) {
        auto & v = iub.second;
        reverse(begin(v), end(v));
        sumub += v.back().second;
        lastub[iub.first] = v.back().second;
        v.pop_back();
        if( !v.empty() )
            nextubchange.insert(make_pair(v.back().first, iub.first));
    }

    normubfs << "0.0 " << sumub/ninst << "\n";
    while( !nextubchange.empty() ) {
        auto first = *nextubchange.begin();
        nextubchange.erase(nextubchange.begin());
        double nexttime = first.first; //awesome
        string& iname = first.second;
        auto & v = allubs[iname];
        //cout << iname << " had " << lastub[iname] << " before, now " << v.back().second << "\n";
        sumub += v.back().second - lastub[iname];
        normubfs << nexttime << " " << sumub/ninst << "\n";
        lastub[iname] = v.back().second;
        v.pop_back();
        if( !v.empty() )
            nextubchange.insert(make_pair(v.back().first, iname));
    }

    normubfs << 1.1 << " " << sumub/ninst << "\n";
}
