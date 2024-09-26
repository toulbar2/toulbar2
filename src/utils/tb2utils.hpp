/** \file tb2utils.hpp
 *  \brief Miscelaneous usefull functions.
 *
 */

#ifndef TB2UTILS_HPP_
#define TB2UTILS_HPP_

// these includes are needed if compiled on new g++ versions (>4.0?)
#include <cstdint>
#include <climits>
#include <cstdlib>
#include <cstring>
#include "libgen.h"

#ifdef OPENMPI
#define BOOST_MPI_HOMOGENEOUS
#include <boost/mpi.hpp>
namespace mpi = boost::mpi;
#endif

#ifdef ILOGLUE
#include <ilsolver/ilosolverint.h>
#else
#include <cassert>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#endif

#include <atomic>
#include <limits>
#include <iterator>
#include <vector>
#include <map>
#include <unordered_map>
#include <sstream>
#include <set>
#include <list>
#include <queue>
#include <stack>
#include <functional>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <memory>
using std::ostream;
using std::pair;
using std::vector;

#ifdef OPENMPI
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
namespace serialization = boost::serialization;
#include <boost/optional.hpp>
#endif

#ifdef NDEBUG
#define DEBONLY(x)
#else
#define DEBONLY(x) x
#endif

// #include <numbers>   // C++-20 std::numbers::pi
// const double PI = boost::math::constants::pi<double>();
const double PI = 3.1415926535897932384626433832795;

template <typename T1, typename T2, typename T3>
struct triplet {
    T1 first;
    T2 second;
    T3 third;
};

template <typename T1, typename T2, typename T3>
triplet<T1, T2, T3> make_triplet(const T1& m1, const T2& m2, const T3& m3)
{
    triplet<T1, T2, T3> ans;
    ans.first = m1;
    ans.second = m2;
    ans.third = m3;
    return ans;
}

template <typename T1, typename T2, typename T3>
std::ostream& operator<<(std::ostream& os, triplet<T1, T2, T3> const& p)
{
    return os << "triplet{" << p.first << "," << p.second << "," << p.third << "}";
}

template <typename U, typename T>
std::ostream& operator<<(std::ostream& os, std::pair<U, T> const& p)
{
    return os << "pair{" << p.first << "," << p.second << "}";
}

template <typename T>
std::ostream& operator<<(std::ostream& os, vector<T> const& v)
{
    os << "v(sz=" << v.size() << ")[";
    bool first = true;
    for (auto&& t : v) {
        if (first)
            first = false;
        else
            os << ",";
        os << t;
    }
    os << "]";
    return os;
}

// template<class T>
// T abs(T x) {
//     if (x < 0) return -(x);
//     else return x;
// }

// Warning! Already defined in STL
// template<class T>
// T min(T x, T y) {
//    if (x < y) return x;
//    else return y;
//}
//
// template<class T>
// T max(T x, T y) {
//    if (x > y) return x;
//    else return y;
//}

// template <class T>
// T min(T* array, int size)
//{
//     assert(size >= 1);
//     T res = array[0];
//     for (int i = 1; i < size; i++) {
//         if (array[i] < res) {
//             res = array[i];
//         }
//     }
//     return res;
// }

// template <class T>
// T max(T* array, int size)
//{
//     assert(size >= 1);
//     T res = array[0];
//     for (int i = 1; i < size; i++) {
//         if (array[i] > res) {
//             res = array[i];
//         }
//     }
//     return res;
// }

template <class T>
inline std::string to_string(const T& t)
{
    std::stringstream ss;
    ss << t;
    return ss.str();
}

static inline void rtrim(std::string& s)
{
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
        return !std::isspace(ch);
    }).base(),
        s.end());
}

template <typename T>
void free_all(T& t)
{
    T tmp;
    t.swap(tmp);
}

// warning! forbidden characters /#[]{}:, and spaces in cfn format for object names
static inline std::string name2cfn(std::string s)
{
    for(auto it = s.begin(); it != s.end(); it ++) {
        if(*it == '[' || *it == '{') {
            *it = '(';
        } else if(*it == ']' || *it == '}') {
            *it = ')';
        } else if(*it == '/' || *it == '#' || *it == ':' || *it == ',' || *it == ' ' || *it == '\t') {
            *it = '_';
        }
    }
    return s;
}

#include "tb2system.hpp"

// Cormen et al, 1990. pages 152, 158, and 184

template <class T>
int partition(T A[], int p, int r)
{
    T x = A[p];
    int i = p - 1;
    int j = r + 1;
    while (true) {
        do {
            j = j - 1;
        } while (A[j] > x);
        do {
            i = i + 1;
        } while (A[i] < x);
        if (i < j) {
            T tmp = A[i];
            A[i] = A[j];
            A[j] = tmp;
        } else
            return j;
    }
}

template <class T>
int stochastic_partition(T A[], int p, int r)
{
    int i = (myrand() % (r - p + 1)) + p;
    T tmp = A[p];
    A[p] = A[i];
    A[i] = tmp;
    return partition(A, p, r);
}

template <class T>
T stochastic_selection(T A[], int p, int r, int i)
{
    if (p == r)
        return A[p];
    int q = stochastic_partition(A, p, r);
    int k = q - p + 1;
    if (i <= k)
        return stochastic_selection(A, p, q, i);
    else
        return stochastic_selection(A, q + 1, r, i - k);
}

#endif /* TB2UTILS_HPP_ */

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
