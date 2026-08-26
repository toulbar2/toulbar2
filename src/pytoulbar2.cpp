/** \file pytoulbar2.cpp
 *  \brief Python wrapper to toulbar2 library
 *
<pre>
    Copyright (c) 2006-2020, toulbar2 team

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

    toulbar2 is currently maintained by Simon de Givry, INRAE - MIAT, Toulouse, France (simon.de-givry@inrae.fr)
</pre>
 */

// How to manually extract class properties to bind in Python:
//  awk '/^class /{ok=1;class=$2} go&&/static/{gsub(";.*","",$0); print "        .def_readwrite_static(\"" $NF "\", &" class "::" $NF ")"} ok&&/public/{go=1}' core/tb2types.hpp
//  awk '/^class/{class=$2} /virtual/{gsub("//.*","",$0);gsub("[(].*[)].*","",$0); print "        .def(\"" $NF "\", &" class "::" $NF ")"}' toulbar2lib.hpp

// How to compile Python3 pytb2 module library on Linux:
//  apt install pybind11-dev (or else pip3 install pybind11)
//  git clone https://github.com/toulbar2/toulbar2.git
//  cd toulbar2; mkdir build; cd build
//  #compile toulbar2 to produce the python C++ library
//  cmake -DPYTB2=ON ..
//  make
//  the module will be in lib/Linux

// Examples using pytb2 module from Python3:
//  NB: pytb2.cpython* must be in your Python3 path or export PYTHONPATH=.
//  python3 -c "import sys; sys.path.append('.'); import pytb2 as tb2; tb2.init(); m = tb2.Solver(); m.read('../validation/default/example.wcsp'); tb2.option.showSolutions = 1; res = m.solve(); print(res); print(m.solutions())"
//  python3 -c "import sys; sys.path.append('.'); import pytb2 as tb2; tb2.init(); m = tb2.Solver(); m.read('../validation/default/1aho.cfn.gz'); res = m.solve(); print(res); print(m.wcsp.getDPrimalBound()); print(m.solution())"
//  python3 -c "import sys; sys.path.append('.'); import random; import pytb2 as tb2; tb2.init(); m = tb2.Solver(); x=m.wcsp.makeEnumeratedVariable('x', 1, 10); y=m.wcsp.makeEnumeratedVariable('y', 1, 10); z=m.wcsp.makeEnumeratedVariable('z', 1, 10); m.wcsp.postUnaryConstraint(x, [random.randint(0,10) for i in range(10)]); m.wcsp.postUnaryConstraint(y, [random.randint(0,10) for i in range(10)]); m.wcsp.postUnaryConstraint(z, [random.randint(0,10) for i in range(10)]); m.wcsp.postBinaryConstraint(x,y, [random.randint(0,10) for i in range(10) for j in range(10)]); m.wcsp.postBinaryConstraint(x,z,[random.randint(0,10) for i in range(10) for j in range(10)]); m.wcsp.postBinaryConstraint(y,z,[random.randint(0,10) for i in range(10) for j in range(10)]); m.wcsp.sortConstraints(); res = m.solve(); print(res); print(m.wcsp.getDPrimalBound()); print(m.solution());"
//  python3 -c "import sys; sys.path.append('.'); import random; import pytb2 as tb2; tb2.init(); m = tb2.Solver(); tb2.option.verbose = 0; tb2.option.elimDegree_preprocessing=1; tb2.check(); x=m.wcsp.makeEnumeratedVariable('x', 1, 10); y=m.wcsp.makeEnumeratedVariable('y', 1, 10); z=m.wcsp.makeEnumeratedVariable('z', 1, 10); w=m.wcsp.makeEnumeratedVariable('w', 1, 10); m.wcsp.postUnaryConstraint(x, [random.randint(0,10) for i in range(10)]); m.wcsp.postUnaryConstraint(y, [random.randint(0,10) for i in range(10)]); m.wcsp.postUnaryConstraint(z, [random.randint(0,10) for i in range(10)]); m.wcsp.postBinaryConstraint(x,y, [random.randint(0,10) for i in range(10) for j in range(10)]); m.wcsp.postBinaryConstraint(x,z,[random.randint(0,10) for i in range(10) for j in range(10)]); m.wcsp.postBinaryConstraint(y,z,[random.randint(0,10) for i in range(10) for j in range(10)]); nary = m.wcsp.postNaryConstraintBegin([x,y,z,w], 10, 1); m.wcsp.postNaryConstraintTuple(nary, [1,1,1,1], 0); m.wcsp.postNaryConstraintEnd(nary); m.wcsp.sortConstraints(); res = m.solve(); print(res); print(m.wcsp.getDPrimalBound()); print(m.solution());"

#include "core/tb2types.hpp"
#include "utils/tb2system.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

// PYBIND11_MAKE_OPAQUE(std::vector<int>);

namespace py = pybind11;

#include "toulbar2lib.hpp"
#include "utils/tb2store.hpp"
#include "utils/tb2btlist.hpp"
#include "search/tb2solver.hpp"
#include "core/tb2wcsp.hpp"
#include "mcriteria/multicfn.hpp"
#include "mcriteria/bicriteria.hpp"

extern void newsolution(int wcspId, void* solver);

// return true if the data type is equivalent to a a signed integer
inline bool is_dtype_sintegers(py::buffer_info& buf_info)
{
    return buf_info.item_type_is_equivalent_to<int8_t>() || buf_info.item_type_is_equivalent_to<int16_t>() || buf_info.item_type_is_equivalent_to<int32_t>() || buf_info.item_type_is_equivalent_to<int64_t>();
}

// return true if the data type is equivalent to a floating point type
inline bool is_dtype_floating_point(py::buffer_info& buf_info)
{
    return buf_info.item_type_is_equivalent_to<float>() || buf_info.item_type_is_equivalent_to<double>();
}

// create n enumerated variables
// var names are formed from base_name plus an index
// return the index of the first variable created
int makeEnumeratedVariableVec(WeightedCSP& s, int n, std::string base_name, Value iinf, Value isup)
{
    int result = s.makeEnumeratedVariable(base_name + "_0", iinf, isup);
    for (size_t ind = 1; ind < static_cast<size_t>(n); ind++) {
        s.makeEnumeratedVariable(base_name + "_" + to_string(ind), iinf, isup);
    }
    return result;
}

template <typename T>
inline void extractUnaryScopes(std::vector<int>& unary_scopes, py::buffer_info& scopes_info)
{
    unary_scopes.resize(scopes_info.shape[0]);
    size_t stride = scopes_info.strides[0] / scopes_info.itemsize;
    T* temp_ptr = static_cast<T*>(scopes_info.ptr);
    for (int i = 0; i < scopes_info.shape[0]; i++) {
        unary_scopes[i] = static_cast<int>(temp_ptr[i * stride]);
    }
}

template <typename T>
inline void extractUnaryCosts(size_t func_ind, vector<Double>& costs, py::buffer_info& costs_info, size_t s1, size_t s2)
{
    T* temp_ptr = static_cast<T*>(costs_info.ptr);
    for (int j = 0; j < costs_info.shape[1]; j++) {
        costs[j] = static_cast<Double>(temp_ptr[func_ind * s1 + j * s2]);
    }
}

// post several unary cost functions at once
// all variables are expected to have the same domain size
// scopes is expected to be a 1-dimensional integer array
// costs is expected to be a 2 dimensional array: n_variables x domain_size
void postUnaryVecConstraints(WeightedCSP& s, py::buffer& scopes, py::buffer& costs, bool incremental)
{

    /* Request a buffer descriptor from Python */
    py::buffer_info scopes_info = scopes.request();
    py::buffer_info costs_info = costs.request();

    // check scope size
    if (scopes_info.ndim != 1 || !is_dtype_sintegers(scopes_info)) {
        std::cerr << "Error, scopes must be provided as a 1-dimensional vector of integer indices!" << std::endl;
        throw BadConfiguration();
    }
    // check costs data types
    if (costs_info.ndim != 2) {
        // error, costs must be floating point values
        std::cerr << "error, costs must be a 2-dimensional vector!" << std::endl;
        throw BadConfiguration();
    }
    // check consistency between costs sizes and scopes sizes
    if (costs_info.shape[0] != scopes_info.shape[0]) {
        // error, costs must be floating point values
        std::cerr << "error, costs and scopes must have same size!" << std::endl;
        throw BadConfiguration();
    }

    // read the scopes
    std::vector<int> unary_scopes;
    if (scopes_info.item_type_is_equivalent_to<int8_t>()) {
        extractUnaryScopes<int8_t>(unary_scopes, scopes_info);
    } else if (scopes_info.item_type_is_equivalent_to<int16_t>()) {
        extractUnaryScopes<int16_t>(unary_scopes, scopes_info);
    } else if (scopes_info.item_type_is_equivalent_to<int32_t>()) {
        extractUnaryScopes<int32_t>(unary_scopes, scopes_info);
    } else if (scopes_info.item_type_is_equivalent_to<int64_t>()) {
        extractUnaryScopes<int64_t>(unary_scopes, scopes_info);
    } else {
        std::cerr << "error, unsupported data types for scopes!" << std::endl;
        throw BadConfiguration();
    }

    // read the costs and create the cost functions
    size_t s1 = costs_info.strides[0] / costs_info.itemsize; // var ind
    size_t s2 = costs_info.strides[1] / costs_info.itemsize; // val ind
    vector<Double> unary_costs(costs_info.shape[1]);
    for (int i = 0; i < costs_info.shape[0]; i++) {
        if (costs_info.item_type_is_equivalent_to<float>()) {
            extractUnaryCosts<float>(i, unary_costs, costs_info, s1, s2);
        } else if (costs_info.item_type_is_equivalent_to<double>()) {
            extractUnaryCosts<double>(i, unary_costs, costs_info, s1, s2);
        } else if (costs_info.item_type_is_equivalent_to<long double>()) {
            extractUnaryCosts<long double>(i, unary_costs, costs_info, s1, s2);
        } else if (costs_info.item_type_is_equivalent_to<int8_t>()) {
            extractUnaryCosts<int8_t>(i, unary_costs, costs_info, s1, s2);
        } else if (costs_info.item_type_is_equivalent_to<int16_t>()) {
            extractUnaryCosts<int16_t>(i, unary_costs, costs_info, s1, s2);
        } else if (costs_info.item_type_is_equivalent_to<int32_t>()) {
            extractUnaryCosts<int32_t>(i, unary_costs, costs_info, s1, s2);
        } else if (costs_info.item_type_is_equivalent_to<int64_t>()) {
            extractUnaryCosts<int64_t>(i, unary_costs, costs_info, s1, s2);
        } else { // unsupported
            std::cerr << "error, unsupported data types for costs!" << std::endl;
            throw BadConfiguration();
        }
        s.postUnaryConstraint(unary_scopes[i], unary_costs, incremental);
    }
}

// extract a binary cost function table
template <typename T>
inline void extractBinaryCosts(std::vector<Double>& binary_costs, py::buffer_info& costs_info)
{
    binary_costs.resize(costs_info.shape[0] * costs_info.shape[1]);
    size_t s1 = costs_info.strides[0] / costs_info.itemsize;
    size_t s2 = costs_info.strides[1] / costs_info.itemsize;
    T* temp_ptr = static_cast<T*>(costs_info.ptr);
    size_t cost_ind = 0;
    for (int i = 0; i < costs_info.shape[0]; i++) {
        for (int j = 0; j < costs_info.shape[1]; j++) {
            binary_costs[cost_ind] = static_cast<Double>(temp_ptr[i * s1 + j * s2]);
            cost_ind++;
        }
    }
}

// extract a single binary scope given a cost function index
template <typename T>
inline void extractBinaryScopes(size_t func_ind, py::buffer_info& scopes_info, size_t ss1, size_t ss2, int& xIndex, int& yIndex)
{
    T* temp_ptr = static_cast<T*>(scopes_info.ptr);
    xIndex = temp_ptr[func_ind * ss1];
    yIndex = temp_ptr[func_ind * ss1 + ss2];
}

inline void extractBinaryScopes(size_t func_ind, py::buffer_info& scopes_info, size_t ss1, size_t ss2, int& xIndex, int& yIndex)
{
    if (scopes_info.item_type_is_equivalent_to<int8_t>()) {
        extractBinaryScopes<int8_t>(func_ind, scopes_info, ss1, ss2, xIndex, yIndex);
    } else if (scopes_info.item_type_is_equivalent_to<int8_t>()) {
        extractBinaryScopes<int16_t>(func_ind, scopes_info, ss1, ss2, xIndex, yIndex);
    } else if (scopes_info.item_type_is_equivalent_to<int32_t>()) {
        extractBinaryScopes<int32_t>(func_ind, scopes_info, ss1, ss2, xIndex, yIndex);
    } else if (scopes_info.item_type_is_equivalent_to<int64_t>()) {
        extractBinaryScopes<int64_t>(func_ind, scopes_info, ss1, ss2, xIndex, yIndex);
    } else {
        std::cerr << "error, unsupported data types for scopes!" << std::endl;
        throw BadConfiguration();
    }
}

// post several binary cost functions with a single costs tensor
// scopes are provided as a n_func*2 matrix
int postBinaryVecConstraints(WeightedCSP& s, py::buffer& scopes, py::buffer& costs, bool incremental)
{

    // scope indices
    int xIndex, yIndex;
    int result = -1;

    /* Request a buffer descriptor from Python */
    py::buffer_info scopes_info = scopes.request();
    py::buffer_info costs_info = costs.request();

    // check scope size
    if (scopes_info.ndim != 2 || scopes_info.shape[1] != 2) {
        std::cerr << "Error, must provide a list of scopes with size 2!" << std::endl;
        throw BadConfiguration();
    }
    // check scopes data type
    if (!is_dtype_sintegers(scopes_info)) {
        std::cerr << "error, scopes must be integers values!" << std::endl;
        throw BadConfiguration();
    }
    if (costs_info.ndim != 2) {
        // error, costs must be floating point values
        std::cerr << "error, costs must be a (dom_size x dom_size x dom_size) tensor!" << std::endl;
        throw BadConfiguration();
    }

    // read the costs
    size_t ss1 = scopes_info.strides[0] / scopes_info.itemsize;
    size_t ss2 = scopes_info.strides[1] / scopes_info.itemsize;
    vector<Double> binary_costs;
    if (costs_info.item_type_is_equivalent_to<int8_t>()) {
        extractBinaryCosts<int8_t>(binary_costs, costs_info);
    } else if (costs_info.item_type_is_equivalent_to<int16_t>()) {
        extractBinaryCosts<int16_t>(binary_costs, costs_info);
    } else if (costs_info.item_type_is_equivalent_to<int32_t>()) {
        extractBinaryCosts<int32_t>(binary_costs, costs_info);
    } else if (costs_info.item_type_is_equivalent_to<int64_t>()) {
        extractBinaryCosts<int64_t>(binary_costs, costs_info);
    } else if (costs_info.item_type_is_equivalent_to<float>()) {
        extractBinaryCosts<float>(binary_costs, costs_info);
    } else if (costs_info.item_type_is_equivalent_to<double>()) {
        extractBinaryCosts<double>(binary_costs, costs_info);
    } else if (costs_info.item_type_is_equivalent_to<long double>()) {
        extractBinaryCosts<long double>(binary_costs, costs_info);
    } else { // unsupported
        std::cerr << "error, costs must be float or double!" << std::endl;
        throw BadConfiguration();
    }
    // read the scopes and create the cost functions
    for (int i = 0; i < scopes_info.shape[0]; i++) {
        extractBinaryScopes(i, scopes_info, ss1, ss2, xIndex, yIndex);
        int temp_result = s.postBinaryConstraint(xIndex, yIndex, binary_costs, incremental);
        if (i == 0) {
            result = temp_result;
        }
    }

    return result; // return index of the first added cost function
}

// extract multiple binary costs function scopes
template <typename T>
void extractBinaryScopes(vector<vector<int>>& binary_scopes, py::buffer_info& scopes_info)
{
    binary_scopes.resize(scopes_info.shape[0], std::vector<int>(2));
    T* temp_ptr = static_cast<T*>(scopes_info.ptr);
    size_t ss1 = scopes_info.strides[0] / scopes_info.itemsize;
    size_t ss2 = scopes_info.strides[1] / scopes_info.itemsize;
    for (int i = 0; i < scopes_info.shape[0]; i++) {
        binary_scopes[i][0] = temp_ptr[i * ss1];
        binary_scopes[i][1] = temp_ptr[i * ss1 + ss2];
    }
}

// extract a binary cost function table
template <typename T>
int extractBinaryCosts(py::buffer_info& costs_info, vector<vector<int>>& binary_scopes, WeightedCSP& s, bool incremental)
{
    int result = -1;
    T* temp_ptr = static_cast<T*>(costs_info.ptr);
    size_t s1 = costs_info.strides[0] / costs_info.itemsize; // n cost functions
    size_t s2 = costs_info.strides[1] / costs_info.itemsize;
    size_t s3 = costs_info.strides[2] / costs_info.itemsize;
    std::vector<Double> binary_costs(costs_info.shape[1] * costs_info.shape[2]);
    for (size_t func_ind = 0; func_ind < binary_scopes.size(); func_ind++) {
        size_t cost_ind = 0;
        for (int j = 0; j < costs_info.shape[1]; j++) {
            for (int k = 0; k < costs_info.shape[2]; k++) {
                binary_costs[cost_ind] = static_cast<Double>(temp_ptr[func_ind * s1 + j * s2 + k * s3]);
                cost_ind++;
            }
        }
        int result_temp = s.postBinaryConstraint(binary_scopes[func_ind][0], binary_scopes[func_ind][1], binary_costs, incremental);
        if (result < 0) {
            result = result_temp;
        }
    }
    return result;
}

// post several binary cost functions from tensors
// scopes are provided as a n_func*2 matrix
// costs are provided as a n_func*dom_size*dom_size tensor
int postMultBinaryVecConstraints(WeightedCSP& s, py::buffer& scopes, py::buffer& costs, bool incremental)
{

    // scope indices
    int result = -1;

    /* Request a buffer descriptor from Python */
    py::buffer_info scopes_info = scopes.request();
    py::buffer_info costs_info = costs.request();

    // check scope size
    if (scopes_info.ndim != 2 || scopes_info.shape[1] != 2) {
        std::cerr << "Error, must provide a list of scopes with size 2!" << std::endl;
        throw BadConfiguration();
    }
    // check that there are as many scopes are cost tables
    if (scopes_info.shape[0] != costs_info.shape[0]) {
        std::cerr << "Error, must provide the same number of scopes and cost tables!" << std::endl;
        throw BadConfiguration();
    }
    // check scopes data type
    if (!is_dtype_sintegers(scopes_info)) {
        std::cerr << "error, scopes must be integers values!" << std::endl;
        throw BadConfiguration();
    }
    // check costs shape
    if (costs_info.ndim != 3) {
        std::cerr << "Error, must provide costs as a n_function x dom_size_1 x dom_size_2 tensor!" << std::endl;
        throw BadConfiguration();
    }

    // read the scopes and create the cost functions
    vector<vector<int>> binary_scopes;
    if (scopes_info.item_type_is_equivalent_to<int8_t>()) {
        extractBinaryScopes<int8_t>(binary_scopes, scopes_info);
    } else if (scopes_info.item_type_is_equivalent_to<int16_t>()) {
        extractBinaryScopes<int16_t>(binary_scopes, scopes_info);
    } else if (scopes_info.item_type_is_equivalent_to<int32_t>()) {
        extractBinaryScopes<int32_t>(binary_scopes, scopes_info);
    } else if (scopes_info.item_type_is_equivalent_to<int64_t>()) {
        extractBinaryScopes<int64_t>(binary_scopes, scopes_info);
    } else {
        std::cerr << "error, unsupported data types for scopes!" << std::endl;
        throw BadConfiguration();
    }

    // read the costs and post the binary functions
    if (costs_info.item_type_is_equivalent_to<float>()) {
        result = extractBinaryCosts<float>(costs_info, binary_scopes, s, incremental);
    } else if (costs_info.item_type_is_equivalent_to<double>()) {
        result = extractBinaryCosts<double>(costs_info, binary_scopes, s, incremental);
    } else if (costs_info.item_type_is_equivalent_to<long double>()) {
        result = extractBinaryCosts<long double>(costs_info, binary_scopes, s, incremental);
    } else if (costs_info.item_type_is_equivalent_to<int8_t>()) {
        result = extractBinaryCosts<int8_t>(costs_info, binary_scopes, s, incremental);
    } else if (costs_info.item_type_is_equivalent_to<int16_t>()) {
        result = extractBinaryCosts<int16_t>(costs_info, binary_scopes, s, incremental);
    } else if (costs_info.item_type_is_equivalent_to<int32_t>()) {
        result = extractBinaryCosts<int32_t>(costs_info, binary_scopes, s, incremental);
    } else if (costs_info.item_type_is_equivalent_to<int64_t>()) {
        result = extractBinaryCosts<int64_t>(costs_info, binary_scopes, s, incremental);
    } else { // unsupported cost types
        std::cerr << "error, unsupported costs data type!" << std::endl;
        throw BadConfiguration();
    }

    return result; // return index of the first added cost function
}

template <typename T>
void extractTernaryCosts(vector<Double>& ternary_costs, py::buffer_info& costs_info)
{
    ternary_costs.resize(costs_info.shape[0] * costs_info.shape[1] * costs_info.shape[2]);
    // read the costs
    size_t s1 = costs_info.strides[0] / costs_info.itemsize;
    size_t s2 = costs_info.strides[1] / costs_info.itemsize;
    size_t s3 = costs_info.strides[2] / costs_info.itemsize;
    T* temp_ptr = static_cast<T*>(costs_info.ptr);
    size_t cost_ind = 0;
    for (int i = 0; i < costs_info.shape[0]; i++) {
        for (int j = 0; j < costs_info.shape[1]; j++) {
            for (int k = 0; k < costs_info.shape[2]; k++) {
                ternary_costs[cost_ind] = static_cast<Double>(temp_ptr[i * s1 + j * s2 + k * s3]);
                cost_ind++;
            }
        }
    }
}

// extract a single binary scope given a cost function index
template <typename T>
inline void extractTernaryScopes(size_t func_ind, py::buffer_info& scopes_info, size_t ss1, size_t ss2, int& xIndex, int& yIndex, int& zIndex)
{
    T* temp_ptr = static_cast<T*>(scopes_info.ptr);
    xIndex = temp_ptr[func_ind * ss1];
    yIndex = temp_ptr[func_ind * ss1 + ss2];
    zIndex = temp_ptr[func_ind * ss1 + ss2 * 2];
}

inline void extractTernaryScopes(size_t func_ind, py::buffer_info& scopes_info, size_t ss1, size_t ss2, int& xIndex, int& yIndex, int& zIndex)
{
    if (scopes_info.item_type_is_equivalent_to<int8_t>()) {
        extractTernaryScopes<uint8_t>(func_ind, scopes_info, ss1, ss2, xIndex, yIndex, zIndex);
    } else if (scopes_info.item_type_is_equivalent_to<int16_t>()) {
        extractTernaryScopes<uint16_t>(func_ind, scopes_info, ss1, ss2, xIndex, yIndex, zIndex);
    } else if (scopes_info.item_type_is_equivalent_to<int32_t>()) {
        extractTernaryScopes<uint32_t>(func_ind, scopes_info, ss1, ss2, xIndex, yIndex, zIndex);
    } else if (scopes_info.item_type_is_equivalent_to<int64_t>()) {
        extractTernaryScopes<uint64_t>(func_ind, scopes_info, ss1, ss2, xIndex, yIndex, zIndex);
    } else {
        std::cerr << "error, unsupported data types for scopes!" << std::endl;
        throw BadConfiguration();
    }
}

// post several ternary cost functions with a single costs tensor
// scopes are provided as a n_func*3 matrix
int postTernaryVecConstraints(WeightedCSP& s, py::buffer& scopes, py::buffer& costs, bool incremental)
{

    // scope indices
    int xIndex, yIndex, zIndex;
    int result = -1;

    /* Request a buffer descriptor from Python */
    py::buffer_info scopes_info = scopes.request();
    py::buffer_info costs_info = costs.request();

    // check scope size
    if (scopes_info.ndim != 2 || scopes_info.shape[1] != 3) {
        std::cerr << "Error, must provide a list of scopes with size 3!" << std::endl;
        throw BadConfiguration();
    }
    // check scopes data type
    if (!is_dtype_sintegers(scopes_info)) {
        std::cerr << "error, scopes must be unsigned integer values!" << std::endl;
        throw BadConfiguration();
    }
    if (costs_info.ndim != 3) {
        // error, costs must be floating point values
        std::cerr << "error, costs must be a 3x3x3 tensor!" << std::endl;
        throw BadConfiguration();
    }

    // read the costs table
    size_t ss1 = scopes_info.strides[0] / scopes_info.itemsize;
    size_t ss2 = scopes_info.strides[1] / scopes_info.itemsize;

    vector<Double> ternary_costs;
    if (costs_info.item_type_is_equivalent_to<double>()) {
        extractTernaryCosts<double>(ternary_costs, costs_info);
    } else if (costs_info.item_type_is_equivalent_to<float>()) {
        extractTernaryCosts<float>(ternary_costs, costs_info);
    } else if (costs_info.item_type_is_equivalent_to<long double>()) {
        extractTernaryCosts<long double>(ternary_costs, costs_info);
    } else if (costs_info.item_type_is_equivalent_to<int8_t>()) {
        extractTernaryCosts<int8_t>(ternary_costs, costs_info);
    } else if (costs_info.item_type_is_equivalent_to<int16_t>()) {
        extractTernaryCosts<int16_t>(ternary_costs, costs_info);
    } else if (costs_info.item_type_is_equivalent_to<int32_t>()) {
        extractTernaryCosts<int32_t>(ternary_costs, costs_info);
    } else if (costs_info.item_type_is_equivalent_to<int64_t>()) {
        extractTernaryCosts<int64_t>(ternary_costs, costs_info);
    } else { // unsupported
        std::cerr << "error, unsupported costs data type!" << std::endl;
        throw BadConfiguration();
    }

    // read the scopes and create the cost functions
    for (int i = 0; i < scopes_info.shape[0]; i++) {
        extractTernaryScopes(i, scopes_info, ss1, ss2, xIndex, yIndex, zIndex);
        int temp_res = s.postTernaryConstraint(xIndex, yIndex, zIndex, ternary_costs, incremental);
        if (i == 0) {
            result = temp_res;
        }
    }

    return result; // return index of the first added cost function
}

// extract multiple binary costs function scopes
template <typename T>
inline void extractTernaryScopes(vector<vector<int>>& ternary_scopes, py::buffer_info& scopes_info)
{
    ternary_scopes.resize(scopes_info.shape[0], std::vector<int>(3));
    T* temp_ptr = static_cast<T*>(scopes_info.ptr);
    size_t ss1 = scopes_info.strides[0] / scopes_info.itemsize;
    size_t ss2 = scopes_info.strides[1] / scopes_info.itemsize;
    for (int i = 0; i < scopes_info.shape[0]; i++) {
        ternary_scopes[i][0] = temp_ptr[i * ss1];
        ternary_scopes[i][1] = temp_ptr[i * ss1 + ss2];
        ternary_scopes[i][2] = temp_ptr[i * ss1 + ss2 * 2];
    }
}

// extract a ternary cost function table
template <typename T>
int extractTernaryCosts(vector<vector<int>>& ternary_scopes, py::buffer_info& costs_info, WeightedCSP& s, bool incremental)
{

    int result = -1;

    size_t s1 = costs_info.strides[0] / costs_info.itemsize; // n cost functions
    size_t s2 = costs_info.strides[1] / costs_info.itemsize;
    size_t s3 = costs_info.strides[2] / costs_info.itemsize;
    size_t s4 = costs_info.strides[3] / costs_info.itemsize;

    std::vector<Double> ternary_costs(costs_info.shape[1] * costs_info.shape[2] * costs_info.shape[3]);

    T* temp_ptr = static_cast<T*>(costs_info.ptr);
    for (int i = 0; i < costs_info.shape[0]; i++) { // cost function loop
        size_t cost_ind = 0;
        for (int j = 0; j < costs_info.shape[1]; j++) {
            for (int k = 0; k < costs_info.shape[2]; k++) {
                for (int l = 0; l < costs_info.shape[3]; l++) {
                    ternary_costs[cost_ind] = static_cast<Double>(temp_ptr[i * s1 + j * s2 + k * s3 + l * s4]);
                    cost_ind++;
                }
            }
        }
        int result_temp = s.postTernaryConstraint(ternary_scopes[i][0], ternary_scopes[i][1], ternary_scopes[i][2], ternary_costs, incremental);
        if (result < 0) {
            result = result_temp;
        }
    }
    return result;
}

// post several ternary cost functions from tensors
// scopes are provided as a n_func*3 matrix
// costs are provided as a n_func*dom_size*dom_size*dom_size tensor
int postMultTernaryVecConstraints(WeightedCSP& s, py::buffer& scopes, py::buffer& costs, bool incremental)
{

    // scope indices
    int result = -1;

    /* Request a buffer descriptor from Python */
    py::buffer_info scopes_info = scopes.request();
    py::buffer_info costs_info = costs.request();

    // check scope size
    if (scopes_info.ndim != 2 || scopes_info.shape[1] != 3) {
        std::cerr << "Error, must provide a list of scopes with size 3!" << std::endl;
        throw BadConfiguration();
    }
    // check that there are as many scopes are cost tables
    if (scopes_info.shape[0] != costs_info.shape[0]) {
        std::cerr << "Error, must provide the same number of scopes and cost tables!" << std::endl;
        throw BadConfiguration();
    }
    // check scopes data type
    if (!is_dtype_sintegers(scopes_info)) {
        std::cerr << "error, scopes must be integers values!" << std::endl;
        throw BadConfiguration();
    }
    // check costs shape
    if (costs_info.ndim != 4) {
        std::cerr << "Error, must provide costs as a n_function x dom_size_1 x dom_size_2 x dom_size_3 tensor!" << std::endl;
        throw BadConfiguration();
    }

    // read the scopes and create the cost functions
    vector<vector<int>> ternary_scopes;
    if (scopes_info.item_type_is_equivalent_to<int8_t>()) {
        extractTernaryScopes<int8_t>(ternary_scopes, scopes_info);
    } else if (scopes_info.item_type_is_equivalent_to<int16_t>()) {
        extractTernaryScopes<int16_t>(ternary_scopes, scopes_info);
    } else if (scopes_info.item_type_is_equivalent_to<int32_t>()) {
        extractTernaryScopes<int32_t>(ternary_scopes, scopes_info);
    } else if (scopes_info.item_type_is_equivalent_to<int64_t>()) {
        extractTernaryScopes<int64_t>(ternary_scopes, scopes_info);
    } else {
        std::cerr << "error, unsupported data types for scopes!" << std::endl;
        throw BadConfiguration();
    }

    if (costs_info.item_type_is_equivalent_to<float>()) {
        result = extractTernaryCosts<float>(ternary_scopes, costs_info, s, incremental);
    } else if (costs_info.item_type_is_equivalent_to<double>()) {
        result = extractTernaryCosts<double>(ternary_scopes, costs_info, s, incremental);
    } else if (costs_info.item_type_is_equivalent_to<long double>()) {
        result = extractTernaryCosts<long double>(ternary_scopes, costs_info, s, incremental);
    } else if (costs_info.item_type_is_equivalent_to<int8_t>()) {
        result = extractTernaryCosts<int8_t>(ternary_scopes, costs_info, s, incremental);
    } else if (costs_info.item_type_is_equivalent_to<int16_t>()) {
        result = extractTernaryCosts<int16_t>(ternary_scopes, costs_info, s, incremental);
    } else if (costs_info.item_type_is_equivalent_to<int32_t>()) {
        result = extractTernaryCosts<int32_t>(ternary_scopes, costs_info, s, incremental);
    } else if (costs_info.item_type_is_equivalent_to<int64_t>()) {
        result = extractTernaryCosts<int64_t>(ternary_scopes, costs_info, s, incremental);
    } else { // unsupported
        std::cerr << "error, unsupported costs data type!" << std::endl;
        throw BadConfiguration();
    }

    return result; // return index of the first added cost function
}

PYBIND11_MODULE(pytb2, m)
{
    m.def("init", []() { tb2init(); }); // must be called at the very beginning
    m.def("reinit", []() { tb2reinit(); }); // must be called before solving again
    m.attr("MAX_COST") = py::int_(MAX_COST);
    m.attr("MIN_COST") = py::int_(MIN_COST);

    py::register_exception<Contradiction>(m, "Contradiction");
    py::register_exception<InternalError>(m, "InternalError");
    py::register_exception<SolverOut>(m, "SolverOut");

    py::class_<ToulBar2, std::unique_ptr<ToulBar2, py::nodelete>>(m, "option")

        .def("writeSolution", [](const string& fileName) {
            if (atoi(fileName.data()) != 0) {
                ToulBar2::writeSolution = atoi(fileName.data());
            } else {
                ToulBar2::solutionFile = fopen(fileName.data(), "w");
            }
        })
        .def("closeSolution", []() {
            fclose(ToulBar2::solutionFile);
        })
        .def("getVarOrder", []() -> string {
            if (reinterpret_cast<uintptr_t>(ToulBar2::varOrder) < ELIM_MAX) {
                return to_string(-reinterpret_cast<uintptr_t>(ToulBar2::varOrder));
            } else {
                return to_string(ToulBar2::varOrder);
            }
        })
        .def("setVarOrder", [](const string& fileName) {
            if (fileName[0] == '-') {
                ToulBar2::varOrder = reinterpret_cast<char*>(-atoi(fileName.data()));
            } else {
                ToulBar2::varOrder = new char[fileName.size() + 1];
                strcpy(ToulBar2::varOrder, fileName.data());
            }
        })

        .def_property_readonly_static("version", [](py::object) { return ToulBar2::version; })
        .def_property_static("verbose", [](py::object) { return ToulBar2::verbose; }, [](py::object, int verbose) { ToulBar2::verbose = verbose; })
        .def_property_static("debug", [](py::object) { return ToulBar2::debug; }, [](py::object, int debug) { ToulBar2::debug = debug; })
        .def_property_static("externalUB", [](py::object) { return ToulBar2::externalUB; }, [](py::object, string externalUB) { ToulBar2::externalUB = externalUB; })
        .def_property_static("showSolutions", [](py::object) { return ToulBar2::showSolutions; }, [](py::object, int showSolutions) { ToulBar2::showSolutions = showSolutions; })

        .def_property_static("showHidden", [](py::object) { return ToulBar2::showHidden; }, [](py::object, bool showHidden) { ToulBar2::showHidden = showHidden; })
        .def_property_static("allSolutions", [](py::object) { return ToulBar2::allSolutions; }, [](py::object, Long allSolutions) { ToulBar2::allSolutions = allSolutions; })
        .def_property_static("dumpWCSP", [](py::object) { return ToulBar2::dumpWCSP; }, [](py::object, int dumpWCSP) { ToulBar2::dumpWCSP = dumpWCSP; })
        .def_property_static("dumpOriginalAfterPreprocessing", [](py::object) { return ToulBar2::dumpOriginalAfterPreprocessing; }, [](py::object, bool dumpOriginalAfterPreprocessing) { ToulBar2::dumpOriginalAfterPreprocessing = dumpOriginalAfterPreprocessing; })
        .def_property_static("approximateCountingBTD", [](py::object) { return ToulBar2::approximateCountingBTD; }, [](py::object, bool approximateCountingBTD) { ToulBar2::approximateCountingBTD = approximateCountingBTD; })
        .def_property_static("binaryBranching", [](py::object) { return ToulBar2::binaryBranching; }, [](py::object, bool binaryBranching) { ToulBar2::binaryBranching = binaryBranching; })
        .def_property_static("dichotomicBranching", [](py::object) { return ToulBar2::dichotomicBranching; }, [](py::object, int dichotomicBranching) { ToulBar2::dichotomicBranching = dichotomicBranching; })
        .def_property_static("dichotomicBranchingSize", [](py::object) { return ToulBar2::dichotomicBranchingSize; }, [](py::object, unsigned int dichotomicBranchingSize) { ToulBar2::dichotomicBranchingSize = dichotomicBranchingSize; })
        .def_property_static("sortDomains", [](py::object) { return ToulBar2::sortDomains; }, [](py::object, bool sortDomains) { ToulBar2::sortDomains = sortDomains; })
        .def_property_static("solutionBasedPhaseSaving", [](py::object) { return ToulBar2::solutionBasedPhaseSaving; }, [](py::object, bool solutionBasedPhaseSaving) { ToulBar2::solutionBasedPhaseSaving = solutionBasedPhaseSaving; })
        .def_property_static("bisupport", [](py::object) { return ToulBar2::bisupport; }, [](py::object, Double bisupport) { ToulBar2::bisupport = bisupport; })
        .def_property_static("elimDegree", [](py::object) { return ToulBar2::elimDegree; }, [](py::object, int elimDegree) { ToulBar2::elimDegree = elimDegree; })
        .def_property_static("elimDegree_preprocessing", [](py::object) { return ToulBar2::elimDegree_preprocessing; }, [](py::object, int elimDegree_preprocessing) { ToulBar2::elimDegree_preprocessing = elimDegree_preprocessing; })
        .def_property_static("elimSpaceMaxMB", [](py::object) { return ToulBar2::elimSpaceMaxMB; }, [](py::object, int elimSpaceMaxMB) { ToulBar2::elimSpaceMaxMB = elimSpaceMaxMB; })
        .def_property_static("minsumDiffusion", [](py::object) { return ToulBar2::minsumDiffusion; }, [](py::object, int minsumDiffusion) { ToulBar2::minsumDiffusion = minsumDiffusion; })
        .def_property_static("preprocessTernaryRPC", [](py::object) { return ToulBar2::preprocessTernaryRPC; }, [](py::object, int preprocessTernaryRPC) { ToulBar2::preprocessTernaryRPC = preprocessTernaryRPC; })
        .def_property_static("hve", [](py::object) { return ToulBar2::hve; }, [](py::object, int hve) { ToulBar2::hve = hve; })
        .def_property_static("pwc", [](py::object) { return ToulBar2::pwc; }, [](py::object, int pwc) { ToulBar2::pwc = pwc; })
        .def_property_static("pwcMinimalDualGraph", [](py::object) { return ToulBar2::pwcMinimalDualGraph; }, [](py::object, bool pwcMinimalDualGraph) { ToulBar2::pwcMinimalDualGraph = pwcMinimalDualGraph; })
        .def_property_static("preprocessFunctional", [](py::object) { return ToulBar2::preprocessFunctional; }, [](py::object, int preprocessFunctional) { ToulBar2::preprocessFunctional = preprocessFunctional; })
        .def_property_static("costfuncSeparate", [](py::object) { return ToulBar2::costfuncSeparate; }, [](py::object, bool costfuncSeparate) { ToulBar2::costfuncSeparate = costfuncSeparate; })
        .def_property_static("preprocessNary", [](py::object) { return ToulBar2::preprocessNary; }, [](py::object, int preprocessNary) { ToulBar2::preprocessNary = preprocessNary; })
        .def_property_static("QueueComplexity", [](py::object) { return ToulBar2::QueueComplexity; }, [](py::object, bool QueueComplexity) { ToulBar2::QueueComplexity = QueueComplexity; })
        .def_property_static("Static_variable_ordering", [](py::object) { return ToulBar2::Static_variable_ordering; }, [](py::object, bool Static_variable_ordering) { ToulBar2::Static_variable_ordering = Static_variable_ordering; })
        .def_property_static("lastConflict", [](py::object) { return ToulBar2::lastConflict; }, [](py::object, bool lastConflict) { ToulBar2::lastConflict = lastConflict; })
        .def_property_static("weightedDegree", [](py::object) { return ToulBar2::weightedDegree; }, [](py::object, int weightedDegree) { ToulBar2::weightedDegree = weightedDegree; })
        .def_property_static("weightedTightness", [](py::object) { return ToulBar2::weightedTightness; }, [](py::object, int weightedTightness) { ToulBar2::weightedTightness = weightedTightness; })
        .def_property_static("constrOrdering", [](py::object) { return ToulBar2::constrOrdering; }, [](py::object, int constrOrdering) { ToulBar2::constrOrdering = constrOrdering; })
        .def_property_static("MSTDAC", [](py::object) { return ToulBar2::MSTDAC; }, [](py::object, bool MSTDAC) { ToulBar2::MSTDAC = MSTDAC; })
        .def_property_static("DEE", [](py::object) { return ToulBar2::DEE; }, [](py::object, int DEE) { ToulBar2::DEE = DEE; })
        .def_property_static("nbDecisionVars", [](py::object) { return ToulBar2::nbDecisionVars; }, [](py::object, int nbDecisionVars) { ToulBar2::nbDecisionVars = nbDecisionVars; })
        .def_property_static("lds", [](py::object) { return ToulBar2::lds; }, [](py::object, int lds) { ToulBar2::lds = lds; })
        .def_property_static("limited", [](py::object) { return ToulBar2::limited; }, [](py::object, bool limited) { ToulBar2::limited = limited; })
        .def_property_static("restart", [](py::object) { return ToulBar2::restart; }, [](py::object, Long restart) { ToulBar2::restart = restart; })
        .def_property_static("backtrackLimit", [](py::object) { return ToulBar2::backtrackLimit; }, [](py::object, Long backtrackLimit) { ToulBar2::backtrackLimit = backtrackLimit; })
        .def_property_static("cfn", [](py::object) { return ToulBar2::cfn; }, [](py::object, bool cfn) { ToulBar2::cfn = cfn; })
        .def_property_static("gz", [](py::object) { return ToulBar2::gz; }, [](py::object, bool gz) { ToulBar2::gz = gz; })
        .def_property_static("bz2", [](py::object) { return ToulBar2::bz2; }, [](py::object, bool bz2) { ToulBar2::bz2 = bz2; })
        .def_property_static("xz", [](py::object) { return ToulBar2::xz; }, [](py::object, bool xz) { ToulBar2::xz = xz; })
        .def_property_static("bayesian", [](py::object) { return ToulBar2::bayesian; }, [](py::object, bool bayesian) { ToulBar2::bayesian = bayesian; })
        .def_property_static("uai", [](py::object) { return ToulBar2::uai; }, [](py::object, int uai) { ToulBar2::uai = uai; })
        .def_property_static("resolution", [](py::object) { return ToulBar2::resolution; }, [](py::object, int resolution) { ToulBar2::resolution = resolution; })
        .def_property_static("resolution_Update", [](py::object) { return ToulBar2::resolution_Update; }, [](py::object, bool resolution_Update) { ToulBar2::resolution_Update = resolution_Update; })
        .def_property_static("errorg", [](py::object) { return ToulBar2::errorg; }, [](py::object, TProb errorg) { ToulBar2::errorg = errorg; })
        .def_property_static("NormFactor", [](py::object) { return ToulBar2::NormFactor; }, [](py::object, TLogProb NormFactor) { ToulBar2::NormFactor = NormFactor; })
        .def_property_static("vac", [](py::object) { return ToulBar2::vac; }, [](py::object, int vac) { ToulBar2::vac = vac; })
        .def_property_static("costThresholdS", [](py::object) { return ToulBar2::costThresholdS; }, [](py::object, string costThresholdS) { ToulBar2::costThresholdS = costThresholdS; })
        .def_property_static("costThresholdPreS", [](py::object) { return ToulBar2::costThresholdPreS; }, [](py::object, string costThresholdPreS) { ToulBar2::costThresholdPreS = costThresholdPreS; })
        .def_property_static("costThreshold", [](py::object) { return ToulBar2::costThreshold; }, [](py::object, Cost costThreshold) { ToulBar2::costThreshold = costThreshold; })
        .def_property_static("costThresholdPre", [](py::object) { return ToulBar2::costThresholdPre; }, [](py::object, Cost costThresholdPre) { ToulBar2::costThresholdPre = costThresholdPre; })
        .def_property_static("FullEAC", [](py::object) { return ToulBar2::FullEAC; }, [](py::object, bool FullEAC) { ToulBar2::FullEAC = FullEAC; })
        .def_property_static("VACthreshold", [](py::object) { return ToulBar2::VACthreshold; }, [](py::object, bool VACthreshold) { ToulBar2::VACthreshold = VACthreshold; })
        .def_property_static("useRASPS", [](py::object) { return ToulBar2::useRASPS; }, [](py::object, int useRASPS) { ToulBar2::useRASPS = useRASPS; })
        .def_property_static("RASPSreset", [](py::object) { return ToulBar2::RASPSreset; }, [](py::object, bool RASPSreset) { ToulBar2::RASPSreset = RASPSreset; })
        .def_property_static("RASPSangle", [](py::object) { return ToulBar2::RASPSangle; }, [](py::object, int RASPSangle) { ToulBar2::RASPSangle = RASPSangle; })
        .def_property_static("RASPSnbBacktracks", [](py::object) { return ToulBar2::RASPSnbBacktracks; }, [](py::object, Long RASPSnbBacktracks) { ToulBar2::RASPSnbBacktracks = RASPSnbBacktracks; })
        .def_property_static("ReducedCostsFiltering", [](py::object) { return ToulBar2::ReducedCostsFiltering; }, [](py::object, int ReducedCostsFiltering) { ToulBar2::ReducedCostsFiltering = ReducedCostsFiltering; })
        .def_property_static("trwsAccuracy", [](py::object) { return ToulBar2::trwsAccuracy; }, [](py::object, Double trwsAccuracy) { ToulBar2::trwsAccuracy = trwsAccuracy; })
        .def_property_static("trwsOrder", [](py::object) { return ToulBar2::trwsOrder; }, [](py::object, bool trwsOrder) { ToulBar2::trwsOrder = trwsOrder; })
        .def_property_static("trwsNIter", [](py::object) { return ToulBar2::trwsNIter; }, [](py::object, unsigned int trwsNIter) { ToulBar2::trwsNIter = trwsNIter; })
        .def_property_static("trwsNIterNoChange", [](py::object) { return ToulBar2::trwsNIterNoChange; }, [](py::object, unsigned int trwsNIterNoChange) { ToulBar2::trwsNIterNoChange = trwsNIterNoChange; })
        .def_property_static("trwsNIterComputeUb", [](py::object) { return ToulBar2::trwsNIterComputeUb; }, [](py::object, unsigned int trwsNIterComputeUb) { ToulBar2::trwsNIterComputeUb = trwsNIterComputeUb; })
        .def_property_static("decimalPoint", [](py::object) { return ToulBar2::decimalPoint; }, [](py::object, unsigned int decimalPoint) { ToulBar2::decimalPoint = decimalPoint; })
        .def_property_static("absgapstr", [](py::object) { return ToulBar2::deltaUbS; }, [](py::object, string deltaUbS) { ToulBar2::deltaUbS = deltaUbS; })
        .def_property_static("deltaUb", [](py::object) { return ToulBar2::deltaUb; }, [](py::object, Cost deltaUb) { ToulBar2::deltaUb = deltaUb; })
        .def_property_static("absgap", [](py::object) { return ToulBar2::deltaUbAbsolute; }, [](py::object, Cost deltaUbAbsolute) { ToulBar2::deltaUbAbsolute = deltaUbAbsolute; })
        .def_property_static("relgap", [](py::object) { return ToulBar2::deltaUbRelativeGap; }, [](py::object, Cost deltaUbRelativeGap) { ToulBar2::deltaUbRelativeGap = deltaUbRelativeGap; })
        .def_property_static("singletonConsistency", [](py::object) { return ToulBar2::singletonConsistency; }, [](py::object, int singletonConsistency) { ToulBar2::singletonConsistency = singletonConsistency; })
        .def_property_static("singletonAccuracy", [](py::object) { return ToulBar2::singletonAccuracy; }, [](py::object, Double singletonAccuracy) { ToulBar2::singletonAccuracy = singletonAccuracy; })
        .def_property_static("GilmoreLawler", [](py::object) { return ToulBar2::GilmoreLawler; }, [](py::object, int GilmoreLawler) { ToulBar2::GilmoreLawler = GilmoreLawler; })
        .def_property_static("vacValueHeuristic", [](py::object) { return ToulBar2::vacValueHeuristic; }, [](py::object, int vacValueHeuristic) { ToulBar2::vacValueHeuristic = vacValueHeuristic; })
        .def_property_static("LcLevel", [](py::object) { return ToulBar2::LcLevel; }, [](py::object, LcLevelType LcLevel) { ToulBar2::LcLevel = LcLevel; })
        .def_property_static("maxEACIter", [](py::object) { return ToulBar2::maxEACIter; }, [](py::object, int maxEACIter) { ToulBar2::maxEACIter = maxEACIter; })
        .def_property_static("wcnf", [](py::object) { return ToulBar2::wcnf; }, [](py::object, bool wcnf) { ToulBar2::wcnf = wcnf; })
        .def_property_static("qpbo", [](py::object) { return ToulBar2::qpbo; }, [](py::object, bool qpbo) { ToulBar2::qpbo = qpbo; })
        .def_property_static("qpboQuadraticCoefMultiplier", [](py::object) { return ToulBar2::qpboQuadraticCoefMultiplier; }, [](py::object, Double qpboQuadraticCoefMultiplier) { ToulBar2::qpboQuadraticCoefMultiplier = qpboQuadraticCoefMultiplier; })
        .def_property_static("opb", [](py::object) { return ToulBar2::opb; }, [](py::object, bool opb) { ToulBar2::opb = opb; })
        .def_property_static("cardinality", [](py::object) { return ToulBar2::cardinality; }, [](py::object, bool cardinality) { ToulBar2::cardinality = cardinality; })
        .def_property_static("lp", [](py::object) { return ToulBar2::lp; }, [](py::object, bool lp) { ToulBar2::lp = lp; })
#ifdef BOOST
        .def_property_static("addAMOConstraints", [](py::object) { return ToulBar2::addAMOConstraints; }, [](py::object, int addAMOConstraints) { ToulBar2::addAMOConstraints = addAMOConstraints; })
#endif
        .def_property_static("knapsackDP", [](py::object) { return ToulBar2::knapsackDP; }, [](py::object, int knapsackDP) { ToulBar2::knapsackDP = knapsackDP; })
        .def_property_static("VAClin", [](py::object) { return ToulBar2::VAClin; }, [](py::object, bool VAClin) { ToulBar2::VAClin = VAClin; })
        .def_property_static("divNbSol", [](py::object) { return ToulBar2::divNbSol; }, [](py::object, unsigned int divNbSol) { ToulBar2::divNbSol = divNbSol; })
        .def_property_static("divBound", [](py::object) { return ToulBar2::divBound; }, [](py::object, unsigned int divBound) { ToulBar2::divBound = divBound; })
        .def_property_static("divWidth", [](py::object) { return ToulBar2::divWidth; }, [](py::object, unsigned int divWidth) { ToulBar2::divWidth = divWidth; })
        .def_property_static("divMethod", [](py::object) { return ToulBar2::divMethod; }, [](py::object, unsigned int divMethod) { ToulBar2::divMethod = divMethod; })
        .def_property_static("divRelax", [](py::object) { return ToulBar2::divRelax; }, [](py::object, unsigned int divRelax) { ToulBar2::divRelax = divRelax; })

        .def_property_static("btdMode", [](py::object) { return ToulBar2::btdMode; }, [](py::object, int btdMode) { ToulBar2::btdMode = btdMode; })
        .def_property_static("btdSubTree", [](py::object) { return ToulBar2::btdSubTree; }, [](py::object, int btdSubTree) { ToulBar2::btdSubTree = btdSubTree; })
        .def_property_static("btdRootCluster", [](py::object) { return ToulBar2::btdRootCluster; }, [](py::object, int btdRootCluster) { ToulBar2::btdRootCluster = btdRootCluster; })
        .def_property_static("rootHeuristic", [](py::object) { return ToulBar2::rootHeuristic; }, [](py::object, int rootHeuristic) { ToulBar2::rootHeuristic = rootHeuristic; })
        .def_property_static("reduceHeight", [](py::object) { return ToulBar2::reduceHeight; }, [](py::object, bool reduceHeight) { ToulBar2::reduceHeight = reduceHeight; })
        // .def_property_static("maxsateval", [](py::object) { return ToulBar2::maxsateval; }, [](py::object, bool maxsateval) { ToulBar2::maxsateval = maxsateval; })
        .def_property_static("xmlflag", [](py::object) { return ToulBar2::xmlflag; }, [](py::object, bool xmlflag) { ToulBar2::xmlflag = xmlflag; })
        .def_property_static("xmlcop", [](py::object) { return ToulBar2::xmlcop; }, [](py::object, bool xmlcop) { ToulBar2::xmlcop = xmlcop; })
        .def_property_static("markov_log", [](py::object) { return ToulBar2::markov_log; }, [](py::object, TLogProb markov_log) { ToulBar2::markov_log = markov_log; })
        .def_property_static("evidence_file", [](py::object) { return ToulBar2::evidence_file; }, [](py::object, string evidence_file) { ToulBar2::evidence_file = evidence_file; })
        .def_property_static("solution_uai_filename", [](py::object) { return ToulBar2::solution_uai_filename; }, [](py::object, string solution_uai_filename) { ToulBar2::solution_uai_filename = solution_uai_filename; })
        .def_property_static("problemsaved_filename", [](py::object) { return ToulBar2::problemsaved_filename; }, [](py::object, string problemsaved_filename) { ToulBar2::problemsaved_filename = problemsaved_filename; })
        .def_property_static("isZ", [](py::object) { return ToulBar2::isZ; }, [](py::object, bool isZ) { ToulBar2::isZ = isZ; })
        .def_property_static("logZ", [](py::object) { return ToulBar2::logZ; }, [](py::object, TLogProb logZ) { ToulBar2::logZ = logZ; })
        .def_property_static("logU", [](py::object) { return ToulBar2::logU; }, [](py::object, TLogProb logU) { ToulBar2::logU = logU; })
        .def_property_static("logepsilon", [](py::object) { return ToulBar2::logepsilon; }, [](py::object, TLogProb logepsilon) { ToulBar2::logepsilon = logepsilon; })
        .def_property_static("epsilon", [](py::object) { return ToulBar2::epsilon; }, [](py::object, Double epsilon) { ToulBar2::epsilon = epsilon; })
        .def_property_static("uaieval", [](py::object) { return ToulBar2::uaieval; }, [](py::object, bool uaieval) { ToulBar2::uaieval = uaieval; })
        .def_property_static("stdin_format", [](py::object) { return ToulBar2::stdin_format; }, [](py::object, string stdin_format) { ToulBar2::stdin_format = stdin_format; })
        .def_property_static("startCpuTime", [](py::object) { return ToulBar2::startCpuTime; }, [](py::object, double startCpuTime) { ToulBar2::startCpuTime = startCpuTime; })
        .def_property_static("startRealTime", [](py::object) { return ToulBar2::startRealTime; }, [](py::object, double startRealTime) { ToulBar2::startRealTime = startRealTime; })
        .def_property_static("startRealTimeAfterPreProcessing", [](py::object) { return ToulBar2::startRealTimeAfterPreProcessing; }, [](py::object, double startRealTimeAfterPreProcessing) { ToulBar2::startRealTimeAfterPreProcessing = startRealTimeAfterPreProcessing; })
        .def_property_static("splitClusterMaxSize", [](py::object) { return ToulBar2::splitClusterMaxSize; }, [](py::object, int splitClusterMaxSize) { ToulBar2::splitClusterMaxSize = splitClusterMaxSize; })
        .def_property_static("boostingBTD", [](py::object) { return ToulBar2::boostingBTD; }, [](py::object, double boostingBTD) { ToulBar2::boostingBTD = boostingBTD; })
        .def_property_static("maxSeparatorSize", [](py::object) { return ToulBar2::maxSeparatorSize; }, [](py::object, int maxSeparatorSize) { ToulBar2::maxSeparatorSize = maxSeparatorSize; })
        .def_property_static("minProperVarSize", [](py::object) { return ToulBar2::minProperVarSize; }, [](py::object, int minProperVarSize) { ToulBar2::minProperVarSize = minProperVarSize; })
        .def_property_static("Berge_Dec", [](py::object) { return ToulBar2::Berge_Dec; }, [](py::object, bool Berge_Dec) { ToulBar2::Berge_Dec = Berge_Dec; })
        .def_property_static("learning", [](py::object) { return ToulBar2::learning; }, [](py::object, bool learning) { ToulBar2::learning = learning; })
        // .def_property_static("interrupted", [](py::object) { return ToulBar2::interrupted; }, [](py::object, std::atomic<bool> interrupted) { ToulBar2::interrupted = interrupted; }) // pybind11 not compatible with type atomic<bool>?
        .def_property_static("seed", [](py::object) { return ToulBar2::seed; }, [](py::object, int seed) { ToulBar2::seed = seed; })
        .def_property_static("sigma", [](py::object) { return ToulBar2::sigma; }, [](py::object, Double sigma) { ToulBar2::sigma = sigma; })
        .def_property_static("incop_cmd", [](py::object) { return ToulBar2::incop_cmd; }, [](py::object, string incop_cmd) { ToulBar2::incop_cmd = incop_cmd; })
        .def_property_static("pils_cmd", [](py::object) { return ToulBar2::pils_cmd; }, [](py::object, string pils_cmd) { ToulBar2::pils_cmd = pils_cmd; })
        .def_property_static("lrBCD_cmd", [](py::object) { return ToulBar2::lrBCD_cmd; }, [](py::object, string lrBCD_cmd) { ToulBar2::lrBCD_cmd = lrBCD_cmd; })
        .def_property_static("searchMethod", [](py::object) { return static_cast<int>(ToulBar2::searchMethod); }, [](py::object, int searchMethod) { ToulBar2::searchMethod = static_cast<SearchMethod>(searchMethod); })
        .def_property_static("clusterFile", [](py::object) { return ToulBar2::clusterFile; }, [](py::object, string clusterFile) { ToulBar2::clusterFile = clusterFile; })
        .def_property_static("vnsInitSol", [](py::object) { return static_cast<int>(ToulBar2::vnsInitSol); }, [](py::object, int vnsInitSol) { ToulBar2::vnsInitSol = static_cast<VNSSolutionInitMethod>(vnsInitSol); })
        .def_property_static("vnsLDSmin", [](py::object) { return ToulBar2::vnsLDSmin; }, [](py::object, int vnsLDSmin) { ToulBar2::vnsLDSmin = vnsLDSmin; })
        .def_property_static("vnsLDSmax", [](py::object) { return ToulBar2::vnsLDSmax; }, [](py::object, int vnsLDSmax) { ToulBar2::vnsLDSmax = vnsLDSmax; })
        .def_property_static("vnsLDSinc", [](py::object) { return static_cast<int>(ToulBar2::vnsLDSinc); }, [](py::object, int vnsLDSinc) { ToulBar2::vnsLDSinc = static_cast<VNSInc>(vnsLDSinc); })
        .def_property_static("vnsKmin", [](py::object) { return ToulBar2::vnsKmin; }, [](py::object, int vnsKmin) { ToulBar2::vnsKmin = vnsKmin; })
        .def_property_static("vnsKmax", [](py::object) { return ToulBar2::vnsKmax; }, [](py::object, int vnsKmax) { ToulBar2::vnsKmax = vnsKmax; })
        .def_property_static("vnsKinc", [](py::object) { return static_cast<int>(ToulBar2::vnsKinc); }, [](py::object, int vnsKinc) { ToulBar2::vnsKinc = static_cast<VNSInc>(vnsKinc); })
        .def_property_static("vnsLDScur", [](py::object) { return ToulBar2::vnsLDScur; }, [](py::object, int vnsLDScur) { ToulBar2::vnsLDScur = vnsLDScur; })
        .def_property_static("vnsKcur", [](py::object) { return ToulBar2::vnsKcur; }, [](py::object, int vnsKcur) { ToulBar2::vnsKcur = vnsKcur; })
        .def_property_static("vnsNeighborVarHeur", [](py::object) { return ToulBar2::vnsNeighborVarHeur; }, [](py::object, VNSVariableHeuristic vnsNeighborVarHeur) { ToulBar2::vnsNeighborVarHeur = vnsNeighborVarHeur; })
        .def_property_static("vnsNeighborChange", [](py::object) { return ToulBar2::vnsNeighborChange; }, [](py::object, bool vnsNeighborChange) { ToulBar2::vnsNeighborChange = vnsNeighborChange; })
        .def_property_static("vnsNeighborSizeSync", [](py::object) { return ToulBar2::vnsNeighborSizeSync; }, [](py::object, bool vnsNeighborSizeSync) { ToulBar2::vnsNeighborSizeSync = vnsNeighborSizeSync; })
        .def_property_static("vnsParallelLimit", [](py::object) { return ToulBar2::vnsParallelLimit; }, [](py::object, bool vnsParallelLimit) { ToulBar2::vnsParallelLimit = vnsParallelLimit; })
        .def_property_static("vnsParallelSync", [](py::object) { return ToulBar2::vnsParallelSync; }, [](py::object, bool vnsParallelSync) { ToulBar2::vnsParallelSync = vnsParallelSync; })
        .def_property_static("vnsOptimumS", [](py::object) { return ToulBar2::vnsOptimumS; }, [](py::object, string bestsol) { ToulBar2::vnsOptimumS = bestsol; ToulBar2::newsolution = newsolution; })
        .def_property_static("vnsOptimum", [](py::object) { return ToulBar2::vnsOptimum; }, [](py::object, const Cost cost) { ToulBar2::vnsOptimum = cost; ToulBar2::newsolution = newsolution; })
        .def_property_static("parallel", [](py::object) { return ToulBar2::parallel; }, [](py::object, bool parallel) { ToulBar2::parallel = parallel; })
        .def_property_static("hbfs", [](py::object) { return ToulBar2::hbfs; }, [](py::object, Long hbfs) { ToulBar2::hbfs = hbfs; })
        .def_property_static("hbfsGlobalLimit", [](py::object) { return ToulBar2::hbfsGlobalLimit; }, [](py::object, Long hbfsGlobalLimit) { ToulBar2::hbfsGlobalLimit = hbfsGlobalLimit; })
        .def_property_static("hbfsAlpha", [](py::object) { return ToulBar2::hbfsAlpha; }, [](py::object, Long hbfsAlpha) { ToulBar2::hbfsAlpha = hbfsAlpha; })
        .def_property_static("hbfsBeta", [](py::object) { return ToulBar2::hbfsBeta; }, [](py::object, Long hbfsBeta) { ToulBar2::hbfsBeta = hbfsBeta; })
        .def_property_static("hbfsCPLimit", [](py::object) { return ToulBar2::hbfsCPLimit; }, [](py::object, ptrdiff_t hbfsCPLimit) { ToulBar2::hbfsCPLimit = hbfsCPLimit; })
        .def_property_static("hbfsOpenNodeLimit", [](py::object) { return ToulBar2::hbfsOpenNodeLimit; }, [](py::object, ptrdiff_t hbfsOpenNodeLimit) { ToulBar2::hbfsOpenNodeLimit = hbfsOpenNodeLimit; })
        .def_property_static("sortBFS", [](py::object) { return ToulBar2::sortBFS; }, [](py::object, Long sortBFS) { ToulBar2::sortBFS = sortBFS; })
        .def_property_static("eps", [](py::object) { return ToulBar2::eps; }, [](py::object, Long eps) { ToulBar2::eps = eps; })
#ifdef OPENMPI
        .def_property_static("burst", [](py::object) { return ToulBar2::burst; }, [](py::object, bool burst) { ToulBar2::burst = burst; })
#endif
        .def_property_static("epsFilename", [](py::object) { return ToulBar2::epsFilename; }, [](py::object, string epsFilename) { ToulBar2::epsFilename = epsFilename; })
        .def_property_static("verifyOpt", [](py::object) { return ToulBar2::verifyOpt; }, [](py::object, bool verifyOpt) { ToulBar2::verifyOpt = verifyOpt; })
        .def_property_static("verifiedOptimum", [](py::object) { return ToulBar2::verifiedOptimum; }, [](py::object, Cost verifiedOptimum) { ToulBar2::verifiedOptimum = verifiedOptimum; })
        .def_property_static("bilevel", [](py::object) { return ToulBar2::bilevel; }, [](py::object, int bilevel) { ToulBar2::bilevel = bilevel; })
        .def_property_static("decimalPointBLP", [](py::object) { return ToulBar2::decimalPointBLP; }, [](py::object, vector<unsigned int> decimalPointBLP) { ToulBar2::decimalPointBLP = decimalPointBLP; })
        .def_property_static("costMultiplierBLP", [](py::object) { return ToulBar2::costMultiplierBLP; }, [](py::object, vector<Double> costMultiplierBLP) { ToulBar2::costMultiplierBLP = costMultiplierBLP; })
        .def_property_static("negCostBLP", [](py::object) { return ToulBar2::negCostBLP; }, [](py::object, vector<Cost> negCostBLP) { ToulBar2::negCostBLP = negCostBLP; })
        .def_property_static("initialLbBLP", [](py::object) { return ToulBar2::initialLbBLP; }, [](py::object, vector<Cost> initialLbBLP) { ToulBar2::initialLbBLP = initialLbBLP; })
        .def_property_static("initialUbBLP", [](py::object) { return ToulBar2::initialUbBLP; }, [](py::object, vector<Cost> initialUbBLP) { ToulBar2::initialUbBLP = initialUbBLP; });

    m.def("check", &tb2checkOptions); // should be called after setting the options (and before reading a problem)

    py::class_<Store, std::unique_ptr<Store, py::nodelete>>(m, "store")
        .def("getDepth", &Store::getDepth)
        .def("store", &Store::store)
        .def("restore", static_cast<void (*)(int)>(&Store::restore));

    py::class_<WeightedObjInt>(m, "WeightedObjInt")
        .def(py::init<int, Cost>())
        .def_readwrite("val", &WeightedObjInt::val)
        .def_readwrite("weight", &WeightedObjInt::weight);

    py::class_<DFATransition>(m, "DFATransition")
        .def(py::init<int, Value, int, Cost>())

        .def_readwrite("start", &DFATransition::start)
        .def_readwrite("end", &DFATransition::end)
        .def_readwrite("symbol", &DFATransition::symbol)
        .def_readwrite("weight", &DFATransition::weight);

    py::class_<BoundedObjValue>(m, "BoundedObjValue")
        .def(py::init<Value, unsigned int, unsigned int>())

        .def_readwrite("val", &BoundedObjValue::val)
        .def_readwrite("upper", &BoundedObjValue::upper)
        .def_readwrite("lower", &BoundedObjValue::lower);

    py::class_<MultiCFN>(m, "MultiCFN")
        .def(py::init())
        .def(
            "push_back", [](MultiCFN& multicfn, WeightedCSP* wcsp, Double weight) { multicfn.push_back(dynamic_cast<WCSP*>(wcsp), weight); }, py::arg("wcsp"), py::arg("weight") = 1.)
        .def("setWeight", &MultiCFN::setWeight)
        .def("nbNetworks", &MultiCFN::nbNetworks)
        .def("nbVariables", &MultiCFN::nbVariables)
        .def("nbCostFunctions", &MultiCFN::nbCostFunctions)
        .def("getNetworkName", &MultiCFN::getNetworkName)
        .def("getVariableIndex", &MultiCFN::getVariableIndex)
        .def("nbValues", &MultiCFN::nbValues)
        .def("getScope", &MultiCFN::getScope)
        .def("getType", &MultiCFN::getType)
        .def("getCost", &MultiCFN::getCost)
        .def("print", [](MultiCFN& multicfn) { multicfn.print(cout); })
        .def(
            "makeWeightedCSP", [](MultiCFN& multicfn, const set<unsigned int>& vars, const vector<set<unsigned int>>& scopes, const vector<unsigned int>& constrs) { return multicfn.makeWeightedCSP(vars, scopes, constrs); }, py::arg("vars") = py::set(), py::arg("scopes") = py::list(), py::arg("constrs") = py::list())
        .def(
            "makeWeightedCSP", [](MultiCFN& multicfn, WeightedCSP* wcsp, const set<unsigned int>& vars, const vector<set<unsigned int>>& scopes, const vector<unsigned int>& constrs) { multicfn.makeWeightedCSP(wcsp, vars, scopes, constrs); }, py::arg("wcsp"), py::arg("vars") = py::set(), py::arg("scopes") = py::list(), py::arg("constrs") = py::list())
        .def("getSolution", &MultiCFN::getSolution)
        .def("getSolutionValues", &MultiCFN::getSolutionValues)
        .def("computeSolutionValues", &MultiCFN::computeSolutionValues)
        .def("convertToSolution", &MultiCFN::convertToSolution)
        .def("setNoiseActivation", &MultiCFN::setNoiseActivation)
        .def("setNoiseLevel", &MultiCFN::setNoiseLevel);

    py::class_<Bicriteria> bcrit(m, "Bicriteria");

    bcrit.def(
             "computeSupportedPoints", [](MultiCFN* multicfn, int first_cfn_index, int second_cfn_index, py::tuple optim_dir, Double delta) { Bicriteria::computeSupportedPoints(multicfn, first_cfn_index, second_cfn_index, std::make_pair(optim_dir[0].cast<Bicriteria::OptimDir>(), optim_dir[1].cast<Bicriteria::OptimDir>()), delta); }, py::arg("first_cfn_index"), py::arg("second_cfn_index"), py::arg("optim_dir"), py::arg("delta") = Bicriteria::Delta)
        .def(
            "computeAdditionalSolutions", [](MultiCFN* multicfn, py::tuple optim_dir, unsigned int solIndex, unsigned int nbLimit, Double pct) { Bicriteria::computeAdditionalSolutions(multicfn, std::make_pair(optim_dir[0].cast<Bicriteria::OptimDir>(), optim_dir[1].cast<Bicriteria::OptimDir>()), solIndex, nbLimit, pct); }, py::arg("optim_dir"), py::arg("solIndex"), py::arg("nbLimit") = 100, py::arg("pct") = 1.)
        .def(
            "computeNonSupported", [](MultiCFN* multicfn, py::tuple optim_dir, unsigned int nbLimit) { Bicriteria::computeNonSupported(multicfn, std::make_pair(optim_dir[0].cast<Bicriteria::OptimDir>(), optim_dir[1].cast<Bicriteria::OptimDir>()), nbLimit); }, py::arg("optim_dir"), py::arg("nbLimit") = 100)
        .def("getSolutions", &Bicriteria::getSolutions)
        .def("getPoints", &Bicriteria::getPoints)
        .def("getWeights", &Bicriteria::getWeights)
        .def("getLowerBounds", &Bicriteria::getLowerBounds)
        .def("setShowSolutions", &Bicriteria::setShowSolutions)
        .def("setVAC", &Bicriteria::setVAC)
        .def("setSeed", &Bicriteria::setSeed)
        .def("setVerbose", &Bicriteria::setVerbose)
        .def("setSolutionTimeout", &Bicriteria::setSolutionTimeout)
        .def("setGlobalTimeout", &Bicriteria::setGlobalTimeout)
        .def("setMaxSolutionCount", &Bicriteria::setMaxSolutionCount);

    py::enum_<Bicriteria::OptimDir>(bcrit, "OptimDir")
        .value("Min", Bicriteria::OptimDir::Optim_Min)
        .value("Max", Bicriteria::OptimDir::Optim_Max)
        .export_values();

    py::class_<WeightedCSP>(m, "WCSP")
        .def(py::init([](Cost ub) { return WeightedCSP::makeWeightedCSP(ub); })) // create this object to interface with the multicfn class
        .def("getIndex", &WeightedCSP::getIndex)
        .def("getName", (string(WeightedCSP::*)() const) & WeightedCSP::getName)
        .def("setName", &WeightedCSP::setName)
        .def("getLb", &WeightedCSP::getLb)
        .def("getUb", &WeightedCSP::getUb)
        .def("getDPrimalBound", &WeightedCSP::getDPrimalBound)
        .def("getDDualBound", &WeightedCSP::getDDualBound)
        .def("getDLb", &WeightedCSP::getDLb)
        .def("getDUb", &WeightedCSP::getDUb)
        .def("updateDUb", &WeightedCSP::updateDUb)
        .def("setLb", &WeightedCSP::setLb)
        .def("setUb", &WeightedCSP::setUb)
        .def("updateUb", &WeightedCSP::updateUb)
        .def("enforceUb", &WeightedCSP::enforceUb)
        .def("increaseLb", &WeightedCSP::increaseLb)
        .def("decreaseLb", &WeightedCSP::decreaseLb)
        .def("getNegativeLb", &WeightedCSP::getNegativeLb)
        .def("finiteUb", &WeightedCSP::finiteUb)
        .def("setInfiniteCost", &WeightedCSP::setInfiniteCost)
        .def("isfinite", &WeightedCSP::isfinite)
        .def("enumerated", &WeightedCSP::enumerated)
        .def("getName", (string(WeightedCSP::*)(int) const) & WeightedCSP::getName)
        .def("getVarIndex", &WeightedCSP::getVarIndex)
        .def("getInf", &WeightedCSP::getInf)
        .def("getSup", &WeightedCSP::getSup)
        .def("getValue", &WeightedCSP::getValue)
        .def("getDomainSize", &WeightedCSP::getDomainSize)
        .def("getEnumDomain", (vector<Value>(WeightedCSP::*)(int varIndex)) & WeightedCSP::getEnumDomain)
        .def("getEnumDomainAndCost", (vector<pair<Value, Cost>>(WeightedCSP::*)(int varIndex)) & WeightedCSP::getEnumDomainAndCost)
        .def("getDomainInitSize", &WeightedCSP::getDomainInitSize)
        .def("toValue", &WeightedCSP::toValue)
        .def("toIndex", (unsigned int (WeightedCSP::*)(int varIndex, Value value)) & WeightedCSP::toIndex)
        .def("toIndex", (unsigned int (WeightedCSP::*)(int varIndex, const string& valueName)) & WeightedCSP::toIndex)
        .def("getDACOrder", &WeightedCSP::getDACOrder)
        .def("assigned", &WeightedCSP::assigned)
        .def("unassigned", &WeightedCSP::unassigned)
        .def("canbe", &WeightedCSP::canbe)
        .def("cannotbe", &WeightedCSP::cannotbe)
        .def("nextValue", &WeightedCSP::nextValue)
        .def("increase", &WeightedCSP::increase)
        .def("decrease", &WeightedCSP::decrease)
        .def("assign", &WeightedCSP::assign)
        .def("remove", &WeightedCSP::remove)
        .def("assignLS", (void(WeightedCSP::*)(vector<int> & varIndexes, vector<Value> & newValues, bool force)) & WeightedCSP::assignLS)
        .def("deconnect", &WeightedCSP::deconnect)
        .def("getUnaryCost", &WeightedCSP::getUnaryCost)
        .def("getMaxUnaryCost", &WeightedCSP::getMaxUnaryCost)
        .def("getMaxUnaryCostValue", &WeightedCSP::getMaxUnaryCostValue)
        .def("getSupport", &WeightedCSP::getSupport)
        .def("getBestValue", &WeightedCSP::getBestValue)
        .def("setBestValue", &WeightedCSP::setBestValue)
        //        .def("getIsPartOfOptimalSolution", &WeightedCSP::getIsPartOfOptimalSolution)
        //        .def("setIsPartOfOptimalSolution", &WeightedCSP::setIsPartOfOptimalSolution)
        .def("getDegree", &WeightedCSP::getDegree)
        .def("getTrueDegree", &WeightedCSP::getTrueDegree)
        .def("getWeightedDegree", &WeightedCSP::getWeightedDegree)
        .def("resetWeightedDegree", &WeightedCSP::resetWeightedDegree)
        .def("resetTightness", &WeightedCSP::resetTightness)
        .def("resetTightnessAndWeightedDegree", &WeightedCSP::resetTightnessAndWeightedDegree)
        .def("preprocessing", &WeightedCSP::preprocessing)
        .def("sortConstraints", &WeightedCSP::sortConstraints) // must be called after creating the model
        .def("getBergeDecElimOrder", &WeightedCSP::getBergeDecElimOrder)
        .def("setDACOrder", &WeightedCSP::setDACOrder)
        .def("whenContradiction", &WeightedCSP::whenContradiction)
        .def("deactivatePropagate", &WeightedCSP::deactivatePropagate)
        .def("isactivatePropagate", &WeightedCSP::isactivatePropagate)
        .def("reactivatePropagate", &WeightedCSP::reactivatePropagate)
        .def("propagate", &WeightedCSP::propagate, py::arg("fromscratch") = false)
        .def("verify", &WeightedCSP::verify)
        .def("propagateConstraint", &WeightedCSP::propagateConstraint)
        .def("numberOfVariables", &WeightedCSP::numberOfVariables)
        .def("numberOfUnassignedVariables", &WeightedCSP::numberOfUnassignedVariables)
        .def("numberOfConstraints", &WeightedCSP::numberOfConstraints)
        .def("numberOfConnectedConstraints", &WeightedCSP::numberOfConnectedConstraints)
        .def("numberOfConnectedBinaryConstraints", &WeightedCSP::numberOfConnectedBinaryConstraints)
        .def("numberOfConnectedKnapsackConstraints", &WeightedCSP::numberOfConnectedKnapsackConstraints)
        .def("medianDomainSize", &WeightedCSP::medianDomainSize)
        .def("medianDegree", &WeightedCSP::medianDegree)
        .def("medianArity", &WeightedCSP::medianArity)
        .def("getMaxDomainSize", &WeightedCSP::getMaxDomainSize)
        .def("getMaxCurrentDomainSize", &WeightedCSP::getMaxCurrentDomainSize)
        .def("getDomainSizeSum", &WeightedCSP::getDomainSizeSum)
        .def("cartProd", &WeightedCSP::cartProd)
        .def("getNbDEE", &WeightedCSP::getNbDEE)
        .def("makeEnumeratedVariable", (int(WeightedCSP::*)(string n, Value iinf, Value isup)) & WeightedCSP::makeEnumeratedVariable)
        .def("addValueName", &WeightedCSP::addValueName)
        .def("getValueName", &WeightedCSP::getValueName)
        .def("makeIntervalVariable", &WeightedCSP::makeIntervalVariable)
        .def("postNullaryConstraint", (void(WeightedCSP::*)(Double cost)) & WeightedCSP::postNullaryConstraint)
        .def(
            "postUnaryConstraint", [](WeightedCSP& s, int xIndex, vector<Double>& costs, bool incremental) {
                s.postUnaryConstraint(xIndex, costs, incremental);
            },
            py::arg("xIndex"), py::arg("costs"), py::arg("incremental") = false)
        .def("postBinaryConstraint", [](WeightedCSP& s, int xIndex, int yIndex, vector<Double>& costs, bool incremental) { return s.postBinaryConstraint(xIndex, yIndex, costs, incremental); }, py::arg("xIndex"), py::arg("yIndex"), py::arg("costs"), py::arg("incremental") = false)
        .def("postTernaryConstraint", [](WeightedCSP& s, int xIndex, int yIndex, int zIndex, vector<Double>& costs, bool incremental) { return s.postTernaryConstraint(xIndex, yIndex, zIndex, costs, incremental); }, py::arg("xIndex"), py::arg("yIndex"), py::arg("zIndex"), py::arg("costs"), py::arg("incremental") = false)

        // vectorized functions, numpy-compatible
        .def("postUnaryVecConstraints", postUnaryVecConstraints, py::arg("scopes"), py::arg("costs"), py::arg("incremental") = false)

        .def("postBinaryVecConstraints", postBinaryVecConstraints, py::arg("scopes"), py::arg("costs"), py::arg("incremental") = false)

        .def("postMultBinaryVecConstraints", postMultBinaryVecConstraints, py::arg("scopes"), py::arg("costs"), py::arg("incremental") = false)

        .def("postTernaryVecConstraints", postTernaryVecConstraints, py::arg("scopes"), py::arg("costs"), py::arg("incremental") = false)

        .def("postMultTernaryVecConstraints", postMultTernaryVecConstraints, py::arg("scopes"), py::arg("costs"), py::arg("incremental") = false)

        .def("makeEnumeratedVariableVec", makeEnumeratedVariableVec, py::arg("n"), py::arg("base_name"), py::arg("iinf"), py::arg("isup"))

        .def("postNaryConstraintBegin", [](WeightedCSP& s, vector<int> scope, Cost defval, Long nbtuples, bool forcenary) { return s.postNaryConstraintBegin(scope, defval, nbtuples, forcenary); }, py::arg("scope"), py::arg("defval"), py::arg("nbtuples"), py::arg("forcenary") = !NARY2CLAUSE)
        .def("postNaryConstraintTuple", (void(WeightedCSP::*)(int ctrindex, vector<Value>& tuple, Cost cost)) & WeightedCSP::postNaryConstraintTuple)
        .def("postNaryConstraintEnd", &WeightedCSP::postNaryConstraintEnd)
        .def("postSupxyc", &WeightedCSP::postSupxyc)
        .def("postDisjunction", &WeightedCSP::postDisjunction)
        .def("postSpecialDisjunction", &WeightedCSP::postSpecialDisjunction)
        .def("postAllDifferentConstraint", (int(WeightedCSP::*)(vector<int> scope, const string& arguments)) & WeightedCSP::postAllDifferentConstraint)
        .def("postGlobalCardinalityConstraint", (int(WeightedCSP::*)(vector<int> scope, const string& arguments)) & WeightedCSP::postGlobalCardinalityConstraint)
        .def("postCliqueConstraint", (int(WeightedCSP::*)(vector<int> scope, const string& arguments)) & WeightedCSP::postCliqueConstraint)
        .def("postKnapsackConstraint", [](WeightedCSP& s, vector<int> scope, const string& arguments, bool isclique, int kp, bool conflict) { return s.postKnapsackConstraint(scope, arguments, isclique, kp, conflict); }, py::arg("scope"), py::arg("arguments"), py::arg("isclique") = false, py::arg("kp") = 0, py::arg("conflict") = false)
        .def("postWeightedCSPConstraint", [](WeightedCSP& s, vector<int> scope, WeightedCSP* problem, WeightedCSP* negproblem, Cost lb, Cost ub, bool duplicateHard, bool strongDuality) { return s.postWeightedCSPConstraint(scope, problem, negproblem, lb, ub, duplicateHard, strongDuality); }, py::arg("scope"), py::arg("problem"), py::arg("negproblem"), py::arg("lb") = MIN_COST, py::arg("ub") = MAX_COST, py::arg("duplicateHard") = false, py::arg("strongDuality") = false)
        .def("postWAmong", (int(WeightedCSP::*)(vector<int> scope, const string& semantics, const string& propagator, Cost baseCost, const vector<Value>& values, int lb, int ub)) & WeightedCSP::postWAmong)
        .def("postWVarAmong", (void(WeightedCSP::*)(vector<int> scope, const string& semantics, Cost baseCost, vector<Value>& values)) & WeightedCSP::postWVarAmong)
        .def("postWRegular", (int(WeightedCSP::*)(vector<int> scope, const string& semantics, const string& propagator, Cost baseCost, int nbStates, const vector<WeightedObjInt>& initial_States, const vector<WeightedObjInt>& accepting_States, const vector<DFATransition>& Wtransitions)) & WeightedCSP::postWRegular)
        .def("postWAllDiff", (int(WeightedCSP::*)(vector<int> scope, const string& semantics, const string& propagator, Cost baseCost)) & WeightedCSP::postWAllDiff)
        .def("postWGcc", (int(WeightedCSP::*)(vector<int> scope, const string& semantics, const string& propagator, Cost baseCost, const vector<BoundedObjValue>& values)) & WeightedCSP::postWGcc)
        //        .def("postWSame", (int (WeightedCSP::*)(int* scopeIndexG1, int arityG1, int* scopeIndexG2, int arityG2, const string& semantics, const string& propagator, Cost baseCost)) &WeightedCSP::postWSame)
        //        .def("postWSameGcc", &WeightedCSP::postWSameGcc)
        //        .def("postWGrammarCNF", &WeightedCSP::postWGrammarCNF)
        //        .def("postMST", &WeightedCSP::postMST)
        //        .def("postMaxWeight", &WeightedCSP::postMaxWeight)
        //        .def("postWSum", &WeightedCSP::postWSum)
        //        .def("postWVarSum", &WeightedCSP::postWVarSum)
        //        .def("postWOverlap", &WeightedCSP::postWOverlap)
        .def("postWDivConstraint", &WeightedCSP::postWDivConstraint)
        .def("initDivVariables", &WeightedCSP::initDivVariables)
        .def("postGlobalFunction", (void(WeightedCSP::*)(vector<int> scope, const string& gcname, const string& arguments)) & WeightedCSP::postGlobalFunction)
        .def("isKnapsack", &WeightedCSP::isKnapsack)
        .def("isGlobal", &WeightedCSP::isGlobal)
        .def("getSolution", (const vector<Value> (WeightedCSP::*)()) & WeightedCSP::getSolution)
        .def("initSolutionCost", &WeightedCSP::initSolutionCost)
        .def("getSolutionValue", &WeightedCSP::getSolutionValue)
        .def("getSolutionCost", &WeightedCSP::getSolutionCost)
        .def("getSolutions", &WeightedCSP::getSolutions)
        //        .def("setSolution", &WeightedCSP::setSolution)
        .def("printSolution", (void(WeightedCSP::*)(ostream&)) & WeightedCSP::printSolution)
        .def("print", [](WeightedCSP& wcsp) { wcsp.print(cout); })
        .def("dump", &WeightedCSP::dump)
        .def("dump_CFN", &WeightedCSP::dump_CFN)
        .def("decimalToCost", &WeightedCSP::decimalToCost)
        .def("DoubletoCost", &WeightedCSP::DoubletoCost)
        .def("Cost2ADCost", &WeightedCSP::Cost2ADCost) // translate internal WCSP cost value to original problem real cost value (CFN)
        .def("Cost2RDCost", &WeightedCSP::Cost2RDCost)
        .def("Prob2Cost", &WeightedCSP::Prob2Cost)
        .def("Cost2Prob", &WeightedCSP::Cost2Prob)
        .def("Cost2LogProb", &WeightedCSP::Cost2LogProb)
        .def("LogProb2Cost", &WeightedCSP::LogProb2Cost)
        .def("LogSumExp", (Cost(WeightedCSP::*)(Cost c1, Cost c2) const) & WeightedCSP::LogSumExp)
        .def("LogSumExp", (TLogProb(WeightedCSP::*)(TLogProb logc1, Cost c2) const) & WeightedCSP::LogSumExp)
        .def("LogSumExp", (TLogProb(WeightedCSP::*)(TLogProb logc1, TLogProb logc2) const) & WeightedCSP::LogSumExp)
        .def("read", [](WeightedCSP& wcsp, const char* fileName) {
            if (strstr(fileName, ".xz") == &fileName[strlen(fileName) - strlen(".xz")])
                ToulBar2::xz = true;
            if (strstr(fileName, ".gz") == &fileName[strlen(fileName) - strlen(".gz")])
                ToulBar2::gz = true;
            if (strstr(fileName, ".bz2") == &fileName[strlen(fileName) - strlen(".bz2")])
                ToulBar2::bz2 = true;
            if (strstr(fileName, ".cfn"))
                ToulBar2::cfn = true;
            if (strstr(fileName, ".wcnf") || strstr(fileName, ".cnf"))
                ToulBar2::wcnf = true;
            if (strstr(fileName, ".qpbo"))
                ToulBar2::qpbo = true;
            if (strstr(fileName, ".opb"))
                ToulBar2::opb = true;
            if (strstr(fileName, ".wbo"))
                ToulBar2::opb = true;
            if (strstr(fileName, ".lp"))
                ToulBar2::lp = true;
            if (strstr(fileName, ".uai")) {
                ToulBar2::uai = 1;
                ToulBar2::bayesian = true;
            }
            if (strstr(fileName, ".LG")) {
                ToulBar2::uai = 2;
                ToulBar2::bayesian = true;
            }
#if defined(XMLFLAG) || defined(XMLFLAG3)
            if (strstr(fileName, ".xml")) {
                ToulBar2::xmlflag = true;
            }
#endif
            tb2checkOptions();
            return wcsp.read_wcsp(fileName); });

    py::class_<WeightedCSPSolver>(m, "Solver")
        .def(py::init([](Cost ub, WeightedCSP* wcsp) {
            ToulBar2::startCpuTime = cpuTime();
            ToulBar2::startRealTime = realTime();
            initCosts();
            if (ToulBar2::seed < 0) { // initialize seed using current time
                ToulBar2::seed = abs((int)time(NULL) * getpid() * ToulBar2::seed);
                if (ToulBar2::verbose >= 0)
                    cout << "Initial random seed is " << ToulBar2::seed << endl;
            }
            mysrand(ToulBar2::seed);
            if (ToulBar2::incop_cmd.size() > 0 && ToulBar2::seed != 1 && ToulBar2::incop_cmd.find("0 1 ") == 0) {
                string sseed = to_string(ToulBar2::seed);
                ToulBar2::incop_cmd.replace(2, 1, sseed);
            }
            return WeightedCSPSolver::makeWeightedCSPSolver(ub, wcsp);
        }),
            py::arg("ub") = MAX_COST,
            py::arg("wcsp") = nullptr)
        .def_property_readonly("wcsp", &WeightedCSPSolver::getWCSP, py::return_value_policy::reference_internal)
        .def("read", [](WeightedCSPSolver& s, const char* fileName) {
            if (strstr(fileName, ".xz") == &fileName[strlen(fileName) - strlen(".xz")])
                ToulBar2::xz = true;
            if (strstr(fileName, ".gz") == &fileName[strlen(fileName) - strlen(".gz")])
                ToulBar2::gz = true;
            if (strstr(fileName, ".bz2") == &fileName[strlen(fileName) - strlen(".bz2")])
                ToulBar2::bz2 = true;
            if (strstr(fileName, ".cfn"))
                ToulBar2::cfn = true;
            if (strstr(fileName, ".wcnf") || strstr(fileName, ".cnf"))
                ToulBar2::wcnf = true;
            if (strstr(fileName, ".qpbo"))
                ToulBar2::qpbo = true;
            if (strstr(fileName, ".opb"))
                ToulBar2::opb = true;
            if (strstr(fileName, ".wbo"))
                ToulBar2::opb = true;
            if (strstr(fileName, ".lp"))
                ToulBar2::lp = true;
            if (strstr(fileName, ".uai")) {
                ToulBar2::uai = 1;
                ToulBar2::bayesian = true;
            }
            if (strstr(fileName, ".LG")) {
                ToulBar2::uai = 2;
                ToulBar2::bayesian = true;
            }
#if defined(XMLFLAG) || defined(XMLFLAG3)
            if (strstr(fileName, ".xml")) {
                ToulBar2::xmlflag = true;
            }
#endif
            tb2checkOptions();
            return s.read_wcsp(fileName);
        })
#ifndef __WIN32__
        .def("timer", [](WeightedCSPSolver& s, int timeout) {
            signal(SIGINT, timeOut);
            if (timeout > 0)
                timer(timeout);
        })
        .def("timerStop", [](WeightedCSPSolver& s) {
            timerStop();
        })
#endif
        .def("solve", [](WeightedCSPSolver& s, bool first) {
            bool res = false;
            try {
                res = s.solve(first);
            } catch (Contradiction) {
                s.getWCSP()->whenContradiction();
                if (ToulBar2::verbose >= 0) cout << "No solution found by initial propagation!" << endl;
                return false;
            }
            return res; }, py::arg("first") = true)
        .def("beginSolve", &WeightedCSPSolver::beginSolve)
        .def("preprocessing", &WeightedCSPSolver::preprocessing)
        .def("recursiveSolve", &WeightedCSPSolver::recursiveSolve)
        .def("recursiveSolveLDS", &WeightedCSPSolver::recursiveSolveLDS)
        .def("hybridSolve", &WeightedCSPSolver::hybridSolve)
        .def("endSolve", &WeightedCSPSolver::endSolve)
        .def("solution", (const vector<Value> (WeightedCSPSolver::*)()) & WeightedCSPSolver::getSolution)
        .def("solutionValue", &WeightedCSPSolver::getSolutionValue)
        .def("solutionCost", &WeightedCSPSolver::getSolutionCost)
        .def("solutions", &WeightedCSPSolver::getSolutions)
        .def("getNbNodes", &WeightedCSPSolver::getNbNodes)
        .def("getNbBacktracks", &WeightedCSPSolver::getNbBacktracks)
        .def("getDDualBound", &WeightedCSPSolver::getDDualBound)
        .def("increase", &WeightedCSPSolver::increase)
        .def("decrease", &WeightedCSPSolver::decrease)
        .def("assign", &WeightedCSPSolver::assign)
        .def("remove", &WeightedCSPSolver::remove)
        .def("generate", &WeightedCSPSolver::read_random)
        //        .def("narycsp", &WeightedCSPSolver::narycsp)
        //        .def("solve_symmax2sat", &WeightedCSPSolver::solve_symmax2sat)
        .def("dump_wcsp", (void(WeightedCSPSolver::*)(const char*, bool, int)) & WeightedCSPSolver::dump_wcsp)
        .def("read_solution", &WeightedCSPSolver::read_solution)
        .def("parse_solution", &WeightedCSPSolver::parse_solution);
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
