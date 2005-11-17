/** \file tb2utils.hpp
 *  \brief Miscelaneous usefull functions.
 * 
 */

#ifndef TB2UTILS_HPP_
#define TB2UTILS_HPP_

template<class T>
T abs(T x) {
    if (x < 0) return -(x);
    else return x;
}

// Warning! Already defined in STL
//template<class T>
//T min(T x, T y) {
//    if (x < y) return x;
//    else return y;
//}
//
//template<class T>
//T max(T x, T y) {
//    if (x > y) return x;
//    else return y;
//}

template<class T>
T min(T *array, int size)
{
    assert(size >= 1);
    T res = array[0];
    for (int i=1; i < size; i++) {
        if (array[i] < res) {
            res = array[i];
        }
    }
    return res;
}

template<class T>
T max(T *array, int size)
{
    assert(size >= 1);
    T res = array[0];
    for (int i=1; i < size; i++) {
        if (array[i] > res) {
            res = array[i];
        }
    }
    return res;
}

#endif /* TB2UTILS_HPP_ */
