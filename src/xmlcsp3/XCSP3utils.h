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
#ifndef UTILS_H
#define    UTILS_H


#include <string>
#include <sstream>
#include <vector>
#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>
#include <cwctype>

namespace XCSP3Core {
    std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);


    std::vector<std::string> split(const std::string &s, char delim);


    template<typename Base, typename T> inline bool instanceof(const T *);


    void ReplaceStringInPlace(std::string &subject, const std::string &search, const std::string &replace);


// trim from start
    std::string &ltrim(std::string &s);



// trim from end
    std::string &rtrim(std::string &s);


    std::string &removeChar(std::string &s, char c);



// trim from both ends
    std::string &trim(std::string &s);


}

#endif	/* UTILS_H */

