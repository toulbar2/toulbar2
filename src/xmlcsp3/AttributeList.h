/*=============================================================================
 * parser for CSP instances represented in XCSP3 Format
 *
 * Copyright (c) 2015 xcsp.org (contact <at> xcsp.org)
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

#ifndef COSOCO_ATTRIBUTELIST_H
#define COSOCO_ATTRIBUTELIST_H

#include <libxml/xmlstring.h>
namespace XCSP3Core {
    /**
 * represents the attribute list of a XML tag
 */
    class AttributeList {
    public:
        typedef unsigned char Byte;


        /**
         * an empty list of attributes
         */
        AttributeList() {
            n = 0;
            list = NULL;
        }


        AttributeList(const Byte **attr) {
            list = attr;

            n = 0;
            if(list == NULL)
                return;

            while(list[n] != NULL)
                n += 2;

            n /= 2;
        }


        inline int size() const {
            return n;
        }


        UTF8String operator[](const char *name) const {
            for(int i = 0; i < n; ++i)
                if(xmlStrEqual(list[2 * i], reinterpret_cast<const Byte *> (name)))
                    return UTF8String(list[2 * i + 1]);

            return UTF8String();
        }


        inline UTF8String getName(int i) const {
            return UTF8String(list[2 * i]);
        }


        inline UTF8String getValue(int i) const {
            return UTF8String(list[2 * i + 1]);
        }


    private:
        int n; // number of attributes
        const Byte **list; // list[2*i] is the name of the i-th attribute,
        // list[2*i+1] is its value
    };

}
#endif //COSOCO_ATTRIBUTELIST_H
