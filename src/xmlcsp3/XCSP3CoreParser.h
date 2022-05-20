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
#ifndef _XMLParser_libxml2_h_
#define _XMLParser_libxml2_h_

#include <iostream>
#include <stdexcept>
#include <cerrno>
#include <climits>
#include <libxml/parser.h>

#include "XMLParser.h"
#include "XCSP3CoreCallbacks.h"
#include "XCSP3Constants.h"
#include "UTF8String.h"
//#define debug

namespace XCSP3Core {
    using namespace std;


    /**
     * @brief the parser using the libxml2 library
     */
    class XCSP3CoreParser {

    protected:
        XMLParser cspParser;

    public:

        XCSP3CoreParser(XCSP3CoreCallbacks *cb) : cspParser(cb) {
            LIBXML_TEST_VERSION
        }


        int parse(istream &in);


        int parse(const char *filename);


    protected:


        /*************************************************************************
         *
         * SAX Handler
         *
         *************************************************************************/

        static void comment(void *parser, const xmlChar *value);


        static void startDocument(void *parser);


        static void endDocument(void *parser);


        static void characters(void *parser, const xmlChar *ch, int len);

        static void startElement(void *parser, const xmlChar *name, const xmlChar **attr);

        static void endElement(void *parser, const xmlChar *name);
    };

}

#endif
