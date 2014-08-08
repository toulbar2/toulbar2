#ifndef  INCOP_HPP_
#define  INCOP_HPP_
#include <cerrno>
#include <stdio.h>
#include <list>
#include <vector>
#include <string>
#include <set>
#include <algorithm>
using namespace std;
#include <iostream>
#include <fstream>
#include "timer.h"
#include "incop.h"
#include "csproblem.h"
#include "incoputil.h"
#include "narycsproblem.h"

#include <string.h>

#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/** le lecteur des fichiers au format wcsp de Simon de Givry
    pour le moment, ne lit que des problèmes avec contraintes  en extension , les valeurs des domaines
doivent être des entiers.

*/


int split (char *str, char c, char ***arr);
int narycspmain( string filename, string outputfile , string cmd, int verbose) ;


#endif /*INCOP_HPP_*/
