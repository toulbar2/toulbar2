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





int split (char *str, char c, char ***arr)
{
	int count = 1;
	int token_len = 1;
	int i = 0;
	char *p;
	char *t;

	p = str;
	while (*p != '\0')
	{
		if (*p == c) count++;
		p++;
	}

	*arr = (char**) malloc(sizeof(char*) * count);
	if (*arr == NULL)
		exit(1);

	p = str;
	while (*p != '\0')
	{
		if (*p == c)
		{
			(*arr)[i] = (char*) malloc( sizeof(char) * token_len );
			if ((*arr)[i] == NULL)
				exit(1);

			token_len = 0;
			i++;
		}
		p++;
		token_len++;
	}
	(*arr)[i] = (char*) malloc( sizeof(char) * token_len );
	if ((*arr)[i] == NULL) exit(1);

	i = 0;
	p = str;
	t = ((*arr)[i]);
	while (*p != '\0')
	{
		if (*p != c && *p != '\0')
		{
			*t = *p;
			t++;
		}
		else
		{
			*t = '\0';
			i++;
			t = ((*arr)[i]);
		}
		p++;
	}

	return count;
}

void removeSpaces(string& str)
{
	/* remove multiple spaces */
	int k=0;
	for (unsigned int j=0; j<str.size(); ++j)
	{
		if ( (str[j] != ' ') || (str[j] == ' ' && str[j+1] != ' ' ))
		{
			str [k] = str [j];
			++k;
		}

	}
	str.resize(k);

	/* remove space at the end */   
	if (str [k-1] == ' ')
		str.erase(str.end()-1);
	/* remove space at the begin */
	if (str [0] == ' ')
		str.erase(str.begin());
}


int narycspmain( string filename, string outputfile , string cmd, int verbose) {

	char line[1024];
	char **tokens= NULL;
	int tuningmode=0;
	int argc=0;

	 // remove leading space  from incop command line
	cmd.erase(cmd.begin(), std::find_if(cmd.begin(), cmd.end(), std::bind1st(std::not_equal_to<char>(), ' ')));

	// remove multiples space in cmd
	removeSpaces(cmd);

	cout << "Narycsp parameter : " << cmd << endl;
	sprintf(line,"bin/Linux/narycsp %s %s %s",  outputfile.c_str(), filename.c_str(), cmd.c_str());

	argc =  split(line, ' ', &tokens);

	if ( verbose > 0 ) {

		cout << "---------------------------" << endl ;
		printf("NARYCSP ARGC found : %d .\n", argc);
		cout << "---------------------------" << endl ;
		for (int i = 0; i <argc; i++) printf("arg #%d --> %s ; ", i, tokens[i]);
	}

	narycsp(argc,tokens,tuningmode);


	return 0;

}

