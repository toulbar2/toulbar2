/**
 * \file tb2files.hpp
 * \brief General purpose class to play with files
 * \author kad
 * \date 24 mai 2019
 * \version 0.1
 */

#ifndef SRC_UTILS_TB2FILES_KAD_HPP_
#define SRC_UTILS_TB2FILES_KAD_HPP_
#include <string>
#include <fstream>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
//#include "toulbar2lib.hpp"
//#include "utils/tb2store.hpp"

class Tb2Files {





public:
static long nbProcess(const std::string nbProcessFic, const int nbCores, const int nbProcPerCore);
static void write_file(const std::string fileName , const std::string textToWrite);
static void read_file(const std::string );
std::vector<long> readDomains(const std::string fic);

	Tb2Files();
	virtual ~Tb2Files();
};

#endif /* SRC_UTILS_TB2FILES_KAD_HPP_ */
