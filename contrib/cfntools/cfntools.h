#ifndef CFNTOOLS
#define CFNTOOLS
#include "./rapidjson/document.h"
#include "./rapidjson/filewritestream.h"
#include "./rapidjson/reader.h"
#include "./rapidjson/filewritestream.h"
#include "./rapidjson/writer.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <regex>

using namespace rapidjson;
using namespace std;


///////////////////////////////////
///   CFNManager Class          ///
/// goal : manage CFN           ///
/// initialized with a cfn file ///
///////////////////////////////////
class CFNManager{
public:
  Document cfn;
  string name_cfn;
  CFNManager(string& cfn_file);

  void check_format(Document& doc, map<string, string> & regexp_map , int verbose); // verify the format
  int check_digits(int digits); // verify the digits

  void filter_value(double threshold); // delete value above a threshold
  void merge(string& cfn_file, bool fast, map<string, string> & regexp_map , int verbose); // merge two cfn
  void updateUB(double new_UB){updateUB(new_UB, 1);}
  void updateUB(double new_UB, int verbose); // update the value of UB
  void extract(string& var_list,int verbose); // extract the subproblem corresponding to a list of variables
  void extract_cost(string& sol_file, int verbose); // write the cost of the solution in a csv
  void add_map(string& map_file, double ps, double pp, bool fast, int verbose  ,  map<string, string> & regexp_map ); // merge with an other cfn corresponding to the map (cryo em)
  void comp(string& cfn_file, bool fast , map<string, string> & regexp_map, int verbose); // compare two cfn
  void multiply(double m); // multiply all the costs by m
  void multiply(double mself, double mpair); // multiply the unary cost by mself and the binary cost by mpair (if values != 0)

  void compute_meandelta(double& meanself, double& meanpair, double& medself, double& medpair, int verbose); // compute mean delta energy (to merge with cfn map)
  //void extract_up(double threshold);
  //void extract_mean(double threshold);
  //void extract_threshold(double threshold, string cond)
  // si "mean" on fait mean de l'array, sinon si "max" on fait le max
  void infos(); // gives information about cfn
  void save(string& output_file, int precision); // save the corresponding cfn
};

#endif


