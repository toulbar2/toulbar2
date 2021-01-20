#include "utils/SimpleOpt.h"
#include "utils/SimpleGlob.h"
#include "utils/SimpleIni.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <regex>
#include "cfntools.h"
#include <regex>
#include <string>


#if defined(_MSC_VER)
# include <windows.h>
# include <tchar.h>
#else
# define TCHAR		char
# define _T(x)		x
# define _tprintf	printf
# define _tmain		main
# define _ttoi      atoi
#endif

#include "./rapidjson/document.h"
#include "./rapidjson/filereadstream.h"
#include "./rapidjson/reader.h"
#include "./rapidjson/filewritestream.h"
#include "./rapidjson/writer.h"

using namespace rapidjson;
using namespace std;

//////////////////////////////////////////
///this is a tool to manage cfn        ///
//////////////////////////////////////////

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////





//////////////////////////////////////////
///this is a tool to manage cfn        ///
//////////////////////////////////////////

void ShowFiles(int argc, TCHAR ** argv) {
	// glob files to catch expand wildcards
	CSimpleGlob glob(SG_GLOB_NODOT|SG_GLOB_NOCHECK);
	if (SG_SUCCESS != glob.Add(argc, argv)) {
		_tprintf(_T("Error while globbing files\n"));
		return;
	}

	for (int n = 0; n < glob.FileCount(); ++n) {
		_tprintf(_T("file %2d: '%s'\n"), n, glob.File(n));

	}
}

////////////////////////////////////////
/// Management of options

// OPTION LIST 

enum{
	OPT_DIGIT,
	OPT_THR,
	OPT_M,
	OPT_F,
	OPT_O,
	OPT_UB,
	OPT_EXT,
	OPT_ADD,
	OPT_PS,
	OPT_PP,
	OPT_COMP,
	OPT_MUL,
	OPT_MULSELF,
	OPT_MULPAIR,
	OPT_SOL,
	OPT_FAST,
	OPT_MEAN,
	OPT_INFO,
	OPT_HELP,
	OPT_verbose,
	OPT_regexp_ext
};

CSimpleOpt::SOption g_rgOptions[] = {
	{ OPT_DIGIT,  "-d", SO_REQ_SEP },
	{ OPT_THR,  "-t", SO_REQ_SEP },
	{ OPT_M,  "-merge", SO_REQ_SEP },
	{ OPT_F,  "-f", SO_REQ_SEP },
	{ OPT_O, "-o", SO_REQ_SEP},
	//{ OPT_UB, "-ub", SO_REQ_SEP},
	//{ OPT_UB, "-UB", SO_REQ_SEP},
	{ OPT_UB, (char*)"-UB", SO_OPT }, // Upperbound setting 
	{ OPT_EXT, "-extract", SO_REQ_SEP},
	{ OPT_ADD, "-add", SO_REQ_SEP},
	{ OPT_PS, "-pctself", SO_REQ_SEP},
	{ OPT_PP, "-pctpair", SO_REQ_SEP},
	{ OPT_PS, "-pctUnary", SO_REQ_SEP}, // renaming of pctse
	{ OPT_PP, "-pctBin", SO_REQ_SEP},
	{ OPT_COMP, "-comp", SO_REQ_SEP},
	{ OPT_MUL, "-mul", SO_REQ_SEP},
	{ OPT_MULSELF, "-mulself", SO_REQ_SEP},
	{ OPT_MULPAIR, "-mulpair", SO_REQ_SEP},
	{ OPT_SOL, "-sol", SO_REQ_SEP},
	{ OPT_FAST, "-fast", SO_NONE},
	{ OPT_MEAN, "-meandeltafunction", SO_NONE},
	{ OPT_INFO, "-info", SO_NONE},
	{ OPT_HELP, "-h", SO_NONE },
	{ OPT_HELP, "--help", SO_NONE },
	{ OPT_regexp_ext, (char*)"--regexp_ext", SO_REQ_SEP },
	{ OPT_verbose, (char*)"-v", SO_OPT }, // verbose level
	{ OPT_verbose, (char*)"--verbose", SO_OPT }, // verbose level
//	{ OPT_verbose, "--verbose", SO_REQ_SEP},
	SO_END_OF_OPTIONS
};

/*----------------*/
	static const TCHAR * 
GetLastErrorText(
		int a_nError
		) 
{
	switch (a_nError) {
		case SO_SUCCESS:            return _T("Success");
		case SO_OPT_INVALID:        return _T("Unrecognized option");
		case SO_OPT_MULTIPLE:       return _T("Option matched multiple strings");
		case SO_ARG_INVALID:        return _T("Option does not accept argument");
		case SO_ARG_INVALID_TYPE:   return _T("Invalid argument format");
		case SO_ARG_MISSING:        return _T("Required argument is missing");
		case SO_ARG_INVALID_DATA:   return _T("Invalid argument data");
		default:                    return _T("Unknown error");
	}
}


// chek if a filename end as ext
/*----------------*/
bool check_file_ext(const string fileName, const string ext)
{
	size_t extLen = ext.length();
	return ((extLen <= fileName.length()) && (fileName.substr(fileName.length() - extLen) == ext));
}

/*----------------*/
/*----------------*/

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

/// Management of options
void Help() {
    cout<<"Usage: cfntools -f file.cfn  [-o output.cfn] [-comp file] [-merge file] [-add file] [-pctself ps] [-pctpair pp] "<<endl;
    cout<<"                             [-extract list] [-mul k] [-mulself ks] [-mulpair kp] [-UB value]  [-t threshold] [-d nb_digit]"<<endl;
    cout<<"                             [--verbose v] [-sol file.sol] [-fast][-h] [-meandeltafunction] [-info]"<<endl;
    cout<< "\n" << endl;
    cout<<"-f file.cfn : cfn file used as reference "<<endl;
    cout<<"-o output.cfn : location of the output"<<endl;

    cout<< "\n" << endl;
    cout<<"-merge file2.cfn : merge  -f file.cfn and file2.cfn "<<endl;
    cout<<"-add file2.cfn : merge with file2.cfn constraints with -f file1.cfn the result will be save into new cfn file : -o foo.fcn "<<endl;
    cout<<"-pctself ps :  unary cost function rescaling when 2 cfn are added (cf. -add option)  , unary cost of the second one  will be rescale by a factor : ps*reader.meandeltaunary*map.meandeltaunary "<<endl;
    cout<<"-pctpair pp : when adding 2 cfn,  binary cost of the second cfn  will be rescale by  a factor : pp*reader.meandeltabinary*map.meandeltabinary"<<endl;

    cout<< "\n" << endl;
    cout<<"-extract var_list.txt : list of variables to extract"<<endl;
    cout<<"-mul k : multiply by k (double) all the values (saved in a new cfn)"<<endl;
    cout<<"-mulself k : multiply by k (double) unary cost functions  (saved in a new cfn cf -o option)"<<endl;
    cout<<"-mulpair k : multiply by k (double) binary cost functions  (saved in -o foo_out.cfn )"<<endl;

    cout<< "\n" << endl;
    cout<<"-UB upperbound to adjust (double) : -UB 0 => compute the UB according to the cost in the cfn"<<endl;
    cout<<"-t threshold (double) : all values above threshold are set equal to upper bound"<<endl;

    cout<<"-d number of digits (int) : -d 0 => adjust to Upperbound format (by default prec = 6 digits)"<<endl;
    cout<<"-meandeltafunction : return the delta mean for unary cost and for binary costs "<<endl;

    cout<< "\n" << endl;
    cout<<"-comp file2.cfn : compare two cfn files specified as : -f file1.cfn -comp file2.cfn "<<endl;
    cout<<"-sol file.sol : solution cost checking and costs corresponding to each function is extract in csv file "<<endl;
    cout<<"-info : format checking and printing of information about thei specified cfn"<<endl;

    cout<<"-fast : improve processing speed json format checking and verbose level are turn to off)"<<endl;

    cout<< "\n" << endl;
    cout<< "-------------\n" << endl;
    cout<<"--verbose v (int) : quantity of print msg (0 or 1 by default 1)"<<endl;
    cout<<"-h : display this help"<<endl;
    cout<< "----\n" << endl;
    cout<<"WARNING : the operations are done in this order, please verify if the steps are corresponding to what you want"<<endl;
};
//int main(int argc, char* argv []){
int _tmain(int argc, TCHAR * argv[]) {
  //CSimpleOpt args(argc, argv, g_rgOptions);
  CSimpleOpt args(argc, argv, g_rgOptions, SO_O_NOERR|SO_O_EXACT);

  int digits = -1;
  double threshold = 99999;
  string merge_file = "";
  string cfn_file = "";
  string output_file = "";
  string sol_file = "";
  string f_extract = "";
  string map= "";
  double ps= -1;
  double pp= -1;
  bool modif = false; // do we need an output file ?
  double UB = -1;
  double mult = 1;
  double mulself = 1;
  double mulpair = 1;
  string comp_file="";
  int verbose=0;
  bool fast=false;
  int precision = 6; // by default

  vector<string> command{};

/////////////////////////////////////////////////////////////////////////////////////////////

/*----------------*/



 /* mpa files extension */
 std::map<std::string, string> file_extension_map;
 file_extension_map["regexp_ext"] = ".regex";

 /* map regexp */

 std::map<std::string, string> regexp_map;

 string var_valid;
 string val_valid;
 string f_valid;
 string ub_valid;
 string var_merged_valid;


/////////////////////////////////////////////////////////////////////////////////////////////


// initialize arguments
while (args.Next()) {
	if (args.LastError() == SO_SUCCESS) {

		if (args.OptionId() == OPT_HELP) {
			Help();
			return 0;
		}

		if (args.OptionId() == OPT_verbose) {
                         if (args.OptionArg() != NULL) {
                                 verbose = atoi(args.OptionArg());
                         } else {
                                 verbose = 1;
                         }

		}

		if (args.OptionId() == OPT_F){
			cfn_file = args.OptionArg();
		}
		else if (args.OptionId() == OPT_INFO) {
			command.push_back("info");
		}
		else if (args.OptionId() == OPT_DIGIT){
			digits = stoi(args.OptionArg());
			modif = true;
			command.push_back("digit");
		}
		else if (args.OptionId() == OPT_THR){
			threshold = stof(args.OptionArg());
			modif = true;
			command.push_back("threshold");
		}
		else if (args.OptionId() == OPT_M){
			merge_file = args.OptionArg();
			modif = true;
			command.push_back("merge");
		}

		else if (args.OptionId() == OPT_UB){
			if (args.OptionArg() != NULL) {
				cout << "UB set to " << UB << endl;
				UB = stof(args.OptionArg());
			} else {
				cout << "UB update with max sum value" << endl;
				UB = 0; 
			}
			if(verbose > 1 ){ cout <<" UB set to " << UB << " in command line " << endl; }

			modif = true;
			command.push_back("ub");
		}
		else if (args.OptionId() == OPT_O){
			if (args.OptionArg() != NULL) {
			output_file = args.OptionArg();
			} else {
			output_file = "output.cfn";
			}
			modif = true;
		}
		else if (args.OptionId() == OPT_EXT){
			f_extract = args.OptionArg();
			modif = true;
			command.push_back("extract");
		}
		else if (args.OptionId() == OPT_ADD){
			map = args.OptionArg();
			modif = true;
			command.push_back("addmap");
		}
		else if (args.OptionId() == OPT_PS){
			ps = stof(args.OptionArg());
		}
		else if (args.OptionId() == OPT_PP){
			pp = stof(args.OptionArg());
		}
		else if (args.OptionId() == OPT_COMP){
			comp_file = args.OptionArg();
			command.push_back("comp");
		}
		else if (args.OptionId() == OPT_MUL){
			mult = atof(args.OptionArg());
			modif = true;
			command.push_back("mulall");
		}
		else if (args.OptionId() == OPT_MULSELF){
			mulself = atof(args.OptionArg());
			modif = true;
			command.push_back("mulself");
		}
		else if (args.OptionId() == OPT_MULPAIR){
			mulpair = atof(args.OptionArg());
			modif = true;
			command.push_back("mulpair");
		}
		else if (args.OptionId() == OPT_SOL){
			sol_file = args.OptionArg();
			command.push_back("sol");
		}
		else if (args.OptionId() == OPT_MEAN){
			command.push_back("mean");
		}
		else if (args.OptionId() == OPT_FAST){
			fast = true;
			verbose = 0;
		}
		else {
			// handle error: see ESOError enums
			cerr<<"unrecognized argument"<<endl;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// manage errors
if(modif && output_file==""){
	cerr<<">> ERROR : output location not defined (please use option : -o ./foo/foo.cfn ) "<<endl;
	return(1);
}
if(cfn_file==""){
	cerr<<"Please precise the input location (-f)"<<endl;
	return(1);
}else { cout << "ref cfn : " << cfn_file << endl; }

// in all cases
CFNManager reader(cfn_file);  // read the file
if(!fast){
	reader.check_format(reader.cfn, regexp_map, verbose);  // check the format
	cout<<"cfn format checked "<<endl;
}
if(modif && verbose != 0){
	cout<<"begin the commands"<<endl;
}


// init regexp value 
regexp_map["ub_valid"]= "<-*[0-9]*\\.[0-9]*";
regexp_map["var_valid"]= "[A-Z]_[0-Z]_[0-9]*";
regexp_map["val_valid"]= "[A-Z]*_[0-9]*";
regexp_map["f_valid"]= "(C*B[0-9]*-[0-9]*)|(C*M*S[0-9]*)|(C*M*Etmp)|(C*M[0-9]*)" ;
regexp_map["var_merged_valid"]= "M*[A-Z]_[0-Z]_[0-9]*";


// default  init value regexp
ub_valid  = regexp_map["ub_valid"] ;
var_valid = regexp_map["var_valid"] ;
val_valid = regexp_map["val_valid"] ;
f_valid   = regexp_map["f_valid"] ;
var_merged_valid = regexp_map["var_merged_valid"] ;


////////////////////////////////////////
// --------------------------simple opt ----------------------

// as well as our array of valid options.

// process any files that were passed to us on the command line.
// send them to the globber so that all wildcards are expanded
// into valid filenames (e.g. *.cpp -> a.cpp, b.cpp, c.cpp, etc)
// See the SimpleGlob.h header file for details of the flags.

//////////////////////BEGIN READ FILE //////////////////////////
CSimpleGlob glob(SG_GLOB_NODOT | SG_GLOB_NOCHECK);
if (SG_SUCCESS != glob.Add(args.FileCount(), args.Files())) {
	_tprintf(_T("Error while globbing files\n"));
	return 1;
}

// dump all of the details, the script that was passed on the
// command line and the expanded file names

if ( glob.FileCount() > 0 ) {
	cout << "cmd line parsing ---> " << glob.FileCount() << " Filename(s) found in command line" << endl;

	for (unsigned int n = 0; n < glob.FileCount() ; ++n) {


		if (check_file_ext(glob.File(n), file_extension_map["regexp_ext"])) {
			std::string myContent ;
			cout << "loading rexpep file: " << glob.File(n) << endl;
			//const bool success = loadFileReg(glob.File(n), myContent,  file_extension_map) ;

			//				const bool success = loadFile(glob.File(n), myContent) ;
			//std::cout << "content is: \n ------\n " << myContent << "\n ------- \n";


			CSimpleIniA ini;
			ini.SetUnicode();
			SI_Error rc = ini.LoadFile(glob.File(n));
			if (rc < 0) { /* handle error */ };

			cout << " INIT FROM ini file : " << glob.File(n) << endl;
			for (std::map<string,string>::iterator it=regexp_map.begin(); it!=regexp_map.end(); ++it)
			{ 	
				string label =  it->first ;
				regexp_map[label.c_str()]= ini.GetValue("VAR_REGEXP" , label.c_str());
				std::cout << it->first << " => " << it->second << '\n';
			}

		}
	}

} // end glob.FileCount() > 0 

// echo regexp used during processus 
if (verbose >0 ) {
	cout << "verbose level = " << verbose << endl;
	cout << "---------INIT regexp  for variable and cost function names----------------" << endl;
	for (std::map<string,string>::iterator it=regexp_map.begin(); it!=regexp_map.end(); ++it) {
		std::cout << it->first << " => " << it->second << '\n';
	}
	cout << "-------------------------------------------------------" << endl;
}
/* ------------------------ */
// -----------------file management ---------------------
for(vector<string>::iterator it = command.begin(); it != command.end(); it++){
	// infos
	if(*it=="info"){
		if(verbose != 0) cout<<"begin giving information about current cfn"<<endl;
		reader.infos();
		if(verbose != 0) cout<<"end giving information about current cfn"<<endl;
	}
	// extract cost solution
	else if(*it=="sol"){
		if(verbose != 0) cout<<"begin creating cost solution"<<endl;
		reader.extract_cost(sol_file, verbose);
		if(verbose != 0) cout<<"end creating cost solution"<<endl;
	}
	// comparing
	else if(*it=="comp"){
		if(verbose != 0) cout<<"begin compare :" << comp_file <<endl;
		reader.comp(comp_file, fast, regexp_map, verbose);
	}
	// merging
	else if(*it=="merge"){
		if(verbose != 0) cout<<"begin merging"<<endl;
		reader.merge(merge_file, fast, regexp_map, verbose);
		if(verbose != 0) cout<<"end merging"<<endl;
	}
	// add map
	else if(*it=="addmap"){
		if(verbose != 0) cout<<"begin merge map"<<endl;
		reader.add_map(map, ps, pp, fast, verbose, regexp_map );
		if(verbose != 0) cout<<"end merge map"<<endl;
	}
	// extract
	else if(*it=="extract"){
		if(verbose != 0) cout<<"begin extracting"<<endl;
		reader.extract(f_extract,verbose);
		if(verbose != 0) cout<<"end extracting"<<endl;
	}
	// multiply
	else if(*it=="mulall"){
		if(mult!=1){
			if(verbose != 0) cout<<"begin multiply all"<<endl;
			reader.multiply(mult);
			if(verbose != 0) cout<<"end multiply all"<<endl;
		}
	}
	else if(*it=="mulself"){
		if(verbose != 0) cout<<"begin multiply self"<<endl;
		reader.multiply(mulself, 1);
		if(verbose != 0) cout<<"end multiply self"<<endl;
	}
	else if(*it=="mulpair"){
		if(verbose != 0) cout<<"begin multiply pair"<<endl;
		reader.multiply(1, mulpair);
		if(verbose != 0) cout<<"end multiply pair"<<endl;
	}
	// update ub
	else if(*it=="ub"){
		if(verbose != 0) cout<<"begin update"<<endl;
		reader.updateUB(UB, verbose);
		if(verbose != 0) cout<<"end update"<<endl;
	}
	// threshold
	else if(*it=="threshold"){
		if(verbose != 0) cout<<"begin value thresholding"<<endl;
		reader.filter_value(threshold);
		if(verbose != 0) cout<<"end thresholding"<<endl;
	}
	// check mean delta energy
	else if(*it=="mean"){
		if(verbose != 0) cout<<"begin compute mean delta"<<endl;
		double meanself, meanpair, medself, medpair;
		reader.compute_meandelta(meanself, meanpair, medself, medpair, verbose);
		cout<<"distribution of Cmax - Cmin accross all functions"<<endl;
		cout<<"meanself: "<<meanself<<"\n"<<"meanpair: "<<meanpair<<endl;
		cout<<"median self: "<<medself<<"\n"<<"median pair: "<<medpair<<endl;
		if(verbose != 0) cout<<"end compute mean delta"<<endl;
	}
	// digits
	else if(*it=="digit"){
		if(verbose != 0) cout<<"begin digits"<<endl;
		precision = reader.check_digits(digits);
		if(verbose != 0) cout<<"digits ok"<<endl;
	}
}
// save
if(output_file != ""){
	if(verbose != 0) cout<<"saving file"<<endl;
	reader.save(output_file, precision);
	if(verbose != 0) cout<<"file saved"<<endl;
}

cout << "----------- end -----------" << endl;
return(0);
}
