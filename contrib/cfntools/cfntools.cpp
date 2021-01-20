#include "cfntools.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <bits/stdc++.h>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <regex>
#include <iomanip>

#include "./rapidjson/document.h"
#include "./rapidjson/filereadstream.h"
#include "./rapidjson/reader.h"
#include "./rapidjson/filewritestream.h"
#include "./rapidjson/writer.h"

using namespace rapidjson;
using namespace std;

/// initialize the CFNManager in loading a cfn file and parsing it into a Document Class
CFNManager::CFNManager(string& cfn_file){  // NP : OK
  name_cfn = cfn_file;
  while(name_cfn.find("/")!=string::npos){
    name_cfn = name_cfn.substr(name_cfn.find("/")+1,name_cfn.size());
  }
  FILE* fp;
  fp = fopen(cfn_file.c_str(), "r");  // use "rb" on Windows
  if (fp!= NULL){
  ifstream in(cfn_file.c_str());
  string contents((istreambuf_iterator<char>(in)), istreambuf_iterator<char>());
  fclose(fp);
  const char* input = contents.c_str();
  cfn.Parse(input);
  }
  else{
    cerr<<"cannot open file "<<cfn_file<<endl;
    exit(1);
  }
}

/// save the Document currently in the class into a cfn file
/// input : string& output_file
void CFNManager::save(string& output_file, int precision){ // NP : OK
  StringBuffer buffer;
  Writer<StringBuffer> writer(buffer);
  writer.SetMaxDecimalPlaces(precision);
  cfn.Accept(writer);
  const char* output = buffer.GetString();
  FILE *out;
  out = fopen(output_file.c_str(), "w");
  if (out != NULL){
  fprintf(out, "%s", output);
  fclose(out);
  }
}

/// check if the Document has a cfn format
/// TO SEE exactly what is a cfn format
void CFNManager::check_format(Document& doc, std::map<std::string, string> & regexp_map , int  verbose){  // NP : OK
  // verify general Organisation
	if( !(doc.IsObject() )) { 
			cerr << "\n ----------\n error: document root error \n json basic rules not respected \n  Please :check syntax  in the json format with : \n    python3 -m json.tool foo.cfn  \n\n" << endl ; 
			exit ;
			} 
	// verify header
	if( !(doc.HasMember("problem") )) {  
				cerr << "Header error problem not defined " << endl;
				exit;
				} else if ( ! (doc["problem"].HasMember("name"))) {
				exit;
				cerr << "Header error problem name not defined " << endl;
				} else if (!(doc["problem"]["name"].IsString())) {
					cerr << "Header error problem name is not a string  " <<   endl;
					exit;
				} else if (  !(doc["problem"].HasMember("mustbe"))){ 
					cerr << "Header error mustbe object not definer => error in upper bound   " << endl;
					exit;
				} else if ( ! ( doc["problem"]["mustbe"].IsString() ) ) {

					cerr << "Header error mustbe definition is empty" << endl;
					exit;

				} 
	// verify header
	regex ub_valid(regexp_map["ub_valid"]);
	smatch m_ub;
	string mustbe = doc["problem"]["mustbe"].GetString();

	if (! (regex_search(mustbe, m_ub, ub_valid))) {
		cerr << "current mustbe value is " << mustbe << endl;
		cerr << "Mustbe regexp error please check  ub_valid regexp and json mustbe value (minimisation, maximisation and reel value of ub ) " << endl;
		exit;


	} 
	if( !  ( doc.HasMember("variables") ) ) {

		cerr << "json format error variables section  not detected " << endl;
		exit;

	}
	if (  ! (doc.HasMember("functions"))) 
	{ 

		cerr << "json format error  function section not detected" << endl;
		exit;

	} 

#ifdef DEBUG
	assert(doc.IsObject());
	assert(doc.HasMember("problem"));
	assert(doc.HasMember("variables"));
	cout "error in functions:" 
	assert(doc.HasMember("functions"));

	// verify header
	assert(doc["problem"].IsObject());
	assert(doc["problem"].HasMember("name"));
	assert(doc["problem"].HasMember("mustbe"));

	assert(doc["problem"]["name"].IsString());
	assert(doc["problem"]["mustbe"].IsString());
	assert(doc["variables"].IsObject());

#endif

	//regex m_valid("<-*[0-9]*\\.[0-9]*");


	// verify variables
	
	

	const Value& variables = doc["variables"];

	//regex n_valid("M*[A-Z]_[0-Z]_[0-9]*");  // M for merge
	//regex r_valid("[A-Z]*_[0-9]*");
	// variable regexp
	regex var_valid(regexp_map["var_merged_valid"]);
	// value regexp
	regex val_valid(regexp_map["val_valid"]);

	string var_name;
	// res are value in general
	string res_name;
	// verify if variable name  and values name are valid 
	for(Value::ConstMemberIterator i = variables.MemberBegin(); i < variables.MemberEnd(); i++){
		var_name = i->name.GetString();
		assert(regex_search(var_name,m_ub,var_valid));
		const Value& var = variables[i->name.GetString()];
		assert(var.IsArray());
		for(Value::ConstValueIterator j = var.Begin(); j < var.End(); j++){
			assert(j->IsString());
			res_name = j->GetString();
			if ( verbose > 3 ) cout << "DEBUG value label :" << res_name << endl;

		    // check variable label grammar => cf val_valid regexp
			if  (!(regex_search(res_name,m_ub,val_valid) )) { 
			cerr << " check format ERROR " << "variable :  "<< res_name << "not compliant with regexp : " <<  regexp_map["val_valid"]  << endl;
			exit(1);
			}
		#ifdef DEGUG
			assert(regex_search(res_name,m_ub,val_valid));
		#endif
		}
	}

	// verify functions
	assert(doc["functions"].IsObject());
	if ( doc["functions"].IsObject() ) {
		const Value& fns = cfn["functions"];
		//regex f_valid("(C*B[0-9]*-[0-9]*)|(C*M*S[0-9]*)|(C*M*Etmp)|(C*M[0-9]*)");
		regex f_valid(regexp_map["f_valid"].c_str());

		string fn_name;
		string t_name;
		for(Value::ConstMemberIterator i = fns.MemberBegin(); i < fns.MemberEnd(); i++){
			fn_name = i->name.GetString();
#ifdef DEBUG
			assert(regex_search(fn_name,m_ub,f_valid));
#endif
			if(!(regex_search(fn_name,m_ub,f_valid))) 
			{
				cerr << "Error in : " << fn_name << " function name not compliante with regular expression : please chech f_valid regexp : " << regexp_map["f_valid"] << endl;
				exit(1);
			}

			const Value& fn = fns[i->name.GetString()];
			assert(fn.IsObject());
			for(Value::ConstMemberIterator t = fn.MemberBegin(); t < fn.MemberEnd(); t++){
				t_name = t->name.GetString();
				if(t_name=="costs"){
					const Value& cost = fn[t->name.GetString()];
					assert(cost.IsArray());
					for(Value::ConstValueIterator j = cost.Begin(); j < cost.End(); j++){
						assert(j->IsFloat()||j->IsInt()||j->IsDouble());
					}
				}
				else if(t_name=="scope"){
					const Value& scope = fn[t->name.GetString()];

					if( verbose > 1 ) {
						cout << "fonction : " << fn_name <<  " => " ;
						for (SizeType i = 0; i < scope.Size(); i++) // Uses SizeType instead of size_t
						{
							//david
							//		  printf("scope[%d] = %d\n", i, scope[i].GetInt());
							if(scope[i].IsInt() ) { cout << "scope are : "<< scope[i].GetString() << " is a var index =>  cfn format not compliant with cfntools scope must be var label" << endl ; exit (1) ; }
							printf("scope[%d] = %s , ", i, scope[i].GetString());
						}
						cout << endl;
					}
#ifdef DEBUG
					cout << " label 1 "<< endl;
					assert(scope.IsArray());
#endif
					if( scope.IsArray() ) {
						for(Value::ConstValueIterator j = scope.Begin(); j < scope.End(); j++){


#ifdef DEBUG
							cout << " label 2 "<< endl;
							assert(j->IsString());
#endif
							if( j-> IsString()) { 
								var_name = j->GetString();

								if( verbose >10 ) { cout << "scope is string " << j->GetString() << endl ; }

							} else if ( j-> IsInt() )  {
								if ( verbose > 10  ) { cout << "scope is index" << j->GetInt() << endl ; }
								cerr << "SCOPE DEF ERROR in " << t_name << "  in scope definition: " << j ;
								cerr <<   " is int => index  in scope  are not compliant with the cfntools events it works with toulbar2 " << endl ;
								exit(1) ; 
							} 

							if( ! ( regex_search(var_name,m_ub,var_valid)  )){

								cerr << "var format ERROR" <<  var_name << " doesn\'t match with espected regexp " << & m_ub << endl ;
								exit(1);


							}
#ifdef DEBUG
							cout << " label 3 "<< endl;
							assert(regex_search(var_name,m_ub,var_valid));
#endif
						}
					} else { cerr << "ERROR in "<< t_name << " scope definition is not an erray " << endl;  exit(1); }
				}
				else {
					cout<<t_name<<endl;
					cerr<<"function contains a variable different from 'scope' and 'costs' : error"<<endl;
					exit(1);
				}
			}
		}
	} else { cerr << "format ERROR : section function not found " << endl ; exit(1) ; }// check function
}

/// check if all costs have the same number of digits, if not, homogeneize
int CFNManager::check_digits(int digits){
	string mustbe = cfn["problem"]["mustbe"].GetString();
	cout << "MUSTBE read : " << cfn["problem"]["mustbe"].GetString() << endl;

	int d_mustbe = mustbe.length() - (mustbe.find(".")+1);
	if(digits==0){
		cout<<"all cost will have "<<d_mustbe<<" digits"<<endl;
		digits = d_mustbe; // initialize the value;
	}
	else{
		if(d_mustbe > digits){
			mustbe = mustbe.substr(0, (mustbe.find(".")+1)+digits);
		}
		else if(d_mustbe < digits){
			for(int k=0; k<digits-d_mustbe; k++){
				mustbe += "0";
			}
		}
		cfn["problem"]["mustbe"].SetString(StringRef(mustbe.c_str()), cfn.GetAllocator());
	}
	return digits;
}

/// filering value above a threshold
void CFNManager::filter_value(double threshold){ // NP OK
	// take the UB value
	string sUB = cfn["problem"]["mustbe"].GetString();
	sUB = sUB.substr(1, sUB.length());
	double UB;
	stringstream(sUB) >> UB; // allow convert even if negative value
	// check all the costs
	Value& fns = cfn["functions"];
	double c;
	unsigned long total_nb_tresh =0;
	cout << "Begin values filering :\n -------\n " << endl;
	for(Value::ConstMemberIterator i = fns.MemberBegin(); i < fns.MemberEnd(); i++){
		Value& fn = fns[i->name.GetString()];
		Value& cost = fn["costs"];
		unsigned int current_tresh =0 ;
		
		for(SizeType j = 0; j < cost.Size(); j++){
			c = cost[j].GetDouble();
			if(c > threshold){
				cost[SizeType(j)].SetFloat(UB);
				total_nb_tresh++;
				current_tresh++;
			}
		}
		if( current_tresh > 0 ) {
		cout << " function: " << i->name.GetString() << ": " << current_tresh << " values > to "  <<  threshold << endl;
		}
	}
		cout << " ------\n Total #value removed : "  <<  total_nb_tresh << endl;
}

/// merging two documents
/// two pdb with a same protein but a different backbone for instance
/// TO DO
void CFNManager::merge(string& cfn_file, bool fast ,  std::map<std::string, string> & regexp_map, int   verbose ){
	// loading the second file :
	Document to_merge;
	FILE* f2;
	f2 = fopen(cfn_file.c_str(), "r");
	if (f2!= NULL){
		ifstream in(cfn_file.c_str());
		string contents((istreambuf_iterator<char>(in)), istreambuf_iterator<char>());
		fclose(f2);
		const char* input = contents.c_str();
		to_merge.Parse(input);
	}
	else{
		cerr<<"cannot open file "<<cfn_file<<endl;
		exit(1);
	}
	if(!fast) check_format(to_merge,regexp_map,verbose);

	// change the header
	string name1 = cfn["problem"]["name"].GetString();
	string name2 = to_merge["problem"]["name"].GetString();
	string new_name = "merge-" + name1 + "-" + name2;
	cfn["problem"]["name"].SetString(StringRef(new_name.c_str()), cfn.GetAllocator());

	string sUB = cfn["problem"]["mustbe"].GetString();
	sUB = sUB.substr(1, sUB.length());
	double UB = stof(sUB);
	string sUB2 = to_merge["problem"]["mustbe"].GetString();
	sUB2 = sUB2.substr(1, sUB2.length());
	double UB2 = stof(sUB2);

	string new_mustbe = "<" + to_string(UB+UB2);
	cfn["problem"]["mustbe"].SetString(StringRef(new_mustbe.c_str()), cfn.GetAllocator());


	// for each variable of merge : change name, add, add cost function
	Value& m_variables = to_merge["variables"];
	string name = "M";

	for(Value::MemberIterator v = m_variables.MemberBegin(); v < m_variables.MemberEnd(); v++){
		name += v->name.GetString();
		Value& m_var = m_variables[v->name.GetString()];

		cfn["variables"].AddMember(Value(name.c_str(), cfn.GetAllocator()).Move(),Value(m_var, cfn.GetAllocator()).Move(), cfn.GetAllocator());
		Value l_costs(kArrayType);
		const Value& mvar = cfn["variables"][name.c_str()];
		// find the variable that correspond if exists and add the cost function
		if(cfn["variables"].HasMember(v->name.GetString())){
			Value& var = cfn["variables"][v->name.GetString()];
			for(SizeType j = 0; j < var.Size(); j++){
				for(SizeType k = 0; k < mvar.Size(); k++){
					string rot1 = var[j].GetString();
					string rot2 = mvar[k].GetString();
					rot1 = rot1.substr(0,1);
					rot2 = rot2.substr(0,1);
					if(rot1==rot2){  // same aa
						l_costs.PushBack(0, cfn.GetAllocator());
					}
					else{
						l_costs.PushBack(UB+UB2, cfn.GetAllocator());
					}
				}
			}
			Value l_scope(kArrayType);
			l_scope.PushBack(Value(name.c_str(), cfn.GetAllocator()).Move(), cfn.GetAllocator());
			l_scope.PushBack(Value(v->name.GetString(), cfn.GetAllocator()).Move(), cfn.GetAllocator());

			Value new_fn(kObjectType);

			new_fn.AddMember("scope", Value(l_scope, cfn.GetAllocator()).Move(), cfn.GetAllocator());
			new_fn.AddMember("costs", Value(l_costs, cfn.GetAllocator()).Move(), cfn.GetAllocator());

			string fn_name = "M" + name.substr(5,name.size()-4); // NP : add this format to check

			cfn["functions"].AddMember(Value(fn_name.c_str(), cfn.GetAllocator()).Move(), Value(new_fn, cfn.GetAllocator()).Move(), cfn.GetAllocator());
			name = "M";
		}
	}

	// add functions
	Value& m_functions = to_merge["functions"];
	name = "M";
	for(Value::MemberIterator f = m_functions.MemberBegin(); f < m_functions.MemberEnd(); f++){
		name += f->name.GetString();
		Value& mfn = m_functions[f->name.GetString()];
		cfn["functions"].AddMember(Value(name.c_str(), cfn.GetAllocator()).Move(), Value(mfn, cfn.GetAllocator()).Move(), cfn.GetAllocator());
		name = "M";
	}
}

void CFNManager::updateUB(double new_UB, int verbose){
	// take the UB value
	string sUBold = cfn["problem"]["mustbe"].GetString();
	sUBold = sUBold.substr(1, sUBold.length());
	double oldUB;
	stringstream(sUBold) >> oldUB; // allow convert even if negative value

	// look for all the values equal to UB+1
	map<string, vector<int>> equalUB;

	double UB = 0;
	Value& functions = cfn["functions"];
	for(Value::ConstMemberIterator i = functions.MemberBegin(); i < functions.MemberEnd(); i++){
		equalUB[i->name.GetString()] = vector<int>{};
		Value& fn = functions[i->name.GetString()];
		Value& cost = fn["costs"];
		double c_max = 0;
		for(SizeType j = 0; j < cost.Size(); j++){
			if(cost[j].GetDouble() > c_max && abs(oldUB+1-cost[j].GetDouble())>=0.01){ // max and not impossible value
				// (we consider that impossible values have been set to UB+1)
				c_max = cost[j].GetDouble();
			}
			else if(abs(oldUB+1-cost[j].GetDouble())<0.01){ // this is an impossible value
				equalUB[i->name.GetString()].push_back(j);
			}
		}
		if( verbose > 2 ) {
			cout<<"DEBUG>> "<<i->name.GetString()<<" : "<<c_max<<endl;
		}
		UB += c_max;
	}

	if(new_UB!=0) UB = new_UB+0.01; // +0.01 to allow to have a solution equal to UB in order to always have a solution

	// replace values equal to oldUB+1 by newUB+1
	/*for(Value::ConstMemberIterator i = functions.MemberBegin(); i < functions.MemberEnd(); i++){
	  Value& fn = functions[i->name.GetString()];
	  Value& cost = fn["costs"];
	  for(unsigned int k=0; k!=equalUB[i->name.GetString()].size(); k++){
	//cout<<"DEBUG>> in "<<i->name.GetString()<<", cost "<<equalUB[i->name.GetString()][k]<<" set to impossible"<<endl;
	cost[equalUB[i->name.GetString()][k]].SetDouble(UB+1);
	}
	}*/

	string s_new;
	if(new_UB != 0){
		cout<<"new UB set to "<<new_UB<<endl;
		s_new = to_string(new_UB);
	}
	else{
		cout<<"new UB set to "<<UB<<" By Max sum of each cost function "<<endl;
		s_new = to_string(UB);
	}
		cout<<"UB>> oldUB: "<<oldUB<<endl;
		cout<<"UB>> newUB: "<<s_new<<endl;
	string sUB = cfn["problem"]["mustbe"].GetString();
	sUB = sUB.replace(1, s_new.length(), s_new);
	if(verbose != 0) cout<<"new ub: "<<sUB<<endl;
	cfn["problem"]["mustbe"].SetString(StringRef(sUB.c_str()), cfn.GetAllocator());
}

void CFNManager::extract(string& var_list, int verbose){ // extract the subproblem corresponding to a list of variables
	// NP :
	// on garde Etmp ?
	// /!\ le numero des fn ne correspond pas au numero des variables impliquées (dans 4dd5.cfn)
	//////

	// initialization of the vector
	vector<string> list_var;
	// loading the file
	ifstream list;
	list.open(var_list);
	if (!list) {
		cerr << "Unable to open file";
		exit(1);   // call system to stop
	}
	// verify format
	regex n_valid("[A-Z]_[0-Z]_[0-9]*");
	smatch m;

	while (true){
		string s_var;
		getline(list, s_var);
		if (!list.good()) { // Break on error or end of file
			if (s_var.empty()) list.clear(); // No error if end of file
			else list >> s_var; // (Dummy) error is line not empty
			break;
		}
		if(regex_search(s_var,m,n_valid)){
			list_var.push_back(s_var);
		}
		else{
			cerr<<s_var<<": not considered (incorrect format)"<<endl;
		}
	}
	list.close();

	// keep only the variables with scope ok
	Value& variables = cfn["variables"];
	for(Value::MemberIterator i = variables.MemberBegin(); i < variables.MemberEnd(); i++){
		if(find(list_var.begin(), list_var.end(), i->name.GetString()) == list_var.end()){
			variables.RemoveMember(i);
			i --;
		}
	}

	// keep only the functions with a scope ok
	Value& functions = cfn["functions"];
	for(Value::MemberIterator i = functions.MemberBegin(); i < functions.MemberEnd(); i++){
		Value& fn = functions[i->name.GetString()];
		Value& scope = fn["scope"];
		for(SizeType j = 0; j < scope.Size(); j++){
			if(find(list_var.begin(), list_var.end(), scope[j].GetString()) == list_var.end()){
				functions.RemoveMember(i);
				i --;
				break;
			}
		}
	}

	// rename header and update bound
	string name = cfn["problem"]["name"].GetString();
	name += "_extracted_";
	name += var_list.substr(var_list.find_last_of("/\\")+1, var_list.size() - var_list.find_last_of("."));
	cfn["problem"]["name"].SetString(StringRef(name.c_str()), cfn.GetAllocator());

	updateUB(0);
}

void CFNManager::extract_cost(string& sol_file, int verbose){
	// prec is choosen according to UB prec
	string mustbe = cfn["problem"]["mustbe"].GetString();
	int d_mustbe = mustbe.length() - (mustbe.find(".")+1);
	cout<<"prec to extract cost: "<<d_mustbe<<endl;

	//open file sol
	if(verbose!=0) cout<<"open solution file"<<endl;
	ifstream sol;
	sol.open(sol_file);
	if (!sol) {
		cerr << "Unable to open file: "<<sol_file<<endl;
		exit(1);   // call system to stop
	}

	queue<string> l_val;
	while (true){
		string line_sol;
		getline(sol, line_sol);
		if (!sol.good()) { // Break on error or end of file
			if (line_sol.empty()) sol.clear(); // No error if end of file
			else sol >> line_sol; // (Dummy) error is line not empty
			break;
		}
		istringstream iss_line(line_sol);
		vector<string> l_sol(istream_iterator<string>{iss_line},
				istream_iterator<string>());
		for(vector<string>::iterator it=l_sol.begin(); it!=l_sol.end(); it++){
			l_val.push((*it));
		}
	}
	sol.close();
	if(verbose!=0) cout<<"solution read"<<endl;
	map<string, int> idx_sol;
	map<pair<string,string>, double> unary_costs;
	map<pair<string, string>, map<string, double>> binary_costs;
	double tmp_cost=0.;

	Value& l_var = cfn["variables"];
	for(Value::ConstMemberIterator v = l_var.MemberBegin(); v < l_var.MemberEnd(); v++){
		Value& var = l_var[v->name.GetString()];
		string name_sol = l_val.front();
		l_val.pop();
		int idx=0;
		while(strcmp(var[idx].GetString(),name_sol.c_str())!=0){
			++idx;
		}
		idx_sol[v->name.GetString()] = idx;
	}
	if(verbose!=0) cout<<"index of solutions found"<<endl;
	Value& func = cfn["functions"];
	for(Value::MemberIterator f=func.MemberBegin(); f!=func.MemberEnd(); f++){
		Value& scope = func[f->name.GetString()]["scope"];
		if(scope.Size() == 1){
			int idx = idx_sol[scope[0].GetString()];
			Value& cost = func[f->name.GetString()]["costs"];
			unary_costs[pair<string, string>{f->name.GetString(),scope[0].GetString()}] = cost[idx].GetDouble();
		}
		else if(scope.Size() == 2){
			int idx = idx_sol[scope[0].GetString()];
			int idx2 = idx_sol[scope[1].GetString()];
			Value& cost = func[f->name.GetString()]["costs"];
			// nb variables of res2 ?
			Value& l_var = cfn["variables"][scope[1].GetString()];
			int nb_rotam2 = l_var.Size();
			int final_idx = idx*nb_rotam2+idx2;
			binary_costs[pair<string, string>{f->name.GetString(),scope[0].GetString()}][scope[1].GetString()] = cost[final_idx].GetDouble();
		}
		else if(scope.Size() == 0){
			Value& cost = func[f->name.GetString()]["costs"];
			tmp_cost += cost[0].GetDouble();
		}
	}
	if(verbose!=0) cout<<"cost of solutions found"<<endl;
	string name_csv = sol_file+"_energies.csv";
	ofstream file;
	file.open(name_csv);
	if (!file) {
		cerr << "Unable to open file "<<name_csv<<endl;
		exit(1);   // call system to stop
	}
	double final_cost = tmp_cost;
	// write tmp
	file<<fixed<<setprecision(d_mustbe)<<"Etmp;;;"<<tmp_cost<<";"<<endl;
	//write unary
	for(map<pair<string, string>, double>::iterator it=unary_costs.begin(); it!=unary_costs.end(); it++){
		file<<fixed<<setprecision(d_mustbe)<<(*it).first.first<<";"<<(*it).first.second<<";;"<<(*it).second<<";"<<idx_sol[(*it).first.second]<<endl;
		final_cost += (*it).second;
	}
	// write binary
	for(map<pair<string, string>, map<string, double>>::iterator it=binary_costs.begin(); it!=binary_costs.end(); it++){
		pair<string, string> r1 = (*it).first;
		for(map<string, double>::iterator it2=(*it).second.begin(); it2!=(*it).second.end(); it2++){
			file<<fixed<<setprecision(d_mustbe)<<r1.first<<";"<<r1.second<<";"<<(*it2).first<<";"<<(*it2).second<<";("<<idx_sol[r1.second]<<","<<idx_sol[(*it2).first]<<")"<<endl;
			final_cost += (*it2).second;
		}
	}
	file.close();
	if(verbose!=0) cout<<"cost written in "<<name_csv<<endl;
	cout<<fixed<<setprecision(d_mustbe)<<"CFNTOOLS>> final cost = "<<final_cost<<endl; // pour le moment 6, on verra après
}

void CFNManager::add_map(string& map_file, double ps, double pp, bool fast, int verbose,  std::map<std::string, string> & regexp_map){
	// load the map
	CFNManager map(map_file);
	// verify format
	if(!fast) check_format(map.cfn,regexp_map,verbose);

	// compute ks and kp
	double ks=1, kp=1;
	if(ps != -1 || pp != -1){
		double meanself, meanpair, medself, medpair;
		//cout<<"DEBUG>> computing mean"<<endl;
		compute_meandelta(meanself, meanpair, medself, medpair, verbose);
		//cout<<"DEBUG>> mean computed"<<endl;
		double meanselfmap, meanpairmap, medselfmap, medpairmap;
		map.compute_meandelta(meanselfmap, meanpairmap, medselfmap, medpairmap, verbose);
		if(ps != -1 && meanselfmap!=0.0){
			ks = ps*medself/medselfmap;
			if(verbose!=0) cout<<"DEBUG>> in add_map, ks="<<ps<<"*"<<medself<<"/"<<medselfmap<<"="<<ks<<endl;
		}
		if(pp != -1 && meanpairmap!=0.0){
			kp = pp*medpair/medpairmap;
			if(verbose!=0) cout<<"DEBUG>> in add_map, kp="<<pp<<"*"<<medpair<<"*/"<<medpairmap<<"="<<kp<<endl;
		}
	}


	// check variable names and value :
	// we consider that if map has a variable,
	// it is also in cfn and it has the same possible values
	// i.e. map can have less variables than cfn
	Value& map_var = map.cfn["variables"];
	for(Value::ConstMemberIterator v = map_var.MemberBegin(); v < map_var.MemberEnd(); v++){
		Value& var = map_var[v->name.GetString()];
		if(verbose> 2) cout<<"DEBUG>> "<<v->name.GetString()<<endl;
		assert(cfn["variables"].HasMember(v->name.GetString()));
		assert(cfn["variables"][v->name.GetString()].Size() == var.Size());
		for(SizeType i = 0; i<var.Size(); i++){
			string rot_name = var[i].GetString();
			if(rot_name.substr(0,2)!="Hd" && rot_name.substr(0,2)!="He" && rot_name.substr(0,2)!="Hp" && rot_name.substr(0,2)!="H_"){ // not his
				assert(strcmp(cfn["variables"][v->name.GetString()][i].GetString(),rot_name.c_str())==0);
			}
			else{ // his, we check if cfn=his too
				string rot_name_cfn = cfn["variables"][v->name.GetString()][i].GetString();
				assert(rot_name_cfn.substr(0,2)=="Hd" || rot_name_cfn.substr(0,2)=="He" || rot_name.substr(0,2)=="Hp" || rot_name.substr(0,2)=="H_");
				string idx = rot_name.substr(rot_name.find('_')+1, rot_name.size());
				string idx_cfn = rot_name_cfn.substr(rot_name_cfn.find('_')+1, rot_name_cfn.size());
				assert(strcmp(idx.c_str(),idx_cfn.c_str())==0);
				// if ok, set final name to rot_name
				cfn["variables"][v->name.GetString()][i] = Value(rot_name.c_str(), cfn.GetAllocator()).Move(), cfn.GetAllocator();
			}

		}
	}

	if(ks!=1 || kp!=1){
		if(verbose != 0) cout<<"multiply map values by "<<ks<<"(unary) and "<<kp<<"(binary)"<<endl;
		map.multiply(ks, kp);
	}

	if(verbose != 0) cout<<"adding functions"<<endl;
	// add functions
	Value& map_fn = map.cfn["functions"];
	for(Value::ConstMemberIterator f = map_fn.MemberBegin(); f < map_fn.MemberEnd(); f++){
		string name_fn = "C";
		name_fn += f->name.GetString();
		Value& fn = map_fn[f->name.GetString()];
		cfn["functions"].AddMember(Value(name_fn.c_str(), cfn.GetAllocator()).Move(), Value(fn, cfn.GetAllocator()).Move(), cfn.GetAllocator());
		// need to do Value(fn, cfn.GetAllocator()).Move() to deepcopy the value, otherwise it doesn't exist when saving file
	}
	// change header
	string new_name = cfn["problem"]["name"].GetString();
	new_name += "_cryo_";
	new_name += map.cfn["problem"]["name"].GetString();
	if(verbose != 0) cout<<new_name<<endl;
	cfn["problem"]["name"].SetString(StringRef(new_name.c_str()), cfn.GetAllocator());
	updateUB(0, verbose);
	// to verify
	if(!fast) check_format(cfn, regexp_map,verbose);

}

void CFNManager::comp(string& cfn_file, bool fast ,  std::map<std::string, string> & regexp_map, int verbose){
	double prec; // prec is choosen according to UB prec
	string mustbe = cfn["problem"]["mustbe"].GetString();
	int d_mustbe = mustbe.length() - (mustbe.find(".")+1);
	prec = 1/pow(10, d_mustbe-1);
	cout<<"prec to compare : "<<prec<<endl;

	// load the second cfn
	Document cfn2;
	FILE* fc;
	fc = fopen(cfn_file.c_str(), "r");  // use "rb" on Windows
	if (fc!= NULL){
		ifstream in(cfn_file.c_str());
		string contents((istreambuf_iterator<char>(in)), istreambuf_iterator<char>());
		fclose(fc);
		const char* input = contents.c_str();
		cfn2.Parse(input);
	}
	else{
		cerr<<"cannot open file "<<cfn_file<<endl;
		exit(1);
	}
	// verify format
	if(!fast) check_format(cfn2,   regexp_map, verbose);

	// variable checking
	// and value checking 

	Value& variables = cfn["variables"];

	Value& m_variables = cfn2["variables"];

	string name;
	if ( verbose > 1 )  { 
		cout << " var number in the first cfn : " << variables.MemberEnd() - variables.MemberBegin()  << endl;
		cout << "Var number in the second cfn :" << m_variables.MemberEnd() - m_variables.MemberBegin()  << endl;
		//	cout << "m_variables.MemberEnd ==> " << &(*m_variables.MemberEnd()) << endl;
		//m_variables.MemberBegin() << endl;
	}
	if(! ( m_variables.MemberEnd()-m_variables.MemberBegin()==variables.MemberEnd()-variables.MemberBegin()  )) {

		cerr << " >>> COMP ERROR #var different between the 2 cfn : first cfn includes " << variables.MemberEnd()-variables.MemberBegin() << " second one " << m_variables.MemberEnd()-m_variables.MemberBegin()  << endl;
		exit;
	} 
	// assert 

#ifdef DEBUG
	assert(m_variables.MemberEnd()-m_variables.MemberBegin()==variables.MemberEnd()-variables.MemberBegin()); // check size
#endif


	for(Value::MemberIterator v = m_variables.MemberBegin(); v < m_variables.MemberEnd(); v++){
		name = v->name.GetString();

#ifdef DEBUG
		assert(variables.HasMember(name.c_str()));
#endif
		Value& l_rot2 = m_variables[name.c_str()];
		Value& l_rot1 = variables[name.c_str()];
		if (!( variables.HasMember(name.c_str()))) { 
			// get  value name
			cerr << " error No var name " << name.c_str() << "in original cfn" << endl;
			exit(1);

		} else {  
			if ( verbose > 5 ) {
				cout << " val 1: " <<  &l_rot1   << endl ;
				cout << " val 2: " <<  &l_rot2   << endl ;
			}
		}

		for(SizeType r = 0; r<l_rot2.Size(); r++){
			assert(l_rot1[r]==l_rot2[r]);
		}
	}

	// add functions
	Value& m_functions = cfn2["functions"];
	Value& functions = cfn["functions"];
	assert(m_functions.MemberEnd()-m_functions.MemberBegin()==functions.MemberEnd()-functions.MemberBegin()); // check size
	for(Value::MemberIterator f = m_functions.MemberBegin(); f < m_functions.MemberEnd(); f++){
		name = f->name.GetString();
		assert(functions.HasMember(name.c_str()));

		Value& f_s = functions[name.c_str()]["scope"];
		Value& f_s2 = m_functions[name.c_str()]["scope"];

		Value& f_c = functions[name.c_str()]["costs"];
		Value& f_c2 = m_functions[name.c_str()]["costs"];

		for(SizeType r = 0; r<f_s2.Size(); r++){
			if (!(f_s[r].GetString(),f_s2[r].GetString()))   { 
				cerr << " var " <<  f_s[r].GetString() << " is diff from " << f_s2[r].GetString() << endl; 
			}
#ifdef DEBUG
			assert(strcmp(f_s[r].GetString(),f_s2[r].GetString())==0);
#endif
		}
		for(SizeType r = 0; r<f_c2.Size(); r++){ // cost aren't okay at this moment

			if ( !(abs(f_c[r].GetDouble()-f_c2[r].GetDouble())<prec) ) {  
				cout <<" >> COMP ERROR  function  tuple different in  "<<  name << "  => index " << r ;
				cout <<  " cost 1 =" << f_c[r].GetDouble() << " cost 2 ="<< f_c2[r].GetDouble() ;
				cout << " delta cost : " ;
				cout << abs(f_c[r].GetDouble()-f_c2[r].GetDouble());
				cout << " ( current precision : " << prec << ")" << endl ;
				exit(1);
			} 
			//	assert(abs(f_c[r].GetDouble()-f_c2[r].GetDouble())<prec);
		}
	}
}

void CFNManager::multiply(double m){
	Value& func = cfn["functions"];
	for(Value::MemberIterator f=func.MemberBegin(); f!=func.MemberEnd(); f++){
		Value& cost = func[f->name.GetString()]["costs"];
		for(SizeType c=0; c!=cost.Size(); c++){
			cost[c].SetDouble(m*cost[c].GetDouble());
		}
	}
	updateUB(0, 0);
}

void CFNManager::multiply(double mself, double mpair){
	Value& func = cfn["functions"];
	for(Value::MemberIterator f=func.MemberBegin(); f!=func.MemberEnd(); f++){
		Value& scope = func[f->name.GetString()]["scope"];
		if(scope.Size() == 1 && mself != 1){
			Value& cost = func[f->name.GetString()]["costs"];
			for(SizeType c=0; c!=cost.Size(); c++){
				cost[c].SetDouble(mself*cost[c].GetDouble());
			}
		}
		else if(scope.Size() == 2 && mpair != 1){
			Value& cost = func[f->name.GetString()]["costs"];
			for(SizeType c=0; c!=cost.Size(); c++){
				cost[c].SetDouble(mpair*cost[c].GetDouble());
			}
		}
	}
	updateUB(0, 0);
}

void CFNManager::compute_meandelta(double& meanself, double& meanpair, double& medself, double& medpair, int verbose){
	/// we don't consider costs = uB+1 -> impossible cost (to verif notation)
	// take the UB value
	string sUB = cfn["problem"]["mustbe"].GetString();
	sUB = sUB.substr(1, sUB.length());
	double UB;
	stringstream(sUB) >> UB; // allow convert even if negative value
	//cout<<fixed<<setprecision(6)<<"DEBUG>> UB is "<<UB<<endl;
	vector<double> deltaself, deltapair;
	vector<double> l_self, l_pair;
	meanself = 0;
	meanpair = 0;
	medself = 0;
	medpair = 0;
	Value& func = cfn["functions"];
	for(Value::MemberIterator f=func.MemberBegin(); f!=func.MemberEnd(); f++){
		Value& scope = func[f->name.GetString()]["scope"];
		if(scope.Size() == 1){
			Value& cost = func[f->name.GetString()]["costs"];
			if(cost.Size()!=1){
				double en = cost[0].GetDouble();
				SizeType i = 1;
				while(UB+1-en<=0.0 && i!=cost.Size()){
					en = cost[i].GetDouble();
					i += 1;
				}
				if(i==cost.Size()&&UB+1-en<=0.0) cout<<"WARNING! in function "<<f->name.GetString()<<" no possible values"<<endl;
				else{
					double min_en=en, max_en=en;
					l_self.push_back(en);
					for(SizeType c=i; c!=cost.Size(); c++){
						en = cost[c].GetDouble();
						if(en<min_en) min_en=en;
						if(en>max_en && UB+1-en>0.0) max_en=en;
						if(UB+1-en>0.0) {
							l_self.push_back(en);
						}
					}
					//cout<<"DEBUG>> "<<f->name.GetString()<<" delta "<<max_en-min_en<<endl;
					deltaself.push_back(max_en-min_en);
					meanself += max_en-min_en;
				}
			}
		}
		else if(scope.Size() == 2){
			Value& cost = func[f->name.GetString()]["costs"];
			if(cost.Size()!=1){
				double en = cost[0].GetDouble();
				SizeType i = 1;
				while(UB+1-en<=0.0 && i!=cost.Size()){
					en = cost[i].GetDouble();
					i += 1;
				}
				if(i==cost.Size()&&UB+1-en<=0.0) cout<<"WARNING! in function "<<f->name.GetString()<<" no possible values"<<endl;
				else{
					double min_en=en, max_en=en;
					l_pair.push_back(en);
					for(SizeType c=i; c!=cost.Size(); c++){
						en = cost[c].GetDouble();
						if(en<min_en) min_en=en;
						if(en>max_en && UB+1-en>0.0) max_en=en;
						if(UB+1-en>0.0) {
							l_pair.push_back(en);
						}
					}
					deltapair.push_back(max_en-min_en);
					meanpair += max_en-min_en;
					if(max_en - min_en > 100000){
						cout<<"DEBUG>> "<<f->name.GetString()<<": "<<max_en<<" - "<<min_en<<endl;
					}
				}
			}
		}
	}
	if(verbose!=0) cout<<"size unary costs "<<deltaself.size()<<", sum "<<meanself<<endl;
	if(verbose!=0) cout<<"size binary costs "<<deltapair.size()<<", sum "<<meanpair<<endl;
	if(verbose!=0){
		if(deltaself.size()!=0){
			ofstream file;
			file.open("deltaself_"+name_cfn+".csv");
			if (!file) {
				cerr << "Unable to open file deltaself_"+name_cfn+".csv"<<endl;
				exit(1);   // call system to stop
			}
			for(unsigned int i=0; i!=deltaself.size(); i++){
				file<<deltaself.at(i)<<";"<<endl;
			}
			file.close();
			cout<<"list of deltaself written in deltaself_"+name_cfn+".csv "<<endl;
		}
		if(deltapair.size()!=0){
			ofstream file;
			file.open("deltapair_"+name_cfn+".csv");
			if (!file) {
				cerr << "Unable to open file deltapair_"+name_cfn+".csv"<<endl;
				exit(1);   // call system to stop
			}
			for(unsigned int i=0; i!=deltapair.size(); i++){
				file<<deltapair.at(i)<<";"<<endl;
			}
			file.close();
			cout<<"list of deltapair written in deltapair_"+name_cfn+".csv "<<endl;
		}

		// write all energies to plot distribution
		if(l_self.size()!=0){
			ofstream file;
			file.open("lself_"+name_cfn+".csv");
			if (!file) {
				cerr << "Unable to open file lself_"+name_cfn+".csv"<<endl;
				exit(1);   // call system to stop
			}
			for(unsigned int i=0; i!=l_self.size(); i++){
				file<<l_self.at(i)<<";"<<endl;
			}
			file.close();
			cout<<"list of l_self written in lself_"+name_cfn+".csv "<<endl;
		}
		if(l_pair.size()!=0){
			ofstream file;
			file.open("lpair_"+name_cfn+".csv");
			if (!file) {
				cerr << "Unable to open file lpair_"+name_cfn+".csv"<<endl;
				exit(1);   // call system to stop
			}
			for(unsigned int i=0; i!=l_pair.size(); i++){
				file<<l_pair.at(i)<<";"<<endl;
			}
			file.close();
			cout<<"list of l_pair written in lpair_"+name_cfn+".csv "<<endl;
		}

	}
	if(deltaself.size()!=0){
		size_t n = deltaself.size() / 2;
		nth_element(deltaself.begin(), deltaself.begin()+n, deltaself.end()); // sort the n first elements (avoid to sort all)
		medself = deltaself[n];
	}
	if(deltapair.size()!=0){
		size_t n = deltapair.size() / 2;
		nth_element(deltapair.begin(), deltapair.begin()+n, deltapair.end()); // sort the n first elements (avoid to sort all)
		medpair = deltapair[n];
	}
	if(deltaself.size()!=0) meanself = meanself/deltaself.size();
	if(deltapair.size()!=0) meanpair = meanpair/deltapair.size();
}

void CFNManager::infos(){
	double sizepb = 0.;
	long function_counter=0;
	long c0_function =0;
	long Binary_function=0;
	long Sup_binary_function=0;
	long nbVar = 0 ;
	long Mindom = +2147483647;
	long Maxdom = 0 ;
	std::vector<int> Scope_Hist; 
	Value& variables = cfn["variables"];

		// add CO constraint Etmp ... with no scope definition
		// putative scope Nbar+1
		Scope_Hist.push_back(-1);

	for(Value::MemberIterator v = variables.MemberBegin(); v < variables.MemberEnd(); v++){
		string name = v->name.GetString();
		Value& l_rot = variables[name.c_str()];
		cout<<  v->name.GetString()<<" has "<<l_rot.Size()<<" values"<<endl;
		sizepb += log10(l_rot.Size());
		Scope_Hist.push_back(-1);
		// min domaine size
		Mindom= ( Mindom > l_rot.Size()  ) ? l_rot.Size(): Mindom;
		// max domain size
		Maxdom= ( Maxdom < l_rot.Size()  ) ? l_rot.Size(): Maxdom;
		nbVar++;
	}

	string sUB = cfn["problem"]["mustbe"].GetString();

	cout << "---------------" << endl;
	cout << "Ub : " << sUB << endl;
	cout<< ">> Nb Var = " << nbVar << endl;
	cout<< ">> Max domaine size  = " << Maxdom << endl;
	cout<< ">> Min domaine size  = " << Mindom << endl;
	Value& func = cfn["functions"];
	for(Value::MemberIterator f=func.MemberBegin(); f<func.MemberEnd(); f++){
		function_counter++;
		Value& scope = func[f->name.GetString()]["scope"];
		int index = 0;
		index = ( scope.Size() < 1 ) ? 0: scope.Size();
		if ( Scope_Hist[index] < 0 )  { Scope_Hist[index] = 1 ; } else { Scope_Hist[index]++; } 


	}
	cout << "-------------- \n Total function number  :"  << function_counter << "\n ----------------"<< endl;
	for (std::vector<int>::iterator it = Scope_Hist.begin() ; it < Scope_Hist.end(); it++) { 

		if ( *it > 0 ) { 
			int arity = std::distance(Scope_Hist.begin(), it);
			std::cout <<" function arity : "<< arity << " = " << *it;
			std::cout << '\n';
		}

	}

	cout << "--------------------\n" << endl;
	cout<<">> Search space size: 10^"<<sizepb<<endl;
}
