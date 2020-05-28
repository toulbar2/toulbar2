/**
 * tb2files.cpp
 *
 *  Created on: 24 mai 2019
 *      Author: kad abeldjilali
 */

#include <utils/tb2files_eps.hpp>
#include <iostream>


using namespace std;
/**
 * @brief If the file exists, it is read otherwise it is created with the value 30 x nb of cores
 * @param nbProcessFic
 * @param nbCores
 * @param nbProcPerCore
 */
long Tb2Files::nbProcess(const string nbProcessFic, const int nbCores, const int nbProcPerCore) {

                  	   ifstream file(nbProcessFic ,ios::in);
                  	   long nbProcess=0;

                  	       if(file.good())  // file exists, it is read
                  	       {
                  	           string line; // variable to memorize a line in the file
                  	           while(getline(file, line))
                  	           {
                  	        	   nbProcess= atoll(line.c_str());

                  	           }
                  	           file.close();

                  	       }
                  	       else  // file does not exists it is created and  written
                  	       {

                  	    	   ofstream file(nbProcessFic);
                  	    	   if(file)  // if ok
                  	    	   {
                  	    		 nbProcess= nbCores * nbProcPerCore;
                  	    		   file << nbProcess; //default value 30 * nbCores for now. todo : machine learning to compute optimal value
                  	    	   }
                  	    	   else
                  	    	   {

                  	    	      cout << "File error "<<  nbProcessFic << endl;
                  	    	   }
                  	    	   file.close();
                  	       }
                  	       return nbProcess;

}


void Tb2Files::write_file(const string fic, const string text)
{
 ofstream file(fic);
 if(file)  // if ok
    {

        file << text;
       // cout << "Text : "<< text << endl;
        cout<< "Text has been written in file "<< fic << endl;
    }
    else
    {
        cout << "File error "<<  fic << endl;
    }
    file.close();
}


void Tb2Files::read_file(const string fic)
{


    ifstream file(fic,ios::in);

    if(file)
    {
        string line; //Une variable pour stocker les lignes lues

        while(getline(file, line))
        {

        }
        cout <<endl;
        cout<< "File "<<fic<<" has been read : "<<endl;

    }
    else
    {
        cout << "ERREUR: no file" << endl;
    }
    return ;

}

vector<long> Tb2Files::readDomains(const string fic){
	ifstream file(fic,ios::in);
    string line;
	    if(file)
	    {

	        getline(file, line);
	        getline(file, line);//get second line with domains of vars

	    }
	    else
	    {
	        cout << "ERREUR: no file" << endl;
	    }

	// Vector of string to save tokens
	    vector <long> domain;

	    // stringstream class check1
	    stringstream check1(line);
	    string token;

	    // Tokenizing w.r.t. space ' '
	    while(getline(check1, token, ' '))
	        domain.push_back(stoi(token));

 return domain;
}

Tb2Files::Tb2Files() {
	// TODO Auto-generated constructor stub

}

Tb2Files::~Tb2Files() {
	// TODO Auto-generated destructor stub
}


std::string basename(const std::string &filename)
{
    if (filename.empty()) {
        return {};
    }

    auto len = filename.length();
    auto index = filename.find_last_of("/\\");

    if (index == std::string::npos) {
        return filename;
    }

    if (index + 1 >= len) {

        len--;
        index = filename.substr(0, len).find_last_of("/\\");

        if (len == 0) {
            return filename;
        }

        if (index == 0) {
            return filename.substr(1, len - 1);
        }

        if (index == std::string::npos) {
            return filename.substr(0, len);
        }

        return filename.substr(index + 1, len - index - 1);
    }

    return filename.substr(index + 1, len - index);
}

