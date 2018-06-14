#include "tb2scpbranch.hpp"
#include "tb2cpd.hpp"
#include "core/tb2types.hpp"

// // Tb2ScpBranch::Tb2ScpBranch(const char *filename)
// {
//   read_solution(filename);
//   vector< vector<char> > & rot2aa = ToulBar2::cpd->getRotamers2AA();
//   vector< SortItem > items;
//   for(size_t i=0; i<rot2aa.size(); i++)
//     {
//       items.clear();
//       for(size_t j=0; j<rot2aa[i].size();j++)
//         {
//           SortItem it;
//           tie(it.h,it.i) = AA2Criterium(rot2aa[i][j]);
//           items.push_back(it);
//         }
//       sort_criterium.push_back(items);
//     }
//   for(size_t i=0;i<sort_criterium.size();i++)
//     {
//       if (initial_sol.size()<=i)
//         {
//           cerr << "Invalid solution file " << filename << endl;
//           exit(EXIT_FAILURE);
//         }
//       ub_reference.push_back(sort_criterium[i][initial_sol[i]]);
//     }
// }

// void Tb2ScpBranch::read_solution(const char *filename)
// {
//   wcsp->propagate();

//   // open the file
//   ifstream file(filename);
//   if (!file) {
//     cerr << "Solution file " << filename << " not found!" << endl;
//     exit(EXIT_FAILURE);
//   }

//   while (!file.eof()) {
//     if ((unsigned int) i >= wcsp->numberOfVariables()) break;
//     Value value = 0;
//     file >> value;
//     if (!file) break;
//     initial_sol.push_back(value);
//   }
// }

// tuple<double,int> Tb2ScpBranch::AA2Criterium(char c)
// {

//   switch(c)
//     {
//     case 'A': return make_tuple(1.6,1); break;
//     case 'R': return make_tuple(-12.3,2); break;
//     case 'N': return make_tuple(-4.8,3); break;
//     case 'D': return make_tuple(-9.2,4); break;
//     case 'C': return make_tuple(2.0,5); break;
//     case 'E': return make_tuple(-8.2,6); break;
//     case 'F': return make_tuple(3.7,7); break;
//     case 'G': return make_tuple(1.0,8); break;
//     case 'H': return make_tuple(-3,9); break;
//     case 'I': return make_tuple(3.1,10); break;
//     case 'K': return make_tuple(-8.8,11); break;
//     case 'L': return make_tuple(2.8,12); break;
//     case 'M': return make_tuple(3.4,13); break;
//     case 'P': return make_tuple(-0.2,14); break;
//     case 'Q': return make_tuple(-4.1,15); break;
//     case 'S': return make_tuple(0.6,16); break;
//     case 'T': return make_tuple(1.2,17); break;
//     case 'V': return make_tuple(2.6,18); break;
//     case 'W': return make_tuple(1.9,19); break;
//     case 'Y': return make_tuple(-0.7,20); break;
//     default: std::cout << "Error, unknown amino acid: " << c << std:endl; exit(1);
//     }
// }

// // WARNING: NE MARCHE PAS POUR LE MOMENT
// tuple<size_t,size_t> Tb2ScpBranch::getBounds(size_t var_index)
// {
//   double min = 100;
//   int AA_index = 0;
//   size_t lb;
//   size_t ub;
//   for(size_t i=0; i<sort_criterium[var_index].size();i++)
//     {
//       double hdiff = sqrt(pow(sort_criterium[var_index][i].h-ub_reference[var_index],2));
//       if(hdiff < min)
//         {
//           AA_index = sort_criterium[var_index][i].i;
//           lb=i;
//           ub=i;
//           min=sort_criterium[var_index][i].h;
//         }
//       else if (hdiff == min && sort_criterium[var_index][i].i==AA_index)
//         {
//           ub++;
//         }
//     }
//   return make_tuple(lb,ub);
// }

tuple<size_t, size_t> Tb2ScpBranch::getBounds(int varIndex, Value value)
{
    return make_tuple(ToulBar2::cpd->getLeft(varIndex, value), ToulBar2::cpd->getRight(varIndex, value));
}

size_t Tb2ScpBranch::moveAAFirst(ValueCost* sorted, int domsize, Value left, Value right)
{
    ValueCost tmp;
    size_t cursor = 0;
    for (int i = 0; i < domsize; i++) {
        if (sorted[i].value >= left && sorted[i].value <= right) {
            tmp = sorted[cursor];
            sorted[cursor] = sorted[i];
            sorted[i] = tmp;
            cursor++;
        }
    }
    return cursor;
}

bool Tb2ScpBranch::multipleAA(int varIndex, Value* values, int domsize)
{
    char type = ToulBar2::cpd->getAA(varIndex, values[0]);
    for (int i = domsize - 1; i > 0; i--) {
        if (ToulBar2::cpd->getAA(varIndex, values[i]) != type)
            return true;
    }
    return false;
}

// void Tb2ScpBranch::keep(size_t begin, size_t end, size_t var_index)
// {
//   size_t domsize = sort_criterium[var_index].size();
//   if (begin)
//     sort_criterium[var_index].erase(0,begin-1);
//   if (domsize>end+1)
//     sort_criterium[var_index].erase(end+1,domsize);
// }

// void Tb2ScpBranch::remove(size_t begin, size_t end, size_t var_index)
// {
//   sort_criterium[var_index].erase(begin, end);
// }

// void Tb2ScpBranch::remove(size_t index, size_t var_index)
// {
//   sort_criterium[var_index].erase(index);
// }
