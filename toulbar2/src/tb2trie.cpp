#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<vector>
#include "tb2trie.hpp"

using namespace std;

size_t TrieCpd::total_sequences = 0;

TrieCpd::TrieCpd()
{
	sons.resize(20, NULL);
	sequence_count = 0;
	maxc = 0;
	minc = 0;
}

TrieCpd::~TrieCpd()
{
	for (size_t i = 0; i < sons.size(); i++)
		delete sons[i];
}

int TrieCpd::aa2int(char aa)
{
	switch (aa) {
	case 'A':
		return 0;
		break;
	case 'R':
		return 1;
		break;
	case 'N':
		return 2;
		break;
	case 'D':
		return 3;
		break;
	case 'C':
		return 4;
		break;
	case 'E':
		return 5;
		break;
	case 'Q':
		return 6;
		break;
	case 'G':
		return 7;
		break;
	case 'H':
		return 8;
		break;
	case 'I':
		return 9;
		break;
	case 'L':
		return 10;
		break;
	case 'K':
		return 11;
		break;
	case 'M':
		return 12;
		break;
	case 'F':
		return 13;
		break;
	case 'P':
		return 14;
		break;
	case 'S':
		return 15;
		break;
	case 'T':
		return 16;
		break;
	case 'W':
		return 17;
		break;
	case 'Y':
		return 18;
		break;
	case 'V':
		return 19;
		break;
	default:
		cout << "Error, unrecognized residue: " << aa << endl;
		exit(1);
	}
}

char TrieCpd::int2aa(int i)
{
	switch (i) {
	case 0:
		return 'A';
		break;
	case 1:
		return 'R';
		break;
	case 2:
		return 'N';
		break;
	case 3:
		return 'D';
		break;
	case 4:
		return 'C';
		break;
	case 5:
		return 'E';
		break;
	case 6:
		return 'Q';
		break;
	case 7:
		return 'G';
		break;
	case 8:
		return 'H';
		break;
	case 9:
		return 'I';
		break;
	case 10:
		return 'L';
		break;
	case 11:
		return 'K';
		break;
	case 12:
		return 'M';
		break;
	case 13:
		return 'F';
		break;
	case 14:
		return 'P';
		break;
	case 15:
		return 'S';
		break;
	case 16:
		return 'T';
		break;
	case 17:
		return 'W';
		break;
	case 18:
		return 'Y';
		break;
	case 19:
		return 'V';
		break;
	default:
		cout << "Error, index out of range: " << i << endl;
		exit(1);
	}
}

bool TrieCpd::present(char aa)
{
	if (sons[aa2int(aa)] != NULL)
		return true;
	else
		return false;
}

void TrieCpd::insert(char aa)
{
	TrieCpd *son = new TrieCpd();
	sons[aa2int(aa)] = son;
}

void TrieCpd::insert_sequence(string seq, Cost _cost)
{
	if (seq.length()) {
		if (!present(seq[0])) {
			insert(seq[0]);
		}
		sons[aa2int(seq[0])]->insert_sequence(seq.substr(1, seq.length() - 1), _cost);
	} else {
		sequence_count++;
		if (!maxc)
			maxc = _cost;
		else if (maxc < _cost)
			maxc = _cost;
		if (!minc)
			minc = _cost;
		else if (minc > _cost)
			minc = _cost;
	}
}

void TrieCpd::print_tree()
{
	print_tree("");
}

void TrieCpd::print_tree(string acc)
{
	bool leaf = true;
	for (size_t i = 0; i < 20; i++) {
		if (sons[i] != NULL) {
			leaf = false;
			acc.push_back(int2aa(i));
			sons[i]->print_tree(acc);
			acc.pop_back();
		}
	}
	if (leaf) {
		total_sequences++;
		cout << acc << " min cost: " << minc << " max cost: " << maxc  << " this sequence occured " << sequence_count << " times" << endl;
	}
}
