#ifndef __POLYMORPHISM_H
#define __POLYMORPHISM_H
// to avoid duplicate definition

/*
 * 2013.09.12 Yu Huang. copyright. header file for polymorphism related classes and functions.
 */
#include <vector>
#include <list>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <utility>
#include <stdlib.h>
#include <assert.h>
#include <stack>
#include <queue>
#include <set>
#include <functional>	//2013.09.11 for customize hash
#include <boost/functional/hash.hpp>	//2013.09.10 yh: for customize boost::hash
#include <boost/bimap.hpp>
#include <boost/assign/list_inserter.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/assign/list_inserter.hpp>
#include <boost/bimap/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/bimap/list_of.hpp>
#include <boost/algorithm/string/join.hpp>

using namespace std;

typedef boost::bimap<string, int > nucleotide2integerBiMapType;
typedef nucleotide2integerBiMapType::value_type nucleotideIntegerTuple;

class GenotypeCoder{
public:
	nucleotide2integerBiMapType nucleotide2integerBiMap;

	GenotypeCoder(){
		initializeNucleotide2integerBiMap();
	}
	~GenotypeCoder(){
		;
	}


	void initializeNucleotide2integerBiMap() {
		/*
		 * 4 smallest prime numbers correspond to 4 nucleotides

		nucleotide2integerBiMap.insert(nucleotideIntegerTuple("A", 2));
		nucleotide2integerBiMap.insert(nucleotideIntegerTuple("C", 3));
		nucleotide2integerBiMap.insert(nucleotideIntegerTuple("G", 5));
		nucleotide2integerBiMap.insert(nucleotideIntegerTuple("T", 7));
		*/
		//nucleotide2integerBiMap = boost::assign::list_of<nucleotide2integerBiMapType ::relation > ("A", 2);
		boost::assign::insert(nucleotide2integerBiMap.left)("A", 2)("C", 3)("G", 5)("T", 7);
	}

	int encode(string genotype){
		/*
		 * 2013.09.11 Encode genotype as multiplication of prime numbers.
		 * genotype is like 'AA', 'A', 'A/G', 'A|G', 'A/A'.
		 * '/' or '|' within genotype is ignored.
		 */
		int genotypeInteger =1;
		if (genotype=="NA"){
			genotypeInteger=0;
		}
		else{
			for(string::iterator gIterator=genotype.begin(); gIterator!=genotype.end(); gIterator++){
				char singleNucleotide = *gIterator;
				if (singleNucleotide!='/' && singleNucleotide!='|'){
					genotypeInteger = genotypeInteger*nucleotide2integerBiMap.left.at(string(gIterator, gIterator+1));
				}
			}
		}
		return genotypeInteger;
	}
	string decode(int genotypeInteger){
		/*
		 * 2013.09.11 divide the genotypeInteger by each number in nucleotide2integerBiMap until it is ==1
		 */
		nucleotide2integerBiMapType::right_const_iterator right_iter = nucleotide2integerBiMap.right.begin();
		vector<string> nucleotideList;
		while (genotypeInteger>1){
			for (nucleotide2integerBiMapType::right_const_iterator rend = nucleotide2integerBiMap.right.end(); \
			right_iter != rend; ++right_iter) {
				// right_iter->first  : key  : int
				// right_iter->second : data : std::string
				if (genotypeInteger%right_iter->first==0){	//found one nucleotide
					nucleotideList.push_back(right_iter->second);
					genotypeInteger = genotypeInteger/right_iter->first;	//peel one prime number off
					break;
				}
			}
		}
		string genotype = boost::algorithm::join(nucleotideList, "");
		return genotype;
	}

};


#endif
