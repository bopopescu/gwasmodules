#ifndef __UTILS_H
#define __UTILS_H
// to avoid duplicate definition
/*
 * 2013.09.12 Yu Huang copyright. file to store various utility functions, classes
 */
#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <string.h>
/*
 * 2013.08.21 convenient function from http://stackoverflow.com/questions/6097927/is-there-a-way-to-implement-analog-of-pythons-separator-join-in-c
 */
template<typename Iter>
std::string joinIteratorToString(Iter begin, Iter end, std::string const& separator) {
	std::ostringstream result;
	if (begin != end)
		result << *begin++;
	while (begin != end)
		result << separator << *begin++;
	return result.str();
}

double sumOfReciprocals(int n){
	/*
	 * 2011-10-21
		for normalized nucleotide diversity
		\pi = no-of-polymorphic-loci/sumOfReciprocals
	 */
	double sum = 0.0;
	for (int i=0; i<n-1; i++){
		sum = sum + 1/(i+1.0);
	}
	return sum;
}

#endif



#ifndef __JOIN_H
#define __JOIN_H
/*
 * 2013.09.11 copied from Erik Garrison's file
 */
// join a vector of elements by a delimiter object.  ostream<< must be defined
// for both class S and T and an ostream, as it is e.g. in the case of strings
// and character arrays
template<class S, class T>
std::string join(std::vector<T>& elems, S& delim) {
    std::stringstream ss;
    typename std::vector<T>::iterator e = elems.begin();
    ss << *e++;
    for (; e != elems.end(); ++e) {
        ss << delim << *e;
    }
    return ss.str();
}

// same for lists
template<class S, class T>
std::string join(std::list<T>& elems, S& delim) {
    std::stringstream ss;
    typename std::list<T>::iterator e = elems.begin();
    ss << *e++;
    for (; e != elems.end(); ++e) {
        ss << delim << *e;
    }
    return ss.str();
}

#endif
