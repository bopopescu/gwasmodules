/*
 * 2013.08.19 a c++ version of AbstractMatrixFileWalker.py
 *
 */
#include <iostream>
#include <string>
#include <cmath>	//sqrt
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <ext/hash_map>	//for hash_map
#include <boost/program_options.hpp>	//for program options
#include <fstream>
#include <exception>
#include <boost/tokenizer.hpp>


#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/format.hpp>
//#include <boost/generator_iterator.hpp>


// This is a typedef for a random number generator.
// Try boost::mt19937 or boost::ecuyer1988 instead of boost::minstd_rand
typedef boost::mt19937 base_generator_type;
typedef boost::tokenizer<boost::char_separator<char> > tokenizerCharType;	//to break a line into list by some delimiter

using namespace std;
using namespace boost;
namespace po = boost::program_options;

class AbstractMatrixFileWalker{
protected:
	int minNoOfTotal;
	int maxNoOfTotal;
	float fractionToSample;
	int whichColumn;
	string whichColumnHeader;
	vector<string> inputFnameList;
	int inputFileFormat;
	string outputFname;
	int outputFileFormat;
	string missingDataNotation;
	int debug;
	int report;

	std::ifstream inputFile;

	std::ofstream outputFile;
	string _currentFilename;
	string _statInStr;
	long _stat;

	po::options_description optionDescription;	//("Allowed options")
	po::positional_options_description positionOptionDescription;
	po::variables_map optionVariableMap;
	void _constructOptionDescriptionStructure();


public:
	AbstractMatrixFileWalker(int argc, char* argv[]);
	~AbstractMatrixFileWalker();
	void setup();
	void reduce();
	void closeFiles();
	void preFileFunction();
	void postFileFunction();
	void fileWalker(string &inputFname);
	int processRow(tokenizerCharType &line_toks);
	void run();
};
