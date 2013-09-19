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
#include <map>	//for hash_map
#include <boost/program_options.hpp>	//for program options
#include <fstream>
#include <list>
#include <exception>
#include <algorithm>    // std::set_intersection, std::sort
#include <boost/tokenizer.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/format.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/assign/list_inserter.hpp>
#include <boost/any.hpp>
#include <boost/variant.hpp>
#include <boost/spirit/home/support/detail/hold_any.hpp>	//hold_any is said to be better than any, AND it can be streamed.
//#include <boost/generator_iterator.hpp>
#include "utils.h"
#include "MatrixFile.h"

// This is a typedef for a random number generator.
// Try boost::mt19937 or boost::ecuyer1988 instead of boost::minstd_rand
typedef boost::mt19937 base_generator_type;
typedef boost::tokenizer<boost::char_separator<char> > tokenizerCharType;	//to break a line into list by some delimiter

typedef boost::variant<long, std::string> boostLongStringVariant;
//do not mix int, long, or float into the type, i.e. variant<int, double, float, long, string>.
// Since one component can be constructed from another component, it'll cause this error :
//  "variant.hpp  error:call of overloaded ‘initialize(void*, const long int&)’ is ambiguous.."
typedef boost::variant<double, std::string> boostDoubleStringVariant;
typedef boost::spirit::hold_any anyType;	//anyType is streamable (outf << anyTypeObject;), while boost::any is not.


using namespace std;
using namespace boost;
namespace po = boost::program_options;


class AbstractMatrixFileWalkerCC{
protected:
	int argc;
	char** argv;		//2013.08.20 somehow, can't use "char* argv[]" even though they are same.
	//otherwise "argv=_argv;" will not work for this error: incompatible types in assignment of ‘char**’ to ‘char* [0]’

	string programName;
	boost::format usageDoc;
	boost::format examplesDoc;

	int minNoOfTotal;
	int maxNoOfTotal;
	float fractionToSample;
	int inputFileSortMode;	//2013.09.03
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
	//MatrixFile<int> outputFile;
	MatrixFilePtr outputFile;
	boost::iostreams::filtering_streambuf<boost::iostreams::output> outputFilterStreamBuffer;
	long noOfOutput;	//2013.09.03 recording how many times output has happened.

	string _currentFilename;
	string _statInStr;
	long _stat;

	po::options_description optionDescription;	//("Allowed options")
	po::positional_options_description positionOptionDescription;
	po::variables_map optionVariableMap;


public:
	AbstractMatrixFileWalkerCC(int _argc, char* _argv[]);
	virtual ~AbstractMatrixFileWalkerCC();
	virtual void constructOptionDescriptionStructure();
	virtual void parseCommandlineOptions();
	virtual void setup();
	virtual void reduce();
	virtual void closeFiles();
	virtual void openOneOutputFile(string &outputFname, \
			boost::iostreams::filtering_streambuf<boost::iostreams::output> &outputFilterStreamBuffer,\
			std::ofstream &outputFile);
	virtual void openOutputFile();
	virtual void preFileFunction();
	virtual void fileWalker(vector<string> &inputFnameList);
	virtual void openOneInputFile(string &inputFname,  boost::iostreams::filtering_streambuf<boost::iostreams::input> &inputFilterStreamBuffer);
	virtual void handleOneFile(string &inputFname);
	virtual int processRow(tokenizerCharType &line_toks);
	virtual int outputRow(tokenizerCharType &line_toks);
	virtual void postFileFunction();
	virtual void run();
};
