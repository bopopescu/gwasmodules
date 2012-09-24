/*
 * 2012.6.14 a faster version of vervet/src/mapper/CalculateMedianModeFromSAMtoolsDepthOutput.py
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

#include "armadillo"	//for mean/median calculation

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
//#include <boost/generator_iterator.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>

// This is a typedef for a random number generator.
// Try boost::mt19937 or boost::ecuyer1988 instead of boost::minstd_rand
typedef boost::mt19937 base_generator_type;

using namespace arma;
using namespace std;
using namespace boost;
namespace ba = boost::accumulators;
namespace po = boost::program_options;

class CalculateMedianMeanOfInputColumn{
	string inputFname;
	std::ifstream inputFile;
	boost::iostreams::filtering_streambuf<boost::iostreams::input> input;
	string outputFname;
	std::ofstream outputFile;
	int alignmentID;
	float fractionToSample;
	int noOfLinesInHeader;
	int whichColumn;
	int maxNumberOfSamplings;
	arma::mat statList;	//armadillo matrix storing all the sampled stats
	__gnu_cxx::hash_map <int, int> stat2Occurrence;
	int modeStat;
	// Define an accumulator set for calculating the mean and the median
	ba::accumulator_set<double, ba::stats<ba::tag::mean, ba::tag::median > > acc;

public:
	CalculateMedianMeanOfInputColumn(string _inputFname, string _outputFname, int _alignmentID,
			float _fractionToSample, int _noOfLinesInHeader, int _whichColumn, int _maxNumberOfSamplings);
	~CalculateMedianMeanOfInputColumn();
	void run();
};
