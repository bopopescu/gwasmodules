/*
 * 2012.6.14
 * a c++ version of vervet/src/mapper/CalculateMedianModeFromSAMtoolsDepthOutput.py
 *
 * Example:
 * 	CalculateMedianMeanOfInputColumn -i /tmp/input.tsv.gz -o /tmp/output -f 0.7
 */

#include "pymodule/include/CalculateMedianMeanOfInputColumn.h"

CalculateMedianMeanOfInputColumn::CalculateMedianMeanOfInputColumn(string _inputFname, string _outputFname, int _alignmentID,
		float _fractionToSample, int _whichColumn, long _maxNumberOfSamplings,
		int _inputFileFormat, string _inputStatName) \
		:inputFname(_inputFname), outputFname(_outputFname), alignmentID(_alignmentID), fractionToSample(_fractionToSample),\
		 whichColumn(_whichColumn), maxNumberOfSamplings(_maxNumberOfSamplings),\
		 inputFileFormat( _inputFileFormat), inputStatName(_inputStatName)
{
	noOfData=0;
	noOfSampledData=0;
	sumOfCoverageAtAllLoci = 0;
	statListExpansionStepSize=100000;
	int inputFnameLength = inputFname.length();
	if (inputFname.substr(inputFnameLength-2, 2)=="gz"){
		inputFilterStreamBuffer.push(boost::iostreams::gzip_decompressor());
		inputFile.open(inputFname.c_str(), std::ios::in | std::ios::binary);
	}
	else{
		inputFile.open(inputFname.c_str(), std::ios::in );
	}
	inputFilterStreamBuffer.push(inputFile);

	outputFile.open(outputFname.c_str(), std::ios::out);
	outputFile<<"alignmentID\ttotal_base_count\tsampled_base_count\tmean" + inputStatName + "\tmedian"
			+ inputStatName + "\tmode"+ inputStatName + "\tsumOf" + inputStatName << endl;

}



CalculateMedianMeanOfInputColumn::~CalculateMedianMeanOfInputColumn(){
	//inputFile.close();
	stat2Frequency.clear();
	outputFile.flush();
	outputFile.close();
}

void CalculateMedianMeanOfInputColumn::skipHeaderInInput(int noOfLinesInHeader){
	std::cerr<<"Skipping " << noOfLinesInHeader << " lines from " << inputFname ;
	std::istream inputStream(&inputFilterStreamBuffer);
	std::string line;
	int noOfLinesRead=0;
	while (noOfLinesRead<noOfLinesInHeader){
		std::getline(inputStream, line);
		noOfLinesRead++;
	}
	std::cerr<< std::endl;
}


int CalculateMedianMeanOfInputColumn::readInRawDataFromSAMtoolsDepthFile(){
	std::cerr<<"Reading data from " << inputFname << std::endl;
	std::string line;
	size_t bytes_read = 0;
	std::istream inputStream(&inputFilterStreamBuffer);
	boost::char_separator<char> sep("\t ,");		//blank or '\t' or ',' is the separator

	std::getline(inputStream, line);
	//while (inputStream.good() && !inputStream.eof())
	statList.set_size(statListExpansionStepSize, 1);
	base_generator_type generator(42);
	//boost::mt19937 generator;
	boost::uniform_real<> uni_dist(0,1);
	boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);
	//boost::uniform_01<> uni();
	map<long, long >::iterator stat2FrequencyIter = stat2Frequency.begin();
	map<long, long >::iterator modeStatIter = stat2Frequency.begin();

	double toss;	//random number
	string statInStr;
	long stat;
	while (!line.empty()){
		noOfData++;
		toss=uni();
		if (toss<=fractionToSample && noOfSampledData<maxNumberOfSamplings){
			tokenizerCharType line_toks(line, sep);
			tokenizerCharType::iterator tokenizer_iter = line_toks.begin();
			for (int i=0;i<whichColumn;i++){
				/*
				if (tokenizer_iter==line_toks.end()){
					//2012.9.15 debug purpose
					std::cerr<< "end of line has been reached before the whichColumn, abort." << std::endl;
					std::cerr<< "Exception line no. " << noOfData << " has content: " << line <<std::endl;
					continue;	//end has been reached before whichColumn, abort.
				}
				*/
				++tokenizer_iter;
			}
			statInStr = *tokenizer_iter;

			noOfSampledData ++;
			if (noOfSampledData>statList.n_elem){
				statList.reshape(statList.n_elem+statListExpansionStepSize,1);
			}

			stat = atol(statInStr.c_str());	// use atol instead of atoi to make sure no overflow
			//cout<<	stat << endl;
			//statList <<  atol(stat.c_str()) << endr;	//endr is not recognized
			statList(noOfSampledData-1,0) = stat;
			accumulatorSet(stat);
			//for the mode statistic
			stat2FrequencyIter = stat2Frequency.find(stat);
			if (stat2FrequencyIter==stat2Frequency.end()){
				stat2Frequency[stat] = 0;
			}
			stat2Frequency[stat]++;
			stat2FrequencyIter = stat2Frequency.find(stat);
			modeStatIter = stat2Frequency.find(modeStat);
			if (modeStatIter==stat2Frequency.end()){	//doesn't exist yet
				modeStat = stat;
			}
			else{
				if (stat2FrequencyIter->second > modeStatIter->second){
					modeStat = stat;
				}
			}
		}
		/*
		if (noOfSampledData>=maxNumberOfSamplings){
			std::cerr<< "Enough data is sampled. noOfSampledData=" << noOfSampledData << " from " << noOfData << "." << std::endl;
			break;
		}
		*/
		//bytes_read += line.size();
		// progress dlg with bytes_read / uncompressed size
		std::getline(inputStream, line);
	}
	statList.reshape(noOfSampledData,1);
	std::cerr << noOfSampledData <<" data points." << std::endl;

	//2013.06.09 compute the mean/median/mode
	//cout<< statList << endl;
	if (noOfSampledData>0){
		sumOfCoverageAtAllLoci = ba::sum(accumulatorSet);

		medianStat = arma::median(statList.col(0));
		meanStat = arma::mean(statList.col(0));

		// Display the results by boost accumulators ...
		std::cout << "Mean:   " << ba::mean(accumulatorSet) << std::endl;
		std::cout << "Median: " << ba::median(accumulatorSet) << std::endl;
		std::cout << "sumOfCoverageAtAllLoci: " << sumOfCoverageAtAllLoci << std::endl;
	}

	return 1;
}

int CalculateMedianMeanOfInputColumn::readInRawDataFromGATKDOCOutput(){
	/*
	 * output of GATK DOC walker (.summary_statistics) looks like this (2 lines):
	 *
	 *	Source_of_reads from_0_to_1)	from_1_to_2)	from_2_to_3)	from_3_to_4)	from_4_to_5)	from_5_to_6) ...
	 * 	sample_4215_768_1997022_GA_vs_3280	20745136	29988560	55808131	90706436	132088981 ...
	 *
	 *
	 */
	skipHeaderInInput(1);

	std::cerr<<"Reading data from " << inputFname << std::endl;
	std::string line;
	std::istream inputStream(&inputFilterStreamBuffer);
	boost::char_separator<char> sep("\t ");		//blank or '\t' is the field separator

	std::getline(inputStream, line);
	string statInStr;
	statList.set_size(statListExpansionStepSize, 1);

	tokenizerCharType line_toks(line, sep);
	//tokenizerCharType::iterator tokenizer_iter = line_toks.begin();
	long inputColumnNumber=0;
	long coverageAtTheseLoci;
	long noOfLociAtModeDepth = -1;
	long noOfLociAtThisCoverage;
	for(tokenizerCharType::iterator it = line_toks.begin(), ite = line_toks.end(); it!=ite; ++it){
		inputColumnNumber++;
		if (inputColumnNumber==1){	//skip the first column (which is sampleID/ReadGroup)
			continue;
		}
		statInStr = *it;
		noOfLociAtThisCoverage = atol(statInStr.c_str());	//atoi is not enough for some cases
		noOfData += noOfLociAtThisCoverage;
		noOfSampledData += noOfLociAtThisCoverage;	//same as noOfData
		coverageAtTheseLoci = inputColumnNumber-1;
		if (coverageAtTheseLoci>statList.n_elem){
			statList.reshape(statList.n_elem+statListExpansionStepSize,1);
		}

		statList(coverageAtTheseLoci-1,0) = noOfLociAtThisCoverage;
		sumOfCoverageAtAllLoci += (noOfLociAtThisCoverage * coverageAtTheseLoci);

		//for the mode statistic, update modeStat if this depth has more loci than that of modeStat
		if (noOfLociAtThisCoverage>noOfLociAtModeDepth){
			modeStat = coverageAtTheseLoci;
			noOfLociAtModeDepth = noOfLociAtThisCoverage;
		}
	}
	statList.reshape(inputColumnNumber-1,1);
	std::cerr << inputColumnNumber << " columns, number of data points= " << noOfData << ", sumOfCoverageAtAllLoci= " << sumOfCoverageAtAllLoci <<", mode = " << modeStat << ", count at modeStat is "
			<< noOfLociAtModeDepth << "." << std::endl;

	// compute the mean
	meanStat = sumOfCoverageAtAllLoci/double(noOfData);
	std::cerr << " mean = "<<meanStat << std::endl;

	std::cerr << " Computing median ... " ;
	// compute median
	long localSumOfCoverage=0;
	for (long i=0; i<statList.n_elem; i++){
		localSumOfCoverage += (statList(i,0)*(i+1));	// i+1 is coverage at these loci
		if (localSumOfCoverage>0.5*sumOfCoverageAtAllLoci){
			//just over half, set medianStat to be the number right before this.
			// and break the loop,
			medianStat = max(long(1),i);	//set minimum to 1
			break;
		}
	}
	std::cerr << " median = "<<medianStat << std::endl;
	// mode is computed on the fly
	return 1;
}

void CalculateMedianMeanOfInputColumn::run(){
	if (inputFileFormat==1){
		readInRawDataFromSAMtoolsDepthFile();
	}
	else if (inputFileFormat==2){
		readInRawDataFromGATKDOCOutput();
	}
	else{
		std::cerr<< "ERROR: inputFileFormat " << inputFileFormat << " is not supported" <<std::endl;
		exit(3);
	}

	//output
	outputFile<< alignmentID << "\t" << noOfData << "\t" << noOfSampledData << "\t"<< meanStat << "\t"
			<< medianStat << "\t" << modeStat << "\t" << sumOfCoverageAtAllLoci << endl;
	//cout<< medianStat << endl;
	//cout<< meanStat << endl;
	//output it
}

int main(int argc, char* argv[]) {
	int alignmentID;
	float fractionToSample;
	int whichColumn;
	string outputFname;
	string inputFname;
	int inputFileFormat;
	string inputStatName;
	long maxNumberOfSamplings;
	po::options_description desc("Allowed options");
	desc.add_options()("help,h", "produce help message")
			("alignmentID,a", po::value<int>(&alignmentID)->default_value(0),
				"ID of this alignment from which all the stats are extracted.")
			("fractionToSample,f", po::value<float>(&fractionToSample)->default_value(0.001),
				"fraction of input data to sample for median/mean calculation. --inputFileFormat=1 only")
			("whichColumn,w", po::value<int >(&whichColumn)->default_value(1), "which column of inputFname is the target stat, --inputFileFormat=1 only.")
			("outputFname,o", po::value<string>(&outputFname), "output filename")
			("inputFname,i", po::value<string>(&inputFname), "input filename, space/tab/coma-delimited, gzipped or not.")
			("inputFileFormat,y", po::value<int >(&inputFileFormat)->default_value(1),
					"input file format. 1: one column of each line contains one data, i.e. samtools depth output; 2: stats about distribution of data points. GATK DOC/DepthOfCoverage walker .summary_statistics output.")
			("inputStatName,S", po::value<string>(&inputStatName)->default_value("Depth"), "name of the input statistics, used only for reporting")
			("maxNumberOfSamplings,m", po::value<long>(&maxNumberOfSamplings)->default_value(10000000),
				"max number of samples to take into memory for median/mode/mean calculation to avoid memory blowup. --inputFileFormat=1 only");

	//po::positional_options_description p;
	//p.add("input-file", -1);

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	if (vm.count("help") || inputFname.empty() || outputFname.empty()){
		cout << desc << endl;
		return 1;
	}



	CalculateMedianMeanOfInputColumn instance(inputFname, outputFname,
			alignmentID, fractionToSample, whichColumn, maxNumberOfSamplings,
			inputFileFormat, inputStatName);
	instance.run();
}
