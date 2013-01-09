/*
 * 2012.6.14
 * a c++ version of vervet/src/mapper/CalculateMedianModeFromSAMtoolsDepthOutput.py
 *
 * Example:
 * 	CalculateMedianMeanOfInputColumn -i /tmp/input.tsv.gz -o /tmp/output -f 0.7
 */

#include "CalculateMedianMeanOfInputColumn.h"

CalculateMedianMeanOfInputColumn::CalculateMedianMeanOfInputColumn(string _inputFname, string _outputFname, int _alignmentID,
		float _fractionToSample, int _noOfLinesInHeader, int _whichColumn, int _maxNumberOfSamplings,
		string _inputStatName) \
		:inputFname(_inputFname), outputFname(_outputFname), alignmentID(_alignmentID), fractionToSample(_fractionToSample),\
		 noOfLinesInHeader(_noOfLinesInHeader), whichColumn(_whichColumn), maxNumberOfSamplings(_maxNumberOfSamplings),\
		 inputStatName(_inputStatName)
{
	int inputFnameLength = inputFname.length();
	if (inputFname.substr(inputFnameLength-2, 2)=="gz"){
		inputFile.open(inputFname.c_str(), std::ios::in | std::ios::binary);
		input.push(boost::iostreams::gzip_decompressor());
	}
	else{
		inputFile.open(inputFname.c_str(), std::ios::in );
	}
	input.push(inputFile);

	outputFile.open(outputFname.c_str(), std::ios::out);
	outputFile<<"alignmentID\ttotal_base_count\tsampled_base_count\tmean" + inputStatName + "\tmedian"
			+ inputStatName + "\tmode"+ inputStatName<<endl;

}



CalculateMedianMeanOfInputColumn::~CalculateMedianMeanOfInputColumn(){
	//inputFile.close();
	stat2Occurrence.clear();
	outputFile.flush();
	outputFile.close();
}


void CalculateMedianMeanOfInputColumn::run(){

	std::string line;
	size_t bytes_read = 0;
	std::istream incoming(&input);
	boost::char_separator<char> sep("\t ,");		//blank or '\t' or ',' is the separator
	typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

	std::getline(incoming, line);
	//while (incoming.good() && !incoming.eof())
	statList.set_size(10000, 1);
	int64_t noOfData = 0;
	int64_t noOfSampledData = 0;
	base_generator_type generator(42);
	//boost::mt19937 generator;
	boost::uniform_real<> uni_dist(0,1);
	boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);
	//boost::uniform_01<> uni();
	__gnu_cxx::hash_map<int, int >::iterator stat2OccurrenceIter = stat2Occurrence.begin();
	__gnu_cxx::hash_map<int, int >::iterator modeStatIter = stat2Occurrence.begin();

	double toss;	//random number
	string statInStr;
	while (!line.empty()){
		noOfData++;
		toss=uni();
		if (toss<=fractionToSample && noOfSampledData<maxNumberOfSamplings){
			tokenizer line_toks(line, sep);
			tokenizer::iterator tokenizer_iter = line_toks.begin();
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
				statList.reshape(statList.n_elem+10000,1);
			}
			
			int stat = atoi(statInStr.c_str());
			//cout<<	stat << endl;
			//statList <<  atoi(stat.c_str()) << endr;	//endr is not recognized
			statList(noOfSampledData-1,0) = stat;
			acc(stat);
			//for the mode statistic
			stat2OccurrenceIter = stat2Occurrence.find(stat);
			if (stat2OccurrenceIter==stat2Occurrence.end()){
				stat2Occurrence[stat] = 0;
			}
			stat2Occurrence[stat]++;
			stat2OccurrenceIter = stat2Occurrence.find(stat);
			modeStatIter = stat2Occurrence.find(modeStat);
			if (modeStatIter==stat2Occurrence.end()){	//doesn't exist yet
				modeStat = stat;
			}
			else{
				if (stat2OccurrenceIter->second > modeStatIter->second){
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
		std::getline(incoming, line);
	}
	statList.reshape(noOfSampledData,1);
	//cout<< statList << endl;
	double medianStat;
	double meanStat;
	if (noOfSampledData>0){
		medianStat = arma::median(statList.col(0));
		meanStat = arma::mean(statList.col(0));

		// Display the results by boost accumulators ...
		std::cout << "Mean:   " << ba::mean(acc) << std::endl;
		std::cout << "Median: " << ba::median(acc) << std::endl;
	}
	outputFile<< alignmentID << "\t" << noOfData << "\t" << noOfSampledData << "\t"<< meanStat << "\t"
			<< medianStat << "\t" << modeStat <<endl;
	//cout<< medianStat << endl;
	//cout<< meanStat << endl;
	//output it
}

int main(int argc, char* argv[]) {
	int alignmentID;
	float fractionToSample;
	int maxNumberOfSamplings;
	int whichColumn;
	int noOfLinesInHeader;
	string inputFname;
	string outputFname;
	string inputStatName;
	po::options_description desc("Allowed options");
	desc.add_options()("help,h", "produce help message")
			("alignmentID,a", po::value<int>(&alignmentID)->default_value(0),
				"ID of this alignment from which all the stats are extracted.")
			("fractionToSample,f", po::value<float>(&fractionToSample)->default_value(0.001),
				"fraction of input data to sample for median/mean calculation.")
			("noOfLinesInHeader,n", po::value<int >(&noOfLinesInHeader)->default_value(1), "how many lines the header contains")
			("whichColumn,w", po::value<int >(&whichColumn)->default_value(1), "which column of inputFname is the target stat.")
			("outputFname,o", po::value<string>(&outputFname), "output filename")
			("inputFname,i", po::value<string>(&inputFname), "input filename, space/tab/coma-delimited, gzipped or not.")
			("inputStatName,S", po::value<string>(&inputStatName)->default_value("Depth"), "name of the input statistics")
			("maxNumberOfSamplings,m", po::value<int>(&maxNumberOfSamplings)->default_value(10000000),
				"max number of samples to take into memory for median/mode/mean calculation to avoid memory blowup.");

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
			alignmentID, fractionToSample, noOfLinesInHeader, whichColumn, maxNumberOfSamplings,
			inputStatName);
	instance.run();
}
