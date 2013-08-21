/*
 * 2013.08.19 a c++ version of AbstractMatrixFileWalker.py
 *
 */
#include "pymodule/include/AbstractMatrixFileWalkerCC.h"

AbstractMatrixFileWalkerCC::AbstractMatrixFileWalkerCC(int _argc, char* _argv[]): argc(_argc), argv(_argv)
{
	debug=0;
	report=0;
	programName = _argv[0];
	//strcpy(_argv[0], programName.c_str());	//2013.08.20 does no work

	cerr << "program name is " << programName << "." <<endl;

	usageDoc = boost::format("%1% -i INPUTFNAME -o OUTPUTFNAME [OPTIONS]\n")% programName;
	examplesDoc = boost::format("%1% -i /tmp/input.tsv.gz -o /tmp/output -w 0 \n")% programName;
}

AbstractMatrixFileWalkerCC::~AbstractMatrixFileWalkerCC()
{
	//clearing
	if (debug){
		std::cerr<<"Deconstructing from AbstractMatrixFileWalkerCC." << endl;
	}
}


void AbstractMatrixFileWalkerCC::constructOptionDescriptionStructure(){

	optionDescription.add_options()("help,h", "produce help message")
				("minNoOfTotal,m", po::value<int>(&minNoOfTotal)->default_value(0),
					"minimum no of data from one file for afterFileFunction() to run")
				("maxNoOfTotal,x", po::value<int>(&maxNoOfTotal)->default_value(-1),
					"maximum no of data to sample from one file. if not set, no limit.")
				("fractionToSample,f", po::value<float>(&fractionToSample)->default_value(1.0),
					"how often you include the data, a probability between 0 and 1.")
				("whichColumn,w", po::value<int >(&whichColumn)->default_value(0),
						"which column of inputFname is the target stat, --inputFileFormat=1 only.")
				("whichColumnHeader,W", po::value<string >(&whichColumnHeader),
						"column header of the data to be used, substitute whichColumn")
				("missingDataNotation", po::value<string>(&missingDataNotation),
						"missing data notation. missing data will be skipped.")
				("debug,b", "toggle debug mode")
				("report,r", "toggle report mode")
				("inputFname,i", po::value<vector<string> >(&inputFnameList),
						"input filename, space/tab/coma-delimited, gzipped or not. could be specified as options or positional arguments")
				("inputFileFormat", po::value<int >(&inputFileFormat)->default_value(1),
						"input file format")
				("outputFname,o", po::value<string>(&outputFname), "output filename")
				("outputFileFormat", po::value<int >(&outputFileFormat)->default_value(1),
						"output file format");

}


void AbstractMatrixFileWalkerCC::parseCommandlineOptions(){
	//all positional arguments are input files.
	positionOptionDescription.add("inputFname", -1);

	po::store(po::command_line_parser(argc, argv).
				  options(optionDescription).positional(positionOptionDescription).run(), optionVariableMap);

	//po::store(po::parse_command_line(argc, argv, optionDescription), optionVariableMap);
	po::notify(optionVariableMap);
	if (optionVariableMap.count("help") || inputFnameList.size()<=0){
		cout << "Usage:" << endl << usageDoc << endl;

		cout << optionDescription << endl << endl;
		cout << "Examples:" << endl << examplesDoc << endl;
		exit(1);
	}
	if (optionVariableMap.count("debug")){
		debug = 1;
	}
	else{
		debug = 0;
	}
	if (optionVariableMap.count("report")){
		report = 1;
	}
	else{
		report = 0;
	}
	if (!outputFname.empty()){
		if (debug){
			std::cerr<<"Open file " << outputFname << " for writing " ;
		}
		outputFile.open(outputFname.c_str(), std::ios::out);
		if (debug){
			std::cerr<< endl ;
		}
	}
}

void AbstractMatrixFileWalkerCC::setup(){

}
void AbstractMatrixFileWalkerCC::reduce(){

}

void AbstractMatrixFileWalkerCC::preFileFunction(){
}

void AbstractMatrixFileWalkerCC::postFileFunction(){
}

int AbstractMatrixFileWalkerCC::processRow(tokenizerCharType &line_toks){
	int returnValue = 0;
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
	_statInStr = *tokenizer_iter;
	_stat = atol(_statInStr.c_str());	// use atol instead of atoi to make sure no overflow
	returnValue = 1;
	return returnValue;
}

int AbstractMatrixFileWalkerCC::outputRow(tokenizerCharType &line_toks){
	int returnValue = 0;
	tokenizerCharType::iterator tokenizer_iter = line_toks.begin();
	for (;tokenizer_iter!=line_toks.end();tokenizer_iter++){
			/*
			if (tokenizer_iter==line_toks.end()){
				//2012.9.15 debug purpose
				std::cerr<< "end of line has been reached before the whichColumn, abort." << std::endl;
				std::cerr<< "Exception line no. " << noOfData << " has content: " << line <<std::endl;
				continue;	//end has been reached before whichColumn, abort.
			}
			*/
		if (tokenizer_iter==line_toks.begin()){
			outputFile << *tokenizer_iter ;
		}
		else{
			outputFile << "\t" << *tokenizer_iter ;
		}

	}
	outputFile << endl;
	returnValue = 1;
	return returnValue;
}

void AbstractMatrixFileWalkerCC::fileWalker(string &inputFname){
	if (debug){
		std::cerr<<"walking through " << inputFname << "..." << std::endl;
	}
	_currentFilename = inputFname;

	int counter = 0;
	int real_counter = 0;
	int noOfSampled = 0;
	preFileFunction();
	boost::char_separator<char> sep("\t ,");		//blank or '\t' or ',' is the separator

	base_generator_type generator(42);
	//boost::mt19937 generator;
	boost::uniform_real<> uni_dist(0,1);
	boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);
	double toss;	//random number
	std::string line;
	int rowReturnValue;
	boost::iostreams::filtering_streambuf<boost::iostreams::input> inputFilterStreamBuffer;
	try{
		int inputFnameLength = inputFname.length();
		if (inputFname.substr(inputFnameLength-2, 2)=="gz"){
			inputFilterStreamBuffer.push(boost::iostreams::gzip_decompressor());
			inputFile.open(inputFname.c_str(), std::ios::in | std::ios::binary);
		}
		else{
			inputFile.open(inputFname.c_str(), std::ios::in );
		}
		inputFilterStreamBuffer.push(inputFile);
		std::istream inputStream(&inputFilterStreamBuffer);

		//col_name2index = reader.constructColName2IndexFromHeader()
		//header = reader.getHeader()

		//#output the header to the output file if necessary
		//self.processHeader(header=header, pdata=pdata) #2012.8.13

		std::getline(inputStream, line);

		while (!line.empty()){
			counter++;
			if (fractionToSample>0 && fractionToSample<1){	//random sampling
				toss=uni();
				if (toss>fractionToSample){
					continue;
				}
			}
			if (maxNoOfTotal>0 && real_counter>maxNoOfTotal){
				break;
			}

			tokenizerCharType line_toks(line, sep);
			rowReturnValue = processRow(line_toks);
			noOfSampled++;
			real_counter += rowReturnValue;
			std::getline(inputStream, line);
		}

		if (noOfSampled>=minNoOfTotal){
			postFileFunction();
		}
		inputFile.close();
	}
	catch (std::exception& e)
	{
		cerr << "An exception occurred: " << e.what() << endl;
		exit(4);
	}
	float fraction;
	if (counter>0){
		fraction = float(noOfSampled)/float(counter);
	}
	else{
		fraction = 0;
	}
	std::cerr<< boost::format("%1% / %2% (%3%) data sampled. real_counter=%4%.") % noOfSampled % counter % fraction % real_counter << endl;

}



void AbstractMatrixFileWalkerCC::closeFiles(){
	//
	if (debug){
		std::cerr<<"closing files " << "..." ;
	}
	if (outputFile){
		outputFile.close();
	}
}
void AbstractMatrixFileWalkerCC::run(){
	constructOptionDescriptionStructure();
	parseCommandlineOptions();

	if (debug){
		std::cerr<<"Entering AbstractMatrixFileWalkerCC.run()..." << std::endl;
	}

	setup();

	for (vector<string>::iterator inputFnameIter = inputFnameList.begin(); inputFnameIter!= inputFnameList.end(); ++inputFnameIter){
		fileWalker(*inputFnameIter);
	}

	reduce();

	closeFiles();
	if (debug){
		std::cerr<<"Exit AbstractMatrixFileWalkerCC.run()." << std::endl;
	}
}


#ifndef __MAIN__
	int main(int argc, char* argv[]) {
		AbstractMatrixFileWalkerCC instance(argc, argv);
		instance.run();
	}
#define __MAIN__
#endif
