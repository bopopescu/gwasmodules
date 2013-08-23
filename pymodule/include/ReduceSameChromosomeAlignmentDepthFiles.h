/*
 * 2013.08.21 a reducer that sums alignment depth of multiple alignments, filling missing bases with 0 depth.
 * Input are one-chromosome output of "samtools depth".
 * Output is single-column file with each value
 */

#include "AbstractMatrixFileWalkerCC.h"

class ReduceSameChromosomeAlignmentDepthFiles : public AbstractMatrixFileWalkerCC {
protected:
	long chromosomeSize;
	int chromosomePositionColumnIndex;


public:
	//must define non-default constructor, the parental non-default constructor won't be inherited by default.
	ReduceSameChromosomeAlignmentDepthFiles(int _argc, char* _argv[]);
	~ReduceSameChromosomeAlignmentDepthFiles();

	void constructOptionDescriptionStructure();
	void fileWalker(vector<string> &inputFnameList);
};



class InputFileDataStructure{
public:
	std::string inputFname;
	boost::iostreams::filtering_streambuf<boost::iostreams::input> inputFilterStreamBuffer;
	std::ifstream inputFile;
	std::string currentLine;
	long currentPosition;
	float currentValue;
	long currentLineNumber;
	int chromosomePositionColumnIndex;
	int valueColumnIndex;
	tokenizerCharType::iterator tokenizer_iter ;
	string _statInStr;



public:
	void openFile(){

		int inputFnameLength = inputFname.length();
		if (inputFname.substr(inputFnameLength-3, 3)==".gz"){
			inputFilterStreamBuffer.push(boost::iostreams::gzip_decompressor());
			inputFile.open(inputFname.c_str(), std::ios::in | std::ios::binary);
		}
		else{
			inputFile.open(inputFname.c_str(), std::ios::in );
		}
		inputFilterStreamBuffer.push(inputFile);
	}
	//InputFileDataStructure(){};
	InputFileDataStructure(std::string _inputFname,  int _chromosomePositionColumnIndex, int _valueColumnIndex):
			inputFname(_inputFname), chromosomePositionColumnIndex(_chromosomePositionColumnIndex), valueColumnIndex(_valueColumnIndex){
		openFile();
		currentLineNumber=0;
		if (!inputFile.eof()){
			getNextLine();
		}
	}

	~InputFileDataStructure(){
		closeFile();
	}

	void closeFile(){
		if (inputFile.is_open()){
			inputFile.close();
		}
	}
	long getCurrentPosition(){
		return currentPosition;
	}

	std::string getCurrentLine(){
		return currentLine;
	}

	float getCurrentValue(){
		return currentValue;
	}
	long getCurrentLineNumber(){
		return currentLineNumber;
	}

	std::string getNextLine(){
		std::istream inputStream(&inputFilterStreamBuffer);
		std::getline(inputStream, currentLine);
		if (!currentLine.empty()){
			parseLine(currentLine);
			currentLineNumber++;
		}
		return currentLine;
	}

	void parseLine(std::string line){
		boost::char_separator<char> sep("\t ,");		//blank or '\t' or ',' is the separator
		tokenizerCharType line_toks(line, sep);
		int i=0;
		for (tokenizer_iter= line_toks.begin();tokenizer_iter!=line_toks.end();tokenizer_iter++){
			_statInStr = *tokenizer_iter;
			if (i==chromosomePositionColumnIndex){
				currentPosition = atol(_statInStr.c_str());
			}
			else if (i==valueColumnIndex){
				currentValue = atof(_statInStr.c_str());
			}
			i++;

		}
	}

};

typedef std::vector< boost::shared_ptr< InputFileDataStructure > >::iterator InputFileDataStructureIterator;
typedef boost::shared_ptr<InputFileDataStructure> InputFileDataStructurePtr;

