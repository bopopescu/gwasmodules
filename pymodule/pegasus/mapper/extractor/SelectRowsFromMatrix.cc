/*
 * 2013.08.20 c++ version of SelectRowsFromMatrix.py
 */
#include "pymodule/include/SelectRowsFromMatrix.h"


SelectRowsFromMatrix::SelectRowsFromMatrix(int _argc, char* _argv[]): AbstractMatrixFileWalker(_argc, _argv) {

}

void SelectRowsFromMatrix::constructOptionDescriptionStructure(){
	AbstractMatrixFileWalker::constructOptionDescriptionStructure();
	if (debug){
		cerr<< "adding more options within SelectRowsFromMatrix::_constructOptionDescriptionStructure() ...";
	}
	optionDescription.add_options()
		("whichColumnValue", po::value<string>(&whichColumnValue),
								"rows whose whichColumn value matches this would be selected.");
	if (debug){
		cerr<< endl;
	}
}

int SelectRowsFromMatrix::processRow(tokenizerCharType &line_toks){
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
	//_stat = atol(_statInStr.c_str());	// use atol instead of atoi to make sure no overflow
	if (_statInStr==whichColumnValue){
		this->outputRow(line_toks);
		returnValue = 1;
	}
	return returnValue;
}

SelectRowsFromMatrix::~SelectRowsFromMatrix(){

}

#ifndef __MAIN__

	int main(int argc, char* argv[]) {
		SelectRowsFromMatrix instance(argc, argv);
		instance.run();
	}
#define __MAIN__
#endif
