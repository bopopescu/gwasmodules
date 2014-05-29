/*
 * 2013.08.20 c++ version of SelectRowsFromMatrix.py
 */
#include "pymodule/include/SelectRowsFromMatrixCC.h"


SelectRowsFromMatrixCC::SelectRowsFromMatrixCC(int _argc, char* _argv[]): AbstractMatrixFileWalkerCC(_argc, _argv) {

	//overwrite these doc variables
	usageDoc = boost::format("%1% -i INPUTFNAME -o OUTPUTFNAME --whichColumnValue WHICHCOLUMNVALUE [OPTIONS]\n")% programName;
	examplesDoc = boost::format("%1% -i /tmp/input.tsv.gz -o /tmp/output -w 0 --whichColumnValue CAE1 \n")% programName;
}

void SelectRowsFromMatrixCC::constructOptionDescriptionStructure(){
	AbstractMatrixFileWalkerCC::constructOptionDescriptionStructure();
	if (debug){
		cerr<< "adding more options within SelectRowsFromMatrixCC::_constructOptionDescriptionStructure() ...";
	}
	optionDescription.add_options()
		("whichColumnValue", po::value<string>(&whichColumnValue),
								"rows whose whichColumn value matches this would be selected.");
	if (debug){
		cerr<< endl;
	}
}

int SelectRowsFromMatrixCC::processRow(tokenizerCharType &line_toks){
	int returnValue = 0;
	tokenizerCharType::iterator tokenizer_iter = line_toks.begin();
	for (int i=0;i<whichColumn;i++){
		++tokenizer_iter;
	}
	_statInStr = *tokenizer_iter;
	//_stat = atol(_statInStr.c_str());	// use atol instead of atoi to make sure no overflow
	if (_statInStr==whichColumnValue){
		this->outputRow(line_toks);
		returnValue = 1;
	}
	else if (inputFileSortMode>0 && noOfOutput>0){	//file is sorted and has passed the section to be selected, early exit.
		returnValue = -1;
	}
	return returnValue;
}

SelectRowsFromMatrixCC::~SelectRowsFromMatrixCC(){

}

#ifndef __MAIN__

	int main(int argc, char* argv[]) {
		SelectRowsFromMatrixCC instance(argc, argv);
		instance.run();
	}
#define __MAIN__
#endif
