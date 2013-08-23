/*
 * 2013.08.21 a reducer that sums alignment depth of multiple alignments, filling missing bases with 0 depth.
 * Input are one-chromosome output of "samtools depth".
 * Output is single-column file with each value
 */

#include "pymodule/include/ReduceSameChromosomeAlignmentDepthFiles.h"


ReduceSameChromosomeAlignmentDepthFiles::ReduceSameChromosomeAlignmentDepthFiles(int _argc, char* _argv[]): AbstractMatrixFileWalkerCC(_argc, _argv) {

	//overwrite these doc variables
	usageDoc = boost::format("%1% -i INPUTFNAME -o OUTPUTFNAME --chromosomeSize CHROMOSOMESIZE [OPTIONS]\n")% programName;
	examplesDoc = boost::format("%1% -i data/ReduceSameChromosomeAlignmentDepthFiles_input1.txt.gz "
			"-i data/ReduceSameChromosomeAlignmentDepthFiles_input2.txt.gz -o data/ReduceSameChromosomeAlignmentDepthFiles_output.txt.gz -w 2 --chromosomePositionColumnIndex 1 --chromosomeSize 15 \n")% programName;
}

ReduceSameChromosomeAlignmentDepthFiles::~ReduceSameChromosomeAlignmentDepthFiles(){

}
void ReduceSameChromosomeAlignmentDepthFiles::constructOptionDescriptionStructure(){
	AbstractMatrixFileWalkerCC::constructOptionDescriptionStructure();
	if (debug){
		cerr<< "adding more options within ReduceSameChromosomeAlignmentDepthFiles::_constructOptionDescriptionStructure() ...";
	}
	optionDescription.add_options()
		("chromosomeSize", po::value<long>(&chromosomeSize),
								"size of the chromosome in input")
		("chromosomePositionColumnIndex", po::value<int>(&chromosomePositionColumnIndex)->default_value(1),
								"column index of the chromosome position column");
	if (debug){
		cerr<< endl;
	}
}



void ReduceSameChromosomeAlignmentDepthFiles::fileWalker(vector<string> &inputFnameList){
	long currentPositionInFile;
	float currentValueInFile;
	//std::vector< InputFileDataStructure> inputFileDataStructureList;
	std::vector< boost::shared_ptr< InputFileDataStructure > > inputFileDataStructureList;	//2013.08.22 have to use boost pointer , otherwise, run into error as InputFileDataStructure contains stream-type members.

	for (vector<string>::iterator inputFnameIter = inputFnameList.begin(); inputFnameIter!= inputFnameList.end(); ++inputFnameIter){
		InputFileDataStructurePtr inputFileDSPtr = InputFileDataStructurePtr( new InputFileDataStructure(*inputFnameIter, chromosomePositionColumnIndex, whichColumn));
		//inputFileDS.openFile();
		inputFileDataStructureList.push_back(inputFileDSPtr);
	}

	std::ostream outputStream(&outputFilterStreamBuffer);

	InputFileDataStructureIterator iter;
	InputFileDataStructurePtr _inputFileDSPtr;
	for (long i=1;i<=chromosomeSize; i++){
		float sumValue = 0.0;
		for (iter=inputFileDataStructureList.begin();iter!=inputFileDataStructureList.end(); iter++){
			_inputFileDSPtr = *iter;
			if (!_inputFileDSPtr->getCurrentLine().empty()){
				currentPositionInFile = _inputFileDSPtr->getCurrentPosition();
				if (currentPositionInFile==i){
					currentValueInFile = _inputFileDSPtr->getCurrentValue();
					sumValue += currentValueInFile;
					_inputFileDSPtr->getNextLine();
				}
				else if (currentPositionInFile>i){
					sumValue += 0;

				}
				else{	//should not happen
					sumValue += 0;
					std::cerr<< boost::format("ERROR: currentPositionInFile (%1%), %2%, (currentLineNumber=%3%) is less than current position across files, %4%.") %
							_inputFileDSPtr->inputFname % currentPositionInFile % _inputFileDSPtr->getCurrentLineNumber() % i
							<< endl;
					//outputStream.flush();
					//outputFile.flush();
					exit(5);
				}
			}
		}
		outputStream << sumValue << endl;
	}
	//close all input files
/*
	for (iter=inputFileDataStructureList.begin();iter!=inputFileDataStructureList.end(); iter++){
		InputFileDataStructure inputFileDS = *iter;
		*iter.closeFile();
	}
	*/
	inputFileDataStructureList.clear();

}


#ifndef __MAIN__
	int main(int argc, char* argv[]) {
		ReduceSameChromosomeAlignmentDepthFiles instance(argc, argv);
		instance.run();
	}
#define __MAIN__
#endif
