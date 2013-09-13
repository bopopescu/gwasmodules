#ifndef __HDF5MatrixFileCC_H
#define __HDF5MatrixFileCC_H
// to avoid duplicate definition

/*
 *2013.09.13 Yu Huang copyright. class for generic matrix-like file.
 *2013.09.13 Usage:
 *	outputFile = MatrixFilePtr(new MatrixFile<int>(outputFname, "w", debug, "\t"));
 *	outputFile->writeRow(line_toks.begin(), line_toks.end());
 *	outputFile->write(sumValue);
 *	outputFile->write("\n");
 *
 *
 */
#include <iostream>
#include <string>
#include <map>
#include <fstream>
#include <cmath>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <functional>	//2013.09.11 for customize hash
#include <boost/functional/hash.hpp>	//2013.09.10 yh: for customize boost::hash
#include <boost/tokenizer.hpp>
#include <boost/assign/list_inserter.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/assign/list_inserter.hpp>
#include <boost/bimap.hpp>
#include <boost/bimap/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/device/file.hpp>	//file_sink ,file_source
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/format.hpp>
#include <boost/shared_ptr.hpp>
//#include <boost/generator_iterator.hpp>
#include "utils.h"

using namespace std;

typedef boost::shared_ptr<std::ifstream> InputFileStreamPtr;
typedef boost::shared_ptr<std::ofstream> OutputFileStreamPtr;

template <typename dataStruct>
class MatrixFile
{
	/*
	 *
	 */
	string filename;
	string openMode;
	string delimiter;
	int debug;
	int isFileGzipped;
	dataStruct *dataMatrix;	//a template class
	string header;

	vector<string> colIDLs;
	vector<string> rowIDLs;
	//boost::bimap rowID2rowIndex;	//rowID be a template class
	//boost::bimap colID2colIndex;	//colID be a template class
	int noOfRows;
	int noOfCols;
	bool isClosed;

public:
	boost::iostreams::filtering_streambuf<boost::iostreams::input> inputFilterStreamBuffer;
	std::ifstream inputFileStream;

	boost::iostreams::filtering_ostream outputStream;
	//std::ofstream outputFileStream;
	//boost::iostreams::file_sink outputFileStream;	//this makes sure gzip_compressor() flushes in the end. while ofstream does not.
	//std::ostream* outputStream;

	MatrixFile():delimiter("\t"), debug(0){
		filename="";
		openMode="";
		isClosed=true;
	}
	MatrixFile(string _filename, string _openMode, int _debug, string _delimiter): filename(_filename),
		openMode(_openMode), delimiter(_delimiter), debug(_debug){
		if (openMode.find("r")==0){
			openOneInputFile(filename, inputFilterStreamBuffer, inputFileStream);
		}
		else if (openMode.find("w")==0){
			openOneOutputFile(filename, outputStream);

		}
		isClosed=false;
	}
	void open(char* _filename, string _openMode, int _debug, string _delimiter){
		filename = _filename;
		openMode = _openMode;
		delimiter = _delimiter;
		debug = _debug;
		if (openMode.find("r")==0){
			openOneInputFile(filename, inputFilterStreamBuffer, inputFileStream);
		}
		else if (openMode.find("w")==0){
			openOneOutputFile(filename, outputStream);

		}
		isClosed = false;
	}
	virtual ~MatrixFile(){
		close();
	}
	virtual void cleanupMemory(){
		delete [] dataMatrix;
		//delete [] colIDLs;
		//delete [] rowIDLs;
	}
	virtual void close(){
		if (inputFileStream.is_open() && !isClosed){
			if (debug){
				std::cerr<<"Closing input file " << filename << " ..." ;
			}
			inputFileStream.close();
			if (debug){
				std::cerr<< endl;
			}
			isClosed=true;
		}
		else if (!isClosed){
			if (debug){
				std::cerr<<"Flushing output file " << filename << " ..." ;
			}
			outputStream.flush();
			//outputFileStream.flush();
			//outputFileStream.close();
			if (debug){
				std::cerr<< endl;
			}
			if (isFileGzipped==0){	//only close on non-gzipped, otherwise, it'll mess up the gzipped output file (EOF)
				if (debug){
					std::cerr<<"Closing output file " << filename << " ..." ;
				}
				//outputFileStream.flush();
				//outputFileStream.close();
				//outputFile.close();	//2013.08.21 if output is a .gz, closing it here would result a truncated .gz file. don't know why.
				//maybe because the buffer or the gzip compressor filter needs to be closed beforehand.
				if (debug){
					std::cerr<< endl;
				}
			}
			//delete outputStream;
			isClosed=true;
		}
	}
	virtual void openOneInputFile(string &filename, boost::iostreams::filtering_streambuf<boost::iostreams::input> &filterStreamBuffer,\
				std::ifstream &fileStream){
		if (debug){
			cerr<< boost::format("Opening file %1% for input ... ")% filename;
		}
		int filenameLength = filename.length();
		if (filename.substr(filenameLength-3, 3)==".gz"){
			filterStreamBuffer.push(boost::iostreams::gzip_decompressor());
			//filterStreamBuffer.push(boost::iostreams::file_source(filename.c_str(), ios_base::in | ios_base::binary));
			fileStream.open(filename.c_str(), std::ios::in | std::ios::binary);
			isFileGzipped=1;
		}
		else{
			fileStream.open(filename.c_str(), std::ios::in );
			isFileGzipped=0;
		}
		filterStreamBuffer.push(fileStream);
		if (debug){
			cerr<< boost::format(" opened and ready for read.") << endl;
		}
		isClosed=false;
	}
	virtual void openOneOutputFile(string &filename, \
				boost::iostreams::filtering_ostream &filterOStream){
		if (debug){
			cerr<< boost::format("Opening file %1% for output ... ")% filename;
		}
		boost::iostreams::filtering_ostream out;
		int filenameLength = filename.length();
		if (filename.substr(filenameLength-3, 3)==".gz"){
			filterOStream.push(boost::iostreams::gzip_compressor());
			//fileStream.open(filename.c_str(), std::ios::out | std::ios::binary);	// doesnot work
			filterOStream.push(boost::iostreams::file_sink(filename.c_str(), std::ios::out | std::ios::binary));
			//file_sink makes sure gzip_compressor() flushes in the end, while ofstream does not.
			isFileGzipped=1;
		}
		else{
			filterOStream.push(boost::iostreams::file_sink(filename.c_str(), std::ios::out));
			//fileStream.open(filename.c_str(), std::ios::out );
			//filterOStream.push(fileStream);
			isFileGzipped=0;
		}
		isClosed=false;
		if (debug){
			cerr<< boost::format(" opened and ready for write.") << endl;
		}

	}

	int getRowIndex(string rowID){
		return 0;
	}
	int getColIndex(string colID){
		return 0;
	}
	virtual int readInDataMatrix(){
		return 0;
	}
	virtual int readInColIDLs(){
		return 0;
	}
	virtual int readInRowIDLs(){
		return 0;
	}

	virtual void run(){
	}

	void smartReadHeader(string headerPattern="", string commentPattern="#"){
		/*
		Note:
			If an input file does not have a header, this function over-reads by one line (stored in self._row)
			so need to process the last self._row before further reading
		2013.08.30 read the header, while ignoring lines fitting the comment pattern
			and construct col_name2index when a line matching headerPattern is encountered

		if headerPattern is None:
			headerPattern = self.headerPattern
		if commentPattern is None:
			commentPattern = self.commentPattern
		row = self.next()
		while commentPattern.search(row[0]):	#passing all comments
			self.comment_row_list.append(row)
			row = self.next()
		if headerPattern.search(row[0]):
			self.header = row
			self.col_name2index = utils.getColName2IndexFromHeader(self.header)
		else:
			self.col_name2index = None
		return self.col_name2index;
		*/
	}

	string getHeader(){
		return header;
	}

	template <class contentType>
	int writeRow(vector<contentType>& contentList){
		int returnValue = 0;
		string combinedString = join(contentList, delimiter);
		outputStream << combinedString << endl;
		returnValue = 1;
		noOfRows++;
		return returnValue;
	}

	template<typename Iter>
	int writeRow(Iter begin, Iter end){
		int returnValue = 0;
		string combinedString = joinIteratorToString(begin, end, delimiter);
		outputStream << combinedString << endl;
		returnValue = 1;
		noOfRows++;
		return returnValue;
	}

	template<typename anyType>
	int write(anyType anything){
		int returnValue = 0;
		outputStream << anything;
		returnValue = 1;
		return returnValue;
	}

	template<typename anyType>
	ostream& operator<<(anyType anything){
		outputStream << anything;
	}

	void* flush(){
		if (openMode.find("r")==0){
			return 0;
		}
		else if (openMode.find("w")==0){
			return outputStream.flush();
		}
		else{
			return 0;
		}
	}
	virtual bool is_open(){
		return !isClosed;
	}
};

typedef boost::shared_ptr<MatrixFile<int> > MatrixFilePtr;	//stream is non-copyable in old c++ standard (2011 corrected it but not widespread).


#endif
