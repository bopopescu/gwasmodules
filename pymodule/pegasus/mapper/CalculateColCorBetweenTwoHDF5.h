/*
 *
//	2012.2.29
/
 * */
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>	//sqrt
#include <getopt.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <ext/hash_map>	//for hash_map

#ifndef H5_NO_NAMESPACE
#ifndef H5_NO_STD
using std::cout;
using std::endl;
#endif  // H5_NO_STD
#endif
using namespace std;

#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

const H5std_string DataMatrixDataSetName("data_matrix");
const H5std_string ColIDDataSetName("col_id_ls");
const H5std_string outputDatasetName("correlation");
const H5std_string outputTypeMember1("input1LocusID");
const H5std_string outputTypeMember2("input2LocusID");
const H5std_string outputTypeMember3("correlation");

/* structure for the output data*/
typedef struct s1_t {
	int a;
	int b;
	float c;
} s1_t;

class CalculateColCorBetweenTwoHDF5
{
	/*
	 *
	 * */
	H5std_string input1Fname;
	H5std_string input2Fname;
	int i1_start;
	int i1_stop;
	int i2_start;
	int i2_stop;
	H5File input1;
	H5File input2;
	DataSet dataset1;
	DataSet dataset2;
	DataSpace dataspace1;
	DataSpace dataspace2;

	int *dataMatrix1InMemory;
	int *input1_col_id_ls;
	int *dataMatrix2InMemory;	//** won't work with dataset.read() as it requires ''void*''
	int *input2_col_id_ls;
	int no_of_rows;
	int input1_no_of_cols;
	int input2_no_of_cols;
	float min_cor;

	H5std_string outputFname;
	H5File out;

	int outputDatasetMaxLength;
	int outputDatasetLength;
	int outputDatasetRank;
	s1_t* outputDataInMemory;

	public:
		CalculateColCorBetweenTwoHDF5(char* _input1Fname, char* _input2Fname, char* outf_name, int _i1_start, int _i1_stop,
			int _i2_start, int _i2_stop, float cor_cut_off_given);
		~CalculateColCorBetweenTwoHDF5();
		int readDataIntoMemory(DataSet &dataset, DataSpace &dataspace, int* &dataMatrixInMemory, int &col_start, int &no_of_cols);
		int readColIDIntoMemory(H5File &h5file, int* &col_id_ls, int &col_start, int &no_of_cols);
		int input();
		double cor(int* &data_matrix1, int* &data_matrix2, int matrix1_col_index, int matrix2_col_index, int &no_of_rows);
		void cleanupMemory();
		void run();
		int output(s1_t* &dataInMemory, H5std_string &outputFname);

};


class FindMaxLDBetweenPeakAndEachLocus
{
	/*
	 *
	 * */
	H5std_string input1Fname;
	int i1_start;
	int i1_stop;
	int noOfSelectedEntries;
	H5File input1;
	DataSet dataset1;
	DataSpace dataspace1;

	H5std_string withinPeakLocusIDFname;
	s1_t *inputCorrelationInMemory;

	__gnu_cxx::hash_map <int, s1_t> fstLocusId2CorStruc;
	__gnu_cxx::hash_map <int, int> sndLocusIdMap;

	H5std_string outputFname;
	H5File out;
	int outputDatasetMaxLength;
	int outputDatasetLength;
	int outputDatasetRank;
	s1_t* outputDataInMemory;

	public:
		FindMaxLDBetweenPeakAndEachLocus(char* _input1Fname, char* _withinPeakLocusIDFname, char* _outputFname, int _i1_start, int _i1_stop);
		~FindMaxLDBetweenPeakAndEachLocus();
		int readDataIntoMemory(DataSet &dataset, DataSpace &dataspace);
		int readLocusIDLs();
		void cleanupMemory();
		void run();
		int output(H5std_string &outputFname);
};
