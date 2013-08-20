/*
 * 2012.3.6
 *
 */

#include "pymodule/include/CalculateColCorBetweenTwoHDF5.h"


FindMaxLDBetweenPeakAndEachLocus::FindMaxLDBetweenPeakAndEachLocus(char* _input1Fname, char* _withinPeakLocusIDFname, char* _outputFname, int _i1_start, int _i1_stop)\
	:outputFname(_outputFname), i1_start(_i1_start), i1_stop(_i1_stop), input1Fname(_input1Fname), withinPeakLocusIDFname(_withinPeakLocusIDFname)
{
	std::cerr<<"FindMaxLDBetweenPeakAndEachLocus starting ...";
	/*
	 * 2012.3.28 check whether input file exists
	 */
	exitIfFileNotGood(_input1Fname, 0);	//exitCode is 0 because in the context of FindGenomeWideLDPatternBetweenSNPsAndPeakWorkflow.py
	// i don't want clustered jobs in the workflow to be disrupted when the input is empty (no SNPs in LD with given peak).
	exitIfFileNotGood(_withinPeakLocusIDFname, 1); //exitCode is 1 because this should not happen.

	input1 = H5File(_input1Fname, H5F_ACC_RDONLY);
	dataset1 = input1.openDataSet(outputDatasetName);
	std::cerr << endl;
	/*
	 * Get the class of the datatype that is used by the dataset.
	 */
	H5T_class_t type_class = dataset1.getTypeClass();
#if defined(DEBUG)
	if (type_class == H5T_INTEGER) {
		cout << "Data set has INTEGER type" << endl;

		/*
		 * Get the integer datatype
		 */
		IntType intype = dataset1.getIntType();

		/*
		 * Get order of datatype and print message if it's a little endian.
		 */H5std_string order_string;
		H5T_order_t order = intype.getOrder(order_string);
		cout << order_string << endl;

		/*
		 * Get size of the data element stored in file and print it.
		 */
		size_t size = intype.getSize();
		cout << "Data size is " << size << endl;
	}
#endif
	/*
	 * Get dataspace of the dataset.
	 */
	dataspace1 = dataset1.getSpace();

	/*
	 * Get the number of dimensions in the dataspace.
	 */
	int rank = dataspace1.getSimpleExtentNdims();

	/*
	 * Get the dimension size of each dimension in the dataspace and
	 * display them.
	 */
	hsize_t dims_out[1];
	int ndims = dataspace1.getSimpleExtentDims(dims_out, NULL);
	if (i1_stop >= dims_out[0] || i1_stop==-1) { //keep the i1_stop within the limit of the number of columns
		i1_stop = dims_out[0] - 1; //could be negative
	}
	noOfSelectedEntries = i1_stop - i1_start + 1;
	std::cerr << "input1 dimension:" << dims_out[0] << ". " << noOfSelectedEntries << " rows to be read." << endl;
	if (noOfSelectedEntries<=0){	//2012.3.28 exit if there is no data in input.
		std:cerr<< "Exit as there is "<< noOfSelectedEntries  << " data."<<endl;
		dataspace1.close();
		dataset1.close();
		input1.close();
		exit(0);
	}
	/*
	 * for the output
	 */
	outputDatasetRank = 1;
}

FindMaxLDBetweenPeakAndEachLocus::~FindMaxLDBetweenPeakAndEachLocus()
{
	fstLocusId2CorStruc.clear();
	sndLocusIdMap.clear();
	delete inputCorrelationInMemory;
	cleanupMemory();

	dataspace1.close();
	dataset1.close();
	input1.close();
}

void FindMaxLDBetweenPeakAndEachLocus::exitIfFileNotGood(char* inputFname, int exitCode){
	/*
	 * 2012.3.28 check whether file exists, it not exits the program
	 */
	ifstream my_file(inputFname);
	if (my_file.good())
	{
		my_file.close();
	}
	else{
		std::cerr<<"Can't read " << inputFname << ". It probably doesn't exist."<<std::endl;
		exit(exitCode);
	}
}

int FindMaxLDBetweenPeakAndEachLocus::readLocusIDLs(){
	std::cerr<<"Reading the locus ID list ...";
	H5File file = H5File(withinPeakLocusIDFname, H5F_ACC_RDONLY);
	DataSet dataset = file.openDataSet("locus_id_ls");

	DataSpace dataspace = dataset.getSpace();
	hsize_t dims_out[1];
	int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);
	int no_of_loci = dims_out[0];
	int* locusIDLs = new int[no_of_loci];
	dataset.read(locusIDLs, PredType::NATIVE_INT);

	for (int i=0; i<no_of_loci; i++){
		sndLocusIdMap[locusIDLs[i]] = i;

	}
	delete [] locusIDLs;
	dataspace.close();
	dataset.close();
	file.close();
	cerr << sndLocusIdMap.size() << " loci." << endl;
	return 0;
}

int FindMaxLDBetweenPeakAndEachLocus::readDataIntoMemory(DataSet &dataset, DataSpace &dataspace)
{
	std::cerr<<"Reading the input correlation data ...";
	/*
	 * Create a datatype
	 */
	CompType mtype1(sizeof(s1_t));
	mtype1.insertMember(outputTypeMember1, HOFFSET(s1_t, a), PredType::NATIVE_INT);
	mtype1.insertMember(outputTypeMember2, HOFFSET(s1_t, b), PredType::NATIVE_INT);
	mtype1.insertMember(outputTypeMember3, HOFFSET(s1_t, c), PredType::NATIVE_FLOAT);


	/*
	 * Define hyperslab in the dataset; implicitly giving strike and
	 * block NULL.
	 */
	noOfSelectedEntries = i1_stop-i1_start +1;
	hsize_t offset[outputDatasetRank]; // hyperslab offset in the file
	hsize_t count[outputDatasetRank]; // size of the hyperslab in the file
	offset[0] = i1_start;
	count[0] = noOfSelectedEntries;
	dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);

	/*
	 * Define the memory dataspace.
	 */
	hsize_t dimsm[outputDatasetRank]; /* memory space dimensions */
	dimsm[0] = noOfSelectedEntries;
	DataSpace memspace(outputDatasetRank, dimsm);

	/*
	 * Define memory hyperslab.
	 */
	hsize_t offset_out[outputDatasetRank]; // hyperslab offset in memory
	hsize_t count_out[outputDatasetRank]; // size of the hyperslab in memory
	offset_out[0] = 0;
	count_out[0] = noOfSelectedEntries;
	memspace.selectHyperslab(H5S_SELECT_SET, count_out, offset_out);

	/*
	 * Read two fields c and a from s1 dataset. Fields in the file
	 * are found by their names "c_name" and "a_name".
	 */
	inputCorrelationInMemory = new s1_t[noOfSelectedEntries];
	dataset.read(inputCorrelationInMemory, mtype1, memspace, dataspace);
	cerr << noOfSelectedEntries << " entries into memory."<< endl;

	// need the sndLocusIdMap
	readLocusIDLs();

	cerr << "Finding 2nd locus with max abs(cor) to 1st locus ...";
	s1_t oneCorrelationStruc;
	__gnu_cxx::hash_map<int, s1_t >::iterator fstLocusIter = fstLocusId2CorStruc.begin();
	__gnu_cxx::hash_map<int, int >::iterator sndLocusIter = sndLocusIdMap.begin();
	bool replaceExistingData = false;

	for (int i =0; i<noOfSelectedEntries; i++){
		oneCorrelationStruc = inputCorrelationInMemory[i];
		sndLocusIter = sndLocusIdMap.find(oneCorrelationStruc.b);	//check if 2nd locus is within the peak
		if (sndLocusIter!=sndLocusIdMap.end()){
			fstLocusIter = fstLocusId2CorStruc.find(oneCorrelationStruc.a);
			if (fstLocusIter!=fstLocusId2CorStruc.end()){
				if (abs(fstLocusIter->second.c)<abs(oneCorrelationStruc.c)){	//use absolute value
					replaceExistingData = true;
				}
			}
			else{
				replaceExistingData = true;
			}
			fstLocusId2CorStruc[oneCorrelationStruc.a] = oneCorrelationStruc;
		}
	}

	cerr << fstLocusId2CorStruc.size() << " 1st loci have max abs(cor) with 2nd loci set." << endl;
	return 0;
}


void FindMaxLDBetweenPeakAndEachLocus::cleanupMemory(){
;
}

void FindMaxLDBetweenPeakAndEachLocus::run(){
	readDataIntoMemory(dataset1, dataspace1);
	outputDatasetLength = fstLocusId2CorStruc.size();
	if (outputDatasetLength>0){
		output(outputFname);
	}
}

int FindMaxLDBetweenPeakAndEachLocus::output(H5std_string &outputFname){
	outputDatasetLength = fstLocusId2CorStruc.size();

	out = H5File(outputFname, H5F_ACC_TRUNC);
	/*
	 * Create the data space.
	 */
	hsize_t dim[] = {outputDatasetLength}; /* Dataspace dimensions */

	DataSpace space(outputDatasetRank, dim);

	/*
	 * Create the memory datatype.
	 */
	CompType mtype1(sizeof(s1_t));
	mtype1.insertMember(outputTypeMember1, HOFFSET(s1_t, a), PredType::NATIVE_INT);
	mtype1.insertMember(outputTypeMember2, HOFFSET(s1_t, b), PredType::NATIVE_INT);
	mtype1.insertMember(outputTypeMember3, HOFFSET(s1_t, c), PredType::NATIVE_FLOAT);

	/*
	 * add chunks and compression to the output data
	 */
	DSetCreatPropList propList = DSetCreatPropList();
	propList.setChunk(outputDatasetRank, dim);
	propList.setDeflate(4);

	/*
	 * Create the dataset.
	 */
	DataSet* dataset;
	dataset = new DataSet(out.createDataSet(outputDatasetName, mtype1, space, propList));

	outputDataInMemory = new s1_t[outputDatasetLength];
	__gnu_cxx::hash_map<int, s1_t >::iterator fstLocusIter = fstLocusId2CorStruc.begin();
	int i =0;
	for (;fstLocusIter!=fstLocusId2CorStruc.end();fstLocusIter++){
		outputDataInMemory[i].a = fstLocusIter->first;
		outputDataInMemory[i].b = fstLocusIter->second.b;
		outputDataInMemory[i].c = fstLocusIter->second.c;
		i++;
	}
	/*
	 * Write data to the dataset;
	 */
	dataset->write(outputDataInMemory, mtype1);

	/*
	 * Release resources
	 */
	delete outputDataInMemory;
	delete dataset;
	out.close();
	std::cerr << outputDatasetLength << " outputted."<<endl;
}


const char* program_name;

void print_usage(FILE* stream,int exit_code)
{
	assert(stream !=NULL);
		fprintf(stream,"Usage: %s options inputfile\n",program_name);
	fprintf(stream,"\t-h  --help	Display the usage infomation.\n"\
		"\t-o ..., --output=...	output HDF5. similar structure as input.\n"\
		"\t-i ..., --input1=...	input HDF5 file which contains correlation between 1st locus and 2nd locus\n"\
		"\t-j ..., --withinPeakLocusIDFname=...	a HDF5 file contains a list of 2nd locus ID\n"\
		"\t-s ..., --i1_start=...,	start index of the input HDF5. default: 0\n"\
		"\t-t ..., --i1_stop=...	stop index for the input HDF5. The column itself will be included. default: -1 (=end)\n"\
		"\tFor long option, = or ' '(blank) is same.\n"\
		"\tLine tokenizer is one space, tab, or \\r\n");
	exit(3);
}


int main(int argc, char* argv[])
{
	int next_option;
	const char* const short_options="ho:i:j:s:t:";
	const struct option long_options[]={
	  {"help",0,NULL,'h'},
	  {"output",1,NULL,'o'},
	  {"input1",1,NULL,'i'},
	  {"withinPeakLocusIDFname", 1, NULL, 'j'},
	  {"i1_start", 1, NULL, 's'},
	  {"i1_stop", 1, NULL, 't'},
	  {NULL,0,NULL,0}
	};

	program_name=argv[0];
	char* output_filename = NULL;
	char* input1Fname = NULL;
	char* withinPeakLocusIDFname = NULL;
	int i1_start = 0;
	int i1_stop = -1;



	do
	{
		next_option=getopt_long(argc,argv,short_options,long_options,NULL);
		switch(next_option)
		{
		case 'h':
			print_usage(stdout,0);
			exit(1);
		case 'o':
			output_filename=optarg;
			break;
		case 'i':
			input1Fname = optarg;
			break;
		case 'j':
			withinPeakLocusIDFname = optarg;
			break;
		case 's':
			i1_start = atoi(optarg);
			break;
		case 't':
			i1_stop = atoi(optarg);
			break;
		case '?':
			print_usage(stderr,-1);
		case -1:
			break;
		default:
			abort();
		}
	}while(next_option!=-1);

	if (optind == argc && argc>1)
	{
		FindMaxLDBetweenPeakAndEachLocus instance(input1Fname, withinPeakLocusIDFname, output_filename, i1_start, i1_stop);
		instance.run();

	}
	else
		print_usage(stderr, 1);
}
