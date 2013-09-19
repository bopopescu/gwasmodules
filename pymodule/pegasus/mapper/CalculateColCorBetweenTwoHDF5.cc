/****
 * 2012.2.29
 * 	a program to calculate correlation between two HDF5 SNP data matrix.
 * 	use ConvertSNPData2HDF5.py to convert StrainXSNP data into HDF5 format
 * ~/script/pymodule/pegasus/mapper/CalculateColCorBetweenTwoHDF5  -o /tmp/out -i /tmp/call_80.hdf5 -j /tmp/call_80.hdf5 -s 0 -t 10 -u 280 -v 300
 *
 */
#include "pymodule/include/CalculateColCorBetweenTwoHDF5.h"


CalculateColCorBetweenTwoHDF5::CalculateColCorBetweenTwoHDF5(char* _input1Fname, char* _input2Fname, char* outf_name, int _i1_start, int _i1_stop,
			int _i2_start, int _i2_stop, float cor_cut_off_given)\
			:outputFname(outf_name), i1_start(_i1_start), i1_stop(_i1_stop), i2_start(_i2_start), i2_stop(_i2_stop), min_cor(cor_cut_off_given)

{
	std::cerr<<"start ...";
	input1 = H5File(_input1Fname, H5F_ACC_RDONLY);
	input2 = H5File(_input2Fname, H5F_ACC_RDONLY);
	dataset1 = input1.openDataSet(DataMatrixDataSetName);
	dataset2 = input2.openDataSet(DataMatrixDataSetName);
	std::cerr<< endl;
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
	dataspace2 = dataset2.getSpace();

	/*
	 * Get the number of dimensions in the dataspace.
	 */
	int rank = dataspace1.getSimpleExtentNdims();

	/*
	 * Get the dimension size of each dimension in the dataspace and
	 * display them.
	 */
	hsize_t dims_out[2];
	int ndims = dataspace1.getSimpleExtentDims(dims_out, NULL);
	no_of_rows = dims_out[0];
	if (i1_stop>=dims_out[1]){	//keep the i1_stop within the limit of the number of columns
		i1_stop=dims_out[1]-1;	//could be negative
	}
	input1_no_of_cols = i1_stop-i1_start+1;
	std::cerr<<"input1 dimension:" << dims_out[0] << " X " << dims_out[1] << ". " << input1_no_of_cols << " columns to be read." << endl;

	ndims = dataspace2.getSimpleExtentDims(dims_out, NULL);
	if (i2_stop>=dims_out[1]){	//keep the i2_stop within the limit of the number of columns
		i2_stop=dims_out[1]-1;
	}
	input2_no_of_cols = i2_stop-i2_start+1;	//could be negative
	std::cerr<<"input2 dimension:" << dims_out[0] << " X " << dims_out[1] << ". " << input2_no_of_cols << " columns to be read."<< endl;


	/*
	 * Create the output file
	 */
	outputDatasetMaxLength = input1_no_of_cols*input2_no_of_cols;
	outputDatasetRank=1;

	//out.open(outf_name);
	//output the header
	//out<<"col_id1"<<"\t"<<"col_id2"<<"\t"<<"correlation"<<endl;


}


CalculateColCorBetweenTwoHDF5::~CalculateColCorBetweenTwoHDF5()
{
	dataspace1.close();
	dataspace2.close();
	dataset1.close();
	dataset2.close();
	input1.close();
	input2.close();
	//out.close();

	/*
	for(int i = 0; i < no_of_rows; ++i) {
		delete [] dataMatrix1InMemory[i];
		delete [] dataMatrix2InMemory[i];
	}
	*/
}

void CalculateColCorBetweenTwoHDF5::cleanupMemory()
{
	delete [] dataMatrix1InMemory;
	delete [] dataMatrix2InMemory;

	delete [] input1_col_id_ls;
	delete [] input2_col_id_ls;
}

double CalculateColCorBetweenTwoHDF5::cor(int* &data_matrix1, int* &data_matrix2, int matrix1_col_index, int matrix2_col_index, int &no_of_rows)
{
	float xx = 0.0;
	float yy = 0.0;
	float xy = 0.0;
	float mean_x = 0.0;
	float mean_y = 0.0;

	for (int i=0; i<no_of_rows; i++)
	{
		mean_x += data_matrix1[i*input1_no_of_cols+matrix1_col_index];
		mean_y += data_matrix2[i*input2_no_of_cols+matrix2_col_index];
	}

	mean_x /= no_of_rows;
	mean_y /= no_of_rows;

	int v1;
	int v2;
	for (int i=0; i<no_of_rows; i++)
	{
			v1 = data_matrix1[i*input1_no_of_cols + matrix1_col_index];
			v2 = data_matrix2[i*input2_no_of_cols + matrix2_col_index];
			xy += (v1-mean_x) * (v2-mean_y);
			xx += (v1-mean_x) * (v1-mean_x) ;
			yy += (v2-mean_y)  * (v2-mean_y) ;
	}

	float correlation;
	if(xx==0.0 || yy==0.0)
		//all NAN
		correlation = 1.1;
	else
		correlation = xy/(sqrt(xx*yy));
	//correlation is modeled as a t distribution of n-2 degree.
	//edge_data.degree = no_of_valids-2;
	//edge_data.significance = 0;
	return correlation;
}


void CalculateColCorBetweenTwoHDF5::run()
{
#if defined(DEBUG)
	cerr << input1_no_of_cols << endl;
	cerr << input2_no_of_cols << endl;
#endif
	if ( (input1_no_of_cols>0) && (input2_no_of_cols>0)){
		readColIDIntoMemory(input1, input1_col_id_ls,  i1_start, input1_no_of_cols);
		#if defined(DEBUG)
			cerr<<input1_col_id_ls[0]<<endl;
		#endif
		readColIDIntoMemory(input2, input2_col_id_ls,  i2_start, input2_no_of_cols);
		#if defined(DEBUG)
			cerr<<input2_col_id_ls[2]<<endl;
		#endif
		readDataIntoMemory(dataset1, dataspace1, dataMatrix1InMemory, i1_start, input1_no_of_cols);
		readDataIntoMemory(dataset2, dataspace2, dataMatrix2InMemory, i2_start, input2_no_of_cols);
		std::cerr<<"Calculating correlation " << input1_no_of_cols <<  " X " << input2_no_of_cols << " ...";
		outputDataInMemory = new outputStruct[outputDatasetMaxLength];	//over-allocate memory. because some pairs have correlation below min_cor.
		int i,j;
		int counter=0;
		int real_counter =0;
		int outputIndex;
		for (i=0;i<input1_no_of_cols;i++)
		{
			for (j=0; j<input2_no_of_cols; j++)
			{
				double correlation = cor(dataMatrix1InMemory, dataMatrix2InMemory, i, j, no_of_rows);
				if (abs(correlation)>=min_cor)
				{
					outputIndex = real_counter;
					outputDataInMemory[outputIndex].a = input1_col_id_ls[i];
					outputDataInMemory[outputIndex].b = input2_col_id_ls[j];
					outputDataInMemory[outputIndex].c = correlation;
					//out <<input1_col_id_ls[i]<<"\t"<<input2_col_id_ls[j]<<"\t"<<correlation<<endl;
					real_counter ++;
				}
				counter ++;
			}
		}
		cerr << counter << " pairs." << endl;
		outputDatasetLength = real_counter;
		output(outputDataInMemory, outputFname);
		cleanupMemory();
	}
	else{
		outputDatasetLength = 0;
	}

}

int CalculateColCorBetweenTwoHDF5::input()
{
	return 0;
}

int CalculateColCorBetweenTwoHDF5::readDataIntoMemory(DataSet &dataset, DataSpace &dataspace, int* &dataMatrixInMemory, int &col_start, int &no_of_cols)
{
	std::cerr<<"Reading the data matrix..." << col_start << "  " << no_of_cols << " ...";
	/*
	 * Define hyperslab in the dataset; implicitly giving strike and
	 * block NULL.
	 */
	hsize_t offset[2]; // hyperslab offset in the file
	hsize_t count[2]; // size of the hyperslab in the file
	offset[0] = 0;
	offset[1] = col_start;
	count[0] = no_of_rows;
	count[1] = no_of_cols;
	dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);

	/*
	 * Define the memory dataspace.
	 */
	hsize_t dimsm[2]; /* memory space dimensions */
	dimsm[0] = no_of_rows;
	dimsm[1] = no_of_cols;
	DataSpace memspace(2, dimsm);

	/*
	 * Define memory hyperslab.
	 */
	hsize_t offset_out[2]; // hyperslab offset in memory
	hsize_t count_out[2]; // size of the hyperslab in memory
	offset_out[0] = 0;
	offset_out[1] = 0;
	count_out[0] = no_of_rows;
	count_out[1] = no_of_cols;
	memspace.selectHyperslab(H5S_SELECT_SET, count_out, offset_out);

	/*
	 * Read data from hyperslab in the file into the hyperslab in
	 * memory and display the data.
	 */
	dataMatrixInMemory= new int[no_of_rows*no_of_cols];
	//dataMatrixInMemory = (int*) malloc(no_of_rows*no_of_cols * sizeof(int));
	/*
	for(int i = 0; i < no_of_rows; ++i) {
		dataMatrixInMemory[i] = new int[no_of_cols];
	}
	*/
	dataset.read(dataMatrixInMemory, PredType::NATIVE_INT, memspace, dataspace);
	/*		2012.2.29 check whether data is right
	int j,i;
	for (j = 0; j < 10; j++) {
		for (i = 0; i < no_of_cols; i++)
			cout << dataMatrixInMemory[j*no_of_cols+i] << " ";
		cout << endl;
	}
	*/
	std::cerr<<no_of_rows<<" X " << no_of_cols <<endl;

	return 0;
}

int CalculateColCorBetweenTwoHDF5::readColIDIntoMemory(H5File &h5file, int* &col_id_ls, int &col_start, int &no_of_cols)
{
	std::cerr<<"Reading the column IDs...";
	DataSet dataset = h5file.openDataSet(ColIDDataSetName);
	/*
	 * Get the class of the datatype that is used by the dataset.
	 */
	H5T_class_t type_class = dataset.getTypeClass();

	/*
	 * Get dataspace of the dataset.
	 */
	DataSpace dataspace = dataset.getSpace();
	/*
	 * Define hyperslab in the dataset; implicitly giving strike and
	 * block NULL.
	 */
	hsize_t offset[1]; // hyperslab offset in the file
	hsize_t count[1]; // size of the hyperslab in the file
	offset[0] = col_start;
	count[0] = no_of_cols;
	dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);

	/*
	 * Define the memory dataspace.
	 */
	hsize_t dimsm[1]; /* memory space dimensions */
	dimsm[0] = no_of_cols;
	DataSpace memspace(1, dimsm);

	/*
	 * Define memory hyperslab.
	 */
	hsize_t offset_out[1]; // hyperslab offset in memory
	hsize_t count_out[1]; // size of the hyperslab in memory
	offset_out[0] = 0;
	count_out[0] = no_of_cols;
	memspace.selectHyperslab(H5S_SELECT_SET, count_out, offset_out);

	/*
	 * Read data from hyperslab in the file into the hyperslab in
	 * memory and display the data.
	 */
	col_id_ls = new int[no_of_cols];
	dataset.read(col_id_ls, PredType::NATIVE_INT, memspace, dataspace);
	std::cerr<< no_of_cols <<endl;

	return 0;
}


int CalculateColCorBetweenTwoHDF5::output(outputStruct* &dataInMemory, H5std_string &outputFname)
{
	out = H5File(outputFname, H5F_ACC_TRUNC);
	/*
	 * Create the data space.
	 */
	hsize_t dim[] = {outputDatasetLength}; /* Dataspace dimensions */

	DataSpace space(outputDatasetRank, dim);

	/*
	 * Create the memory datatype.
	 */
	CompType mtype1(sizeof(outputStruct));
	mtype1.insertMember(outputTypeMember1, HOFFSET(outputStruct, a), PredType::NATIVE_INT);
	mtype1.insertMember(outputTypeMember2, HOFFSET(outputStruct, b), PredType::NATIVE_INT);
	mtype1.insertMember(outputTypeMember3, HOFFSET(outputStruct, c), PredType::NATIVE_FLOAT);

	/*
	 * add chunks and compression to the output data
	 */
	DSetCreatPropList propList = DSetCreatPropList();
	propList.setChunk(outputDatasetRank, dim );
	propList.setDeflate(4);
	/*
	 * Create the dataset.
	 */
	DataSet* dataset;
	dataset = new DataSet(out.createDataSet(outputDatasetName, mtype1, space, propList));

	/*
	 * Write data to the dataset;
	 */
	if (outputDatasetLength>0){
		dataset->write(outputDataInMemory, mtype1);
	}

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
	fprintf(stream,"This program calculates correlation between columns of the data matrices from two input. Their rows must match.\n"\
			"\t-h  --help	Display the usage infomation.\n"\
		"\t-o ..., --output=...	output HDF5 to store the pairwise correlation\n"\
		"\t-i ..., --input1Fname=...	1st polymorphism dataset in HDF5\n"\
		"\t-j ..., --input2Fname=...	2nd polymorphism dataset in HDF5\n"\
		"\t-c ..., --min_cor=...	minimum abs(correlation) for output, -1.0(default)\n"\
		"\t-s ..., --i1_start=...,	starting locus index for the 1st polymorphism dataset, 0(default)\n"\
		"\t-t ..., --i1_stop=...	stop locus index for the 1st polymorphism dataset (inclusive), 0(default)\n"\
		"\t-u ..., --i2_start=...	starting locus index for the 2nd polymorphism dataset, 0(default)\n"\
		"\t-v ..., --i2_stop=...	stop locus index for the 2nd polymorphism dataset (inclusive), 0(default)\n"\
		"\tFor long option, = or ' '(blank) is same.\n"\
		"\tLine tokenizer is one space, tab, or \\r\n");
	exit(3);
}


int main(int argc, char* argv[])
{
	int next_option;
	const char* const short_options="ho:i:j:s:t:u:v:c:";
	const struct option long_options[]={
	  {"help",0,NULL,'h'},
	  {"output",1,NULL,'o'},
	  {"input1",1,NULL,'i'},
	  {"input2", 1, NULL, 'j'},
	  {"i1_start", 1, NULL, 's'},
	  {"i1_stop", 1, NULL, 't'},
	  {"i2_start", 1, NULL, 'u'},
	  {"i2_stop", 1, NULL, 'v'},
	  {"min_cor", 0, NULL, 'c'},
	  {NULL,0,NULL,0}
	};

	program_name=argv[0];
	char* output_filename = NULL;
	char* input1Fname = NULL;
	char* input2Fname = NULL;
	int i1_start = 0;
	int i1_stop = 0;

	int i2_start = 0;
	int i2_stop = 0;
	float min_cor = -1.0;


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
			input2Fname = optarg;
			break;
		case 's':
			i1_start = atoi(optarg);
			break;
		case 't':
			i1_stop = atoi(optarg);
			break;
		case 'u':
			i2_start = atoi(optarg);
			break;
		case 'v':
			i2_stop = atoi(optarg);
			break;
		case 'c':
			min_cor= atof(optarg);
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
		CalculateColCorBetweenTwoHDF5 instance(input1Fname, input2Fname, output_filename, i1_start, i1_stop, i2_start, i2_stop,
				min_cor);
		instance.run();

	}
	else
		print_usage(stderr, 1);
}
