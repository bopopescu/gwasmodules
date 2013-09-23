/*
 GADA v1.0 Genome Alteration Detection Algorithm
 Copyright (C) 2008  Childrens Hospital of Los Angeles

 GADA is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 any later version.

 GADA is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with GADA.  If not, see <http://www.gnu.org/licenses/>.

 Author:
 Roger Pique-Regi    piquereg@usc.edu

 */
#include "BaseGADA.h"

#include <boost/program_options.hpp>	//for program options
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>

using namespace std;
using namespace boost;
namespace po = boost::program_options;


class RogerGADA{
	// Global variables with algorithm parameters

	double T; //Backward elimination threshold
	//double T2=5.0; //Segment collapse to base Non-Alteration level threshold
	double BaseAmp; //Base-level
	double a; //SBL hyperprior parameter
	double sigma2; //Variance observed, if negative value, it will be estimated by the mean of the differences
	// I would recommend to be estimated on all the chromosomes and as a trimmed mean.
	long MinLen; //Lenght in number of probes for a CNA segment to be called significan.
	long SelectClassifySegments; //Classify segment into altered state (1), otherwise 0
	long SelectEstimateBaseAmp; //Estimate Neutral hybridization amplitude.
	char *InputFile;
	char *OutputFile;
public:

	RogerGADA(){
		T = 5.0; //Backward elimination threshold
		//double T2=5.0; //Segment collapse to base Non-Alteration level threshold
		BaseAmp = 0.0; //Base-level
		a = 0.2; //SBL hyperprior parameter
		sigma2 = -1; //Variance observed, if negative value, it will be estimated by the mean of the differences
		// I would recommend to be estimated on all the chromosomes and as a trimmed mean.
		MinLen = 0; //Lenght in number of probes for a CNA segment to be called significan.
		SelectClassifySegments = 0; //Classify segment into altered state (1), otherwise 0
		SelectEstimateBaseAmp = 1; //Estimate Neutral hybridization amplitude.
	}
	~RogerGADA(){};
	void help_message(FILE *fd) {
		//fprintf(fd,"# Welcome to GADA 1.0 \n");
		fprintf(fd, "# Usage:\n");
		fprintf(fd,
				"# \t GADA [-T 5] [-a 0.8] [-s -1] [-M 3] -i input.txt -o output.txt\n");
		fprintf(fd, "# \t input.txt is a single column text file with no header\n");
		fprintf(fd, "# Possible options:\n");
		fprintf(fd, "# \t -i\t Input file. Otherwise standard input assumed\n");
		fprintf(fd, "# \t -o\t Output file. Otherwise standard output assumed\n");
		fprintf(fd,
				"# \t -a\t is the SBL hyperprior parameter for a breakpoint. It is the \n");
		fprintf(fd,
				"# \t\t shape parameter of the Gamma distribution. Default value %g.\n",
				a);
		fprintf(fd,
				"# \t\t Higher (lower) value less (more) breakpoints expected a priori\n");
		fprintf(fd,
				"# \t -T\t is the backward elimination critical value for a breakpoint. \n");
		fprintf(fd,
				"# \t\t i.e. the minimum difference between the (mean1-mean2)/segment divided\n");
		fprintf(fd, "# \t\t by sigma. The default value for T is %g\n", T);
		fprintf(fd,
				"# \t -M\t is the minimum size in number of probes for a segment to be \n");
		fprintf(fd, "# \t\t deemed significant. Default value is %ld.\n", MinLen);
		fprintf(fd,
				"# \t -s\t The variance estimate for the noise. If not provided or \n");
		fprintf(fd,
				"# \t\t negative it will be estimated from the provided data. \n");
		fprintf(fd,
				"# \t\t We recomend to estimate this on the entire data and not \n");
		fprintf(fd, "# \t\t separately for each chromosome.\n");
		fprintf(fd,
				"# \t -c\t Classify segments into altered state (L)oss -- (N)eutral \n");
		fprintf(fd,
				"# \t\t (G)ain). If c option is not specified, segments are returned \n");
		fprintf(fd, "# \t\t with their mean.\n");
		fprintf(fd,
				"# \t -b\t Mean amplitude associated to the Neutral state. If not \n");
		fprintf(fd,
				"# \t\t provided, and c option is used, then it is estimated as the\n");
		fprintf(fd,
				"# \t\t median value of all probes hybridization values after running \n");
		fprintf(fd,
				"# \t\t the algorithm. We recomend to estimate this on chromosomes\n");
		fprintf(fd,
				"# \t\t that are known to have a Neutral state on most areas. In some\n");
		fprintf(fd,
				"# \t\t cases this value may be known if we have been applied some \n");
		fprintf(fd,
				"# \t\t normalization, preprocessing or using another sample as ref.\n");
		fprintf(fd, "# \t -h \t Prints this help message.\n");
	}

	void help_and_exit(FILE *fd, int code) {
		fprintf(fd, "Invalid syntax. Use GADA -h for help\n");
		exit(code);
	}

	void Configure(int ac, char *av[], string &input_fname, string &output_fname) {
	#if defined(DEBUG)
		cerr<< boost::format("# Getting commandline arguments ...");
	#endif
		long CLcount, i;
		FILE *fd;

		// Parse the command line
		CLcount = 1;

		fd = stdout;
		fprintf(fd, "# NumArgs = %d \n", ac);
		fprintf(fd, "# CallStr = ");

		for (i = 0; i < ac; i++)
			fprintf(fd, "%s ", av[i]);
		fprintf(fd, "\n# Parsing Arguments: \n");

		while (CLcount < ac) {
			if (0 == strncmp(av[CLcount], "-h", 2)) {
				help_message(stderr);
				exit(0);
			}

			else if (0 == strncmp(av[CLcount], "-c", 2)) {
				SelectClassifySegments = 1;
				CLcount += 1;
				fprintf(fd,
						"# -c option activated to classify segments into altered states\n");
			} else if (0 == strncmp(av[CLcount], "-i", 2)) //! Input file
					{
				if (0 == strncmp(av[CLcount + 1], "-", 1))
					help_and_exit(stderr, 1);
				InputFile = (char*) malloc(strlen(av[CLcount + 1]));
				//strcpy(InputFile, av[CLcount+1]);
				input_fname = av[CLcount + 1];
				cout << boost::format("# Input file: %1%") % input_fname
						<< std::endl;
				InputFile = av[CLcount + 1];
				//fprintf(fd,"# Input file: %s \n",InputFile);
				CLcount += 2;
			} else if (0 == strncmp(av[CLcount], "-o", 2)) //! Output File
					{
				if (0 == strncmp(av[CLcount + 1], "-", 1))
					help_and_exit(stderr, 1);
				OutputFile = (char*) malloc(strlen(av[CLcount + 1]));
				//strcpy(OutputFile, av[CLcount+1]);	// 2009-11-10 strcpy is a dangerous function.
				// One weird bug is if it's /tmp/GADA/GADA_output_48 (48 could be any 2-digit), it'll cause "malloc: memory corruption" in fin = fopen(input_fname.c_str(), "r");.
				output_fname = av[CLcount + 1];
				OutputFile = av[CLcount + 1];
				CLcount += 2;
				fprintf(fd, "# Output file: %s \n", OutputFile);
			} else if (0 == strncmp(av[CLcount], "-a", 2)) //! a parameter
					{
				if (0 == strncmp(av[CLcount + 1], "-", 1))
					help_and_exit(stderr, 1);
				sscanf(av[CLcount + 1], "%lf", &a);
				CLcount += 2;
				fprintf(fd, "# a= %g \n", a);
			} else if (0 == strncmp(av[CLcount], "-T", 2)) //! T parameter
					{
				if (0 == strncmp(av[CLcount + 1], "-", 1))
					help_and_exit(stderr, 1);
				sscanf(av[CLcount + 1], "%lf", &T);
				CLcount += 2;
				fprintf(fd, "# T= %g \n", T);
			} else if (0 == strncmp(av[CLcount], "-M", 2)) //! MinLen parameter
					{
				if (0 == strncmp(av[CLcount + 1], "-", 1))
					help_and_exit(stderr, 1);
				sscanf(av[CLcount + 1], "%ld", &MinLen);
				CLcount += 2;
				fprintf(fd, "# MinLen= %ld \n", MinLen);
			} else if (0 == strncmp(av[CLcount], "-b", 2)) //! BaseAmp parameter
					{
				//if(0 == strncmp (av[CLcount+1], "-", 1))help_and_exit(stderr,1);
				sscanf(av[CLcount + 1], "%lf", &BaseAmp);
				CLcount += 2;
				fprintf(fd, "# BaseAmp= %g \n", BaseAmp);
				SelectEstimateBaseAmp = 0;
			} else if (0 == strncmp(av[CLcount], "-s", 3)) //! sigma2 parameter
					{
				//if(0 == strncmp (av[CLcount+1], "-", 1))help_and_exit(stderr,1);
				sscanf(av[CLcount + 1], "%lf", &sigma2);
				CLcount += 2;
				fprintf(fd, "# sigma2= %g \n", sigma2);
			} else {
				help_and_exit(stderr, 1);
			}
		}
	#if defined(DEBUG)
		cerr<< "Done.\n";
	#endif
	}

	void run(int argc, char *argv[]) {
	#if defined(DEBUG)
		cerr<< boost::format("# Starting ....\n");
	#endif
		long M = 1000;
		long i;
		long *Iext;
		long K;
		double *tn;
		long *SegLen;
		double *SegAmp;
		double *SegState;
		double *Wext;
		long debug = 0;

		FILE *fin, *fout;

		string input_fname;
		string output_fname;

		tn = (double *) calloc(M, sizeof(double));

		Configure(argc, argv, input_fname, output_fname);

		// fin = stdin;
		fin = fopen(input_fname.c_str(), "r");
		// ifstream in;
		// in.open(InputFile);
		// fout = stdout;
		fout = fopen(OutputFile, "w");
		//ofstream out;
		//out.open(output_fname.c_str());

		fprintf(fout, "# GADA v1.0 Genome Alteration Detection Algorithm\n");
		fprintf(fout, "# Copyright (C) 2008  Childrens Hospital of Los Angeles\n");
		fprintf(fout, "# author: Roger Pique-Regi piquereg@usc.edu\n");

		// cout << boost::format("# Parameter setting: a=%1%,T=%2%,MinLen=%3%,sigma2=%4%,BaseAmp=%5%.") % a % T % MinLen % sigma2 % BaseAmp << std::endl;
		fprintf(fout,
				"# Parameter setting: a=%g,T=%g,MinLen=%ld,sigma2=%g,BaseAmp=%g\n",
				a, T, MinLen, sigma2, BaseAmp);

		/*
		 #if defined(DEBUG)
		 std::cerr<<"Read in the data...";
		 #endif
		 string line;
		 while(getline(in, line))
		 {
		 cout << line[0];
		 }
		 */
		i = 0;
		while (!feof(fin)) {
			fscanf(fin, "%lf", &tn[i++]);
			if (i >= M) {
				M = M + 1000;
				tn = (double *) realloc(tn, M * sizeof(double));
			}
		}
		M = i - 1;
		tn = (double*) realloc(tn, M * sizeof(double));

		fprintf(fout, "# Reading M=%ld probes in input file\n", M);

		double delta;
		long numEMsteps;
		long noOfBreakpointsAfterSBL;
		K = SBLandBE(tn, M, &sigma2, a, 0.0, 0, &Iext, &Wext, debug, delta, numEMsteps, noOfBreakpointsAfterSBL, 1E-10, 50000, 1E8,1E-20);

		fprintf(fout, "# Overall mean %g\n", Wext[0]);
		fprintf(fout, "# Sigma^2=%g\n", sigma2);
		fprintf(fout, "# Found %ld breakpoints after SBL\n", K);

		BEwTandMinLen(Wext, Iext, &K, sigma2, T, MinLen, 0);
		fprintf(fout, "# Kept %ld breakpoints after BE\n", K);

		SegLen = (long*) calloc(K + 1, sizeof(long));
		SegAmp = (double *) calloc(K + 1, sizeof(double));
		IextToSegLen(Iext, SegLen, K);
		IextWextToSegAmp(Iext, Wext, SegAmp, K);
		fprintf(fout, "# Making segments\n");

		//Collapse Segments
		if (SelectClassifySegments == 1) {
			if (SelectEstimateBaseAmp == 1) {
				BaseAmp = CompBaseAmpMedianMethod(SegLen, SegAmp, K);
				fprintf(fout, "# Estimating BaseAmp\n");
			}
			fprintf(fout, "# BaseAmp=%g \n", BaseAmp);
			fprintf(fout, "# Classify Segments \n", BaseAmp);

			SegState = (double *) calloc(K + 1, sizeof(double));
			for (i = 0; i <= K; i++)
				SegState[i] = SegAmp[i];
			CollapseAmpTtest(SegState, SegLen, K, BaseAmp, sigma2, T);
		}

		if (SelectClassifySegments == 0) {
			fprintf(fout, "Start\tStop\tLength\tAmpl\n");
			for (i = 0; i < K + 1; i++)
				fprintf(fout, "%ld\t%ld\t%ld\t%g\n", Iext[i] + 1, Iext[i + 1],
						SegLen[i], SegAmp[i]);
		} else if (SelectClassifySegments == 1) {
			fprintf(fout, "Start\tStop\tLenght\tAmpl\tState\n");
			for (i = 0; i < K + 1; i++) {
				fprintf(fout, "%ld\t%ld\t%ld\t%g\t", Iext[i] + 1, Iext[i + 1],
						SegLen[i], SegAmp[i]);
				if (SegState[i] > BaseAmp)
					fprintf(fout, "G");
				else if (SegState[i] < BaseAmp)
					fprintf(fout, "L");
				else
					fprintf(fout, "N");
				fprintf(fout, "\n");
			}

		}

	}

};

class GADA {
public:
	int argc;
	char** argv;		//2013.08.20 somehow, can't use "char* argv[]" even though they are same.
	//otherwise "argv=_argv;" will not work for this error: incompatible types in assignment of ‘char**’ to ‘char* [0]’
	string programName;
	boost::format usageDoc;
	boost::format examplesDoc;

	long M;
	long i;
	long *Iext;
	long K;
	double *tn;
	long *SegLen;
	double *SegAmp;
	double *SegState;
	double *Wext;
	double delta;
	long numEMsteps;
	long noOfBreakpointsAfterSBL;

	int report;
	long debug; //verbosity... set equal to 1 to see messages of SBLandBE(). 0 to not see them
	double T; //Backward elimination threshold
	//double T2=5.0; //Segment collapse to base Non-Alteration level threshold
	double BaseAmp; //Base-level
	double a; //SBL hyperprior parameter
	double sigma2; //Variance observed, if negative value, it will be estimated by the mean of the differences
				   // I would recommend to be estimated on all the chromosomes and as a trimmed mean.
	long MinLen; //Length in number of probes for a CNA segment to be called significan.
	long SelectClassifySegments; //Classify segment into altered state (1), otherwise 0
	long SelectEstimateBaseAmp; //Estimate Neutral hybridization amplitude.
	double convergenceDelta;	//1E-10 or 1E-8 seems to work well for this parameter. -- => ++ conv time
			//1E8 better than 1E10 seems to work well for this parameter. -- => -- conv time
	long maxNoOfIterations;	//=50000, //10000 is enough usually
	double convergenceMaxAlpha;	// 1E8 Maximum number of iterations to reach convergence...
	double convergenceB;	// a number related to convergence = 1E-20

	string inputFname;
	string outputFname;

	GADA(int _argc, char* _argv[]);	//2013.08.28 commandline version
	GADA();
	GADA(long _debug);
	virtual ~GADA(){
		//cleanupMemory();	//comment it out to avoid repetitive free, in python module cleanupMemory is called by the end of run().
	}

	void initParameters();
//#ifdef GADABIN

#ifndef GADABIN	// 2009-11-21 boost python module code included under if macro GADABIN (GADA standalone) is not defined.
	void readInIntensity(boost::python::list intensity_list);
	boost::python::list run(boost::python::list intensity_list, double aAlpha,
			double TBackElim, long MinSegLen);
#else
	// 2013.08.28 stream causes error " note: synthesized method ... required here" because stream is noncopyable.
	po::options_description optionDescription;	//("Allowed options")
	po::positional_options_description positionOptionDescription;
	po::variables_map optionVariableMap;

	std::ifstream inputFile;
	std::ofstream outputFile;
	boost::iostreams::filtering_streambuf<boost::iostreams::input> inputFilterStreamBuffer;
	boost::iostreams::filtering_streambuf<boost::iostreams::output> outputFilterStreamBuffer;

	// 2013.08.28 stream causes error " note: synthesized method ... required here" because stream is noncopyable.
	virtual void constructOptionDescriptionStructure();
	virtual void parseCommandlineOptions();
	virtual void openOutputFile();
	virtual void openOneInputFile(string &inputFname,  boost::iostreams::filtering_streambuf<boost::iostreams::input> &inputFilterStreamBuffer);
	virtual void closeFiles();
	void commandlineRun();
#endif
	void cleanupMemory(){
		free(tn);
		free(Iext);
		free(SegLen);
		free(SegAmp);
		free(Wext);
		//free(SegState);	//2013.08.30 SegState is not always allocated with extra memory
	}
};

GADA::GADA() :
		debug(0) {
	initParameters();
}

GADA::GADA(int _argc, char* _argv[]) : argc(_argc), argv(_argv)
	{
	debug=0;
	report=0;
	programName = _argv[0];
	//strcpy(_argv[0], programName.c_str());	//2013.08.20 does no work

	cerr << "program name is " << programName << "." <<endl;

	usageDoc = boost::format("%1% -i INPUTFNAME -o OUTPUTFNAME [OPTIONS]\n")% programName;
	examplesDoc = boost::format("%1% -i /tmp/input.tsv.gz -o /tmp/output.gz -M 1000 --convergenceDelta 0.01 \n")% programName;
}

GADA::GADA(long _debug) :
		debug(_debug) {
	initParameters();
}


void GADA::initParameters() {
	T = 5.0; //Backward elimination threshold
	//double T2=5.0; //Segment collapse to base Non-Alteration level threshold
	BaseAmp = 0.0; //Base-level
	a = 0.2; //SBL hyperprior parameter
	sigma2 = -1; //Variance observed, if negative value, it will be estimated by the mean of the differences
				 // I would recommend to be estimated on all the chromosomes and as a trimmed mean.
	MinLen = 0; //Lenght in number of probes for a CNA segment to be called significan.
	SelectClassifySegments = 0; //Classify segment into altered state (1), otherwise 0
	SelectEstimateBaseAmp = 1; //Estimate Neutral hybridization amplitude.

	convergenceDelta=1E-8;	// or 1E-8 seems to work well for this parameter. -- => ++ conv time
				//1E8 better than 1E10 seems to work well for this parameter. -- => -- conv time
	maxNoOfIterations=50000; //10000 is enough usually
	convergenceMaxAlpha=1E8; // 1E8 Maximum number of iterations to reach convergence...
	convergenceB=1E-20;	// a number related to convergence = 1E-20
}


#ifndef GADABIN	// 2009-11-21 boost python module code included under if macro GADABIN (GADA standalone) is not defined.

void GADA::readInIntensity(boost::python::list intensity_list)
/*
 * 2010-6-8 fix a bug in moving intensity from intensity_list to tn: the last intensity was forgotten.
 */
{
	#if defined(DEBUG)
		cerr<< boost::format("# Start reading ... \n");
	#endif
	long no_of_probes = boost::python::extract<long>(
			intensity_list.attr("__len__")());
	M = no_of_probes;
	tn = (double *) calloc(no_of_probes, sizeof(double));

	for (long i = 0; i < no_of_probes; i++) {
		tn[i] = boost::python::extract<double>(intensity_list[i]);
	}

	#if defined(DEBUG)
		cerr<< boost::format("# M=%1% probes in input file\n") % M;
	#endif
}

boost::python::list GADA::run(boost::python::list intensity_list, double aAlpha,
		double TBackElim, long MinSegLen) {
	readInIntensity(intensity_list);
	// 2010-6-8 sigma2 has to be set here. sigma2 set in GADA::initParameters() refers to the GADA::sigma2;
	// here i suspect is the global one. Without the sentence below, the results seem to be weird.
	// It used to be fine without it.
	sigma2 = -1; //Variance observed, if negative value, it will be estimated by the mean of the differences
	K = SBLandBE(tn, M, &sigma2, aAlpha, TBackElim, MinSegLen, &Iext, &Wext, debug, delta, numEMsteps, noOfBreakpointsAfterSBL, convergenceDelta, maxNoOfIterations, convergenceMaxAlpha, convergenceB);

#if defined(DEBUG)
	cerr<< boost::format("# Overall mean %1%.\n")%Wext[0];
	cerr<< boost::format("# Sigma^2=%1%.\n")%sigma2;
	cerr<< boost::format("# Found %1% breakpoints after SBL\n")%noOfBreakpointsAfterSBL;
#endif

	//BEwTandMinLen(Wext, Iext, &K, sigma2, TBackElim, MinSegLen, debug);

#if defined(DEBUG)
	cerr<< boost::format("# Kept %1% breakpoints after BE\n")%K;
#endif

	SegLen = (long*) calloc(K + 1, sizeof(long));
	SegAmp = (double *) calloc(K + 1, sizeof(double));
	IextToSegLen(Iext, SegLen, K);
	IextWextToSegAmp(Iext, Wext, SegAmp, K);
#if defined(DEBUG)
	cerr<< boost::format("# Making segments\n");
#endif

	boost::python::list return_ls;

	//fprintf(fout,"Start\tStop\tLength\tAmpl\n");
	for (i = 0; i < K + 1; i++) {
		boost::python::list d_row;
		d_row.append(Iext[i] + 1);
		d_row.append(Iext[i + 1]);
		d_row.append(SegLen[i]);
		d_row.append(SegAmp[i]);
		//cerr<< boost::format("%1% \t %2% \t %3% %4% \n")%Iext[i] % Iext[i+1] % SegLen[i] % SegAmp[i];
		return_ls.append(d_row);
	}
	cleanupMemory();	//release tn, SegLen, SegAmp, and others
	return return_ls;

}

BOOST_PYTHON_MODULE(GADA)
{
	using namespace boost::python;
	class_<GADA>("GADA").def(init<long>()).def("run", &GADA::run);

}
#else
// 2013.08.28 stream causes error " note: synthesized method ... required here" because stream is noncopyable.

void GADA::constructOptionDescriptionStructure(){
	optionDescription.add_options()("help,h", "produce help message")
				("TBackElim,T", po::value<double>(&T)->default_value(5.0),
					" is the backward elimination critical value for a breakpoint. i.e. minimum (mean1-mean2)/stddev difference between two adjacent segments.")
				("aAlpha,a", po::value<double>(&a)->default_value(0.5),
					"is the SBL hyper prior parameter for a breakpoint. It is the  shape parameter of the Gamma distribution. Higher (lower) value means less (more) breakpoints expected a priori.")
				("MinSegLen,M", po::value<long>(&MinLen)->default_value(0),
					"is the minimum size in number of probes for a segment to be deemed significant.")
				("BaseAmp", po::value<double>(&BaseAmp)->default_value(0.0),
					"Mean amplitude associated to the Neutral state. If not "
					"provided, and c option is used, then it is estimated as the "
					"median value of all probes hybridization values after running "
					"the algorithm. We recomend to estimate this on chromosomes "
					"that are known to have a Neutral state on most areas. In some "
					"cases this value may be known if we have been applied some "
					"normalization, preprocessing or using another sample as ref.")
				("sigma2,s", po::value<double >(&sigma2)->default_value(-1),
					"Variance observed, if negative value, it will be estimated by the mean of the differences. "
					"I would recommend to be estimated on all the chromosomes and as a trimmed mean.")
				("SelectClassifySegments,c", po::value<long>(&SelectClassifySegments)->default_value(0),
						"Classify segment into altered state (1), otherwise 0")
				("SelectEstimateBaseAmp", po::value<long>(&SelectEstimateBaseAmp)->default_value(1),
						"missing data notation. missing data will be skipped.")
				("convergenceDelta", po::value<double>(&convergenceDelta)->default_value(1E-8),
							"a delta number controlling convergence")
				("maxNoOfIterations", po::value<long>(&maxNoOfIterations)->default_value(50000),
							"maximum number of iterations for EM convergence algorithm to run before being stopped")
				("convergenceMaxAlpha", po::value<double>(&convergenceMaxAlpha)->default_value(1E8),
							"one convergence related number, not sure what it does.")
				("convergenceB", po::value<double>(&convergenceB)->default_value(1E-20),
							"one convergence related number, not sure what it does")
				("debug,b", "toggle debug mode")
				("report,r", "toggle report mode")
				("inputFname,i", po::value<string >(&inputFname),
						"input filename, gzipped or not. could be specified as option or positional argument."
						"It is a single column text file with no header.")
				("outputFname,o", po::value<string>(&outputFname), "output filename");

}

void GADA::parseCommandlineOptions(){
	//all positional arguments are input files.
	positionOptionDescription.add("inputFname", -1);

	po::store(po::command_line_parser(argc, argv).
				  options(optionDescription).positional(positionOptionDescription).run(), optionVariableMap);

	//po::store(po::parse_command_line(argc, argv, optionDescription), optionVariableMap);
	po::notify(optionVariableMap);
	if (optionVariableMap.count("help") || inputFname.empty() || outputFname.empty()){
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

}


void  GADA::openOneInputFile(string &inputFname, boost::iostreams::filtering_streambuf<boost::iostreams::input> &inputFilterStreamBuffer){

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


void GADA::openOutputFile(){
	if (!outputFname.empty()){
			if (debug){
				std::cerr<<"Open file " << outputFname << " for writing " ;
			}
			int outputFnameLength = outputFname.length();
			if (outputFname.substr(outputFnameLength-3, 3)==".gz"){
				//boost::iostreams::gzip_compressor gzipCompressor;
				outputFilterStreamBuffer.push(boost::iostreams::gzip_compressor());
				//outputFilterStreamBuffer.push(boost::iostreams::base64_encoder());
				//outputFile.open(outputFname.c_str(), std::ios::out | std::ios::binary);
				outputFile.open(outputFname.c_str(), std::ios::out | std::ios::binary);
			}
			else{
				outputFile.open(outputFname.c_str(), std::ios::out);
			}
			outputFilterStreamBuffer.push(outputFile);
			//outputStream.rdbuf(&outputFilterStreamBuffer);
			//(&outputFilterStreamBuffer);
			if (debug){
				std::cerr<< endl ;
			}
		}
	else{
		if (debug){
			std::cerr << "Warning: Output file, " << outputFname << ", is an empty string." << endl;
		}

	}
}


void GADA::closeFiles(){
	//
	inputFile.close();
	if (outputFile!=NULL){
		if (debug){
			std::cerr<<"closing files " << "..." ;
		}
		//delete outputFilterStreamBuffer;
		outputFile.flush();
		//outputFile.close();	//2013.08.21 if output is a .gz, closing it here would result a truncated .gz file. don't know why.
		//maybe because the buffer or the gzip compressor filter needs to be closed beforehand.
	}

}

void GADA::commandlineRun(){
	if (debug){
		std::cerr<<"Entering GADA.commandlineRun()..." << std::endl;
	}

	constructOptionDescriptionStructure();
	parseCommandlineOptions();

	if (debug){
		std::cerr<<"Reading data from " << inputFname << " ... ";
	}
	openOneInputFile(inputFname, inputFilterStreamBuffer);
	std::istream inputStream(&inputFilterStreamBuffer);

	i = 0;
	M=1000;
	tn = (double *) calloc(M, sizeof(double));
	std::string line;
	std::getline(inputStream, line);
	while (!line.empty()){
		tn[i++] = atof(line.c_str());
		if (i >= M) {
			M = M + 1000;
			tn = (double *) realloc(tn, M * sizeof(double));
		}
		std::getline(inputStream, line);
	}
	M = i;
	tn = (double*) realloc(tn, M * sizeof(double));
	if (debug){
		std::cerr<< M << " data points." << endl;
	}

	if (debug){
			std::cerr<< "Running SBLandBE ... " << endl;
	}
	K = SBLandBE(tn, M, &sigma2, a, T, MinLen, &Iext, &Wext, debug , delta, numEMsteps, noOfBreakpointsAfterSBL, convergenceDelta, maxNoOfIterations, convergenceMaxAlpha, convergenceB);
	if (debug){
		std::cerr<< boost::format(" %1% breakpoints after SBL, %2% breakpoints after BE.\n")%noOfBreakpointsAfterSBL % K;
	}
	if (debug){
		//std::cerr<< boost::format("Backward elimination (T=%2%) and remove segments that are shorter than %1% ... ") % MinLen % T;
		std::cerr<< boost::format("Starting IextToSegLen() & IextWextToSegAmp() ... ");
	}

	//BEwTandMinLen(Wext, Iext, &K, sigma2, T, MinLen, debug);

	SegLen = (long*) calloc(K + 1, sizeof(long));
	SegAmp = (double *) calloc(K + 1, sizeof(double));
	IextToSegLen(Iext, SegLen, K);
	IextWextToSegAmp(Iext, Wext, SegAmp, K);

	if (debug){
		std::cerr<< " IextToSegLen() & IextWextToSegAmp() done." << endl;
	}

	if (debug){
		std::cerr<< "Outputting final result ... ";
	}
	openOutputFile();	//2013.08.30 open this file here. Do not open it way before the main writing starts.
	//it will leave a long period of zero-writing-activity (due to computation), which could hang the program sometimes on panfs system

	std::ostream outputStream(&outputFilterStreamBuffer);
	outputStream << "# GADA v1.0 Genome Alteration Detection Algorithm\n";
	outputStream << "# Copyright (C) 2008  Childrens Hospital of Los Angeles\n";
	outputStream << "# author: Roger Pique-Regi piquereg@usc.edu, Yu Huang polyactis@gmail.com\n";
	outputStream << boost::format("# Parameters: a=%1%,T=%2%,MinLen=%3%,sigma2=%4%,BaseAmp=%5%, convergenceDelta=%6%, maxNoOfIterations=%7%, convergenceMaxAlpha=%8%, convergenceB=%9%.")%
			a % T % MinLen % sigma2 % BaseAmp % convergenceDelta % maxNoOfIterations % convergenceMaxAlpha % convergenceB << std::endl;
	outputStream << boost::format("# Reading M=%1% probes in input file\n")% M;
	outputStream << boost::format("# Overall mean %1%\n")%Wext[0];
	outputStream << boost::format("# Sigma^2=%1%\n")%sigma2;
	outputStream << boost::format("# Convergence: delta=%1% after %2% EM iterations.\n")% delta % numEMsteps;
	outputStream << boost::format("# Found %1% breakpoints after SBL\n")% noOfBreakpointsAfterSBL;
	outputStream<< boost::format("# Kept %1% breakpoints after BE\n")%K;

	//Collapse Segments
	if (SelectClassifySegments == 1) {
		if (debug){
			std::cerr<<" Select and Classify Segments ...";
		}
		if (SelectEstimateBaseAmp == 1) {
			outputStream<< boost::format("# Estimating BaseAmp\n");
			BaseAmp = CompBaseAmpMedianMethod(SegLen, SegAmp, K);
		}
		outputStream<< boost::format("# BaseAmp=%1% \n")%BaseAmp;
		//outputStream<< boost::format("# Classify Segments \n", BaseAmp);

		SegState = (double *) calloc(K + 1, sizeof(double));
		for (i = 0; i <= K; i++)
			SegState[i] = SegAmp[i];
		CollapseAmpTtest(SegState, SegLen, K, BaseAmp, sigma2, T);
		if (debug){
			std::cerr<<" SelectClassifySegments done.\n";
		}
	}

	if (SelectClassifySegments == 0) {
		outputStream<< boost::format("Start\tStop\tLength\tAmpl\n");
		for (i = 0; i < K + 1; i++)
			outputStream<< Iext[i] + 1 << "\t" << Iext[i + 1] << "\t" << SegLen[i] << "\t"<< SegAmp[i] <<std::endl;
	} else if (SelectClassifySegments == 1) {
		outputStream<< boost::format("Start\tStop\tLenght\tAmpl\tState\n");
		for (i = 0; i < K + 1; i++) {
			outputStream<< Iext[i] + 1 << "\t" << Iext[i + 1] << "\t" << SegLen[i] << "\t"<< SegAmp[i];
			if (SegState[i] > BaseAmp)
				outputStream<< "G";
			else if (SegState[i] < BaseAmp)
				outputStream<< "L";
			else
				outputStream<< "N";
			outputStream<< endl;
		}

	}
	if (debug){
		std::cerr<< " output done." << endl;
	}
	closeFiles();
	if (debug){
		std::cerr<<"Exit GADA.commandlineRun()." << std::endl;
	}
}
#endif	// 2009-11-21 end of the whole boost python module code




#ifdef GADABIN
	int main(int argc, char* argv[]) {
		GADA instance(argc, argv);
		instance.commandlineRun();
	}
#endif
