/*
 * 2013.09.10 a reducer that sums alignment depth of multiple alignments, filling missing bases with 0 depth.
 * Input are two VCF files.
 * Output is
 */

#include "AbstractMatrixFileWalkerCC.h"
#include "vcflib/src/Variant.h"
#include "MatrixFile.h"
using namespace vcf;

class CheckTwoVCFOverlapCC : public AbstractMatrixFileWalkerCC {
protected:
	std::string perSampleConcordanceOutputFname;
	MatrixFilePtr perSampleConcordanceOutputFile;

	std::string overlappingSitesOutputFname;
	MatrixFilePtr overlappingSitesOutputFile;

	std::string chromosome;
	long chrLength;
	set<Locus> overlapLocusSet;
	set<Locus> unionLocusSet;
	vector<string> overlapping_sample_id_list;


public:
	//must define non-default constructor, the parental non-default constructor won't be inherited by default.
	CheckTwoVCFOverlapCC(int _argc, char* _argv[]);
	~CheckTwoVCFOverlapCC();


	void constructOptionDescriptionStructure();
	void fileWalker(vector<string> &inputFnameList);
	void closeFiles();
	void checkSiteOverlap();
	void findOverlapSamples(VariantCallFile& vcfFile1, VariantCallFile& vcfFile2);
	void checkGenotypeConcordance();
};
