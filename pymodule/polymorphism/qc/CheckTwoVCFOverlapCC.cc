/*
 * 2013.09.10 a c++ version of CheckTwoVCFOverlap.py
 *
 */
#include "pymodule/include/CheckTwoVCFOverlapCC.h"

CheckTwoVCFOverlapCC::CheckTwoVCFOverlapCC(int _argc, char* _argv[]): AbstractMatrixFileWalkerCC(_argc, _argv) {

	//overwrite these doc variables
	usageDoc = boost::format("%S -i INPUT1.vcf -i INPUT2.vcf -o OUTPUTFNAME [OPTIONS]\n")% programName;
	examplesDoc = boost::format("%S -i ~/script/vcflib/samples/sample.vcf -i ~/script/vcflib/samples/sample.vcf"
			" -o /tmp/vcfOverlap.tsv.gz --overlappingSitesOutputFname /tmp/overlapSites.tsv.gz --perSampleConcordanceOutputFname /tmp/sampleConcordance.tsv.gz \n")% programName;
}

CheckTwoVCFOverlapCC::~CheckTwoVCFOverlapCC(){

}
void CheckTwoVCFOverlapCC::constructOptionDescriptionStructure(){
	AbstractMatrixFileWalkerCC::constructOptionDescriptionStructure();
	if (debug){
		cerr<< "adding more options within CheckTwoVCFOverlapCC::_constructOptionDescriptionStructure() ...";
	}
	optionDescription.add_options()
		("chromosome", po::value<string>(&chromosome)->default_value("Chr0"), "which chromosome the input VCF is from."
				"Used only to mark output, not used for filtering ")
		("chrLength", po::value<long>(&chrLength)->default_value(100000), "length of chromosome in the input VCF."
				"Used only to mark output, not used for filtering ")
		("perSampleConcordanceOutputFname", po::value<string>(&perSampleConcordanceOutputFname), "output filename for the per-sample matching result."
				"Per-sample concordance will not be run if this is not given.")
		("overlappingSitesOutputFname", po::value<string>(&overlappingSitesOutputFname), "output filename that would contain locations of all overlap sites between two VCF. no output if not given.");
	if (debug){
		cerr<< endl;
	}
}


void CheckTwoVCFOverlapCC::checkSiteOverlap(){
	if (debug){
		std::cerr << boost::format("CheckTwoVCFOverlapCC::checkSiteOverlap() %1% and %2% ...")%inputFnameList[0] % inputFnameList[1];
	}
	VariantCallFile vcfFile1;
	VariantCallFile vcfFile2;

	if (inputFnameList.size()>=2){
		vcfFile1.open(inputFnameList[0]);
		vcfFile2.open(inputFnameList[1]);
	}
	else{
		std::cerr << boost::format("Only %1% files <2, not enough to do comparison.\n")%inputFnameList.size();
		return;
	}
	if (!vcfFile1.is_open() || !vcfFile2.is_open()) {
		std::cerr << boost::format("Either %1% or %2% is not opened properly.\n")%inputFnameList[0] % inputFnameList[1];
		return;
	}

	/*
	Variant var(vcfFile1);
	while (vcfFile1.getNextVariant(var)) {
		cerr<<var.sequenceName << ":" << var.position << endl;
	}
	*/

	vcfFile1.readInLocusList();
	vcfFile2.readInLocusList();
	int noOfOverlapLoci = 0;
	vector<Locus>::iterator locusListIter;
	for (locusListIter=vcfFile1.locusList.begin();locusListIter!=vcfFile1.locusList.end(); locusListIter++){
		if (vcfFile2.locus2index.left.find(*locusListIter)!=vcfFile2.locus2index.left.end()){
			overlapLocusSet.insert(*locusListIter);
		}
		unionLocusSet.insert(*locusListIter);
		//cerr << *locusListIter << ":" << vcfFile2.locus2index.left.at(*locusListIter) << endl;
	}
	//add loci from 2nd file into union set
	for (locusListIter=vcfFile2.locusList.begin();locusListIter!=vcfFile2.locusList.end(); locusListIter++){
		unionLocusSet.insert(*locusListIter);
	}

	if (debug){
		cerr << boost::format(" %1% overlap loci and %2% union loci.\n")% overlapLocusSet.size() % unionLocusSet.size();
	}
	// output overlapping site summary statistics
	if (outputFile->is_open()){
		vector<string> summaryOutputHeader {"#chromosome", "length", "#sitesInInput1", "#sitesInInput2", "#overlapping", "overlappingOverTotal", \
						"overlappingOverInput1", "overlappingOverInput2", "#segregatingSitesNormalized"};
		outputFile->writeRow(summaryOutputHeader);
		double no_of_overlapping_sites = double(overlapLocusSet.size());
		double overlapping_fraction = no_of_overlapping_sites/unionLocusSet.size();
		double overlappingOverInput1 = no_of_overlapping_sites/vcfFile1.locusList.size();
		double overlappingOverInput2 = no_of_overlapping_sites/vcfFile2.locusList.size();
		float normalizingConstant;
		int noOfSamples = vcfFile1.sampleNames.size();
		if (noOfSamples>1)
			normalizingConstant = float(sumOfReciprocals(noOfSamples*2-1));
		else
			normalizingConstant = 1;
		double noOfSegregatesSitesNormalized = no_of_overlapping_sites/(normalizingConstant*chrLength);
		vector< anyType > dataRow {anyType(chromosome), anyType(chrLength),
			anyType(vcfFile1.locusList.size()), anyType(vcfFile2.locusList.size()),
			anyType(no_of_overlapping_sites),
			anyType(overlapping_fraction), anyType(overlappingOverInput1), anyType(overlappingOverInput2),
			anyType(noOfSegregatesSitesNormalized)};
		outputFile->writeRow(dataRow);
		//outputFile->outputStream << chromosome ;
	}
	if (debug){
		cerr<< endl;
	}
}
void CheckTwoVCFOverlapCC::findOverlapSamples(VariantCallFile& vcfFile1, VariantCallFile& vcfFile2){
	/*
	 * 2013.09.17
	 */
	if (debug){
			cerr<<boost::format("findOverlapSamples(): %2% samples in vcfFile1, %3% samples in vcfFile2, start with %1% overlap samples ") % overlapping_sample_id_list.size() %
					vcfFile1.sampleName2index.size() % vcfFile2.sampleName2index.size() ;
		}

	vector<string>::iterator sampleNameVectorIter;
	for (sampleNameVectorIter=vcfFile1.sampleNames.begin();sampleNameVectorIter!=vcfFile1.sampleNames.end(); sampleNameVectorIter++){
		if (vcfFile2.sampleName2index.left.find(*sampleNameVectorIter)!=vcfFile2.sampleName2index.left.end()){
			overlapping_sample_id_list.push_back(*sampleNameVectorIter);
		}
	}
	/*
	set_intersection(vcfFile1.sampleNames.begin(), vcfFile1.sampleNames.end(),
				vcfFile2.sampleNames.begin(), vcfFile2.sampleNames.end(),
				std::inserter(overlapping_sample_id_list, overlapping_sample_id_list.begin()) );
	*/
	sort(overlapping_sample_id_list.begin(), overlapping_sample_id_list.end());
	if (debug){
		cerr<<boost::format(" %1% overlap samples ") % overlapping_sample_id_list.size() ;
	}
}

void CheckTwoVCFOverlapCC::checkGenotypeConcordance() {
	if (debug){
		std::cerr << boost::format("CheckTwoVCFOverlapCC::checkGenotypeConcordance() %1% and %2% ...")%inputFnameList[0] % inputFnameList[1];
	}
	// do per-sample concordance
	// close and re-open the variant files first , in order to reset the file handler
	// construct locus-id2index and sample-id2index maps for both files
	// read-in whole data matrix

	VariantCallFile vcfFile1;
	VariantCallFile vcfFile2;
	vcfFile1.open(inputFnameList[0]);
	vcfFile2.open(inputFnameList[1]);

	vcfFile1.readInDataMatrix();
	vcfFile2.readInDataMatrix();
	//check concordance


	findOverlapSamples(vcfFile1, vcfFile2);


	int no_of_samples_to_compare = overlapping_sample_id_list.size();
	int* no_of_non_NA_pairs_per_sample_ls = new int [no_of_samples_to_compare]();
	int* no_of_matches_per_sample_ls = new int[no_of_samples_to_compare]();	// add () to initialize everything to 0
	double* matchFractionLs = new double [no_of_samples_to_compare]();// add () to initialize everything to 0

	set<Locus>::iterator locusIter;
	string sample_id;
	int row_index1, row_index2, col_index1, col_index2, call1, call2;
	for (int j=0; j<no_of_samples_to_compare; j++){
		sample_id = overlapping_sample_id_list[j];
		//cerr << endl << sample_id;
		col_index1 = vcfFile1.sampleName2index.left.find(sample_id)->second;
		col_index2 = vcfFile2.sampleName2index.left.find(sample_id)->second;	//.at() does not work because it's multiset
		no_of_non_NA_pairs_per_sample_ls[j]=0;
		no_of_matches_per_sample_ls[j] = 0;
		for (locusIter=overlapLocusSet.begin(); locusIter!=overlapLocusSet.end(); locusIter++){
			row_index1 = vcfFile1.locus2index.left.at(*locusIter);
			row_index2 = vcfFile2.locus2index.left.at(*locusIter);
			call1 = vcfFile1.dataMatrix[row_index1][col_index1];
			call2 = vcfFile2.dataMatrix[row_index2][col_index2];
			//cerr << boost::format("\t %1% \t %2% \t %3%")% locusIter->repr % call1 % call2 ;
			if ((call1 >=2) && (call2 >=2)){	//exclude missing data
				no_of_non_NA_pairs_per_sample_ls[j] += 1;
				if (call1==call2){
					no_of_matches_per_sample_ls[j] += 1;
				}
			}
		}
		//cerr << endl;
	}
	for (int j=0; j<no_of_samples_to_compare; j++){
		if (no_of_non_NA_pairs_per_sample_ls[j]>0){
			matchFractionLs[j] = no_of_matches_per_sample_ls[j]/double(no_of_non_NA_pairs_per_sample_ls[j]);
		}
		else{
			matchFractionLs[j]=-1;
		}
	}
	//output concordance result
	perSampleConcordanceOutputFile = MatrixFilePtr(
			new MatrixFile<int>(perSampleConcordanceOutputFname, "w", debug,
						"\t"));
	vector<string> header { "sample_id", "no_of_matches", "no_of_non_NA_pairs",
			"matchFraction" };
	perSampleConcordanceOutputFile->writeRow(header);
	for (int j=0; j<no_of_samples_to_compare; j++){
		vector<anyType> data_row {anyType(overlapping_sample_id_list[j]), anyType(no_of_matches_per_sample_ls[j]),
			anyType(no_of_non_NA_pairs_per_sample_ls[j]), anyType(matchFractionLs[j])};
		perSampleConcordanceOutputFile->writeRow(data_row);
	}

	overlapping_sample_id_list.clear();
	delete no_of_matches_per_sample_ls;
	delete no_of_non_NA_pairs_per_sample_ls;
	delete matchFractionLs;
	if (debug) {
		cerr << endl;
	}
}

void CheckTwoVCFOverlapCC::fileWalker(vector<string> &inputFnameList){
	/*
	 * 2013.09.10
	 */
	if (debug){
		std::cerr << boost::format("CheckTwoVCFOverlapCC::fileWalker() Walking through %1% files ...")%inputFnameList.size();
	}

	//outputFname is opened elsewhere
	//openOneOutputFile(outputFname, outputFilterStreamBuffer, outputFile);


	checkSiteOverlap();

	//output detailed info on overlapping sites
	if (!overlappingSitesOutputFname.empty()){
		overlappingSitesOutputFile = MatrixFilePtr(new MatrixFile<int>(overlappingSitesOutputFname, "w", debug, "\t"));
		vector<string> outputHeader = { "chromosome", "position", "random"};
		overlappingSitesOutputFile->writeRow(outputHeader);
		for (set<Locus>::iterator i=overlapLocusSet.begin(); i!=overlapLocusSet.end(); i++){
			vector< boost::variant<std::string, long> > dataRow {i->sequenceName, i->position, 0};
			//boost::any causes trouble in writeRow and couldn't be << to outputStream
			overlappingSitesOutputFile->writeRow(dataRow);
		}


	}

	if (!perSampleConcordanceOutputFname.empty()){
		checkGenotypeConcordance();
	}

	if (debug){
		std::cerr << boost::format(" CheckTwoVCFOverlapCC::fileWalker() done.\n");
	}

}

void CheckTwoVCFOverlapCC::closeFiles(){
	//
	AbstractMatrixFileWalkerCC::closeFiles();

	if (!perSampleConcordanceOutputFname.empty()){
		perSampleConcordanceOutputFile->close();
		//outputFile.close();	//2013.08.21 if output is a .gz, closing it here would result a truncated .gz file. don't know why.
		//maybe because the buffer or the gzip compressor filter needs to be closed beforehand.
	}
	if (!overlappingSitesOutputFname.empty()){
		overlappingSitesOutputFile->close();
	}
}

#ifndef __MAIN__
	int main(int argc, char* argv[]) {
		CheckTwoVCFOverlapCC instance(argc, argv);
		instance.run();
	}
#define __MAIN__
#endif
