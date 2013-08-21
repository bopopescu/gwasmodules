/*
 * 2013.08.20 c++ version of SelectRowsFromMatrix.py, less function though.
 * selecting rows whose whichColumn value matches user input.
 */

#include "AbstractMatrixFileWalkerCC.h"


class SelectRowsFromMatrixCC : public AbstractMatrixFileWalkerCC {
protected:
	string whichColumnValue;


public:
	//must define non-default constructor, the parental non-default constructor won't be inherited by default.
	SelectRowsFromMatrixCC(int _argc, char* _argv[]);
	~SelectRowsFromMatrixCC();

	void constructOptionDescriptionStructure();
	int processRow(tokenizerCharType &line_toks);
};
