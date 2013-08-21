/*
 * 2013.08.20 c++ version of SelectRowsFromMatrix.py, less function though.
 * selecting rows whose whichColumn value matches user input.
 */

#include "AbstractMatrixFileWalker.h"


class SelectRowsFromMatrix : public AbstractMatrixFileWalker {
protected:
	string whichColumnValue;


public:
	//must define non-default constructor, the parental non-default constructor won't be inherited by default.
	SelectRowsFromMatrix(int _argc, char** _argv);
	~SelectRowsFromMatrix();

	void constructOptionDescriptionStructure();
	int processRow(tokenizerCharType &line_toks);
};
