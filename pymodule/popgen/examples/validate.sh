#!/bin/sh
# 2013.3.11 works like a unittest. this program runs some programs with example input and check whether they match the old output.
# good for testing after program is changed.

for exID in 1 2; do
	tmpOutFname=/tmp/sfs_code.new.ex$exID.polymorphism.h5
	echo -n "Running SFS_CODE_Output2PolymorphismTableFile.py on example $exID to temporary output file $tmpOutFname ..."
	../SFS_CODE_Output2PolymorphismTableFile.py -i ./sfs_code.ex$exID.output -u sfs_code.ex$exID.fasta -o $tmpOutFname
	echo " Done."
	
	# have to run ptdump. because binary files seem to differ when actual content is same.
	
	tmpOutTextFname=$tmpOutFname.txt
	echo -n "ptdump $tmpOutFname to $tmpOutTextFname ..."
	ptdump -a -v -d $tmpOutFname > $tmpOutTextFname
	echo " Done."
	
	standardOutputFname=sfs_code.ex$exID.polymorphism.h5
	standardOutTextFname=/tmp/$standardOutputFname.txt
	echo -n "ptdump $standardOutputFname to $standardOutTextFname ..."
	ptdump -a -v -d $standardOutputFname > $standardOutTextFname
	echo " Done."
	
	echo "Diff result:"
	diff $standardOutTextFname $tmpOutTextFname
	echo " (no difference if empty)."

done
