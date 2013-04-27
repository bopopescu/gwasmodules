#!/bin/bash

source ~/.bash_profile
deleteStageOutFolderAsWellDefault=0
relativePathToScratchDefault=scratch/
relativePathToPegasusFolderDefault=work/
if test $# -lt 1; then
	echo "Usage:"
	echo "  $0 relativePathToStageOutFolder [deleteStageOutFolderAsWell] [relativePathToScratch] [relativePathToPegasusFolder]"
	echo ""
	echo "Note:"
	echo "	#. This program deletes the workflow pegasus folder (job scripts, stdout/stderr, dag etc.= rm -rf \$relativePathToPegasusFolder/\$relativePathToStageOutFolder), its content in scratch folder (\$relativePathToScratch/\$relativePathToStageOutFolder) and the staged-out folder ( \$relativePathToStageOutFolder, if deleteStageOutFolderAsWell is non-zero)."
	echo "	#. You run this program while your current path is something like /Network/Data/vervet/workflow/. It assumes everything is in relative path."
	echo " 	#. Default relativePathToScratch is $relativePathToScratchDefault"
	echo " 	#. Default relativePathToPegasusFolder is $relativePathToPegasusFolderDefault"
	echo
	echo "Examples:"
	echo "	#"
	echo "	$0 BaseQualityRecalibration/LocalRealignmentBQSR_AlnID2828_2847_vsMethod87.2013.Apr.17T230821/ 1"
	echo
	echo
	exit 1
fi

relativePathToStageOutFolder=$1
deleteStageOutFolderAsWell=$2
relativePathToScratch=$3
relativePathToPegasusFolder=$4


if test -z "$deleteStageOutFolderAsWell"; then
	deleteStageOutFolderAsWell=$deleteStageOutFolderAsWellDefault
fi
if test -z "$relativePathToScratch"; then
	relativePathToScratch=$relativePathToScratchDefault
fi
if test -z "$relativePathToPegasusFolder"; then
	relativePathToPegasusFolder=$relativePathToPegasusFolderDefault
fi

deleteFolderFunc(){
	folderPath=$1
	if test -d "$folderPath"; then
		commandLine="rm -rf $folderPath"
		echo commandLine is $commandLine.
		$commandLine
	else
		echo $folderPath does not exist or not a folder
	fi
}
deleteFolderFunc $relativePathToPegasusFolder/$relativePathToStageOutFolder

deleteFolderFunc $relativePathToScratch/$relativePathToStageOutFolder

if test $deleteStageOutFolderAsWell != "0"; then
	echo "Deleting the stage out folder as well...";
	deleteFolderFunc $relativePathToStageOutFolder
fi
