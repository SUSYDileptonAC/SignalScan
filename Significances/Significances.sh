#!/bin/bash

EXE=combine


COMBINATIONS="$(< ../massXSectionCombinations.txt)" #names from names.txt file
for COMBINATION in $COMBINATIONS; do
	SBOTTOM=$(echo $COMBINATION| cut -d'_' -f 1)
	NEUTRALINO=$(echo $COMBINATION| cut -d'_' -f 2)
	XSECTION=$(echo $COMBINATION| cut -d'_' -f 3)
	echo $SBOTTOM
	echo $NEUTRALINO
	MODEL="T6bbllslepton_${SBOTTOM}_${NEUTRALINO}"
	DATACARD="../DataCards/${MODEL}.txt"
	RESULT="$MODEL.result.txt"
	LOG_ASYMPTOTIC='asymptotic.log'

	rm -f $RESULT
	touch $RESULT

	pwd
	echo "# $MODEL"
	echo ""
	echo "sbottom = $SBOTTOM"
	echo "neutralino = $NEUTRALINO"
	echo "Xsection = $XSECTION"
	
	echo "# $MODEL" >> $RESULT
	echo "" >> $RESULT
	echo " sbottom = $SBOTTOM" >> $RESULT
	echo " neutralino 2 = $NEUTRALINO" >> $RESULT
	echo " Xsection = $XSECTION" >> $RESULT


	#Calculate asymptotic CLs x-section limits
	echo "$EXE -M ProfileLikelihood --uncapped 1 --significance $DATACARD > $LOG_ASYMPTOTIC"
	$EXE -M ProfileLikelihood --uncapped 1 --significance $DATACARD > $LOG_ASYMPTOTIC
	Signif=`grep "Significance:" $LOG_ASYMPTOTIC | cut -b 15-21`
	PValue=`grep "p-value =" $LOG_ASYMPTOTIC | cut -b 19-27`

	echo "" >> $RESULT
	echo "CLs observed significance = $Signif" >> $RESULT
	echo "CLs observed p-value = $PValue" >> $RESULT

	cat $RESULT

	rm -f *.root 

done
