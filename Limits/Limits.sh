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
	#~ DATACARD="../combinedDataCardsFlatSystematics/${MODEL}.txt"
	DATACARD="../combinedDataCardsFlatSystematics/${MODEL}.txt"
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
	echo "$EXE -M Asymptotic $DATACARD > $LOG_ASYMPTOTIC"
	$EXE -M Asymptotic $DATACARD > $LOG_ASYMPTOTIC
	OBSas=`grep "Observed Limit: r <" $LOG_ASYMPTOTIC | cut -b 21-`
	EXPas=`grep "Expected 50.0%: r <" $LOG_ASYMPTOTIC | cut -b 21-`
	EXPm2as=`grep "Expected  2.5%: r <" $LOG_ASYMPTOTIC | cut -b 21-`
	EXPm1as=`grep "Expected 16.0%: r <" $LOG_ASYMPTOTIC | cut -b 21-`
	EXPp1as=`grep "Expected 84.0%: r <" $LOG_ASYMPTOTIC | cut -b 21-`
	EXPp2as=`grep "Expected 97.5%: r <" $LOG_ASYMPTOTIC | cut -b 21-`
	
	echo "CLs observed asymptotic = $OBSas" >> $RESULT
	echo "CLs expected asymptotic = $EXPas" >> $RESULT
	echo "CLs expected m2sigma asymptotic = $EXPm2as" >> $RESULT
	echo "CLs expected m1sigma asymptotic = $EXPm1as" >> $RESULT
	echo "CLs expected p1sigma asymptotic = $EXPp1as" >> $RESULT
	echo "CLs expected p2sigma asymptotic = $EXPp2as" >> $RESULT

	cat $RESULT

	rm -f *.root 

done
