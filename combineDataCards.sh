#!/bin/bash

MASSCOMBINATIONS="$(< massCombinations.txt)" #names from names.txt file
for MASSCOMBINATION in $MASSCOMBINATIONS; do
	echo $MASSCOMBINATION
	combineCards.py lowMass=DataCards/T6bbllslepton_${MASSCOMBINATION}_lowMass.txt belowZ=DataCards/T6bbllslepton_${MASSCOMBINATION}_belowZ.txt onZ=DataCards/T6bbllslepton_${MASSCOMBINATION}_onZ.txt aboveZ=DataCards/T6bbllslepton_${MASSCOMBINATION}_aboveZ.txt highMass=DataCards/T6bbllslepton_${MASSCOMBINATION}_highMass.txt > combinedDataCards/T6bbllslepton_${MASSCOMBINATION}.txt
	#~ combineCards.py lowMass=DataCardsFlatSystematics/T6bbllslepton_${MASSCOMBINATION}_lowMass.txt belowZ=DataCardsFlatSystematics/T6bbllslepton_${MASSCOMBINATION}_belowZ.txt onZ=DataCardsFlatSystematics/T6bbllslepton_${MASSCOMBINATION}_onZ.txt aboveZ=DataCardsFlatSystematics/T6bbllslepton_${MASSCOMBINATION}_aboveZ.txt highMass=DataCardsFlatSystematics/T6bbllslepton_${MASSCOMBINATION}_highMass.txt > combinedDataCardsFlatSystematics/T6bbllslepton_${MASSCOMBINATION}.txt
done
