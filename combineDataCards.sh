#!/bin/bash

MASSCOMBINATIONS="$(< massCombinations.txt)" #names from names.txt file
for MASSCOMBINATION in $MASSCOMBINATIONS; do
	echo $MASSCOMBINATION
	combineCards.py Mass20To60=DataCards/T6bbllslepton_${MASSCOMBINATION}_20To60.txt Mass60To86=DataCards/T6bbllslepton_${MASSCOMBINATION}_60To86.txt Mass96To150=DataCards/T6bbllslepton_${MASSCOMBINATION}_96To150.txt Mass150To200=DataCards/T6bbllslepton_${MASSCOMBINATION}_150To200.txt Mass200To300=DataCards/T6bbllslepton_${MASSCOMBINATION}_200To300.txt Mass300To400=DataCards/T6bbllslepton_${MASSCOMBINATION}_300To400.txt MassAbove400=DataCards/T6bbllslepton_${MASSCOMBINATION}_Above400.txt  > combinedDataCards/T6bbllslepton_${MASSCOMBINATION}.txt
done
