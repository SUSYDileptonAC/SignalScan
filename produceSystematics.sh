#!/bin/bash

MASSES="$(< masses.txt)" #names from names.txt file 

for MASS in $MASSES; do
	echo $MASS
	python produceSystematics.py $MASS 
done
