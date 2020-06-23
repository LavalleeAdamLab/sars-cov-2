#!/bin/bash 

#script to run ontologizer on mcl determined clusters

clustersFile="Krogan_bait_interactions_UniprotAC.txt"
count=1

while IFS=$'\t' read -r -a line
do
	for j in "${line[@]}"
	do
		echo $j >> "starsCluster${count}.txt"
	done
	
	java -Xmx300M -jar Ontologizer.jar -a goa_human.gaf -g go.obo -s "starsCluster${count}.txt" -p kroganProteins.txt -r 1000 -n -m "Benjamini-Hochberg"
	
	count=$(( ${count}+1 ))

done < "$clustersFile"

