#!/bin/bash 

#script to run ontologizer on mcl determined clusters

clustersFile="kroganMCLclusters_i2.txt"
count=1

while IFS=$'\t' read -r -a line
do
	for j in "${line[@]}"
	do
		echo $j >> "mcl2_cluster${count}.txt"
	done
	
	java -Xmx300M -jar Ontologizer.jar -a goa_human.gaf -g go.obo -s "mcl2_cluster${count}.txt" -p kroganProteins.txt -r 1000 -n -m "Benjamini-Hochberg"
	
	count=$(( ${count}+1 ))

done < "$clustersFile"

