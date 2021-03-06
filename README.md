# SARS-CoV-2
Understanding the effects of the SARS-Cov-2 virus on cells is currently limited. Investigating protein-protein interactions between viral and host proteins can provide a better understanding of the mechanisms exploited by the virus and enable the identification of potential drug targets. Here we present our computational analysis of the SARS-CoV-2 and human protein interactome in HEK293 cells published by Gordon *et al.* to reveal processes that are affected by the virus and putative protein binding sites. 

Below instructions demonstrate the use of the various tools (i.e. GoNet, LESMoN-Pro, MCL, Ontologizer) with their input files to produce obtained results.

## GoNet: Gene Ontology enrichment analysis of String supplemented Krogan network
This tool enables statistical analysis of proteins annotated by Gene Ontology annotations in the provided network. For a given Gene Ontology Term, it computes a clustering measure for it's associated proteins and assesses their significance from a Monte Carlo derived normal distribution. A false discovery rate is approximated by comparing the set of Gene Ontology terms to a randomized set. Significantly clustered Gene Ontology annotations are identified at an FDR < 0.01. 

### Dependencies
* Java Version 8+
* Required library: The Apache Commons Mathematics Library (*commons-math3-3.6.1.jar*) 

### Files
#### Required input
Found under: `GoNet > files > input_files`
* String supplemented Krogan Network (*krogan_prey_interaction_network.tsv*)
* Propagated Gene Ontology terms (*GO_annotations-9606-inferred-allev-2.tsv*)

#### Generated intermediate files
Found under: `GoNet > files > IO_files`
* initial distance matrix (*krogan_DistanceMatrix.txt*)
* final distance matrix of fully connected component (*krogan_DistanceMatrix_FullConnected.txt*)
* Monte Carlo distributions (*kroganDistribution_s10X5_n3_186.txt*)
* shuffle Gene Ontology set (*shuffled-go.txt*)

#### Output files
Found under : `GoNet > files > output_files`
* list of calculated false discovery rates at significant thresholds (*2020.04.08.MonoTransf_kroganFDR_10X5.txt*)
* list of GO terms and their calculated p-values (*2020.04.08.MonoTransf_kroganGoTerms_s10X5.txt*)

### Running GoNet
This tool was developped using java. All required scripts are found under GoNet > java. File paths must be modified under *Main* and scripts must be compiled prior to running. GoNet requires command line arguments, the array of proteins to sample for Monte Carlo Sampling and the number of times to sample the network during Monte Carlo Sampling.

## LESMoN-Pro: Sequence Motif enrichment analysis of String supplemented Krogan Network 
This tool enables statistical analysis of proteins with shared sequence motifs in the provided network. For a given motif, it computes a clustering measure for it's associated proteins and assesses their significance from a Monte Carlo derived normal distribution. A false discovery rate is approximated by comparing the motif set to randomized sequence motifs. Motifs who's proteins are significantly clustered are identified at an FDR < 0.01. 

*Prior to running LESMoN pro, sequence randomization and sequence enumeration must be performed.*

### Dependencies
* Java Version 8+
* Required library: The Apache Commons Mathematics Library (*commons-math3-3.6.1.jar*)

* Python Version 3+ 
* Required package: Panda

### Files
#### Required input
Found under: `LESMoN-Pro > files > inout_files`
* String supplemented Krogan Network (*krogan_prey_interaction_network.tsv*)
* Example list of sequence motifs to test and associated proteins (*KroganMotifsAnnotations_real_exampleInput.tsv*)
* Examples list of shuffled sequence motifs to test and associated proteins (*KroganMotifAnnotations_shuffled_exampleInput.tsv*)

#### Generated intermediate files
Found under: `LESMoN-Pro > files > IO_files`
* initial distance matrix (*krogan_DistanceMatrix.txt*)
* final distance matrix of fully connected component (*krogan_DistanceMatrix_FullConnected.txt*)
* list of network proteins and their associated number of real sequence motifs (*annotatedProteinsToOccurence_real.txt*)
* list of network proteins and their associated number of shuffled sequence motifs (*annotatedProteinsToOccurence_shuffled.txt*) 
* Monte Carlo distributions corresponding to real sequence motifs (*kroganDistribution_REAL_s10X5_n3_195.txt*)
* Monte Carlo distributions corresponding to shuffled sequence motifs (*kroganDistribution_SHUFFLED_s10X5_n3_195.txt*)
* Normal distribution parameters derived from Monte Carlo distribution for real sequence motifs (*krogan_REAL_normalDistributionParams.txt*) 
* Normal distribution parameters derived from Monte Carlo distribution for shuffled sequence motifs (*krogan_SHUFFLED_normalDistributionParams*) 

#### Output files
Found under : `LESMoN-Pro > files > output_files`
* example list of real sequence motifs and their calculated p-values (*KroganMotifs_real_exampleOutput.txt*)
* example list of shuffled sequence motifs and their calculated p-values (*KroganMotifs_shuffled_exampleOutput*)
* list of calculated false discovery rates at significant thresholds (*covFDRs2.tsv*)

### Running LESMoN-Pro
This tool was developped using java and python. All required scripts are found under `LESMoN-Pro > java` and `LESMON-Pro > python`. File paths must be modified under *Main* and *FDRcalc* and java scripts must be compiled prior to running. The java implementation must be run prior to python FDRcalc and requires command line arguments, the array of proteins to sample for Monte Carlo Sampling and the number of times to sample the network during Monte Carlo Sampling.

## Shuffling sequence
Protein sequences were randomized using non overlapping sliding windows of size 10. This tool was implemented in python. It requires the input of a Fasta file containing all protein sequences in the String supplement Krogan Network (*Krogan_Protein_database.fasta*) and returns a Fasta file of randomized sequences (*Krogan_Protein_database_SHUFFLED.fasta*). Files are found under `Shuffling Fasta Sequences` folder. 

## Markov Cluster (MCL) Algorithm of the String supplemented Krogan network
The MCL algorithm is a tool developped by Enright A.J., *et al.* to cluster graphs. Here, it is used to identify innate clusters in the String supplemented Krogan network prior to functional enrichment analysis using Ontologizer and sequence motif enrichment analysis using MEME. Files are found under `MCL Algorithm` folder. 

### Dependencies
The MCL algorithm is a commandline tool that was downloaded from [https://micans.org/mcl/src/]. 

### Running MCL
The MCL algorithm was provided the graph input file (*nodes.txt*), with parameters `--abc` for specified input format, `-I` for an inflation variable of 2, and `-o` specified output file name (*kroganMCLclusters_i2.txt*). Commandline output was directed to file (*kroganOutput_i2.txt*). The bash script to execute MCL algorithm is provided (*runMCL.sh*).

## Ontologizer: Functional enrichment of MCL algorithm clusters.
Ontologizer is a command line tool developped by Bauer, S., *et al* to perform Gene Ontology enrichment analysis. Here it is used to assess enriched Gene Ontologies in the protein clusters generated by the MCL algorithm of the String supplemented Krogan network.

### Ontologizer File Requirements
Found under: `Ontologizer-MCL > java`
* Ontologizer application as jar executable (*Ontologizer.jar*)
* Gene Ontologies (*go.obo*)
* Gene association file (*goa_human.gaf*)
* Protein background file (*kroganProteins.txt*)
* Ontologizer input example MCL cluster file (*mcl_i2_cluster1.txt*)
* Ontologizer output table file (*table-mcl2_cluster1-Parent-Child-Union-Benjamini-Hochberg.txt*) 
* bash script to run Ontologizer (*runOntologizer.sh*) 

### BH correction File Requirements
Found under: `Ontologizer-MCL > R`
* Combined all go terms input file (*2020.04.06.kroganMCLi2_allGOs.txt*)
* GO terms with corrected p-values output file (*2020.04.08.kroganMCLi2_allSignificantGOs.txt*)
* R script to perform BH correction (*ontologizer_pvalCorrection.R*) 

### Running Ontologizer
Ontologizer was run on individual MCL clusters. Obtained results were combined and Benjamini-Hochberg correction was performed in R on Ontologizer p-value. Smallest corrected p-value associated to a given GO-term was kept.

## Ontologizer: Functional enrichment of Krogan paper star clusters
Ontologizer is a command line tool developped by Bauer, S., *et al* to perform Gene Ontology enrichment analysis. Here it is used to assess enriched Gene Ontologies in the protein clusters generated by the MCL algorithm of the String supplemented Krogan network.

### Ontologizer File Requirements
Found under: `Ontologizer-Star`
* Ontologizer application as jar executable (*Ontologizer.jar*)
* Gene Ontologies (*go.obo*)
* Gene association file (*goa_human.gaf*)
* Protein background file (*kroganProteins.txt*)
* Ontologizer input example Star cluster file (*starCluster3.txt*)
* Ontologizer output table file (*table-starsCluster3-Parent-Child-Union-Benjamini-Hochberg.txt*) 
* bash script to run Ontologizer (*runOntologizer.sh*) 

### Running Ontologizer
Ontologizer was run on individual star clusters. 
