# SARS-CoV-2
Understanding the effects of the SARS-Cov-2 virus on cells is currently limited. Investigating protein-protein interactions between viral and host proteins can provide a better understanding of the mechanisms exploited by the virus and enable the identification of potential drug targets. Here we present our computational analysis of the SARS-CoV-2 and human protein interactome in HEK293 cells published by Gordon *et al.* to reveal processes that are affected by the virus and putative protein binding sites. 

Below instructions demonstrate the use of the various tools (i.e. GoNet, LESMoN-Pro, MCL, Ontologizer) with their input files to produce obtained results.

## GoNet: Gene Ontology enrichment analysis of String supplemented Krogan network
This tool enables statistical analysis of proteins annotated by Gene Ontology annotations in the provided network. For a given Gene Ontology Term, it computes a clustering measure for it's associated proteins and assesses their significance from a Monte Carlo derived normal distribution. A false discovery rate is approximated by comparing the set of Gene Ontology terms to a randomized set. Significantly clustered Gene Ontology annotations are identified at an FDR < 0.01. 

### Dependencies
* Java Version 8+
* The Apache Commons Mathematics Library (*commons-math3-3.6.1.jar*) 

### Files
#### Required input
Found under: GoNet > files > inout_files
* String supplemented Krogan Network (*krogan_prey_interaction_network.tsv*)
* Propagated Gene Ontology terms (*GO_annotations-9606-inferred-allev-2.tsv*)

#### Generated intermediate files
Found under: GoNet > files > IO_files
* initial distance matrix (*krogan_DistanceMatrix.txt*)
* final distance matrix of fully connected component (*krogan_DistanceMatrix_FullConnected.txt*)
* Monte Carlo distributions (*kroganDistribution_s10X5_n3_186.txt*)
* shuffle Gene Ontology set (*shuffled-go.txt*)

#### Output files
Found under : GoNet > files > output_files
* list of calculated false discovery rates at significant thresholds (*2020.04.08.MonoTransf_kroganFDR_10X5.txt*)
* list of GO terms and their calculated p-values (*2020.04.08.MonoTransf_kroganGoTerms_s10X5.txt*)

### Running GoNet
This tool was developped using java. All required scripts are found under GoNet > java. File paths must be modified under *Main* and scripts must be compiled prior to running. GoNet requires command line arguments, the array of proteins to sample for Monte Carlo Sampling and the number of times to sample the network during Monte Carlo Sampling.

# LESMoN-Pro: Sequence Motif enrichment analysis of String supplemented Krogan Network 
This tool enables statistical analysis of proteins with shared sequence motifs in the provided network. For a given motif, it computes a clustering measure for it's associated proteins and assesses their significance from a Monte Carlo derived normal distribution. A false discovery rate is approximated by comparing the motif set to randomized sequence motifs. Motifs who's proteins are significantly clustered are identified at an FDR < 0.01. 

### Dependencies
* Java Version 8+
* The Apache Commons Mathematics Library (*commons-math3-3.6.1.jar*)

* Python Version **XX 
* Packages **XX

### Files
#### Required input
Found under: LESMoN-Pro > files > inout_files
* String supplemented Krogan Network (*krogan_prey_interaction_network.tsv*)
* List of motifs and associated proteins (**)
* List of randomized motifs and associated proteins (**)

#### Generated intermediate files
Found under: LESMoN-Pro > files > IO_files
* initial distance matrix (*krogan_DistanceMatrix.txt*)
* final distance matrix of fully connected component (*krogan_DistanceMatrix_FullConnected.txt*)
* Monte Carlo distributions (**)

#### Output files
Found under : GoNet > files > output_files
* list of calculated false discovery rates at significant thresholds (**)
* list of motifs and their calculated p-values (**)

### Running LESMoN-Pro

