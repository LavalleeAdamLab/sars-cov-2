# SARS-CoV-2
Understanding the effects of the SARS-Cov-2 virus on cells is currently limited. Investigating protein-protein interactions between viral and host proteins can provide a better understanding of the mechanisms exploited by the virus and enable the identification of potential drug targets. Here we present our computational analysis of the SARS-CoV-2 and human protein interactome in HEK293 cells published by Gordon *et al.* to reveal processes that are affected by the virus and putative protein binding sites. 

Below instructions demonstrate the use of the various tools (i.e. GoNet, LESMoN-Pro, MCL, Ontologizer) with their input files to produce obtained results.

## GoNet: Gene Ontology enrichment analysis of String supplemented Krogan network
This tool enables statistical analysis of proteins annotated by Gene Ontology annotations in the provided network. For a given Gene Ontology Term, it computes a clustering measure for it's associated proteins and assesses their significance from a Monte Carlo derived normal distribution. A false discovery rate is approximated by comparing the set of Gene Ontology terms to a randomized set. Significantly clustered Gene Ontology annotations are identified at an FDR < 0.01. 

### Dependencies
* Java Version 8+
* The Apache Commons Mathematics Library (*commons-math3-3.6.1.jar*) 

### Required Input
* String supplemented Krogan Network (*krogan_prey_interaction_network.tsv*)
* Propagated Gene Ontology terms (*GO_annotations-9606-inferred-allev-2.tsv*)

### Running tool
1. File paths listed in the *Main* will need to be updated
2. Java scripts need to be compiled
