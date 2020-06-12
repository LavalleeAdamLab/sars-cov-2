import java.util.ArrayList;

import graph.Annotation;
import graph.Interaction;
import graph.Protein;
import utils.*;


public class Main {

	/*** Constants method variables ***/
	public static final boolean toComputeDistanceMatrix = false;
	public static final boolean toUpdateDistanceMatrix = false;
	
	public static final boolean toMonteCarlo = false;
		public static final boolean weightedSampling = true;
		public static final boolean testNormalApprox  = true;
	
	public static final boolean testGoTerms = false;
		public static final boolean printFDRresults = true;
		public static final boolean printSignificantGOterms = true;

	/*** Constant running variables ***/
	public static final boolean runComputeCanada = false;

	/**** Constant output variables ***/
	public static final boolean printNumberOfInteractionsPerProtein = false;
	public static final boolean printPPInetwork = false;
	
	public static final boolean toShuffle = true;
	
	/**** Constant  distribution variables ****/ 
	public static final boolean computeDistributionParams = false; //if false don't read MC file and re-calculate parameters
	public static final boolean useNormalDistribution = true; //if false uses monte carlo values
	
	public static void main(String[] args) {
		
		/*****************************************************************************************************
		 * This algorithm aims to identify differentially expressed pathways in two conditions of interest,
		 * here breast cancer subtypes (Her2+ and TN), using a Protein protein interaction network guided-
		 * functional enrichement analysis.
		 * 
		 * Input
		 *  - BioGRID (Homo sapiens). version 3.4.158 release March 1st 2018.
		 *  - Protein expression file (Tyanova et al. study)
		 *  - GO annotations (Homo sapiens). Gene Ontology Consortium obtained June 2018.
		 *****************************************************************************************************/

		int nProtToSampleLowerBound = Integer.parseInt(args[0]);		// Lower limit of number of proteins sampled
		int nProtToSampleUpperBound = Integer.parseInt(args[1]);		// Upper limit of number of proteins sampled
		int numOfTimesNetworkIsSampled = Integer.parseInt(args[2]);		// Number of times the network is sampled during Monte Carlo Method

		String working_directory = "C:/Users/Rachel/Documents/COVID/NetworkMotif-Krogan/";
		
		String PPI_inputFile = "krogan_prey_interaction_network.tsv.txt";
		String annotationInputFile = "RNtest.tsv";
		String annotatedProteinOccurenceFile = "";
		String normalDistributionParametersFile = "";
		String distanceMatrixFile = "unweighted_DistanceMatrix.txt";
		String distanceMatrix2File = "unweighted_DistanceMatrix_FullConnected.txt";

		String distributionFile = "unweightedNetwork3.0_s10^7_n1_1000.txt";
		if(toShuffle) {
			annotationInputFile = "shuffledMotifs";
			distributionFile = "";
		}

		if (!runComputeCanada) {
			PPI_inputFile = working_directory + "input_files/krogan_prey_interaction_network.tsv.txt";
			distanceMatrixFile = working_directory + "IO_files/krogan_DistanceMatrix.txt";
			distanceMatrix2File = working_directory + "IO_files/krogan_DistanceMatrix_FullConnected.txt";
			
			annotationInputFile = working_directory + "input_files/Krogan_motif_annotations_real_FILTERED.tsv";
			annotatedProteinOccurenceFile = working_directory + "IO_files/annotatedProteinsToOccurence_real.txt";
			distributionFile = working_directory + "IO_files/kroganDistribution_REAL_s10X5_n3_195.txt";
			normalDistributionParametersFile = working_directory + "IO_files/krogan_REAL_normalDistributionParams.txt";
			
		
			if(toShuffle) {
				annotationInputFile = working_directory + "Krogan_motif_annotations_shuffled_FILTERED.tsv";
				annotatedProteinOccurenceFile = working_directory + "IO_files/annotatedProteinsToOccurence_shuffled.txt";
				distributionFile = working_directory + "IO_files/kroganDistribution_SHUFFLED_s10X5_n3_195.txt";
				normalDistributionParametersFile = working_directory + "IO_files/krogan_SHUFFLED_normalDistributionParams.txt";
			}
		}

		/* Output files */
		//String proteinListFile = working_directory + "output_files/2020.04.07.ListOfKroganProteinsInConnectedComponent.txt";
		String motifSignificanceFile = working_directory + "output_files/2020.05.27.MonoTransf_kroganFDR_representativeTEST_s10X5.txt";

		/* Read BioGRID file. Extract all possible interactions between proteins and store in Interaction object */
		ArrayList<Interaction> networkInteractionsList = NetworkInteractionsLoader.importInteractionNetwork(PPI_inputFile, runComputeCanada);
//		System.out.println("Number of initial interactions: " + networkInteractionsList.size());
		
		/* Extract all proteins in the network (ie. all proteins within the interaction file) */
		ArrayList<Protein> networkProteinList = NetworkProteins.getProteinsInNetwork(networkInteractionsList);
//		System.out.println("Number of initial proteins: " + networkProteinList.size());

		if(toComputeDistanceMatrix) {
			/* Compute distance matrix using the list of proteins in the network and the list of known interactions in the network. */
			DistanceMatrix.computeDistanceMatrix(networkInteractionsList, networkProteinList, distanceMatrixFile);
		}
		
		/* Load distance matrix (Note: it contains disconnected components) */ 
		double[][] distanceMatrix = Loader.loadDistanceMatrix(distanceMatrixFile, networkProteinList, runComputeCanada);

		/* Determine which proteins are disconnected*/ 
		boolean[] proteinsToKeep = Calculator.determineConnectedProteins(distanceMatrix);
		
		/* Update proteins to keep in the network (ie. those that contribute to the most connected component) */
		ArrayList<Protein> networkProteinsList3 = NetworkProteins.modifyNetworkProteinsList(networkProteinList, proteinsToKeep);
		
		
		if(toUpdateDistanceMatrix) {
			/* Update distance matrix to only include proteins that contribute to the most connected component of the graph */
			DistanceMatrix.updateDistanceMatrix(proteinsToKeep, distanceMatrix, distanceMatrix2File);
		}
		
		/* Load distance matrix representing fully connected component */
		double[][] distanceMatrix2 = Loader.loadDistanceMatrix(distanceMatrix2File, networkProteinsList3, runComputeCanada);
		
		//NetworkInformation.outputGraphForMCL(networkInteractionsList, networkProteinsList3, mclGraphFile);
		//NetworkInformation.outputProteinList(networkProteinsList3, proteinListFile);
	
			
		System.out.println("break");
		if(toMonteCarlo) {
			/* Measure distributions for given number of proteins within the bounds for a given number of sampling */
			//Sampling sampling = new Sampling(annotationGoList, distanceMatrix2, weightedSampling);
			//annotatedProteinOccurenceFile.sampling.computeMultipleDistributions(nProtToSampleLowerBound, nProtToSampleUpperBound, numOfTimesNetworkIsSampled, distributionFile);
			MotifSampling sampling = new MotifSampling(annotatedProteinOccurenceFile, networkProteinsList3, distanceMatrix2);
			sampling.computeMultipleDistributions2(nProtToSampleLowerBound, nProtToSampleUpperBound, numOfTimesNetworkIsSampled, distributionFile);
		}
		
		if(computeDistributionParams) {
    		System.out.println("computing normal distribution params");
			SignificanceCalculator.computeNormalDistributionParameters(distributionFile, nProtToSampleLowerBound, nProtToSampleUpperBound, normalDistributionParametersFile);
		}
		
		if(testGoTerms) {		
			
			System.out.println("Testing go terms");
			
			/* Load annotations */ 
			ArrayList<Annotation> annotationGoList = Loader.importAnnotationGo(annotationInputFile, networkProteinsList3, runComputeCanada);
			
			SignificanceCalculator fdrCalculator = new SignificanceCalculator(annotationGoList);
        	fdrCalculator.modifyGoAnnotationsWithTPD(distanceMatrix2);
			
        	
        	
        	double minimum_pval = 1; 
        	if(useNormalDistribution) {
        		minimum_pval = fdrCalculator.modifyGoAnnotationsWithPvalueFromNormalApproximation(normalDistributionParametersFile);
        	}
        
        	/* Output motifs and their significance score */
        	Exporter.exportAnnotationInfoForFDRcalc(annotationGoList, motifSignificanceFile);
        	
        	
		}
		
   } // close main

} // end Network
