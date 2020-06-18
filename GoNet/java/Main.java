import java.util.ArrayList;

import graph.Annotation;
import graph.FalseDiscoveryRate;
import graph.Interaction;
import graph.Protein;
//import opt.NetworkInformation;
import utils.*;


public class Main {

	/*** Constants method variables ***/
	
	public static final boolean toComputeDistanceMatrix = true;
	public static final boolean toUpdateDistanceMatrix = true;
	
	public static final boolean toMonteCarlo = true;
		public static final boolean weightedSampling = true;
		public static final boolean testNormalApprox  = true;
	
	public static final boolean testGoTerms = true;
		public static final boolean printFDRresults = true;
		public static final boolean printSignificantGOterms = true;

	/**** Constant output variables ***/
	public static final boolean printNumberOfInteractionsPerProtein = false;
	public static final boolean printPPInetwork = false;
	
	public static final boolean toShuffle = true;
		public static final boolean printShuffledGOList = true;
	
	/**** Constant  distribution variables ****/ 
	public static final boolean computeDistributionParams = true; //if true don't recalculate parameters
	public static final boolean useNormalDistribution = true; //if false uses monte carlo values
	
	public static void main(String[] args) {
		
		/*****************************************************************************************************
		 * This algorithm assess the clustering of proteins in the String augmented krogan network. 
		 * We assess the clustering significance of a given set of proteins against a Monte Carlo approximated
		 * Normal Distribution. A false discovery rate is estimated at various thresholds.
		 *****************************************************************************************************/

		int nProtToSampleLowerBound = Integer.parseInt(args[0]);		// Lower limit of number of proteins sampled
		int nProtToSampleUpperBound = Integer.parseInt(args[1]);		// Upper limit of number of proteins sampled
		int numOfTimesNetworkIsSampled = Integer.parseInt(args[2]);		// Number of times the network is sampled during Monte Carlo Method

		String working_directory = "C:/Users/Rachel/Documents/COVID/NetworkGO-Krogan/";

		String PPI_inputFile = working_directory + "input_files/KroganNetwork_EntrezGeneID_GeneName_Organism.txt";
		String annotationGO_inputFile = working_directory + "input_files/GO_annotations-9606-inferred-allev-2.tsv";
			
		String distanceMatrixFile = working_directory + "IO_files/krogan_DistanceMatrix.txt";
		String distanceMatrix2File = working_directory + "IO_files/krogan_DistanceMatrix_FullConnected.txt";

		String distributionFile = working_directory + "IO_files/kroganDistribution_s10X5_n3_193.txt";

		/* Output files */
		String fdrExportFile = working_directory + "output_files/2020.04.08.MonoTransf_kroganFDR_10X5.txt";
		String goExportFile = working_directory + "output_files/2020.04.08.MonoTransf_kroganGoTerms_s10X5.txt";

		String shuffledGOFile = working_directory + "IO_files/shuffled-go.txt";
		String normalDistributionParametersFile = working_directory + "IO_files/krogan_normalDistributionParams_n3_193.txt";
		
		//String mclGraphFile = working_directory + "output_files/unweightedMCLgraph.txt";
//		String proteinListFile = working_directory + "output_files/proteinsInNetworkList.txt";
//		String nInteractionsPerProtFile = "/Users/Rachel/eclipse-files/network2.0/output_files/proteins.txt";


		/* Load Krogan Interactions. Extract all possible interactions between proteins and store in Interaction object */
		ArrayList<Interaction> networkInteractionsList = NetworkInteractionsLoader.importInteractionNetwork(PPI_inputFile);
		
		/* Extract all proteins in the network (ie. all proteins within the interaction file) */
		ArrayList<Protein> networkProteinList = NetworkProteins.getProteinsInNetwork(networkInteractionsList);

		if(toComputeDistanceMatrix) {
			/* Compute distance matrix using the list of proteins in the network and the list of known interactions in the network. */
			DistanceMatrix.computeDistanceMatrix(networkInteractionsList, networkProteinList, distanceMatrixFile);
		}
		
		/* Load distance matrix (Note: it contains disconnected components) */ 
		double[][] distanceMatrix = Loader.loadDistanceMatrix(distanceMatrixFile, networkProteinList);

		/* Determine which proteins are disconnected*/ 
		boolean[] proteinsToKeep = Calculator.determineConnectedProteins(distanceMatrix);
		
		/* Update proteins to keep in the network (ie. those that contribute to the most connected component) */
		ArrayList<Protein> networkProteinsList3 = NetworkProteins.modifyNetworkProteinsList(networkProteinList, proteinsToKeep);

		if(toUpdateDistanceMatrix) {
			/* Update distance matrix to only include proteins that contribute to the most connected component of the graph */
			DistanceMatrix.updateDistanceMatrix(proteinsToKeep, distanceMatrix, distanceMatrix2File);
		}
		
		/* Load distance matrix representing fully connected component */
		double[][] distanceMatrix2 = Loader.loadDistanceMatrix(distanceMatrix2File, networkProteinsList3);
		
		//NetworkInformation.outputGraphForMCL(networkInteractionsList, networkProteinsList3, mclGraphFile);
		//NetworkInformation.outputProteinList(networkProteinsList3, proteinListFile);
		/* Load annotations */ 
		ArrayList<Annotation> annotationGoList = Loader.importAnnotationGo(annotationGO_inputFile, networkProteinsList3);
		if(toMonteCarlo) {
			/* Measure distributions for given number of proteins within the bounds for a given number of sampling */
			Sampling sampling = new Sampling(annotationGoList, distanceMatrix2, weightedSampling);
	        sampling.computeMultipleDistributions(nProtToSampleLowerBound, nProtToSampleUpperBound, numOfTimesNetworkIsSampled);
		}
			
		ArrayList<Annotation> shuffled_goAnnotationList = new ArrayList<Annotation>();
		if(toShuffle) {	
			/**
			 * Could modify to check if shuffled file already exists? Removes the need of an if statement. 
			 */
			/* Shuffle proteins associated to annotations */
			shuffled_goAnnotationList = Calculator.shuffleGoProtAssociations2(annotationGoList);
			if (printShuffledGOList) {
				Exporter.printShuffledNetwork(shuffled_goAnnotationList, shuffledGOFile);
			}
		} else {
			/* Load shuffled annotations */ 
            Loader.loadShuffledGoAnnotations(shuffled_goAnnotationList, shuffledGOFile);
		}
		
		if(testGoTerms) {		
		
			System.out.println("Testing go terms");

			FdrCalculator fdrCalculator = new FdrCalculator(annotationGoList, shuffled_goAnnotationList);
        	fdrCalculator.modifyGoAnnotationsWithTPD(distanceMatrix2);
			
        	if(computeDistributionParams) {
        		System.out.println("computing normal distribution params");
				fdrCalculator.computeNormalDistributionParameters(distributionFile, nProtToSampleLowerBound, nProtToSampleUpperBound, normalDistributionParametersFile);
			}
        	
        	double minimum_pval = 1; 
        	if(useNormalDistribution) {
        		minimum_pval = fdrCalculator.modifyGoAnnotationsWithPvalueFromNormalApproximation(normalDistributionParametersFile);
        	} else {
        		minimum_pval = Loader.setAnnotationPvaluesFromMonteCarloDistribution(distributionFile, nProtToSampleLowerBound, nProtToSampleUpperBound, annotationGoList, shuffled_goAnnotationList);
        	}
        
        	
  			ArrayList<FalseDiscoveryRate> fdrs = fdrCalculator.computeFdr(minimum_pval);
  			fdrs = fdrCalculator.monotonicTransformationForFdr(fdrs);
  			Modifier.modifyClustersFdr(annotationGoList, fdrs);
  			
			if (printSignificantGOterms) {
				System.out.println("Print GO results"); 
				Exporter.printGO_results(annotationGoList, shuffled_goAnnotationList, goExportFile);
			}
				
			if (printFDRresults) {
				System.out.println("Print FDR"); 
				Exporter.testGoFDR(fdrs, fdrExportFile);
			}
		}
		
   } // close main

} // end Network
