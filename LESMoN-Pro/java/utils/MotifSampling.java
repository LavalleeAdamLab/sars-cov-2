package utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

import graph.Protein;

public class MotifSampling {
	/**
	 * Weighted motif sampling inspired by: https://stackoverflow.com/questions/1761626/weighted-random-numbers
	 * Here, we generate a list of cumulative weights to sample from. 
	 * We generate a random number between 0 and the sum of cumulative weights. 
	 * The selected protein will be the one where its assigned cumulative weight is greater or equal to the random number.
	 */

	double[][] distanceMatrix;
	ArrayList<Protein> proteinsInNetworkList;
	HashMap<String, Integer> proteinToOccurrenceMap; 	// Map containing annotated proteins and their occurrence of annotation
	int[] cumulativeSumOfWeights;						// array containing the cumulative sum of weights for random weighted selection

	public MotifSampling(String inputFile, ArrayList<Protein> protList, double[][] dm) {
		distanceMatrix = dm;
		proteinsInNetworkList = protList;
		proteinToOccurrenceMap = loadProteinOccurrenceList(inputFile);
		cumulativeSumOfWeights = computeCumulativeSumOfWeights();
	}

	/** 
	 * Load annotated proteins and their occurrence as a HashMap. 
	 * @param inputFile			Text File containing annotated proteins and their occurrence
	 * @return weightedList 	HashMap<String, Integer> mapping Annotated protein: occurrence of annotation
	 */
	private HashMap<String, Integer> loadProteinOccurrenceList(String inputFile){
		HashMap<String, Integer> weightedList = new HashMap<String, Integer>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));

			String line = in.readLine();
			while(line != null) {
				String[] col = line.split("\t");
				weightedList.put(col[0], Integer.parseInt(col[1])); // col[0] = protein ID, col[1] = occurrence
				line = in.readLine();
			}
			in.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
		return weightedList;
	}

	/**
	 * Generate the cumulative sum of weights list in the same order as the proteins in network list,
	 * ensuring that the index of the weight corresponds to it's network proteins index for ease of look up. 
	 * 
	 * @return cumulativeWeightList	list of cumulative weights.
	 */
	private int[] computeCumulativeSumOfWeights() {
		int[] cumulativeWeightList = new int[this.proteinsInNetworkList.size()]; // initialize list the size of network list
		int cumulativeWeight = 0; // initialize cumulative weight

		/* iterate all proteins in the order of the network protein list; 
		 * obtain weight of protein as it's annotation occurrence,
		 * update it's weight as a cumulative weight and set in list */
		for(int i=0; i<proteinsInNetworkList.size(); i++) {
			int currentWeight = this.proteinToOccurrenceMap.get(this.proteinsInNetworkList.get(i).getProteinId());
			cumulativeWeight += currentWeight; 
			cumulativeWeightList[i] = cumulativeWeight;
		}
		return cumulativeWeightList;
	}

	/**
	 * Computes distributions for multiple number of proteins in range from go_start to go_stop. 
	 *
	 * @param nProtToSampleLowerBound		beginning of the range of proteins that will be sampled
	 * @param nProtToSampleUpperBound		end of the range of proteins that will be sampled
	 * @param numOfTimesNetworkIsSampled 	number of times to calculate the distribution for each amount of proteins
	 */
	public void computeMultipleDistributions2(int nProtToSampleLowerBound, int nProtToSampleUpperBound, int numOfTimesNetworkIsSampled, String mcFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(mcFile)));

			for (int n = nProtToSampleLowerBound; n <= nProtToSampleUpperBound; n++) { // range of proteins to sample
				/* Given a GO Term is associated with n number of proteins, we perform Monte Carlo Sampling
				 * and output the resulting distribution to the console */

				HashMap<Double, Double> distribution = compute_distributionTPD(n, numOfTimesNetworkIsSampled);
				System.out.println("Computing TPD for " + n + " proteins");
				out.write("TPD (n = " + n + ")" + "\t" + "Frequency" + "\n");
				for (double dist : distribution.keySet()) {
					out.write(dist + "\t" + distribution.get(dist) + "\n");
					out.flush();

				} 
				out.write("\n");
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Computes the distribution of TPD for certain amount of randomly selected proteins.
	 *
	 * @param distanceMatrix 			distance matrix of proteins shortest pairwise distances
	 * @param timesToSampleNetwork   	number of iterations to sample Network
	 * @return the distribution of TPD for certain amount of randomly selected proteins
	 */
	private HashMap<Double, Double> compute_distributionTPD(int numProteinsToSample, int timesToSampleNetwork) {

		HashMap<Double, Double> distribution;

		/* Sample network for x amounts of proteins, y amount of times */
		double[] tpdSampleList = sampleNetwork(numProteinsToSample, timesToSampleNetwork);

		/* Measure the frequencies of sampled total pairwise distances from x amount of proteins */
		distribution = Sampling.computeFrequenciesOfSampledTPDs(tpdSampleList, timesToSampleNetwork);

		/* Sanity check to ensure the sum of TPD frequencies equals 1 */
		Sampling.checkFrequencyTotal(distribution, numProteinsToSample);

		return distribution;
	}
	
	/**
	 * Sample network X amount of times for Y number of proteins. 
	 * 
	 * @param numProteinsToSample		Proteins to sample in network
	 * @param timesToSampleNetwork		Number of iterations to sample network
	 * @return	list of TPDs measured in network
	 */
	private double[] sampleNetwork(int numProteinsToSample, int timesToSampleNetwork) {
		double[] tpdSampleList = new double[timesToSampleNetwork]; // array to store the results of sampling {list of TPDs}

		/* Randomly select proteins from the network and compute their total pairwise distance
		 * as many times as the network is to be sampled (e.g. 10X7) */
		for (int i = 0; i < timesToSampleNetwork; i++) {

			/* select proteins from the weighted list (ie. proteins are proportional to their occurrence in annotation list */
			ArrayList<Integer> randomProteins = getRandomWeightedProteins(numProteinsToSample);
			
			/* compute the total pairwise distance from the proteins selected above */
			tpdSampleList[i] = Calculator.computeTPD(distanceMatrix, randomProteins);
		}

		return tpdSampleList;
	}
	
	/**
	 * Randomly selects proteins by identifying a random weight in the cumulative distribution
	 *
	 * @param numProteinsToSample number of proteins to select from network
	 * @return array of selected protein indexes
	 */
	private ArrayList<Integer> getRandomWeightedProteins(int numProteinsToSample) {
		Random r = new Random(); 
		ArrayList<Integer> randomProteins = new ArrayList<>();

		/* Selection process occurs until the number of selected proteins equals number of proteins to sample from */ 
		while (randomProteins.size() < numProteinsToSample) {
			int sumCumulativeWeight = this.cumulativeSumOfWeights[this.cumulativeSumOfWeights.length -1]; 
			/* select a random weight between 1 and sum of cumulative weight */
			int selectedWeight = r.nextInt(sumCumulativeWeight); 
			/* protein corresponding to selected random weight will be the first protein to be greater or equal to the weight */
			for(int protIndex=0; protIndex<this.cumulativeSumOfWeights.length; protIndex++) {
				if(this.cumulativeSumOfWeights[protIndex] >= selectedWeight && !randomProteins.contains(protIndex)) {
					randomProteins.add(protIndex);
					break; // when protein index is identified break out of the loop
				}
			}
		}
		return randomProteins;
	}

}

