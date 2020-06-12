package utils;

import graph.Annotation;

import java.util.ArrayList;

public class SignificanceCalculator {

    /* original goTerms list */
    private ArrayList<Annotation> goAnnotations;

    /**
     * Constructor of the class.
     *
     * @param goAnnotations list of goTerms
     * @param shuffledGoAnnotations list of shuffled goTerms
     */
    public SignificanceCalculator(ArrayList<Annotation> goAnnotations) {
        this.goAnnotations = goAnnotations;
    }

    /**
     * Modifies every goTerm in both shuffled and original lists with TPD.
     *
     * @param distanceMatrix distance matrix which is used to calculate TPD.
     */
    public void modifyGoAnnotationsWithTPD(double[][] distanceMatrix) {
        Modifier.setClusterTPD(goAnnotations, distanceMatrix); // Annotation Go List
    }

    /**
     * Modifies every goTerm in both shuffled and original lists with a p-value. P-value is calculated using normal
     * approximation (which is built according to previously computed parameters), so it is obligatory to have a file
     * with normal distribution parameters before calling this function. If you don't have this file, you can call
     * computeNormalDistributionParameters and then pass the name of the new file to this function.
     *
     * @param distributionParametersFilePath path to the file with distribution parameter
     */
    public double modifyGoAnnotationsWithPvalueFromNormalApproximation(String distributionParametersFilePath) {
    	
    	double[] minimum_pvals = new double[2];
        minimum_pvals[0] = NormalApproximation.importNormalDistributionParameters(goAnnotations, distributionParametersFilePath);
        
        double min_pval = Math.min(minimum_pvals[0], minimum_pvals[1]);
        
        return min_pval;
    }

    /**
     * Computes all normal distributions parameters and stores them in the file.
     *
     * @param distributionFilePath 				path to the file with distributions
     * @param nProtToSampleUpperBound 			largest amount of protein for which we will compute distribution parameters
     * @param distributionParametersFilePath 	path to the file with output (computed parameters)
     */
    public static void computeNormalDistributionParameters(String distributionFilePath, int nProtToSampleLowerBound, int nProtToSampleUpperBound, String distributionParametersFilePath) {
        Loader.loadMonteCarloDistributions(distributionFilePath, nProtToSampleLowerBound, nProtToSampleUpperBound, distributionParametersFilePath);
    	//Loader.loadDistributions2(distributionFilePath, nProtToSampleLowerBound, nProtToSampleUpperBound, distributionParametersFilePath);
    }
}
