import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

public class CompareConditions {
	
	public static final boolean compareNetworkConditions = true;
	public static final boolean compareNetworkToMCL = false;
	
	public static void main(String[] args) {

		String wd = "/Users/Rachel/eclipse-files/network2020/output_files/"; // working directory
		
		/* Input files */
		String weightedNetworkFDRfile = wd + "2020.02.19.MonoTransf_FDRvPval_her2vTNweightedSampling10X7_newGOFile.txt";
		String unweightedNetworkFDRfile = wd + "2020.02.19.MonoTransf_FDRvPval_unweightedNetwork_weightedSampling10X7.txt";
		String unweightedNetworkGoSummaryFile = wd + "2020.02.19.MonoTransf_GoTerms_unweightedNetwork_weightedSampling10X7_newGOFile.txt";
		String weightedNetworkGoSummaryFile = wd + "2020.02.19.MonoTransf_GoTerms_her2vTNweightedSampling10X7_newGOFile.txt";
		
		String mclWeightedNetwork = "/Users/Rachel/eclipse-files/network2020/mcl_analysis/ontologizer_files/2020.03.24.weightedMCLi2_SignificantCorrectedGO_p0.001.txt";
		
		/* Output files */
		String overlapInfoFile = wd + "2020.02.27.overlapInfo_her2vTN_ws.txt";
		
		String unique_wnFile = wd + "2020.03.27.uniqueHer2vTNGOs_0.0005_ws.txt";
		String significantGO_wnFile = wd + "2020.03.27.significantGOHer2vTN_0.0005_ws.txt";
		String significantGO_unFile = wd + "2020.03.27.significantGOunweight_0.0005_ws.txt";
		
		String unique_networkFile = wd + "2020.03.27.uniqueNetworkGOs_FDR0.0005_ws.txt";
		String unique_mclFile = wd + "2020.03.27.uniqueMCLGOs_p0.001.txt";
		
		double fdr = 0.0005; // FDR cut off for significant GOs 
		
		if(compareNetworkConditions) {
			/* Print info for plot */ 
			HashMap<Double, Double> overlapInfoMap = getOverlapInfo(unweightedNetworkFDRfile, weightedNetworkFDRfile, unweightedNetworkGoSummaryFile, weightedNetworkGoSummaryFile);
			printOverlapInfo(overlapInfoFile, overlapInfoMap);
			
			/* Find optimal p-values based on FDR for weighted network and unweighted network */
			double un_pval = getPvalAtFdrCutOff(unweightedNetworkFDRfile, fdr);
			double wn_pval = getPvalAtFdrCutOff(weightedNetworkFDRfile, fdr);
	
			/* Load significant GO Terms for weighted network and unweighted network*/
			HashMap<String, Double> un_GoMap = loadSignificantGoTerms(unweightedNetworkGoSummaryFile, un_pval);
			HashMap<String, Double> wn_GoMap = loadSignificantGoTerms(weightedNetworkGoSummaryFile, wn_pval);
			
			printSignificantGOs(significantGO_unFile, un_GoMap);
			printSignificantGOs(significantGO_wnFile, wn_GoMap);
			
			printUniqueSignifantGOs(unique_wnFile, un_GoMap, wn_GoMap);
		}
		
		if(compareNetworkToMCL) {
			double p_val = getPvalAtFdrCutOff(weightedNetworkFDRfile, fdr);
			HashMap<String, Double> networkMap = loadSignificantGoTerms(weightedNetworkGoSummaryFile, p_val);
			
			HashMap<String, Double> mclMap = loadMCLSignificantGOs(mclWeightedNetwork);
			
			printUniqueSignifantGOs(unique_networkFile, networkMap, mclMap);
			printUniqueSignifantGOs(unique_mclFile, mclMap, networkMap);
		}
	}


	/***
	 * Reads file containing {FDR; p-value; # Go terms} and return p-value at requested FDR cut off
	 * Note: the file is ordered such that smaller FDRs are first
	 * 
	 * @param inputFile	Result metric file
	 * @param fdr		FDR cut off
	 * @return p_val 	significant p-value cut off
	 */
	private static double getPvalAtFdrCutOff(String inputFile, double fdr) {
		double p_val = 0;
		try {
			BufferedReader in = new BufferedReader(new FileReader(new File(inputFile)));

			String line = in.readLine(); // header
			line = in.readLine();

			while(line != null) {
				String[] col = line.split("\t");

				/* Store p-vals as long as it's smaller than the fdr cut off
				 * When p-val > fdr, break out of loop without storing read p-val */
				if(Double.parseDouble(col[0]) <= fdr) {
					p_val = Double.parseDouble(col[1]);
				} else {
					break;
				}
				line = in.readLine();
			}
			in.close();	
		} catch (IOException e) {
			e.printStackTrace();
		}	
		return p_val;
	}
	/***
	 * Load significant GO Terms as determined by the p-value cut-off
	 *
	 * @param inputFile		GO term summary file
	 * @param pval			p-value cut off
	 * @return goTermMap	Map of {significant Go Term; p-value}
	 */
	private static HashMap<String, Double> loadSignificantGoTerms(String inputFile, double pval){
		HashMap<String, Double> goTermMap = new HashMap<>();

		try {
			BufferedReader in = new BufferedReader(new FileReader(new File(inputFile)));

			String line = in.readLine(); // header
			line = in.readLine();

			while(line!=null) {
				String[] col = line.split("\t");
				/* store GO Term and its p-value if the p-value is smaller than the cut off */
				if(Double.parseDouble(col[3]) < pval) {
					goTermMap.put(col[0], Double.parseDouble(col[3]));
				}

				line = in.readLine();
			}
			in.close();	
		} catch (IOException e) {
			e.printStackTrace();
		}	
		return goTermMap;
	}
	/***
	 * Compares un-weighted and weighted network significant GO terms and counts the overlap
	 * 
	 * @param un_GoMap	significant GO terms in un-weighted Network map 
	 * @param wn_GoMap	significant GO terms in weighted Network map
	 * @return overlap	count 
	 */
	private static int compareGoListOverlap(HashMap<String, Double> un_GoMap, HashMap<String, Double> wn_GoMap) {
		int overlap = 0; 
		/* iterate over significant GO Terms in weighted Network, count those found also in unweighted network*/ 
		for(String GO: wn_GoMap.keySet()) {
			if(un_GoMap.containsKey(GO)) {
				overlap += 1;
			}
		}
		return overlap;
	}
	/***
	 * For given fdr thresholds, identify significant p-value and consequent GO terms. Compute unique to overlap ratio.
	 * @param unFDRfile			FDR info file for unweighted network
	 * @param wnFDRfile			FDR info file for weighted network
	 * @param unGOfile			GO summary file for unweighted network
	 * @param wnGOfile			GO summary file for weighted network
	 * @return overlapInfoMap	{fdr ; ratio} 
	 */
	private static HashMap<Double, Double> getOverlapInfo(String unFDRfile, String wnFDRfile, String unGOfile, String wnGOfile) {
		HashMap<Double, Double> overlapInfoMap = new HashMap<Double, Double>();
		
		for(double fdr=0.5; fdr>=0.01; fdr=fdr-0.01) {
			/* Find optimal p-values based on FDR for weighted network and unweighted network */
			double un_pval = getPvalAtFdrCutOff(unFDRfile, fdr);
			double wn_pval = getPvalAtFdrCutOff(wnFDRfile, fdr);

			/* Load significant GO Terms for weighted network and unweighted network*/
			HashMap<String, Double> un_GoMap = loadSignificantGoTerms(unGOfile, un_pval);
			HashMap<String, Double> wn_GoMap = loadSignificantGoTerms(wnGOfile, wn_pval);
			
			/* Compute overlap info (#overlap, #unique in weighted & ratio)*/
			int overlap = compareGoListOverlap(un_GoMap, wn_GoMap);
			int unique = wn_GoMap.size() - overlap;
			double ratio = unique / (double) overlap;
			
			overlapInfoMap.put(fdr, ratio); // store fdr and ratio in map
		}
		
		double[] fdrToIterate = new double[] {0.005, 0.001, 0.0005, 0.0001};
		for(int fdr=0; fdr<fdrToIterate.length; fdr++) {
			/* Find optimal p-values based on FDR for weighted network and unweighted network */
			double un_pval = getPvalAtFdrCutOff(unFDRfile, fdrToIterate[fdr]);
			double wn_pval = getPvalAtFdrCutOff(wnFDRfile, fdrToIterate[fdr]);

			/* Load significant GO Terms for weighted network and unweighted network*/
			HashMap<String, Double> un_GoMap = loadSignificantGoTerms(unGOfile, un_pval);
			HashMap<String, Double> wn_GoMap = loadSignificantGoTerms(wnGOfile, wn_pval);

			/* Compute overlap info (#overlap, #unique in weighted & ratio)*/
			int overlap = compareGoListOverlap(un_GoMap, wn_GoMap);
			int unique = wn_GoMap.size() - overlap;
			double ratio = unique / (double) overlap;
			
			overlapInfoMap.put(fdrToIterate[fdr], ratio); // store fdr and ratio in map
		}
	return overlapInfoMap;	
	}
	/***
	 * Print the overlap info comparing the unweighted and weighted network significant GO terms.
	 * @param outputFile		file to countain overlap info
	 * @param overlapInfoMap	map {fdr: ratio(unique weighted/overlap)}
	 */
	private static void printOverlapInfo(String outputFile, HashMap<Double, Double> overlapInfoMap) {
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			out.write("FDR\tRatio\n"); //header
			for(double fdr: overlapInfoMap.keySet()) {
				out.write(fdr + "\t" + overlapInfoMap.get(fdr) + "\n"); // body {fdr \t ratio}
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	/***
	 * Print the significant GO terms and it's associated p-value
	 * @param outputFile	file to contain significant go terms and p-values
	 * @param goMap			map {GO: p-value}
	 */
	private static void printSignificantGOs(String outputFile, HashMap<String, Double> goMap) {
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			out.write("GO\tp-value\n"); //header
			for(String GO: goMap.keySet()) {
				out.write(GO + "\t" + goMap.get(GO) + "\n"); // body {GO \t p-value}
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	/**
	 * Print unique GO terms deemed significant in the weighted network.
	 * @param outputFile	file to contain unique go terms & p-value
	 * @param un_goMap		unweighted network map {GO: p-value}
	 * @param wn_goMap		weighted network map {GO: p-value}
	 */
	private static void printUniqueSignifantGOs(String outputFile, HashMap<String, Double> un_goMap, HashMap<String, Double> wn_goMap) {
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			out.write("GO\tp-value\n"); //header
			for(String GO: wn_goMap.keySet()) {
				/* print Go terms that are not in the unweighted list*/
				if(!un_goMap.containsKey(GO)) {
					out.write(GO + "\t" + wn_goMap.get(GO) + "\n");
					out.flush();
				}
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	 
	/***
	 * Load significant GO terms identified in mcl/ontologizer approach
	 * @param mclInput_file		String input file containing GO terms and adjusted p-values of mcl/ontologizer approach
	 * @return significantGos	HashMap<String, Double> where String = GO, Double = adj_pval
	 */
	private static HashMap<String, Double> loadMCLSignificantGOs(String mclInput_file){
		HashMap<String, Double> significantGos = new HashMap<>();
		
		try {
			BufferedReader in = new BufferedReader(new FileReader(new File(mclInput_file)));
			
			String line = in.readLine(); // header
			line = in.readLine();
			
			/* Load file containing GO terms that were significant in mcl-ontologizer approach */
			while(line!=null) {
				String[] col = line.split("\t"); // col[0] = GO term, col[1] = adj_pval
				significantGos.put(col[0], Double.parseDouble(col[1]));
				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return significantGos;
	}
}
