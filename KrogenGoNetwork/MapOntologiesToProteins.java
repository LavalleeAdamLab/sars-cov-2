import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class MapOntologiesToProteins {

	public static void main(String[] args) {

		/* File paths*/ 
		String significantGOsFile = "C:/Users/Rachel/Documents/COVID/NetworkGO-Krogan/output_files/2020.04.08.MonoTransf_kroganGoTerms_s10X5.txt";
		String fdrThresholdFile = "C:/Users/Rachel/Documents/COVID/NetworkGO-Krogan/output_files/2020.04.08.MonoTransf_kroganFDR_10X5.txt";
		String annotationFile = "C:/Users/Rachel/Documents/COVID/NetworkGO-Krogan/input_files/GO_annotations-9606-inferred-allev-2.tsv";
		String proteinsInNetworkFile = "C:/Users/Rachel/Documents/COVID/NetworkMotif-Krogan/output_files/2020.04.07.ListOfKroganProteinsInConnectedComponent.txt";
		String significantGOwDetailsFile = "C:/Users/Rachel/Documents/COVID/NetworkGO-Krogan/output_files/2020.06.11.kroganGoTermsMappedToProteins_wDetails_s10X5.txt";
		
		/* Load tested GO terms {GO term : p-value} */
		HashMap<String, Double> significantGOMap = loadSignificantAnnotations(significantGOsFile);
		HashMap<String, String> proteinsInNetworkMap = loadProteinsInNetwork(proteinsInNetworkFile);
		
		ArrayList<ArrayList<Double>> fdrThresholdList = loadFdrThresholds(fdrThresholdFile);
		
		/* Read in GO annotations & output to folder extra info to significant GO terms */
		printAnnotationDetails(annotationFile, significantGOwDetailsFile, significantGOMap, proteinsInNetworkMap, fdrThresholdList);
		
	}

	/**
	 * Load significant GO annotations deemed by statistical approach into a HashMap {GO-term: p-value} 
	 * @param inputFile					String - file path to filtered file containing significant GO-terms
	 * @return significantAnnotations	HashMap<String, Double> - map {Go-term: p-values}
	 */
	private static HashMap<String, Double> loadSignificantAnnotations(String inputFile){
		
		HashMap<String, Double> significantAnnotations = new HashMap<>();
	
		try {
			BufferedReader in = new BufferedReader(new FileReader(new File(inputFile)));
			String line = in.readLine(); // HEADER
			line = in.readLine();
			while(line!=null) {
				String[] col = line.split("\t");
				
				/* col[0] = go annotation, col [3] = p-value*/
				significantAnnotations.put(col[0], Double.parseDouble(col[3]));
				line = in.readLine();
			}
			in.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
		return significantAnnotations;
	}
	
	/**
	 * Load list of proteins that are tested in our network approach
	 * @param inputFile					String - file path to list of proteins in network
	 * @return proteinsInNetworkMap 	HashMap<String,String> - maps {accession : protein} 
	 */
	private static HashMap<String, String> loadProteinsInNetwork(String inputFile){
		
		HashMap<String, String> proteinsInNetworkMap = new HashMap<>();
		
		try {
			BufferedReader in = new BufferedReader(new FileReader(new File(inputFile)));
			String line = in.readLine(); // no HEADER
			while(line!=null) {
				String[] col = line.split("\t");
				
				/* col[0] = accession, col [1] = protein*/
				proteinsInNetworkMap.put(col[1], col[0]);
				line = in.readLine();
			}
			in.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
		return proteinsInNetworkMap;
	}
	
	private static ArrayList<ArrayList<Double>> loadFdrThresholds(String inputFile){
		
		ArrayList<ArrayList<Double>> fdrThresholdList = new ArrayList<ArrayList<Double>>();
		
		try {
			BufferedReader in = new BufferedReader(new FileReader(new File(inputFile)));
			String line = in.readLine(); // HEADER
			line = in.readLine();
			
			while(line!=null) {
				String[] col = line.split("\t");
					ArrayList<Double> fdrThreshold = new ArrayList<>();
					fdrThreshold.add(Double.parseDouble(col[0]));
					fdrThreshold.add(Double.parseDouble(col[1]));
					
					fdrThresholdList.add(fdrThreshold);
				line = in.readLine();
			}
			in.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
		return fdrThresholdList;
	}
	
	/**
	 * Output details of significant GO annotations deemed by statistical approach:
	 * Load initially used annotations; output information only for GO annotations deemed significant
	 * Consideration: insure only protein information that was studied in our approach
	 * 
	 * @param annotationFile		String - file path for all annotations tested
	 * @param outputFile			String - file path for output
	 * @param significantGOMap		HashMap<String, Double> - maps {significant GO-term: p-value}
	 * @param proteinsInNetworkMap	HashMap<String, String> - maps {accession: protein}
	 */
	private static void printAnnotationDetails(String annotationFile, String outputFile, HashMap<String, Double> significantGOMap, HashMap<String, String> proteinsInNetworkMap, ArrayList<ArrayList<Double>> fdrThresholdsList) {
		try {
			BufferedReader in = new BufferedReader(new FileReader(new File(annotationFile)));
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			String line = in.readLine(); //HEADER
			line = in.readLine();
			
			out.write("GO-term\tName\tCategory\tNumberOfProteinsInNetwork\tP-value\tFDR\tProteinAccessions\tProteinNames\n"); // header output file
			
			while(line!=null) {	
				String[] col = line.split("\t");
								
				if(significantGOMap.containsKey(col[0])) {
					/* iterate through proteins of given annotation to only output those part of the network */ 
					
					double fdr = 1; 
					for(int i=0; i<fdrThresholdsList.size(); i++) {
						ArrayList<Double> thresholds = fdrThresholdsList.get(i);
						
						if(thresholds.get(1)>significantGOMap.get(col[0])) {
							fdr = thresholds.get(0);
							if(fdr > 0) {
								break;	
							}
						}
					}
					
					String[] proteins = col[7].split("\\|");
					
					String updatedAccessions = "";
					String updatedProteins = "";
					
					int protCount = 0;
					
					for(int i=0; i<proteins.length; i++) {
						if(proteinsInNetworkMap.containsKey(proteins[i])) {
							updatedAccessions += ( proteinsInNetworkMap.get(proteins[i]) + "|" );
							updatedProteins += ( proteins[i] + "|" );
							protCount++;
						}
					}
					out.write(col[0] + "\t" + col[1] + "\t" + col[2]+ "\t" + protCount + "\t" + 
					+ significantGOMap.get(col[0]) + "\t" + fdr + "\t" +  updatedAccessions + "\t" + updatedProteins + "\n");
					out.flush();
				}
				line = in.readLine();
			}
			in.close();
			out.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
	}
}
