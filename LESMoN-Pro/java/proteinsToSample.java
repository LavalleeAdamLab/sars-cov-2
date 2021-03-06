import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.HashMap;

public class proteinsToSample {

	public static void main(String[] args) {

		/***
		 * Given an annotation list of type 
		 * {annotation; protein id list; protein symbol list}
		 * 
		 * We want to compute the frequency of every protein in the list, and
		 * we want to output these frequencies in a text file
		 * to be used for a weighted Monte Carlo Sampling approach.
		 */
		
		String input_annotationList = "C:/Users/Rachel/Documents/COVID/NetworkMotif-Krogan/input_files/Krogan_motif_annotations_shuffled_FILTERED.tsv";
		String output_frequencyList = "C:/Users/Rachel/Documents/COVID/NetworkMotif-Krogan/IO_files/annotatedProteinsToOccurence_shuffled.txt";
		
		HashMap<String, Integer> proteinToFrequencyMap = computeFrequencyOfProteins(input_annotationList);
		printProteinsToOccurence(output_frequencyList, proteinToFrequencyMap);
		
	}
	
	/**
	 * Load annotation list and store in HashMap proteins and their number of occurrence in annotation list.
	 * @param inputFile					text file containing annotation list {annotation; protein id list; protein symbol list}
	 * @return proteinToFrequencyMap	HashMap<String, Integer> {protein id: number of occurrence}
	 */
	private static HashMap<String, Integer> computeFrequencyOfProteins(String inputFile){
		
		HashMap<String, Integer> proteinToFrequencyMap = new HashMap<String,Integer>(); // Map to contain {Protein ID; occurrence} 
		
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));
			
			//int count = 0;
			String line = in.readLine();
			line = in.readLine();
			while(line!=null) {
				String[] col = line.split("\t"); 
				String[] protein_ids = col[7].split("\\|"); // 7th column holds gene IDs (proteins IDs) (index = 6)
				
				/* For all proteins of a given annotation; if in list update number of occurrence, otherwise initialize */
				for(int i=0; i<protein_ids.length; i++) {
					if(proteinToFrequencyMap.containsKey(protein_ids[i])) {
						proteinToFrequencyMap.replace(protein_ids[i], proteinToFrequencyMap.get(protein_ids[i]) + 1);
					} else {
						proteinToFrequencyMap.put(protein_ids[i], 1);
					}
				}
				line = in.readLine();
				//count++;
			}			
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return proteinToFrequencyMap;
	}
	
	/**
	 * Output annotated proteins and their occurrences in a text file
	 * @param outputFile				text file to contain Protein : Occurrence
	 * @param proteinToOccurrenceMap	map of {protein: occurrence}
	 */
	private static void printProteinsToOccurence(String outputFile, HashMap<String, Integer> proteinToOccurrenceMap) {
		try {
			OutputStreamWriter out = new OutputStreamWriter(new FileOutputStream(new File(outputFile)));
			/* iterate through map */
			for(String protein: proteinToOccurrenceMap.keySet()) {
				out.write(protein + "\t" + proteinToOccurrenceMap.get(protein) + "\n");
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
}