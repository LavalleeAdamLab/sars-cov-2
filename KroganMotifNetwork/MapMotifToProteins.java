import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

public class MapMotifToProteins {

	public static void main(String[] args) {
		
		String inputMotifFile = "C:\\Users\\Rachel\\Documents\\COVID\\NetworkMotif-Krogan\\krogan_MotifAnnotations_details_s10X5.txt";
		String inputAnnotationFile = "C:\\Users\\Rachel\\Documents\\COVID\\NetworkMotif-Krogan\\input_files\\Krogan_motif_annotations_real_FILTERED.tsv";
		String outputMotifProteinMappingFile = "C:\\Users\\Rachel\\Documents\\COVID\\NetworkMotif-Krogan\\krogan_significantMotifs_details_FDR0.5_proteinsMapped.txt";
		String outputProteinFrequencyFile = "C:\\Users\\Rachel\\Documents\\COVID\\NetworkMotif-Krogan\\krogan_significantProteinFrequencies_FDR0.1.txt";
		double threshold = 0.000218000000000065;
		
		HashMap<String, String> significantMotifs = loadMotifs(inputMotifFile, threshold);
		HashMap<String, Integer> protFrequencyMap = mapSignificantMotifsToProteins(significantMotifs, inputAnnotationFile, outputMotifProteinMappingFile);
		//printProteinFrequencies(protFrequencyMap, outputProteinFrequencyFile);
	}
	
	/**
	 * Load file containing motifs and their obtained p-values, 
	 * store significant motifs (p-value smaller than threshold) in hash map
	 * 
	 * @param inputFile				String - File path for motifs and obtained p-values
	 * @param threshold				double - p-value threshold determined by FDR cut-off
	 * @return significantMotifs 	HashMap - {motif : p-value}
	 */
	private static HashMap<String, String> loadMotifs(String inputFile, double threshold) {
		
		HashMap<String, String> significantMotifs = new HashMap<>();
		
		try {
			BufferedReader in = new BufferedReader(new FileReader(new File(inputFile)));
			
			String line = in.readLine(); // header
			line = in.readLine();
			
			while(line != null) {
				String[] col = line.split("\t");
				
				if(Double.parseDouble(col[3]) < threshold) {
					String info = col[1] + "\t" + col[2] + "\t" + col[3];
					significantMotifs.put(col[0], info);
				}
				
				line = in.readLine();
			}
			
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return significantMotifs;
	}
	
	private static HashMap<String, Integer> mapSignificantMotifsToProteins(HashMap<String, String> significantMotifs, String annotationFile, String outputFile) {
		HashMap<String, Integer> protFreqMap = new HashMap<>();
		
		try {
			BufferedReader in = new BufferedReader(new FileReader(new File(annotationFile)));
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			String line = in.readLine(); // header
			line = in.readLine();

			out.write("Motif\tNumOfProteins\tTPD\tP-value\tProtAccessions\tProteinNames\n");
			
			while(line != null) {
		
				String[] col = line.split("\t");

				if(significantMotifs.containsKey(col[0])) {
					String[] proteins = col[8].split("\\|");
					out.write(col[0] + "\t" + significantMotifs.get(col[0]) + "\t" + col[7] +"\t" + col[8] + "\n");	
					
					for(int i=0; i<proteins.length; i++) {
						if(protFreqMap.containsKey(proteins[i])) {
							int count = protFreqMap.get(proteins[i]) + 1;
							protFreqMap.put(proteins[i], count);
						} else { 
							protFreqMap.put(proteins[i], 1);
						}
					}
				}
				line = in.readLine();
			}

			out.close();
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return protFreqMap;
	} 
	
	private static void printProteinFrequencies(HashMap<String, Integer> protFreqMap, String outputFile) {
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			out.write("Protein\tFrequency\n");
			for(String protein: protFreqMap.keySet()) {
				out.write(protein + "\t" + protFreqMap.get(protein) + "\n");
				out.flush();
			}
			
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
}
