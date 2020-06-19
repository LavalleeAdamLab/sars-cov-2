import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class FilterMotifsByFDR {

	public static void main(String[] args) {
		
		String inputMotifs = "C:\\Users\\Rachel\\Documents\\COVID\\NetworkMotif-Krogan\\krogan_MotifAnnotations_details_s10X5.txt";
		String significantMotifs = "C:\\\\Users\\\\Rachel\\\\Documents\\\\COVID\\\\NetworkMotif-Krogan\\significantMotifswDetails_FDR0.05.txt";
		
		double threshold = 0.000218000000000065;
		
		filterMotifs(inputMotifs, significantMotifs, threshold);
	}
	
	private static void filterMotifs(String inputFile, String outputFile, double threshold) {
		
		try {
			BufferedReader in = new BufferedReader(new FileReader(new File(inputFile)));
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			String line = in.readLine(); // header
			line = in.readLine();
			
			while(line != null) {
				String[] col = line.split("\t");
				
				if(Double.parseDouble(col[1]) < threshold) {
					out.write(col[0] + "\t" + col[1] + "\n");
					out.flush();
				}
				line = in.readLine();
			}
			
			out.close();
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
}
