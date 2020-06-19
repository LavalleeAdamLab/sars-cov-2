import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;

public class Ontologizer {

	public static void main(String[] args) {
		String working_directory = "/Users/Rachel/eclipse-files/network2020/mcl_analysis/ontologizer_files/";
		
		String testedGOs_file = working_directory + "TestedGOs.txt";
		String gaf_file = working_directory + "goa_human.gaf";
		String updatedGAF_file = working_directory + "goa_humanFiltered.gaf";
		
		HashSet<String> testedGOSet = loadTestedGOs(testedGOs_file);
		updateGAF(gaf_file, updatedGAF_file, testedGOSet);
		
	}
	
	/***
	 * Load list of GO terms that were tested in the network.
	 * @param testedGOs_file	list of GO terms; every line contains (1) GO:######
	 * @return testedGOSet		HashSet<String> format of list
	 */
	private static HashSet<String> loadTestedGOs(String testedGOs_file){
		/* Read tested GOs into HashSet */
		HashSet<String> testedGOSet = new HashSet<String>();
		try {
			BufferedReader in = new BufferedReader(new FileReader(new File(testedGOs_file)));
			/* Load all tested GO terms in hash set*/
			String line = in.readLine(); //line = {<String> GOterm}
			while(line!=null) {
				testedGOSet.add(line); //add go term to set
				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return testedGOSet;
	}
	
	/***
	 * Update GAF file to contain only GOs that were tested in the network approach
	 * @param gaf_file			tab delimited file with information about proteins and associated GO:####
	 * @param updatedGAF_file	file to contain only lines from GAF file only if associated GO:#### was studied in network approach
	 * @param testedGOSet		HashSet format list of tested GO terms in the network approach
	 */
	private static void updateGAF(String gaf_file, String updatedGAF_file, HashSet<String> testedGOSet) {
		try {
			BufferedReader in = new BufferedReader(new FileReader(new File(gaf_file)));
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(updatedGAF_file)));
			/* Print lines from GAF file in updated file if GO term of a given line was tested in our approach */
			String line = in.readLine();
			while(line != null) {
				if(line.startsWith("!")) { // print header lines (that start w/ !)
					out.write(line + "\n");
				} else { // filter content of GAF
					String[] col = line.split("\t"); 
					if(testedGOSet.contains(col[4])) { // col[4] = {GO:#####}
						out.write(line + "\n");
						out.flush();
					}
				}
				
				line = in.readLine();
			}
			in.close();
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
}
