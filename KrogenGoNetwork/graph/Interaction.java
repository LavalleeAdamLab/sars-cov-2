package graph;

public class Interaction {
		
	private final static double DEFAULT_WEIGHT = 1.0;
	
		// fields
		private String proteinName1;
		private String proteinName2;
		private int entrezID1;
		private int entrezID2;
		private double weight;

		// constructors
		public Interaction(String _Protein1, String _Protein2, int _ID1, int _ID2) {
			this(_Protein1, _Protein2, _ID1, _ID2, DEFAULT_WEIGHT);
		}
		
		public Interaction(String _Protein1, String _Protein2, int _ID1, int _ID2, double _w) {
			proteinName1 = _Protein1;
			proteinName2 = _Protein2;
			entrezID1 = _ID1;
			entrezID2 = _ID2;
			weight = _w;
		}

		// get
		public String getProtein1() {
			return proteinName1;
		}

		public String getProtein2() {
			return proteinName2;
		}

		public int getID1() {
			return entrezID1;
		}

		public int getID2() {
			return entrezID2;
		}

		public double getWeight() {
			return this.weight;
		}
		
		public void setWeight(double w) {
			this.weight = w;
		}
}

