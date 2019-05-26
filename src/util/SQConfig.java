package util;



public class SQConfig {
	/** ============================ basic threshold ================ */
	/** number of K */
	public static int K =15;
	/** N for Top-N */
	public static int TOPN =5;
	/** input file path */

	public static String dataset = "E://data/smtp.txt";    

	/** number of dimensions */
	public static int dims =3;
	
	/*------------KDDcup------------------------------*
//	/** domain range */
	public static float[] domains = { -50.0f, 50.0f };
//	/** domain Space */
	public static float[] domainSpace = { -50.0f, 50.0f, -50.0f, 50.0f,-50.0f, 50.0f};

	/** metric space */
	public static final String strMetricSpace = "vector";
	/** metric */
	public static final String strMetric = "L2Metric";

	/** fixed bandwidth */
	public static float h = 1.5f;
	/** current threshold of Top-N LOF */
	public static float thresholdLof = 0.0f;
	/** Number of window*/
	public static int numWin =600;
	/** Number of slide*/
	public static int slide =100;	
	/** Number of Number of reads*/
	public static int numOnce =50000;    //numOnce   >= numWin

	
	/** start time of window*/
	public static int start = 1;
	/** end time of window*/
	public static int end = 1;
	
	/** output file path */
	public static String outputFile = "result.csv";
	public static int countPointBasedPruned = 0;
	/** Number of threads */
	public static int numThreads = 50;
	
	/**
	 * ============================ parameters for Multi-dimensional Pruning
	 * ==================
	 */

	/**
	 * index of independent dims, have to change when running multidimensional
	 * data
	 */
	public static final int[] indexOfIndependentDims = { 0, 1,2};
	/**
	 * number of independent dims, have to change when running multidimensional
	 * data
	 */
	public static final int num_independentDim = 3;
	/**
	 * number of correlated dims, have to change when running multidimensional
	 * data
	 */
	public static final int num_correlatedDim = 0;
	/**
	 * index of correlated dims, have to change when running multidimensional
	 * data
	 */
	public static final int[] indexOfCorrelatedDims = {};

	/**
	 * ============================ parameters for CF-Tree ==================
	 */
	/**
	 * number of max entries in each node, set to be large since we don't limit
	 * the max entry size
	 */
	public static final int maxNodeEntries = 400000000;
	/** initial distance threshold (= sqrt(radius)) */
	public static double clusterRadiusRate = 0.1;

	/** domain range of the given dataset [0,domainRange] */
	public static double domainRange = 10000;

	public static double distThreshold = domainRange * clusterRadiusRate;
	/** distance function of CF-Tree */
	//public static final int distFunction = CFTree.D0_DIST;

	/** ============================= seperators ================ */
	/** seperator for items of every record in the index */
	public static final String sepStrForRecord = ",";
	public static final String sepStrForKeyValue = "\t";
	public static final String sepStrForIDDist = "|";
	public static final String sepSplitForIDDist = "\\|";

}
