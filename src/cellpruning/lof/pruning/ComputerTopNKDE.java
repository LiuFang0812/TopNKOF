package cellpruning.lof.pruning;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import metricspace.IMetric;
import metricspace.IMetricSpace;
import metricspace.MetricObject;
import metricspace.MetricSpaceUtility;
import metricspace.Record;
import util.LinkQueue;
import util.PriorityQueue;
import util.SQConfig;
import util.LinkQueue.ListNode;

public class ComputerTopNKDE {
	private ArrayList<MetricObject> pointList;
	private ArrayList<MetricObject> readOnceList=new ArrayList<>();
	//HashMap<Integer, LargeCellStore> hashWin = new HashMap<Integer,LargeCellStore>();
	public LinkQueue<LargeCellStore> queueWin=new LinkQueue<>();
	HashMap<Long, MetricObject> KnnPoints = new HashMap<Long, MetricObject>();//当前窗口的对象集合
//	ArrayList<LargeCellStore> hashWin = new ArrayList<LargeCellStore>();
	private IMetricSpace metricSpace = null;
	private IMetric metric = null;
	private float thresholdLof = 0.0f;
	private PriorityQueue topnKDE = new PriorityQueue(PriorityQueue.SORT_ORDER_ASCENDING);
	private HashMap<Long,Float> topnList=new HashMap<>();
	private HashMap<Long,Float> c_topnKDE=new HashMap<>();
	private HashMap<Long,Float> outlierList=new HashMap<>();
	int currentNum=0;//
	int N=0;//
	int end = 0;//end time of window
	public ComputerTopNKDE() throws IOException {
		
		this.readMetricAndMetricSpace();
		//this.thresholdLof = SQConfig.thresholdLof;
		this.ReadInputFile();
		//System.out.println("SQConfig.numWin/SQConfig.slide"+SQConfig.numWin/SQConfig.slide);
		for(int i=0;i<(SQConfig.numWin/SQConfig.slide);i++){
			
			readSlide();
		}		
	}
/**
 * 2019.01.01
 * Single slide
 */
	public void readSlide(){
		end++;
		float[] independentCoordinates = new float[SQConfig.num_independentDim * 2];
		ArrayList<MetricObject> slideList=new ArrayList<>();
		//int a=end%(SQConfig.numWin/SQConfig.slide);
		int endNumSlide=Math.min(N+SQConfig.slide, readOnceList.size());
			for(int i=N;i<endNumSlide;i++){
				MetricObject point=readOnceList.get(i);
				point.setCellIndexPoint(end);
				KnnPoints.put(((Record)point.getObj()).getRId(), point);
				slideList.add(point);
			}
			N+=slideList.size();
			//System.out.println("vbhbhv "+slideList.size());
			independentCoordinates=ClosestPair.slideRange(slideList, independentCoordinates, SQConfig.indexOfIndependentDims, metric, metricSpace);
			LargeCellStore slide=new LargeCellStore(independentCoordinates, slideList, 0,metric, metricSpace);
			slide.setCellIndex(end);
			//hashWin.put(a, slide);//findKnnsForRtreeSlide
			//slide.calExtendDistanceSli(hashWin, curPoint, kdist, independentDims)
			//slide.createRtree(SQConfig.K,0, SQConfig.indexOfIndependentDims, SQConfig.indexOfCorrelatedDims);
			//slide.seperateLargeNoPrune(SQConfig.K, a, SQConfig.indexOfIndependentDims, SQConfig.indexOfCorrelatedDims);
			queueWin.add(slide);
	}

	public void ReadInputFile() {
		//System.out.println("Start Reading Points");
		//pointList = new ArrayList<MetricObject>();
		BufferedReader in;
		int i=0;
	//	System.out.println();
		// long count = 0;
		try {
			in = new BufferedReader(new FileReader(SQConfig.dataset));
			String line = null;
			while ((line = in.readLine()) != null&&readOnceList.size()<=SQConfig.numOnce) {
				i++;
				//int aa=Math.min(20, currentNum+SQConfig.numOnce);
				if(i>currentNum&&i<=currentNum+SQConfig.numOnce){
				//	System.out.println("i "+i+"currentNum "+currentNum);
					Object currentPoint = metricSpace.readObject(line, SQConfig.dims);
					MetricObject mo = new MetricObject(currentPoint);
					this.readOnceList.add(mo);
					
				}	
				// count++;
				// if (count % 10000000 == 0)
				// System.out.println(readOnceList.size() + " Points Added!");
			}
			in.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		currentNum+=SQConfig.numOnce;
		System.out.println(currentNum + " Points Added!");
	}

	/**
	 * get MetricSpace and metric from configuration
	 * 
	 * @param conf
	 * @throws IOException
	 */
	public void readMetricAndMetricSpace() throws IOException {
		try {
			metricSpace = MetricSpaceUtility.getMetricSpace(SQConfig.strMetricSpace);
			metric = MetricSpaceUtility.getMetric(SQConfig.strMetric);
			metricSpace.setMetric(metric);
			
		} catch (InstantiationException e) {
			throw new IOException("InstantiationException");
		} catch (IllegalAccessException e) {
			e.printStackTrace();
			throw new IOException("IllegalAccessException");
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
			throw new IOException("ClassNotFoundException");
		}
	}

	public static float[] maxOfTwoFloatArray(float[] x, float[] y) {
		float[] newArray = new float[x.length];
		for (int i = 0; i < x.length; i++) {
			newArray[i] = Math.max(x[i], y[i]);
		}
		return newArray;
	}
	
	
	public void ComputerTopN(){
		long start=System.currentTimeMillis();
		long readTime=0;
		
		float[] independentCoordinates = new float[SQConfig.num_independentDim * 2];
		for (int i = 0; i < SQConfig.num_independentDim; i++) {
			independentCoordinates[i * 2] = SQConfig.domainSpace[SQConfig.indexOfIndependentDims[i] * 2];
			independentCoordinates[i * 2 + 1] = SQConfig.domainSpace[SQConfig.indexOfIndependentDims[i] * 2 + 1];
		}
		//初始化窗口信息
		
		ListNode a2=queueWin.getFront();
		LargeCellStore InitCell=(LargeCellStore)queueWin.getFront().getE();
		while(a2!=null){
			LargeCellStore curSli=(LargeCellStore)a2.getE();
			//curSli.findRnnsInit(queueWin);
			curSli.findKnnsForLargeCellInsideSlide_add(curSli, queueWin, SQConfig.K, SQConfig.indexOfIndependentDims, independentCoordinates,SQConfig.dims, KnnPoints);
			a2=a2.getNext();
			
		}
		
	
		thresholdLof=InitCell.computerTopNandCt(queueWin,thresholdLof,KnnPoints,SQConfig.dims ,SQConfig.h,SQConfig.K,SQConfig.TOPN,topnKDE,c_topnKDE);
		long []tempId=topnKDE.getValueSet();
		float []kof=topnKDE.getPrioritySet();
		for(int i=0;i<tempId.length;i++){
			topnList.put(tempId[i], kof[i]);
		}
		//System.out.println("topn       "+topnList.size());
		
		//System.out.println("id  "+((Record)KnnPoints.get((long)1).getObj()).getRId()+"   LOf   "+KnnPoints.get((long)1).getLofValue());
		
		int num_l=0;
		for(Map.Entry<Long, Float> entry:topnList.entrySet()){
			int label=((Record)KnnPoints.get(entry.getKey()).getObj()).getLabel();
			if(label==0){
				num_l++;
				if(!outlierList.containsKey(entry.getKey())){
					outlierList.put(entry.getKey(), entry.getValue());
				}
				
			}
		}
		
			
			while(queueWin.length()==(SQConfig.numWin/SQConfig.slide)){
			LargeCellStore needDel=(LargeCellStore)queueWin.getFront().getE();
			
			queueWin.poll();
			thresholdLof=needDel.deleteNewSlide(needDel, queueWin, SQConfig.K,SQConfig.indexOfIndependentDims, independentCoordinates, SQConfig.dims, KnnPoints,SQConfig.h,topnKDE,c_topnKDE,topnList,thresholdLof,SQConfig.TOPN);
		
			ArrayList<MetricObject> li=needDel.getListOfPoints();
			for(int i=0;i<li.size();i++){
				KnnPoints.remove(((Record)li.get(i).getObj()).getRId());
			}
			
			long tempRead=System.currentTimeMillis();//read data for start time
			if(N<readOnceList.size()){
				readSlide();
			}else{
				if(readOnceList.size()>=SQConfig.numOnce){
					N=0;
					readOnceList.clear();
					//System.out.println("read "+readOnceList.size());
					this.ReadInputFile();
					readSlide();
					//System.out.println("currentNum  "+currentNum);
				}else{
					break;
				}
			}
			readTime+=(System.currentTimeMillis()-tempRead);
			
			
			LargeCellStore ls=(LargeCellStore) queueWin.getRear().getE();
			//new slide for kNN
			ls.findKnnsForLargeCellInsideSlide_add(ls, queueWin, SQConfig.K, SQConfig.indexOfIndependentDims, independentCoordinates,SQConfig.dims, KnnPoints);
			
			//update ub or kof or ct
			thresholdLof=ls.enterNewSlide(ls, queueWin, SQConfig.K, SQConfig.indexOfIndependentDims, independentCoordinates, SQConfig.dims, KnnPoints, SQConfig.h,topnKDE,c_topnKDE,topnList,thresholdLof,SQConfig.TOPN);
			int num_label=0;
			for(Map.Entry<Long, Float> entry:topnList.entrySet()){
				int label=((Record)KnnPoints.get(entry.getKey()).getObj()).getLabel();
				if(label==0){
					num_label++;
					//System.out.println("iddd   "+entry.getKey());
					if(!outlierList.containsKey(entry.getKey())){
						outlierList.put(entry.getKey(), entry.getValue());
					}
				}
				
			}
		}//end while queueWin.length()==(SQConfig.numWin/SQConfig.slide)
		System.out.println("Time   "+(System.currentTimeMillis()-start-readTime)+"    "+(System.currentTimeMillis()-start-readTime)/1000+"seconds");
				
		BufferedWriter out;
		try {
			out = new BufferedWriter(new FileWriter(new File(SQConfig.outputFile)));
			while (topnKDE.size() > 0) {
				// System.out.println(topnLOF.getValue() + "," +
				// topnLOF.getPriority());
				out.write(topnKDE.getValue() + "," + topnKDE.getPriority());
				out.newLine();
				topnKDE.pop();
			}
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		
		
	}//end function
	

	public static void main(String[] args) throws InterruptedException {
		System.out.println("Set Arguments.......");

		Options options = new Options();

		Option paraK = new Option("k", true, "K for KNN Search");
		paraK.setRequired(true);
		options.addOption(paraK);

		Option paraTopN = new Option("n", true, "Top-N Number");
		paraTopN.setRequired(true);
		options.addOption(paraTopN);

		Option inputFilePath = new Option("i", true, "Input File Path");
		inputFilePath.setRequired(true);
		options.addOption(inputFilePath);

		Option outputFilePath = new Option("o", true, "Output File Path");
		outputFilePath.setRequired(true);
		options.addOption(outputFilePath);

		Option paraDim = new Option("d", true, "Dimensions");
		paraDim.setRequired(true);
		options.addOption(paraDim);

		Option paraDomain = new Option("r", true, "Domain Range");
		paraDomain.setRequired(true);
		options.addOption(paraDomain);

		Option paraThreshold = new Option("h", true, "Start Threshold");
		paraThreshold.setRequired(true);
		options.addOption(paraThreshold);

		CommandLineParser parser = new DefaultParser();
		HelpFormatter formatter = new HelpFormatter();
		CommandLine cmd;

		try {
			cmd = parser.parse(options, args);
		} catch (ParseException e) {
			System.out.println(e.getMessage());
			formatter.printHelp("utility-name", options);

			System.exit(1);
			return;
		}
		SQConfig.K = Integer.parseInt(cmd.getOptionValue("k"));
		SQConfig.TOPN = Integer.parseInt(cmd.getOptionValue("n"));
		SQConfig.dataset = cmd.getOptionValue("i");
		SQConfig.dims = Integer.parseInt(cmd.getOptionValue("d"));
		float tempDomainRange = Float.parseFloat(cmd.getOptionValue("r"));
		SQConfig.domains = new float[2];
		SQConfig.domains[0] = 0.0f;
		SQConfig.domains[1] = tempDomainRange;
		SQConfig.domainSpace = new float[2 * SQConfig.dims];
		for (int i = 0; i < SQConfig.dims; i++) {
			SQConfig.domainSpace[i * 2] = 0.0f;
			SQConfig.domainSpace[i * 2 + 1] = tempDomainRange;
		}

		SQConfig.thresholdLof = Float.parseFloat(cmd.getOptionValue("h"));
		SQConfig.outputFile = cmd.getOptionValue("o");

		System.out.println("K = " + SQConfig.K);
		System.out.println("Top-N = " + SQConfig.TOPN);
		System.out.println("Input File Path =  " + SQConfig.dataset);
		System.out.println("Output File Path = " + SQConfig.outputFile);
		System.out.println("Dim =  " + SQConfig.dims);
		System.out.println("Domain Range =  " + tempDomainRange);
		System.out.println("Start threshold = " + SQConfig.thresholdLof);
		
		long begin = System.currentTimeMillis();

		/*---------------start   KDE-------------------*/
		try {
			ComputerTopNKDE computeTopNLOFWithPruning = new ComputerTopNKDE();


		computeTopNLOFWithPruning.ComputerTopN();
		
	} catch (IOException e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
	} 
		
		
		
		
	}
}
