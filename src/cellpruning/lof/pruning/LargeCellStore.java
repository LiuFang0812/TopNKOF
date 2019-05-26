package cellpruning.lof.pruning;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Stack;

import cellpruning.lof.pruning.firstknn.prQuadTree.*;
import metricspace.*;
import util.LinkQueue;
import util.PriorityQueue;
import util.SQConfig;
import util.LinkQueue.ListNode;

public class LargeCellStore extends partitionTreeNode {

	/**
	 * Coordinates of the large cell, save independent dims, for example,
	 * dims[0,1,2,3], independentDims [2,3] then save coordinates for dims [2,3]
	 */
	private float[] coordinates;

	/** list of points in the Large grid cell */
	private ArrayList<MetricObject> listOfPoints;

	/** Closest pair distance inside the large cell */
	private float cpDist;

	/** Number of points inside the Large cell */
	private int numOfPoints;

	/** Break up into small cells, only break up those independent dims */
	private float smallCellSize;

	/**
	 * Keep track of number of small cells per dim, only for independent dims
	 */
	private int[] numSmallCells;

	private IMetric metric = null;

	private IMetricSpace metricSpace = null;

	private boolean breakIntoSmallCells = false;

	private prQuadInternal rootForPRTree;
	// save leaves that cannot be pruned
	ArrayList<prQuadLeaf> prLeaves;

	/** priority of the bucket, used for sorting */
	private double bucketPriority = 1;
	
	private int cellIndex=-1;

	public LargeCellStore(float[] coordinates, ArrayList<MetricObject> listOfPoints, float cpDist, IMetric metric,
			IMetricSpace metricspace) {
		this.coordinates = coordinates.clone();
		this.listOfPoints = listOfPoints;
		this.numOfPoints = this.listOfPoints.size();
		this.cpDist = cpDist;
		// this.setIndexForLeaveNodesList(indexForLeaveNodesList);
		this.metric = metric;
		this.metricSpace = metricspace;
	}
	
	

	public int getCellIndex() {
		return cellIndex;
	}



	public void setCellIndex(int cellIndex) {
		this.cellIndex = cellIndex;
	}



	public int[] getNumSmallCells() {
		return numSmallCells;
	}

	public void setNumSmallCells(int[] numSmallCells) {
		this.numSmallCells = numSmallCells.clone();
	}

	public void addPoints(MetricObject newpoint) {
		listOfPoints.add(newpoint);
		numOfPoints++;
	}

	public ArrayList<MetricObject> getListOfPoints() {
		return listOfPoints;
	}

	public int getTotalNumOfDim() {
		return coordinates.length / 2;
	}

	/**
	 * compute priority for each large cell
	 * 
	 * @param isOnBoundary
	 */
	public void computePriorityForLargeCell(boolean isOnBoundary) {
		if (isOnBoundary) {
			this.bucketPriority = 0;
			return;
		} else
			this.bucketPriority = this.bucketPriority * Math.log(this.numOfPoints);
	}


	public  float computerTopNandCt(LinkQueue<LargeCellStore> queueWin, float thresholdLof,HashMap<Long, MetricObject> KnnPoints,int dim,float fixh,int k,int topNNumber,PriorityQueue topnKDE,HashMap<Long,Float> c_topnKDE){
		ListNode a2=queueWin.getFront();
		while(a2!=null){
			LargeCellStore curSli=(LargeCellStore)a2.getE();
			//curSli.findRnnsInit(queueWin);
			ArrayList<MetricObject> curL=curSli.getListOfPoints();
			for(int i=0;i<curL.size();i++){
				MetricObject curP=curL.get(i);
				
				if(topnKDE.size() < topNNumber){
					ComputerKOF.calKOF(curP, KnnPoints, dim, fixh, k);
					float tempKofValue=curP.getKof();
					topnKDE.insert(metricSpace.getID(curP.getObj()), tempKofValue);
				}else{
					thresholdLof=topnKDE.getPriority();
					ComputerUb.calUB(curP, KnnPoints, dim, fixh, k);
					if(curP.getUb()>thresholdLof){
						ComputerKOF.calKOF(curP, KnnPoints, dim, fixh, k);
						float tempKofValue=curP.getKof();
						if(tempKofValue>thresholdLof){
							if (topnKDE.size() < topNNumber) {
								topnKDE.insert(metricSpace.getID(curP.getObj()), tempKofValue);
								
							} else if (tempKofValue > topnKDE.getPriority()) {
								topnKDE.pop();
								topnKDE.insert(metricSpace.getID(curP.getObj()), tempKofValue);
								if (thresholdLof < topnKDE.getPriority())
									thresholdLof = topnKDE.getPriority();
								
							}else
								c_topnKDE.put(metricSpace.getID(curP.getObj()), tempKofValue);
						}else
							c_topnKDE.put(metricSpace.getID(curP.getObj()), tempKofValue);
					}
				}
			
			}
			a2=a2.getNext();
		}//end while
		return thresholdLof;
	}//end function
	
	public float multiplyArray(float[] array) {
		float result = 1;
		for (float element : array) {
			result *= element;
		}
		return result;
	}

	// partitionTreeNode ptn, float[] partition_store,
	public void seperateToSmallCells(HashMap<Long, MetricObject> CanPrunePoints, int indexOfLeaveNodes, float threshold,
			int K, int[] independentDims, int[] correlatedDims, int totalDim) {
		// save independent dim range
		float[] LargeCellRange = new float[independentDims.length];
		for (int i = 0; i < independentDims.length; i++) {
			LargeCellRange[i] = coordinates[i * 2 + 1] - coordinates[i * 2];
			//System.out.println("LargeCellRange[i]  "+LargeCellRange[i]);
		}

		// smallCellSize using closest pair and threshold
		smallCellSize = (float) (threshold * this.cpDist / (2 * Math.sqrt(totalDim)));

		// cell calculated as average cell size
		float smallCellSize_predict = (float) Math.pow(multiplyArray(LargeCellRange) * K / numOfPoints,
				1.0 / LargeCellRange.length);

		// if the small cell size too small, then don't use this size to build
		// PRQuadTree
		if (smallCellSize < smallCellSize_predict / 10) {
			return;
		}
		if (smallCellSize > smallCellSize_predict / 5) {
			smallCellSize = smallCellSize_predict / 5;
		}

		for (float LargeCellRangeTemp : LargeCellRange) {
			if (smallCellSize >= LargeCellRangeTemp)
				return;
		}
		numSmallCells = new int[independentDims.length];
		// calculate how many small cells for each partition per dimension
		for (int i = 0; i < LargeCellRange.length; i++) {
			numSmallCells[i] = (int) Math.ceil(LargeCellRange[i] / smallCellSize);
		}

		for (int i = 0; i < LargeCellRange.length; i++) {
			if (numSmallCells[i] < 10)
				return;
		}

		breakIntoSmallCells = true;

		for (MetricObject mo : listOfPoints) {
			Record record = (Record) mo.getObj();
			int[] indexInSmall = new int[independentDims.length];
			for (int i = 0; i < independentDims.length; i++) {
				float tempValue = record.getValue()[independentDims[i]];
				indexInSmall[i] = (int) (Math.floor(tempValue - coordinates[2 * i])
						/ (smallCellSize + Float.MIN_VALUE));
			}
			mo.setIndexForSmallCell(indexInSmall);
		}
		// build up PR quadtree
		float[] largeCellCoor = coordinates.clone();
		buildPRQuadTree(CanPrunePoints, numSmallCells, smallCellSize, listOfPoints, numOfPoints, largeCellCoor,
				indexOfLeaveNodes, independentDims, correlatedDims, K, true);
	}

	public boolean QuerySurroundingBucketsForCP(partitionTreeNode ptn, float[] expectedRange) {
		Stack<partitionTreeInternal> partitionTree = new Stack<partitionTreeInternal>();
		partitionTree.push((partitionTreeInternal) ptn);
		while (!partitionTree.isEmpty()) {
			partitionTreeInternal tempInternal = partitionTree.pop();
			ArrayList<partitionTreeNode> tempChildNodes = tempInternal.getChildNodes();
			// check children
			for (int i = 0; i < tempChildNodes.size(); i++) {
				if (tempChildNodes.get(i).getClass().getName().endsWith("partitionTreeInternal")) {
					float[] tempCoordinates = ((partitionTreeInternal) tempChildNodes.get(i)).getCoordinates();
					if (checkRange(expectedRange, tempCoordinates)) {
						partitionTree.push((partitionTreeInternal) tempChildNodes.get(i));
					}
				} else if (tempChildNodes.get(i).getClass().getName().endsWith("LargeCellStore")) {
					float[] tempCoordinates = ((LargeCellStore) tempChildNodes.get(i)).getCoordinates();
					if (checkRange(expectedRange, tempCoordinates)
							&& ((LargeCellStore) tempChildNodes.get(i)).getCpDist() < this.cpDist) {
						return false;
					}
				} else {
					System.out.println("Unknown Bucket Node Type!");
				}
			}
		}
		return true;
	}
/**
 * 2019.01.03 LF
 * create rtree
 * indexOfLeaveNodes not use
 * @param K
 * @param indexOfLeaveNodes
 * @param independentDims
 * @param correlatedDims
 */
	public void createRtree(int K, int indexOfLeaveNodes, int[] independentDims, int[] correlatedDims) {
		float[] LargeCellRange = new float[independentDims.length];
		for (int i = 0; i < independentDims.length; i++) {
			LargeCellRange[i] = coordinates[i * 2 + 1] - coordinates[i * 2];
		}

		smallCellSize = (float) Math.pow(multiplyArray(LargeCellRange) * 4 * K / numOfPoints,
				1.0 / LargeCellRange.length);
		for (float LargeCellRangeTemp : LargeCellRange) {
			if (smallCellSize >= LargeCellRangeTemp)
				return;
		}
		
		numSmallCells = new int[independentDims.length];
		// calculate how many small cells for each partition per dimension
		for (int i = 0; i < LargeCellRange.length; i++) {
			numSmallCells[i] = (int) Math.ceil(LargeCellRange[i] / smallCellSize);
		}
		breakIntoSmallCells = true;
		for (MetricObject mo : listOfPoints) {
			Record record = (Record) mo.getObj();
			int[] indexInSmall = new int[independentDims.length];
			for (int i = 0; i < independentDims.length; i++) {
				float tempValue = record.getValue()[independentDims[i]];
				indexInSmall[i] = (int) (Math.floor(tempValue - coordinates[2 * i])
						/ (smallCellSize + Float.MIN_VALUE));
			}
			mo.setIndexForSmallCell(indexInSmall);
		}
		// build up PR quadtree
		float[] largeCellCoor = coordinates.clone();
		buildPRQuadTree(null, numSmallCells, smallCellSize, listOfPoints, numOfPoints, largeCellCoor, indexOfLeaveNodes,
				independentDims, correlatedDims, K, false);
	}
	public void seperateLargeNoPrune(int K, int indexOfLeaveNodes, int[] independentDims, int[] correlatedDims) {
		float[] LargeCellRange = new float[independentDims.length];
		for (int i = 0; i < independentDims.length; i++) {
			LargeCellRange[i] = coordinates[i * 2 + 1] - coordinates[i * 2];	
		}
	
		smallCellSize = (float) Math.pow(multiplyArray(LargeCellRange) * 4 * K / numOfPoints,
				1.0 / LargeCellRange.length);
		for (float LargeCellRangeTemp : LargeCellRange) {
			if (smallCellSize >= LargeCellRangeTemp)
				return;
		}
		numSmallCells = new int[independentDims.length];
		// calculate how many small cells for each partition per dimension
		for (int i = 0; i < LargeCellRange.length; i++) {
			numSmallCells[i] = (int) Math.ceil(LargeCellRange[i] / smallCellSize);
		}
		//System.out.println("gbh"+numSmallCells[0]+"  "+numSmallCells[1]);
		breakIntoSmallCells = true;
		for (MetricObject mo : listOfPoints) {
			Record record = (Record) mo.getObj();
			int[] indexInSmall = new int[independentDims.length];
			for (int i = 0; i < independentDims.length; i++) {
				float tempValue = record.getValue()[independentDims[i]];
				indexInSmall[i] = (int) (Math.floor(tempValue - coordinates[2 * i])
						/ (smallCellSize + Float.MIN_VALUE));
			}
			mo.setIndexForSmallCell(indexInSmall);
		}
		//System.out.println("锟竭筹拷 "+smallCellSize);
		// build up PR quadtree
		float[] largeCellCoor = coordinates.clone();
		
		//System.out.println("hj "+coordinates[0]+"  "+coordinates[1]+"  "+coordinates[2]+"  "+coordinates[3]);
		buildPRQuadTree(null, numSmallCells, smallCellSize, listOfPoints, numOfPoints, largeCellCoor, indexOfLeaveNodes,
				independentDims, correlatedDims, K, false);
	}

	/**
	 * build up PR quad tree with information from the large cell
	 * 
	 * @param numSmallCells
	 * @param smallCellSize
	 * @param listOfPoints
	 * @param numOfPoints
	 * @param largeCellCoor
	 * @return root of PR QuadTree
	 */
	public void buildPRQuadTree(HashMap<Long, MetricObject> CanPrunePoints, int[] numSmallCells, float smallCellSize,
			ArrayList<MetricObject> listOfPoints, int numOfPoints, float[] largeCellCoor, int indexOfLeaveNodes,
			int[] independentDims, int[] correlatedDims, int K, boolean withPrune) {
		
		// init root
		int[] indexRangeInSmallCell = new int[numSmallCells.length * 2];
		for (int i = 0; i < numSmallCells.length; i++) {
			indexRangeInSmallCell[2 * i] = 0;
			indexRangeInSmallCell[2 * i + 1] = numSmallCells[i] - 1;
		}

		rootForPRTree = new prQuadInternal(largeCellCoor, indexRangeInSmallCell, null, numSmallCells, smallCellSize);
		Stack<prQuadInternal> prQuadTree = new Stack<prQuadInternal>();
		HashMap<prQuadInternal, ArrayList<MetricObject>> mapQuadInternalWithPoints = new HashMap<prQuadInternal, ArrayList<MetricObject>>();
		mapQuadInternalWithPoints.put(rootForPRTree, listOfPoints);
		// save leaves in the pr tree
		prLeaves = new ArrayList<prQuadLeaf>();
		prQuadTree.push(rootForPRTree);
		int count = 0;
		while (!prQuadTree.empty()) {
			/**
			 * pop up the quad node and divide into 4 parts check if each part
			 * contains enough points if contains K+1 points, create a
			 * prQuadInternal and push to stack if contains less than K points,
			 * create a prQuadLeaf and save the pointer if can not divide (reach
			 * minimum size), create prQuadLeaf and save the pointer
			 */
		
			prQuadInternal curPRNode = prQuadTree.pop();
		//System.out.println("curPRNode  "+curPRNode.getCoordinates()[0]+"  "+curPRNode.getCoordinates()[1]+"  "+curPRNode.getCoordinates()[2]+"  "+curPRNode.getCoordinates()[3]);
			curPRNode.generateChilden(CanPrunePoints, curPRNode, prQuadTree, mapQuadInternalWithPoints, prLeaves,
					largeCellCoor, numSmallCells, indexOfLeaveNodes, independentDims, correlatedDims, K, withPrune);
		}
	}

	public void traverseLargeCell(MetricObject curPoint, LargeCellStore large_cell_store, int K) {

	
		float dist = 0.0f;
		float theta;
		if (curPoint.pointPQ.size() > 0)
			theta = curPoint.pointPQ.getPriority();
		else
			theta = Float.POSITIVE_INFINITY;
		for (int i = 0; i < large_cell_store.getNumOfPoints(); i++) {
			MetricObject o_S = large_cell_store.getListOfPoints().get(i);
			//System.out.println("((Record) o_S.getObj()).getRId()  "+((Record) o_S.getObj()).getRId());
			if (((Record) o_S.getObj()).getRId() == ((Record) curPoint.getObj()).getRId()) {
				continue;
			} else if (o_S.getType() == 'C')
				continue;
			else {
				try {
					dist = metric.dist(curPoint.getObj(), o_S.getObj());
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				if (curPoint.pointPQ.size() < K) {
					curPoint.insertKNN(metricSpace.getID(o_S.getObj()), dist);
					curPoint.pointPQ.insert(metricSpace.getID(o_S.getObj()), dist);
					theta = curPoint.pointPQ.getPriority();
				} else if (dist < theta) {
//					if(((Record)curPoint.getObj()).getRId()==11){
//						System.out.println("metricSpace.getID(o_S.getObj()) "+metricSpace.getID(o_S.getObj())+" "+dist);
//					}
					curPoint.deleteKNN(curPoint.pointPQ.getValue());
					curPoint.pointPQ.pop();
					curPoint.insertKNN(metricSpace.getID(o_S.getObj()), dist);
					curPoint.pointPQ.insert(metricSpace.getID(o_S.getObj()), dist);
					theta = curPoint.pointPQ.getPriority();
				}
			}
		}
	}

	

	
	public void traverseLargeCell_del(MetricObject curPoint, LargeCellStore large_cell_store, int K,HashMap<Long, MetricObject> KnnPoints) {

	
		float dist = 0.0f;
		float theta;
		if (curPoint.pointPQ.size() > 0)
			theta = curPoint.pointPQ.getPriority();
		else
			theta = Float.POSITIVE_INFINITY;
		for (int i = 0; i < large_cell_store.getNumOfPoints(); i++) {
			MetricObject o_S = large_cell_store.getListOfPoints().get(i);
			//System.out.println("((Record) o_S.getObj()).getRId()  "+((Record) o_S.getObj()).getRId());
			if (((Record) o_S.getObj()).getRId() == ((Record) curPoint.getObj()).getRId()||curPoint.getKnnList().containsKey(metricSpace.getID(o_S.getObj()))) {
				continue;
			} else if (o_S.getType() == 'C')
				continue;
			else {
				try {
					dist = metric.dist(curPoint.getObj(), o_S.getObj());
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				if (curPoint.pointPQ.size() < K) {
					if(curPoint.getKnnList().containsKey(metricSpace.getID(o_S.getObj()))) {
						
					}else{
						o_S.insert(curPoint);
						curPoint.insertKNN(metricSpace.getID(o_S.getObj()), dist);
						curPoint.pointPQ.insert(metricSpace.getID(o_S.getObj()), dist);
						theta = curPoint.pointPQ.getPriority();
					}
					
				} else if (dist < theta) {
					if(curPoint.getKnnList().containsKey(metricSpace.getID(o_S.getObj()))) {
						
					}else{
						KnnPoints.get(curPoint.pointPQ.getValue()).delete(curPoint);
						curPoint.deleteKNN(curPoint.pointPQ.getValue());
						curPoint.pointPQ.pop();	
						curPoint.insertKNN(metricSpace.getID(o_S.getObj()), dist);
						curPoint.pointPQ.insert(metricSpace.getID(o_S.getObj()), dist);
						theta = curPoint.pointPQ.getPriority();
						o_S.insert(curPoint);
					}
					
				}
			}
		}
	}
	public void findKnns(MetricObject curPoint, prQuadLeaf curCheckLeaf, int K) {
		// traverse points
		
		float dist = 0.0f;
		float theta;
		if (curPoint.pointPQ.size() > 0)
			theta = curPoint.pointPQ.getPriority();
		else
			theta = Float.POSITIVE_INFINITY;
		for (int i = 0; i < curCheckLeaf.getNumPoints(); i++) {
			MetricObject o_S = curCheckLeaf.getListOfPoints().get(i);
			if (((Record) o_S.getObj()).getRId() == ((Record) curPoint.getObj()).getRId()) {
				continue;
			} else if (o_S.getType() == 'C')
				continue;
			else {
				try {
					dist = metric.dist(curPoint.getObj(), o_S.getObj());
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				if (curPoint.pointPQ.size() < K) {
					curPoint.insertKNN(metricSpace.getID(o_S.getObj()), dist);
					curPoint.pointPQ.insert(metricSpace.getID(o_S.getObj()), dist);
					theta = curPoint.pointPQ.getPriority();
				} else if (dist < theta) {
					curPoint.deleteKNN(curPoint.pointPQ.getValue());
					curPoint.pointPQ.pop();
					curPoint.insertKNN(metricSpace.getID(o_S.getObj()), dist);
					curPoint.pointPQ.insert(metricSpace.getID(o_S.getObj()), dist);
					theta = curPoint.pointPQ.getPriority();
				}
			}
		}
	}
	/**
	 * 2019.01.08LF
	 * @param queueWin
	 */
	public void findRnnsInit(LinkQueue<LargeCellStore> queueWin){
			for(int i=0;i<this.getListOfPoints().size();i++){
				MetricObject b=this.getListOfPoints().get(i);
				findRnns(b,queueWin, SQConfig.K);
			}
		
	}
	public void findRnns(MetricObject curPoint, LinkQueue<LargeCellStore> queueWin, int K) {
		// traverse points
		
		float dist = 0.0f;
		
		ListNode aa1=queueWin.getFront();
		while(aa1!=null){
			ArrayList<MetricObject> sli=((LargeCellStore)aa1.getE()).getListOfPoints();
			
			for (int i = 0; i <sli.size(); i++) {
				MetricObject o_S = sli.get(i);
				if (((Record) o_S.getObj()).getRId() == ((Record) curPoint.getObj()).getRId()) {
					continue;
				} else if (o_S.getType() == 'C')
					continue;
				else {
					try {
						dist = metric.dist(curPoint.getObj(), o_S.getObj());
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					if(dist<=o_S.getKdist()){
						curPoint.insert(o_S);
					}
					
				}
			}
			
			aa1=aa1.getNext();
		}//end while
		
	}
	public void delUpdateKnnForPoint(MetricObject point){
		HashMap<Long,Float> kNN=point.getKnnList();
		point.pointPQ.clear();
		for(Map.Entry<Long, Float> entry:kNN.entrySet()){
			point.pointPQ.insert(entry.getKey(), entry.getValue());
		}
	}
	/**
	 * 2019.01.02LF
	 * update kNN and KDE and KOF or UB
	 *  slide expired
	 */
	public float deleteNewSlide(LargeCellStore needDelSlide,LinkQueue<LargeCellStore> queueWin,int K, int[] independentDims,
			float []partition_store,int num_dims,HashMap<Long, MetricObject> KnnPoints,float fixh,PriorityQueue topnKDE,
		    	HashMap<Long,Float> c_topn,HashMap<Long,Float> topnList,float thresholdLof,int topNNumber){
		
		ArrayList<MetricObject> delSlide=needDelSlide.getListOfPoints();
		
		for(int i=0;i<delSlide.size();i++){
			
			ArrayList<MetricObject> updateK=new ArrayList<>();
			MetricObject curDel=delSlide.get(i);				
			long []curDelKnn=curDel.pointPQ.getValueSet();
			for(int a=0;a<curDelKnn.length;a++){
				if(!delSlide.contains(KnnPoints.get(curDelKnn[a]))){			
					KnnPoints.get(curDelKnn[a]).delete(curDel);					
				}
			}
	
			HashMap<Long,MetricObject> rNN=curDel.getRnnList();
			for(Map.Entry<Long, MetricObject> entry:rNN.entrySet()){
				MetricObject curP=entry.getValue();				
				curP.deleteKNN(((Record)curDel.getObj()).getRId());
				if((!delSlide.contains(curP))&&(!updateK.contains(curP))){
					updateK.add(curP);
				}	
			}
			
			if(topnList.containsKey(((Record)curDel.getObj()).getRId())){		
				topnList.remove(((Record)curDel.getObj()).getRId());
			}
			if(c_topn.containsKey(((Record)curDel.getObj()).getRId())){
				c_topn.remove(((Record)curDel.getObj()).getRId());
			}
			 
			
			if(topnList.size()<SQConfig.TOPN){
				topnKDE.clear();
				for(Map.Entry<Long, Float> entry:topnList.entrySet()){
					topnKDE.insert(entry.getKey(), entry.getValue());
				}
				
				while(topnKDE.size()<SQConfig.TOPN){
					if(c_topn.size()>0){
						long id=publicFunction.maxValueForTopn(c_topn);
						topnKDE.insert(id,c_topn.get(id));
						topnList.put(id, c_topn.get(id));
						c_topn.remove(id);
					}else{
						break;
					}
				}/////
				if(topnKDE.size()>0){
					thresholdLof = topnKDE.getPriority();
				}
				
			}
			
			for(int j=0;j<updateK.size();j++){
				
				MetricObject p=updateK.get(j);
				//System.out.println("eee  "+((Record)p.getObj()).getRId());
				delUpdateKnnForPoint(p);
				//p.pointPQ.clear();
				findKnnsForOnePointOutsideSlide_del(p, needDelSlide, queueWin,partition_store, SQConfig.K,num_dims,independentDims,KnnPoints);
			
			}
			
			for(int j=0;j<updateK.size();j++){
				
				ArrayList<MetricObject> kofSet=new ArrayList<>();
				ArrayList<MetricObject> ubSet=new ArrayList<>();
				MetricObject p=updateK.get(j);
				if(p.getKof()!=-1){
					ComputerKOF.calKDE(p, KnnPoints, num_dims, fixh, K);
					kofSet.add(p);
					for(Map.Entry<Long, MetricObject> entry:p.getRnnList().entrySet()){
						MetricObject Rnn_p=entry.getValue();
						//System.out.println(((Record)Rnn_p.getObj()).getRId());
						if(!delSlide.contains(Rnn_p)){
							if(Rnn_p.getKof()!=-1){
								kofSet.add(Rnn_p);
							}else{
								ubSet.add(Rnn_p);
							}
						}
						
					}	
				}else{			
					ArrayList<MetricObject> ub_temp=new ArrayList<>();
					if(p.getKde()==-1){
						//System.out.println("else  ");
						ComputerUb.calKDEMax(p, KnnPoints, num_dims, fixh);
						ComputerUb.calKDEMin(p, KnnPoints, num_dims, fixh);
						ubSet.add(p);
						for(Map.Entry<Long, MetricObject> entry:p.getRnnList().entrySet()){
							//System.out.println("else12  ");
							//System.out.println(((Record)entry.getValue().getObj()).getRId());
						
							if(!delSlide.contains(entry.getValue())){
								if(!updateK.contains(entry.getValue())){
									ub_temp.add(entry.getValue());
								}else{
									ubSet.add(entry.getValue());
								}
									
							}
								
						}
					}else{//P.KDE!=-1
						//System.out.println("P.KDE!=-1");
						ComputerKOF.calKDE(p, KnnPoints, num_dims, fixh, K);
						ubSet.add(p);
						for(Map.Entry<Long, MetricObject> entry:p.getRnnList().entrySet()){
							//System.out.println(((Record)entry.getValue().getObj()).getRId());
							if(!delSlide.contains(entry.getValue())){
								if(entry.getValue().getKof()!=-1){
									kofSet.add(entry.getValue());
								}else{
									if(updateK.contains(entry.getValue())){
										ubSet.add(entry.getValue());
									}else{
										ub_temp.add(entry.getValue());
									}
											
								}
							}
						}
					}
					for(int m=0;m<ub_temp.size();m++){
						boolean flag=false;
						//'System.out.println("ub_temp  "+((Record)ub_temp.get(m).getObj()).getRId());
						long []id=ub_temp.get(m).pointPQ.getValueSet();
						for(int m1=0;m1<id.length;m1++){
						if(updateK.contains(KnnPoints.get(id[m1]))){
							flag=true;
							break;
							
						   }
						
						}
						if(!flag){	
							ubSet.add(ub_temp.get(m));
						}
						
					}
				}
				for(int m=0;m<kofSet.size();m++){
					MetricObject pp=kofSet.get(m);
					ComputerKOF.calKOF(pp, KnnPoints, num_dims, fixh, K);
					float tempKofValue=pp.getKof();	
					thresholdLof=determineTopn( pp,topnKDE,	c_topn,topnList, thresholdLof, topNNumber);
				}
				for(int mm=0;mm<ubSet.size();mm++){
					ComputerUb.calUB(ubSet.get(mm), KnnPoints, SQConfig.dims, fixh, K);
					//System.out.println("ub  "+((Record)ubSet.get(mm).getObj()).getRId());
					MetricObject p1=ubSet.get(mm);
					if(p1.getUb()>thresholdLof){
						ComputerKOF.calKOF(p1, KnnPoints, num_dims, fixh, K);
						float tempKofValue=p1.getKof();
						thresholdLof=determineTopn( p1,topnKDE,
						    	c_topn,topnList, thresholdLof, topNNumber);
					}
				}
			}
			
			
		}
		return thresholdLof;
		
	}
	
	/**
	 * 2019.01.15LF
	 * 
	 * @param needDelSlide
	 * @param queueWin
	 * @param K
	 * @param independentDims
	 * @param partition_store
	 * @param num_dims
	 * @param KnnPoints
	 * @param fixh
	 * @param topnKDE
	 * @param c_topn
	 * @param topnList
	 * @param thresholdLof
	 * @param topNNumber
	 * @return
	 */
	public float deleteNewSlide_LazyUpdate(LargeCellStore needDelSlide,LinkQueue<LargeCellStore> queueWin,int K, int[] independentDims,
			float []partition_store,int num_dims,HashMap<Long, MetricObject> KnnPoints,float fixh,PriorityQueue topnKDE,
		    	HashMap<Long,Float> c_topn,HashMap<Long,Float> topnList,float thresholdLof,int topNNumber){
		
		ArrayList<MetricObject> updateK=new ArrayList<>();
		ArrayList<MetricObject> delSlide=needDelSlide.getListOfPoints();
		ArrayList<MetricObject> kofSet=new ArrayList<>();
		ArrayList<MetricObject> ubSet=new ArrayList<>();
		ArrayList<MetricObject> ub_temp=new ArrayList<>();
	//	System.out.println("edede    "+((Record)needDelSlide.getListOfPoints().get(0).getObj()).getRId());
		for(int i=0;i<delSlide.size();i++){
			MetricObject curDel=delSlide.get(i);				
			long []curDelKnn=curDel.pointPQ.getValueSet();
			for(int a=0;a<curDelKnn.length;a++){
				if(!delSlide.contains(KnnPoints.get(curDelKnn[a]))){			
					KnnPoints.get(curDelKnn[a]).delete(curDel);					
				}
			}
			HashMap<Long,MetricObject> rNN=curDel.getRnnList();
			for(Map.Entry<Long, MetricObject> entry:rNN.entrySet()){
				MetricObject curP=entry.getValue();				
				curP.deleteKNN(((Record)curDel.getObj()).getRId());
				if((!delSlide.contains(curP))&&(!updateK.contains(curP))){
					updateK.add(curP);
				}	
			}
			if(topnList.containsKey(((Record)curDel.getObj()).getRId())){		
				topnList.remove(((Record)curDel.getObj()).getRId());
			}
			 if(c_topn.containsKey(((Record)curDel.getObj()).getRId())){
				c_topn.remove(((Record)curDel.getObj()).getRId());
			}
		}// end for delSlide
		
		if(topnList.size()<SQConfig.TOPN){
			topnKDE.clear();
			for(Map.Entry<Long, Float> entry:topnList.entrySet()){
				topnKDE.insert(entry.getKey(), entry.getValue());
			}
			
			while(topnKDE.size()<SQConfig.TOPN){
				if(c_topn.size()>0){
					long id=publicFunction.maxValueForTopn(c_topn);
					topnKDE.insert(id,c_topn.get(id));
					topnList.put(id, c_topn.get(id));
					c_topn.remove(id);
				}else{
					break;
				}
			}
			if(topnKDE.size()>0){
				thresholdLof = topnKDE.getPriority();
			}
			
		}
		
		for(int j=0;j<updateK.size();j++){
			
			MetricObject p=updateK.get(j);
			//System.out.println("eee  "+((Record)p.getObj()).getRId());
			delUpdateKnnForPoint(p);
			//p.pointPQ.clear();
			findKnnsForOnePointOutsideSlide_del(p, needDelSlide, queueWin,partition_store, SQConfig.K,num_dims,independentDims,KnnPoints);
		
		}
		///////////////////////////////////////////////////
	//	System.out.println("bv bvhj   "+updateK.size());
		
		for(int j=0;j<updateK.size();j++){
			
			MetricObject p=updateK.get(j);
			if(p.getKof()!=-1){
				ComputerKOF.calKDE(p, KnnPoints, num_dims, fixh, K);
				kofSet.add(p);
				for(Map.Entry<Long, MetricObject> entry:p.getRnnList().entrySet()){
					MetricObject Rnn_p=entry.getValue();
					//System.out.println(((Record)Rnn_p.getObj()).getRId());
					if(!delSlide.contains(Rnn_p)){
						if(Rnn_p.getKof()!=-1){
							if(!kofSet.contains(Rnn_p)){
								kofSet.add(Rnn_p);
							}
						}else{
							if(!ubSet.contains(Rnn_p)){
								ubSet.add(Rnn_p);
							}
						}
					}
					
				}	
			}else{			
				
				if(p.getKde()==-1){
					//System.out.println("else  ");
					ComputerUb.calKDEMax(p, KnnPoints, num_dims, fixh);
					ComputerUb.calKDEMin(p, KnnPoints, num_dims, fixh);
					ubSet.add(p);
					for(Map.Entry<Long, MetricObject> entry:p.getRnnList().entrySet()){
						//System.out.println("else12  ");
						//System.out.println(((Record)entry.getValue().getObj()).getRId());
					
						if(!delSlide.contains(entry.getValue())){
							if(!updateK.contains(entry.getValue())){
								if(!ub_temp.contains(entry.getValue())){
									ub_temp.add(entry.getValue());
								}
							}else{
								if(!ubSet.contains(entry.getValue())){
									ubSet.add(entry.getValue());
								}
							}
								
						}
							
					}
				}else{//P.KDE!=-1
					//System.out.println("P.KDE!=-1");
					ComputerKOF.calKDE(p, KnnPoints, num_dims, fixh, K);
					ubSet.add(p);
					for(Map.Entry<Long, MetricObject> entry:p.getRnnList().entrySet()){
						//System.out.println(((Record)entry.getValue().getObj()).getRId());
						if(!delSlide.contains(entry.getValue())){
							if(entry.getValue().getKof()!=-1){
								if(!kofSet.contains(entry.getValue())){
									kofSet.add(entry.getValue());
								}
							}else{
								if(updateK.contains(entry.getValue())){
									if(!ubSet.contains(entry.getValue())){
										ubSet.add(entry.getValue());
									}
								}else{
									if(!ub_temp.contains(entry.getValue())){
										ub_temp.add(entry.getValue());
									}
								}
										
							}
						}
					}
				}
			
			}
		}
		
		for(int m=0;m<ub_temp.size();m++){
			boolean flag=false;
			//'System.out.println("ub_temp  "+((Record)ub_temp.get(m).getObj()).getRId());
			long []id=ub_temp.get(m).pointPQ.getValueSet();
			for(int m1=0;m1<id.length;m1++){
			if(updateK.contains(KnnPoints.get(id[m1]))){
				flag=true;
				break;
				
			   }
			
			}
			if(!flag){	
				ubSet.add(ub_temp.get(m));
			}
			
		}
		for(int m=0;m<kofSet.size();m++){
			MetricObject pp=kofSet.get(m);
			ComputerKOF.calKOF(pp, KnnPoints, num_dims, fixh, K);
			//float tempKofValue=pp.getKof();	
			thresholdLof=determineTopn( pp,topnKDE,	c_topn,topnList, thresholdLof, topNNumber);
		}
		for(int mm=0;mm<ubSet.size();mm++){
			ComputerUb.calUB(ubSet.get(mm), KnnPoints, SQConfig.dims, fixh, K);
			//System.out.println("ub  "+((Record)ubSet.get(mm).getObj()).getRId());
			MetricObject p1=ubSet.get(mm);
			if(p1.getUb()>thresholdLof){
				ComputerKOF.calKOF(p1, KnnPoints, num_dims, fixh, K);
				float tempKofValue=p1.getKof();
				thresholdLof=determineTopn( p1,topnKDE,
				    	c_topn,topnList, thresholdLof, topNNumber);
			}
		}
		return thresholdLof;
		
	}
	
	
	
	
	/**
	 * 2019.01.15LF
	 * Determine if the object belongs to top_n or not
	 * @param pp
	 * @param topnKDE
	 * @param c_topn
	 * @param topnList
	 * @param thresholdLof
	 * @param topNNumber
	 */
	public float determineTopn(MetricObject pp,PriorityQueue topnKDE,
	    	HashMap<Long,Float> c_topn,HashMap<Long,Float> topnList,float thresholdLof,int topNNumber){
		float tempKofValue=pp.getKof();
		
		if(tempKofValue>thresholdLof){
			if(topnList.containsKey(metricSpace.getID(pp.getObj()))){
				
				topnList.put(metricSpace.getID(pp.getObj()), tempKofValue);
				topnKDE.clear();
				for(Map.Entry<Long, Float> entry:topnList.entrySet()){
					topnKDE.insert(entry.getKey(), entry.getValue());
				}	
				thresholdLof = topnKDE.getPriority();
				
				if(c_topn.size()>0){
					long id=publicFunction.maxValueForTopn(c_topn);
					//System.out.println("id  "+id);
					float c_dist=c_topn.get(id);
					if(thresholdLof<c_dist){
						topnList.remove(topnKDE.getValue());
						topnKDE.pop();
						topnKDE.insert(id, c_dist);
						topnList.put(id, c_dist);
					}
					
					thresholdLof = topnKDE.getPriority();
				}
				
			}else{
				if (topnKDE.size() < topNNumber) {
					topnList.put(metricSpace.getID(pp.getObj()), tempKofValue);
					topnKDE.insert(metricSpace.getID(pp.getObj()), tempKofValue);
					
				} else if (tempKofValue > topnKDE.getPriority()) {
					topnList.remove(topnKDE.getValue());
					topnKDE.pop();
					topnKDE.insert(metricSpace.getID(pp.getObj()), tempKofValue);
					topnList.put(metricSpace.getID(pp.getObj()), tempKofValue);
					
					if (thresholdLof < topnKDE.getPriority())
						thresholdLof = topnKDE.getPriority();
					
				}else{
					c_topn.put(metricSpace.getID(pp.getObj()), tempKofValue);
				}
					
					
			}
			
		}else
			c_topn.put(metricSpace.getID(pp.getObj()), tempKofValue);
		return thresholdLof;
	}
	/**
	 * 2019.01.02
	 * New slide arrives
	 */
	public float enterNewSlide(LargeCellStore slide,LinkQueue<LargeCellStore> queueWin,int K, int[] independentDims,
			float []partition_store,int num_dims,HashMap<Long, MetricObject> KnnPoints,float fixh,PriorityQueue topnKDE,
	    	HashMap<Long,Float> c_topn,HashMap<Long,Float> topnList,float thresholdLof,int topNNumber){
		ArrayList<MetricObject> newSlide=slide.getListOfPoints();
		float dist=0.0f;
		for(int i=0;i<newSlide.size();i++){
			MetricObject newP=newSlide.get(i);
			ArrayList<MetricObject> updateK=new ArrayList<>();
			ArrayList<MetricObject> kof=new ArrayList<>();
			ArrayList<MetricObject> tempUb=new ArrayList<>();
			ArrayList<MetricObject> ub=new ArrayList<>();
			ListNode a2=queueWin.getFront();
			while(a2!=null){
				ArrayList<MetricObject> curLi=new ArrayList<>();
				//curSli.findRnnsInit(queueWin)
				if((LargeCellStore)a2.getE()==slide){	
					
				}else{
					curLi=((LargeCellStore)a2.getE()).getListOfPoints();
					for(int j=0;j<curLi.size();j++){
						MetricObject po=curLi.get(j);
			
						try {
							dist = metric.dist(newP.getObj(), po.getObj());
						} catch (IOException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
						if(po.getKdist()>dist){
						//	System.out.println("po.pointPQ.getValue()   "+po.pointPQ.getValue()+"  po  "+((Record)po.getObj()).getRId());
							KnnPoints.get(po.pointPQ.getValue()).delete(po);
							po.deleteKNN(po.pointPQ.getValue());
							po.pointPQ.pop();
							po.insertKNN(metricSpace.getID(newP.getObj()), dist);
							po.pointPQ.insert(metricSpace.getID(newP.getObj()), dist);
//							theta = curPoint.pointPQ.getPriority();
							newP.insert(po);
							updateK.add(po);
						}
						
						
					}	
				}
				
				a2=a2.getNext();
			}
			for(int m=0;m<updateK.size();m++){
				
				MetricObject updatePo=updateK.get(m);
				if(updatePo.getKof()!=-1){
					//System.out.println("kof!=-1    "+((Record)updatePo.getObj()).getRId());
					ComputerKOF.calKDE(updatePo, KnnPoints, num_dims, fixh, K);
					kof.add(updatePo);
					for(Map.Entry<Long, MetricObject> entry:updatePo.getRnnList().entrySet()){
						MetricObject Rnn_p=entry.getValue();
						if(Rnn_p.getKof()!=-1){
							if(!kof.contains(Rnn_p)){
								kof.add(Rnn_p);
							}
							
						}else{
							if(updateK.contains(Rnn_p)&&!tempUb.contains(Rnn_p)){
								tempUb.add(Rnn_p);
							}else{
								if(!ub.contains(Rnn_p)){
									ub.add(Rnn_p);
								}
							}
							
						}
					}
					
				}else{
					if(updatePo.getKde()!=-1){
						ComputerKOF.calKDE(updatePo, KnnPoints, num_dims, fixh, K);
						tempUb.add(updatePo);
						for(Map.Entry<Long, MetricObject> entry:updatePo.getRnnList().entrySet()){
							MetricObject Rnn_p=entry.getValue();
							if(Rnn_p.getKof()!=-1){
								if(!kof.contains(Rnn_p)){
									kof.add(Rnn_p);
								}
								
							}else{
								if(updateK.contains(Rnn_p)&&!tempUb.contains(Rnn_p)){
									tempUb.add(Rnn_p);
								}else{
									if(!ub.contains(Rnn_p)){
										ub.add(Rnn_p);
									}
								}
							}
						}
						
					}else{
						ComputerUb.calKDEMax(updatePo, KnnPoints, num_dims, fixh);
						ComputerUb.calKDEMin(updatePo, KnnPoints, num_dims, fixh);
						tempUb.add(updatePo);
						for(Map.Entry<Long, MetricObject> entry:updatePo.getRnnList().entrySet()){
							MetricObject Rnn_p=entry.getValue();
							if(updateK.contains(Rnn_p)&&!tempUb.contains(Rnn_p)){
								tempUb.add(Rnn_p);
							}else{
								if(!ub.contains(Rnn_p)){
									ub.add(Rnn_p);
								}
							}
						}
						
					}
				}
			}//end updateK seach
//			
			for(int m=0;m<tempUb.size();m++){
				boolean flag=false;
				long []id=tempUb.get(m).pointPQ.getValueSet();
				for(int mm=0;mm<id.length;mm++){
					if(updateK.contains(KnnPoints.get(id[mm]))){
						flag=true;
						break;
					}
				}
				if(flag){
					ub.add(tempUb.get(m));
				}
			}
			
			for(int m=0;m<kof.size();m++){
				MetricObject pp=kof.get(m);
				ComputerKOF.calKOF(pp, KnnPoints, num_dims, fixh, K);
				//System.out.println("kof  "+((Record)kofSet.get(m).getObj()).getRId());
			
				thresholdLof=determineTopn( pp,topnKDE,
				    	c_topn,topnList, thresholdLof, topNNumber);
			}
			
			for(int mm=0;mm<ub.size();mm++){
				ComputerUb.calUB(ub.get(mm), KnnPoints, SQConfig.dims, fixh, K);
				MetricObject p1=ub.get(mm);
				if(p1.getUb()>thresholdLof){
					ComputerKOF.calKOF(p1, KnnPoints, num_dims, fixh, K);
					thresholdLof=determineTopn( p1,topnKDE,c_topn,topnList, thresholdLof, topNNumber);
				}
			}
			
			ComputerUb.calUB(newP, KnnPoints, SQConfig.dims, fixh, K);
			if(newP.getUb()>thresholdLof){
				ComputerKOF.calKOF(newP, KnnPoints, num_dims, fixh, K);
				float tempKofValue=newP.getKof();
				if(tempKofValue>thresholdLof){
					if (topnKDE.size() < topNNumber) {
						topnList.put(metricSpace.getID(newP.getObj()), tempKofValue);
						topnKDE.insert(metricSpace.getID(newP.getObj()), tempKofValue);
						
					} else if (tempKofValue > topnKDE.getPriority()) {
						topnList.remove(topnKDE.getValue());
						topnKDE.pop();
						topnKDE.insert(metricSpace.getID(newP.getObj()), tempKofValue);
						topnList.put(metricSpace.getID(newP.getObj()), tempKofValue);
						if (thresholdLof < topnKDE.getPriority())
							thresholdLof = topnKDE.getPriority();
						
					}else
						c_topn.put(metricSpace.getID(newP.getObj()), tempKofValue);
				}else
					c_topn.put(metricSpace.getID(newP.getObj()), tempKofValue);
			}
		}//end newSlide
		return thresholdLof;
		
	}
	
	/**
	 * 2019.01.05LF
	 * @param slide
	 * @param queueWin
	 * @param K
	 * @param independentDims
	 * @param partition_store
	 * @param num_dims
	 * @param KnnPoints
	 * @param fixh
	 * @param topnKDE
	 * @param c_topn
	 * @param topnList
	 * @param thresholdLof
	 * @param topNNumber
	 * @return
	 */
	public float enterNewSlide_LazyUpdate(LargeCellStore slide,LinkQueue<LargeCellStore> queueWin,int K, int[] independentDims,
			float []partition_store,int num_dims,HashMap<Long, MetricObject> KnnPoints,float fixh,PriorityQueue topnKDE,
	    	HashMap<Long,Float> c_topn,HashMap<Long,Float> topnList,float thresholdLof,int topNNumber){
	//	System.out.println("Lazy  ");
		ArrayList<MetricObject> newSlide=slide.getListOfPoints();
		float dist=0.0f;
		ArrayList<MetricObject> updateK=new ArrayList<>();
		ArrayList<MetricObject> kof=new ArrayList<>();
		ArrayList<MetricObject> tempUb=new ArrayList<>();
		ArrayList<MetricObject> ub=new ArrayList<>();
		for(int i=0;i<newSlide.size();i++){
			MetricObject newP=newSlide.get(i);
			ListNode a2=queueWin.getFront();
			while(a2!=null){
				ArrayList<MetricObject> curLi=new ArrayList<>();
				if((LargeCellStore)a2.getE()==slide){	
					
				}else{
					curLi=((LargeCellStore)a2.getE()).getListOfPoints();
					for(int j=0;j<curLi.size();j++){
						MetricObject po=curLi.get(j);
						try {
							dist = metric.dist(newP.getObj(), po.getObj());
						} catch (IOException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
						if(po.getKdist()>dist){
						//	System.out.println("po.pointPQ.getValue()   "+po.pointPQ.getValue()+"  po  "+((Record)po.getObj()).getRId());
							KnnPoints.get(po.pointPQ.getValue()).delete(po);
							po.deleteKNN(po.pointPQ.getValue());
							po.pointPQ.pop();
							po.insertKNN(metricSpace.getID(newP.getObj()), dist);
							po.pointPQ.insert(metricSpace.getID(newP.getObj()), dist);
//							theta = curPoint.pointPQ.getPriority();
							newP.insert(po);
							if(!updateK.contains(po)){
								updateK.add(po);
							}		
						}		
					}	
				}	
				a2=a2.getNext();
			}	
		}
		
		for(int m=0;m<updateK.size();m++){		
			MetricObject updatePo=updateK.get(m);
			if(updatePo.getKof()!=-1){
				//System.out.println("kof!=-1    "+((Record)updatePo.getObj()).getRId());
				ComputerKOF.calKDE(updatePo, KnnPoints, num_dims, fixh, K);
				kof.add(updatePo);
				for(Map.Entry<Long, MetricObject> entry:updatePo.getRnnList().entrySet()){
					MetricObject Rnn_p=entry.getValue();
					if(Rnn_p.getKof()!=-1){
						if(!kof.contains(Rnn_p)){
							kof.add(Rnn_p);
						}
					}else{
						if(updateK.contains(Rnn_p)){
							if(!tempUb.contains(Rnn_p)){
								tempUb.add(Rnn_p);
							}
						}else{
							if(!ub.contains(Rnn_p)){
								ub.add(Rnn_p);
							}
						}
						
					}
				}
				
			}else{
				if(updatePo.getKde()!=-1){
					ComputerKOF.calKDE(updatePo, KnnPoints, num_dims, fixh, K);
					tempUb.add(updatePo);
					for(Map.Entry<Long, MetricObject> entry:updatePo.getRnnList().entrySet()){
						MetricObject Rnn_p=entry.getValue();
						if(Rnn_p.getKof()!=-1){
							if(!kof.contains(Rnn_p)){
								kof.add(Rnn_p);
							}
							
						}else{
							if(updateK.contains(Rnn_p)){
								if(!tempUb.contains(Rnn_p)){
									tempUb.add(Rnn_p);
								}
							}else{
								if(!ub.contains(Rnn_p)){
									ub.add(Rnn_p);
								}
							}
						}
					}
					
				}else{
					ComputerUb.calKDEMax(updatePo, KnnPoints, num_dims, fixh);
					ComputerUb.calKDEMin(updatePo, KnnPoints, num_dims, fixh);
					tempUb.add(updatePo);
					for(Map.Entry<Long, MetricObject> entry:updatePo.getRnnList().entrySet()){
						MetricObject Rnn_p=entry.getValue();
						if(updateK.contains(Rnn_p)){
							if(!tempUb.contains(Rnn_p)){
								tempUb.add(Rnn_p);
							}
						}else{
							if(!ub.contains(Rnn_p)){
								ub.add(Rnn_p);
							}
						}
					}
					
				}
			}
		}//end updateK seach
		
		for(int m=0;m<tempUb.size();m++){
			boolean flag=false;
			long []id=tempUb.get(m).pointPQ.getValueSet();
			for(int mm=0;mm<id.length;mm++){
				if(updateK.contains(KnnPoints.get(id[mm]))){
					flag=true;
					break;
				}
			}
			if(flag){
				ub.add(tempUb.get(m));
			}
		}
		
		for(int m=0;m<kof.size();m++){
			MetricObject pp=kof.get(m);
			ComputerKOF.calKOF(pp, KnnPoints, num_dims, fixh, K);
			//System.out.println("kof  "+((Record)kofSet.get(m).getObj()).getRId());
		
			thresholdLof=determineTopn( pp,topnKDE,
			    	c_topn,topnList, thresholdLof, topNNumber);
		}
		
		for(int mm=0;mm<ub.size();mm++){
			ComputerUb.calUB(ub.get(mm), KnnPoints, SQConfig.dims, fixh, K);
			MetricObject p1=ub.get(mm);
			if(p1.getUb()>thresholdLof){
				ComputerKOF.calKOF(p1, KnnPoints, num_dims, fixh, K);
				thresholdLof=determineTopn( p1,topnKDE,c_topn,topnList, thresholdLof, topNNumber);
			}
		}
		
		for(int i=0;i<newSlide.size();i++){
			MetricObject newP=newSlide.get(i);
			ComputerUb.calUB(newP, KnnPoints, SQConfig.dims, fixh, K);
			if(newP.getUb()>thresholdLof){
				ComputerKOF.calKOF(newP, KnnPoints, num_dims, fixh, K);
				float tempKofValue=newP.getKof();
				if(tempKofValue>thresholdLof){
					if (topnKDE.size() < topNNumber) {
						topnList.put(metricSpace.getID(newP.getObj()), tempKofValue);
						topnKDE.insert(metricSpace.getID(newP.getObj()), tempKofValue);
						
					} else if (tempKofValue > topnKDE.getPriority()) {
						topnList.remove(topnKDE.getValue());
						topnKDE.pop();
						topnKDE.insert(metricSpace.getID(newP.getObj()), tempKofValue);
						topnList.put(metricSpace.getID(newP.getObj()), tempKofValue);
						if (thresholdLof < topnKDE.getPriority())
							thresholdLof = topnKDE.getPriority();
						
					}else
						c_topn.put(metricSpace.getID(newP.getObj()), tempKofValue);
				}else
					c_topn.put(metricSpace.getID(newP.getObj()), tempKofValue);
			}
		}//end newSlide
		return thresholdLof;
		
	}

	
	
	static boolean checkRange(float[] expectedRange, float[] checkedRange) {
		
		for (int i = 0; i < expectedRange.length / 2; i++) {
			if (expectedRange[2 * i] > checkedRange[2 * i + 1] || expectedRange[2 * i + 1]<checkedRange[2 * i]  ){
				//System.out.println("   "+expectedRange[2 * i]+"    "+expectedRange[2 * i + 1]+"   "+checkedRange[2 * i + 1]+"   "+checkedRange[2 * i]);
				return false;
				}
		}
		return true;
	}

	/**
	 * check if the extended area is in checked area
	 * 
	 * @param checkedArea
	 * @param extendedArea
	 * @return
	 */
	boolean insideCheckedArea(float[] checkedArea, float[] extendedArea) {
		for (int i = 0; i < checkedArea.length / 2; i++) {
			if (extendedArea[2 * i] < checkedArea[2 * i] || extendedArea[2 * i + 1] > checkedArea[2 * i + 1])
				return false;
		}
		return true;
	}

	/**
	 * @param curPoint
	 * @param curLeaf
	 * @param large_cell_store
	 * @param K
	 * @param independentDims
	 */
	public void findKnnsWithinOneCell(MetricObject curPoint, prQuadLeaf curLeaf, LargeCellStore large_cell_store, int K,
			int[] independentDims) {//LF
		
		float[] curPointCoor = ((Record) curPoint.getObj()).getValue();
		float kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
		Stack<prQuadNode> checkWithinOneTree = new Stack<prQuadNode>();
		prQuadNode tempCheckBreakNode = curLeaf;
		float[] extendArea = new float[independentDims.length * 2];
		float[] largeCellCoor = large_cell_store.coordinates;
		for (int i = 0; i < independentDims.length; i++) {
			extendArea[2 * i] = Math.max(largeCellCoor[2 * i], curPointCoor[independentDims[i]] - kdist);
			extendArea[2 * i + 1] = Math.min(largeCellCoor[2 * i + 1], curPointCoor[independentDims[i]] + kdist);
		}
		float[] checkedCoordinates = tempCheckBreakNode.getCoordinates();
		if (insideCheckedArea(checkedCoordinates, extendArea)) {
			return;
		}
		while (tempCheckBreakNode.getParentNode() != null) {
			// first add the parent node
			checkWithinOneTree.push(tempCheckBreakNode.getParentNode());
			// then add brothers
			for (prQuadNode brother : ((prQuadInternal) tempCheckBreakNode.getParentNode()).getChildNodes()) {
				if (brother.equals(tempCheckBreakNode)) /////////
					continue;
				else if (checkRange(extendArea, brother.getCoordinates())) {
					checkWithinOneTree.push(brother);
				} else
					continue;
			}
			// traverse the stack until has only one element(the parent) left,
			// all brother traversed
			while (checkWithinOneTree.size() > 1) {
				prQuadNode tempNode = checkWithinOneTree.pop();
				if (tempNode.getClass().getName().endsWith("prQuadLeaf")) {
					if (tempNode.equals(curLeaf)) {
						continue;
					}
					findKnns(curPoint, (prQuadLeaf) tempNode, K);
					kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
					// update extend area
					float[] newExtendArea = new float[independentDims.length * 2];
					for (int i = 0; i < independentDims.length; i++) {
						newExtendArea[2 * i] = Math.max(largeCellCoor[2 * i], curPointCoor[independentDims[i]] - kdist);
						newExtendArea[2 * i + 1] = Math.min(largeCellCoor[2 * i + 1],
								curPointCoor[independentDims[i]] + kdist);
					}
					extendArea = newExtendArea;
					if (insideCheckedArea(checkedCoordinates, extendArea)) {
						return;
					}
				} else {
					// add children
					for (prQuadNode children : ((prQuadInternal) tempNode).getChildNodes()) {
						if (children.equals(curLeaf))
							continue;
						else if (checkRange(extendArea, children.getCoordinates())) {
							checkWithinOneTree.push(children);
						} else
							continue;
					}
				} // end else
			} // end while
				// stack size == 1 only one parent left, all brothers traversed
			tempCheckBreakNode = checkWithinOneTree.pop(); // let parent be the
															// new node
			checkedCoordinates = tempCheckBreakNode.getCoordinates();
			if (insideCheckedArea(checkedCoordinates, extendArea)) {
				return;
			}
		} // end while
	}

	public void savePriorityQueueToKNN(MetricObject curPoint, int K) {
		if (curPoint.pointPQ.size() != K) {
			System.out.println("Less than K points in Priority Queue");
			// curPoint.setCanPrune(true);
			return;
		}
		curPoint.setKdist(curPoint.pointPQ.getPriority());
		

	}

	/**
	 * search the partitionTreeNode and find supporting large cells
	 * 
	 * @param ExtendArea
	 *            the coordinate of the point's extended area
	 * @param ptn
	 *            the root node of the partition tree (binary tree)
	 * @param currentCell
	 *            current LargeCellStore, supportingLargeCells will not include
	 *            current cell
	 * @return the supporting cells for the current point
	 */
	public ArrayList<LargeCellStore> searchSupportingLargeCells(float[] ExtendArea, partitionTreeNode ptn,
			LargeCellStore currentCell) {
		ArrayList<LargeCellStore> supportingLargeCells = new ArrayList<LargeCellStore>();
		Stack<partitionTreeInternal> stackOfInternals = new Stack<partitionTreeInternal>();
		if (ptn.getClass().getName().endsWith("LargeCellStore")) {
			return supportingLargeCells;
		} else if (ptn.getClass().getName().endsWith("partitionTreeInternal")) {
			stackOfInternals.push((partitionTreeInternal) ptn);
		}
		while (!stackOfInternals.isEmpty()) {
			// check the coordinates of each child node
			ArrayList<partitionTreeNode> tempChildNodes = stackOfInternals.pop().getChildNodes();
			for (int i = 0; i < tempChildNodes.size(); i++) {
				partitionTreeNode tempPTN = tempChildNodes.get(i);
				if (tempPTN.getClass().getName().endsWith("LargeCellStore")
						&& checkRange(ExtendArea, ((LargeCellStore) tempPTN).getCoordinates())
						&& !(tempPTN == currentCell)) {
					// if(tempPTN == currentCell){
					// System.out.println("This is exactly the same cell~ not a
					// support");
					// }
					supportingLargeCells.add((LargeCellStore) tempPTN);
				} else if (tempPTN.getClass().getName().endsWith("partitionTreeInternal")
						&& checkRange(ExtendArea, ((partitionTreeInternal) tempPTN).getCoordinates())) {
					stackOfInternals.push((partitionTreeInternal) tempPTN);
				}
			}
		}
		return supportingLargeCells;
	}

	public void findKnnsForOnePointInsideBucket(HashMap<Long, MetricObject> TrueKnnPoints, MetricObject curPoint,
			prQuadLeaf curLeaf, LargeCellStore currentLeafNode, int K, int[] independentDims) {
		curPoint.setInsideKNNfind(true);
		float kdist = Float.POSITIVE_INFINITY;

		// first find kNNs within the large cell and bound a partition area for
		// largeCell
		// first find kNNs within the leaf
		findKnns(curPoint, curLeaf, K);
		kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;

		// then find KNNs within the large cell
		findKnnsWithinOneCell(curPoint, curLeaf, currentLeafNode, K, independentDims);
		kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
		// System.out.println("old kdistance: " + kdist);
		// check if kNNs exceeds the large cell
		curPoint.setLargeCellExpand(calExtendDistance(currentLeafNode, curPoint, kdist, independentDims));
		// if not exceed the large cell, don't need to traverse other large
		// cells
		if (curPoint.getLargeCellExpand() <= 1e-9) {
			savePriorityQueueToKNN(curPoint, K);
			TrueKnnPoints.put(((Record) curPoint.getObj()).getRId(), curPoint);
		}
	}
	/**
	 * 2019.01.06LF
	 * Judging kNN in coincident slide
	 * @param TrueKnnPoints
	 * @param curPoint
	 * @param large_cell_store
	 * @param currentLeafNode
	 * @param ptn
	 * @param partition_store
	 * @param K
	 * @param num_dims
	 * @param independentDims
	 */
	public void findKnnsForOnePointOutsideSlide( MetricObject curPoint, LargeCellStore currentLeafNode, LinkQueue<LargeCellStore> queueWin,
			float[] partition_store, int K, int num_dims, int[] independentDims) {//LF12019.01.06
		
		float kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;

		float[] curPointCoor = ((Record) curPoint.getObj()).getValue();

		float[] ExtendArea = new float[independentDims.length * 2];
		for (int i = 0; i < independentDims.length; i++) {
			
			ExtendArea[2 * i] = (float) Math.max(partition_store[2 * i], curPointCoor[independentDims[i]] - kdist);
			ExtendArea[2 * i + 1] = (float) Math.min(partition_store[2 * i + 1],
					curPointCoor[independentDims[i]] + kdist);
		}
			ArrayList<LargeCellStore> supportingLargeCells =calExtendDistanceSli(currentLeafNode,queueWin, curPoint, kdist,independentDims);
			
		for (LargeCellStore supportingCell : supportingLargeCells) {
	
			if (checkRange(ExtendArea, supportingCell.getCoordinates())) {
				
				// find a leaf to start
				if (supportingCell.breakIntoSmallCells) {
					
					findOverLapTreeNodesAndSearch(supportingCell.getRootForPRTree(), ExtendArea, curPoint, K);
					
					kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
					
					float[] newExtendArea = new float[independentDims.length * 2];
					for (int i = 0; i < independentDims.length; i++) {
						newExtendArea[2 * i] = (float) Math.max(partition_store[2 * i],
								curPointCoor[independentDims[i]] - kdist);
						newExtendArea[2 * i + 1] = (float) Math.min(partition_store[2 * i + 1],
								curPointCoor[independentDims[i]] + kdist);
					}
					ExtendArea = newExtendArea;
				} // end if
				else { // traverse Large cell
					
					traverseLargeCell(curPoint, supportingCell, K);
					kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
					float[] newExtendArea = new float[independentDims.length * 2];
					for (int i = 0; i < independentDims.length; i++) {
						newExtendArea[2 * i] = (float) Math.max(partition_store[2 * i],
								curPointCoor[independentDims[i]] - kdist);
						newExtendArea[2 * i + 1] = (float) Math.min(partition_store[2 * i + 1],
								curPointCoor[independentDims[i]] + kdist);
					}
					ExtendArea = newExtendArea;
				}
			} // end if(checkRange(ExtendArea, supportingCell.getCoordinates()))
		} // end for
			// bound supporting area for the partition
		savePriorityQueueToKNN(curPoint, K);
		//TrueKnnPoints.put(((Record) curPoint.getObj()).getRId(), curPoint);
	}

	public void findKnnsForOnePointOutsideSlide_add( MetricObject curPoint, LargeCellStore currentLeafNode, LinkQueue<LargeCellStore> queueWin,
			float[] partition_store, int K, int num_dims, int[] independentDims,HashMap<Long, MetricObject> KnnPoints) {//LF12019.01.06
		
		float kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;

        
		// if exceed the large cell, traverse nearby large cells
		// include more supporting cells
		float[] curPointCoor = ((Record) curPoint.getObj()).getValue();

		float[] ExtendArea = new float[independentDims.length * 2];
		for (int i = 0; i < independentDims.length; i++) { 
			ExtendArea[2 * i] = (float) Math.max(partition_store[2 * i], curPointCoor[independentDims[i]] - kdist);
			ExtendArea[2 * i + 1] = (float) Math.min(partition_store[2 * i + 1],
					curPointCoor[independentDims[i]] + kdist);
		}
			ArrayList<LargeCellStore> supportingLargeCells =calExtendDistanceSli(currentLeafNode,queueWin, curPoint, kdist,independentDims);

	
		
		for (LargeCellStore supportingCell : supportingLargeCells) {
			//System.out.println("888");
//			if(((Record)curPoint.getObj()).getRId()==4001){
//				System.out.println("gbhj"+checkRange(ExtendArea, supportingCell.getCoordinates()));
//			}
//			System.out.println("supportingCell"+supportingCell.getCoordinates()[0]+"  "+supportingCell.getCoordinates()[1]+" "
//					+supportingCell.getCoordinates()[2]+"  "+supportingCell.getCoordinates()[3]);
			//System.out.println("torf "+checkRange(ExtendArea, supportingCell.getCoordinates()));
			if ((curPoint.getCellIndexPoint()!=supportingCell.getCellIndex())&&checkRange(ExtendArea, supportingCell.getCoordinates())) {
//		/ System.out.println("addKnn");
				// find a leaf to start
				if (supportingCell.breakIntoSmallCells) {
					//System.out.println("bhj");
					findOverLapTreeNodesAndSearch_del(supportingCell.getRootForPRTree(), ExtendArea, curPoint, K,KnnPoints);
					
					kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
					
					float[] newExtendArea = new float[independentDims.length * 2];
					for (int i = 0; i < independentDims.length; i++) {
						newExtendArea[2 * i] = (float) Math.max(partition_store[2 * i],
								curPointCoor[independentDims[i]] - kdist);
						newExtendArea[2 * i + 1] = (float) Math.min(partition_store[2 * i + 1],
								curPointCoor[independentDims[i]] + kdist);
					}
					ExtendArea = newExtendArea;
				} // end if
				else { // traverse Large cell
					//System.out.println("fbvh");
					traverseLargeCell_del(curPoint, supportingCell, K,KnnPoints);
					kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
					float[] newExtendArea = new float[independentDims.length * 2];
					for (int i = 0; i < independentDims.length; i++) {
						newExtendArea[2 * i] = (float) Math.max(partition_store[2 * i],
								curPointCoor[independentDims[i]] - kdist);
						newExtendArea[2 * i + 1] = (float) Math.min(partition_store[2 * i + 1],
								curPointCoor[independentDims[i]] + kdist);
					}
					ExtendArea = newExtendArea;
				}
			} // end if(checkRange(ExtendArea, supportingCell.getCoordinates()))
		} // end for
			// bound supporting area for the partition
		savePriorityQueueToKNN(curPoint, K);
		//TrueKnnPoints.put(((Record) curPoint.getObj()).getRId(), curPoint);
	}
	
	public void findKnnsForOnePointOutsideSlide_del( MetricObject curPoint, LargeCellStore currentLeafNode, LinkQueue<LargeCellStore> queueWin,
			float[] partition_store, int K, int num_dims, int[] independentDims,HashMap<Long, MetricObject> KnnPoints) {//LF12019.01.06
		//System.out.println("14568锟斤拷锟斤拷"+"  "+partition_store[0]+"  "+partition_store[1]+"  "+partition_store[2]+" "+partition_store[3]);
		float kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
//       
//        	if(((Record)curPoint.getObj()).getRId()==11){
//        		System.out.println("kdist "+kdist+"  "+"qu  "+curPoint.getKdist());
//        	}
        
		// if exceed the large cell, traverse nearby large cells
		// include more supporting cells
		float[] curPointCoor = ((Record) curPoint.getObj()).getValue();
		float[] ExtendArea = new float[independentDims.length * 2];
		for (int i = 0; i < independentDims.length; i++) {
			
			ExtendArea[2 * i] = (float) Math.max(partition_store[2 * i], curPointCoor[independentDims[i]] - kdist);
			ExtendArea[2 * i + 1] = (float) Math.min(partition_store[2 * i + 1],
					curPointCoor[independentDims[i]] + kdist);
		}
			ArrayList<LargeCellStore> supportingLargeCells =calExtendDistanceSli(currentLeafNode,queueWin, curPoint, kdist,independentDims);

		for (LargeCellStore supportingCell : supportingLargeCells) {

			if (checkRange(ExtendArea, supportingCell.getCoordinates())) {
				
				// find a leaf to start
				if (supportingCell.breakIntoSmallCells) {
					findOverLapTreeNodesAndSearch_del(supportingCell.getRootForPRTree(), ExtendArea, curPoint, K, KnnPoints);
					
					kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
					
					float[] newExtendArea = new float[independentDims.length * 2];
					for (int i = 0; i < independentDims.length; i++) {
						newExtendArea[2 * i] = (float) Math.max(partition_store[2 * i],
								curPointCoor[independentDims[i]] - kdist);
						newExtendArea[2 * i + 1] = (float) Math.min(partition_store[2 * i + 1],
								curPointCoor[independentDims[i]] + kdist);
					}
					ExtendArea = newExtendArea;
				} // end if
				else { // traverse Large cell
					//System.out.println("2019");
					traverseLargeCell_del(curPoint, supportingCell, K,KnnPoints);
					kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
					float[] newExtendArea = new float[independentDims.length * 2];
					for (int i = 0; i < independentDims.length; i++) {
						newExtendArea[2 * i] = (float) Math.max(partition_store[2 * i],
								curPointCoor[independentDims[i]] - kdist);
						newExtendArea[2 * i + 1] = (float) Math.min(partition_store[2 * i + 1],
								curPointCoor[independentDims[i]] + kdist);
					}
					ExtendArea = newExtendArea;
				}
			} // end if(checkRange(ExtendArea, supportingCell.getCoordinates()))
		} // end for
			// bound supporting area for the partition
		savePriorityQueueToKNN(curPoint, K);
		//TrueKnnPoints.put(((Record) curPoint.getObj()).getRId(), curPoint);
	}

	
	
	
	
	
	// prQuadLeaf curLeaf
	public void findKnnsForOnePointOutsideBucket(HashMap<Long, MetricObject> TrueKnnPoints, MetricObject curPoint,
			ArrayList<LargeCellStore> large_cell_store, LargeCellStore currentLeafNode, partitionTreeNode ptn,
			float[] partition_store, int K, int num_dims, int[] independentDims) {//LF12019.01.06
		
		float kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;

		// if exceed the large cell, traverse nearby large cells
		// include more supporting cells
		float[] curPointCoor = ((Record) curPoint.getObj()).getValue();
		float[] ExtendArea = new float[independentDims.length * 2];
		for (int i = 0; i < independentDims.length; i++) {
			
			ExtendArea[2 * i] = (float) Math.max(partition_store[2 * i], curPointCoor[independentDims[i]] - kdist);
			ExtendArea[2 * i + 1] = (float) Math.min(partition_store[2 * i + 1],
					curPointCoor[independentDims[i]] + kdist);
		}
		
		ArrayList<LargeCellStore> supportingLargeCells = searchSupportingLargeCells(ExtendArea, ptn, currentLeafNode);
		// System.out.println("Size of Support:" + supportingLargeCells.size());
		// for each supporting cell, traverse until not exceed the checked area
		for (LargeCellStore supportingCell : supportingLargeCells) {
			if (checkRange(ExtendArea, supportingCell.getCoordinates())) {
				// find a leaf to start
				if (supportingCell.breakIntoSmallCells) {
					findOverLapTreeNodesAndSearch(supportingCell.getRootForPRTree(), ExtendArea, curPoint, K);
					// prQuadLeaf tempLeaf =
					// RangeQuery(supportingCell.getRootForPRTree(),
					// ExtendArea);
					// if (tempLeaf != null) {
					// // then find KNNs within the large cell
					// findKnnsWithinOneCell(curPoint, tempLeaf, supportingCell,
					// K, independentDims);
					// kdist = curPoint.pointPQ.size() == K ?
					// curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
					// float[] newExtendArea = new float[independentDims.length
					// * 2];
					// for (int i = 0; i < independentDims.length; i++) {
					// newExtendArea[2 * i] = (float) Math.max(partition_store[2
					// * i],
					// curPointCoor[independentDims[i]] - kdist);
					// newExtendArea[2 * i + 1] = (float)
					// Math.min(partition_store[2 * i + 1],
					// curPointCoor[independentDims[i]] + kdist);
					// }
					// ExtendArea = newExtendArea;
					//
					// }
//					traverseLargeCell(curPoint, supportingCell, K);
					kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
					float[] newExtendArea = new float[independentDims.length * 2];
					for (int i = 0; i < independentDims.length; i++) {
						newExtendArea[2 * i] = (float) Math.max(partition_store[2 * i],
								curPointCoor[independentDims[i]] - kdist);
						newExtendArea[2 * i + 1] = (float) Math.min(partition_store[2 * i + 1],
								curPointCoor[independentDims[i]] + kdist);
					}
					ExtendArea = newExtendArea;
				} // end if
				else { // traverse Large cell
					traverseLargeCell(curPoint, supportingCell, K);
					kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
					float[] newExtendArea = new float[independentDims.length * 2];
					for (int i = 0; i < independentDims.length; i++) {
						newExtendArea[2 * i] = (float) Math.max(partition_store[2 * i],
								curPointCoor[independentDims[i]] - kdist);
						newExtendArea[2 * i + 1] = (float) Math.min(partition_store[2 * i + 1],
								curPointCoor[independentDims[i]] + kdist);
					}
					ExtendArea = newExtendArea;
				}
			} // end if(checkRange(ExtendArea, supportingCell.getCoordinates()))
		} // end for
			// bound supporting area for the partition
		savePriorityQueueToKNN(curPoint, K);
		TrueKnnPoints.put(((Record) curPoint.getObj()).getRId(), curPoint);
	}

	public void findOverLapTreeNodesAndSearch(prQuadInternal prRoot, float[] expectedRange, MetricObject curPoint,
			int K) {
//		if(((Record)curPoint.getObj()).getRId()==4001){
//			System.out.println("expectedRange  "+expectedRange[0]+"  "+expectedRange[1]+"  "+expectedRange[2]+"  "+expectedRange[3]);
//		}
		ArrayList<prQuadLeaf> prQuadLeafList = new ArrayList<prQuadLeaf>();
		Stack<prQuadInternal> prQuadTreeInternal = new Stack<prQuadInternal>();
		prQuadTreeInternal.push(prRoot);
		
		//if(((Record)curPoint.getObj()).getRId()==4001){
			while (!prQuadTreeInternal.empty()) {
				prQuadInternal curPRNode = prQuadTreeInternal.pop();
				
				// traverse 4 childs and save to the stack if inside the expecting
				// range
//				if(((Record)curPoint.getObj()).getRId()==4001){
//					System.out.println("leaf12锟斤拷锟斤拷   "+curPRNode.getChildNodes().size());
//				}
				for (prQuadNode tempNode : curPRNode.getChildNodes()) {
					
					// check range
					if (!checkRange(expectedRange, tempNode.getCoordinates())){
						//System.out.println("cheack");
//						if(((Record)curPoint.getObj()).getRId()==4001){
//							System.out.println("leaf12   "+tempNode.getCoordinates()[0]+"   "+tempNode.getCoordinates()[1]+"  "+tempNode.getCoordinates()[2]+"  "+tempNode.getCoordinates()[3]);
//						}
						continue;
					}
						
					if (tempNode.getClass().getName().endsWith("prQuadInternal")) {
						prQuadTreeInternal.push((prQuadInternal) tempNode);
					} else { // leaf
						
						prQuadLeafList.add((prQuadLeaf) tempNode);
					}
				}
			}
		

		// traverse each leaf nodes to find kNN
		float dist = 0.0f;
		float theta;
		if (curPoint.pointPQ.size() == K){
			theta = curPoint.pointPQ.getPriority();
		}	
		else
			theta = Float.POSITIVE_INFINITY;
		for (prQuadLeaf currentLeaf : prQuadLeafList) {
			for (MetricObject o_S : currentLeaf.getListOfPoints()) {
				if (((Record) o_S.getObj()).getRId() == ((Record) curPoint.getObj()).getRId()) {
					continue;
				} else if (o_S.getType() == 'C')
					continue;
				else {
					try {
						dist = metric.dist(curPoint.getObj(), o_S.getObj());
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					if (curPoint.pointPQ.size() < K) {
						curPoint.insertKNN(metricSpace.getID(o_S.getObj()), dist);
						curPoint.pointPQ.insert(metricSpace.getID(o_S.getObj()), dist);
						theta = curPoint.pointPQ.getPriority();
					} else if (dist < theta) {
						//System.out.println(((Record)curPoint.getObj()).getRId()+"锟斤拷去k锟斤拷锟斤拷锟斤拷锟斤拷"+ "dist= "+dist+" theta "+theta);
						curPoint.deleteKNN(curPoint.pointPQ.getValue());
						curPoint.pointPQ.pop();
						curPoint.insertKNN(metricSpace.getID(o_S.getObj()), dist);
						curPoint.pointPQ.insert(metricSpace.getID(o_S.getObj()), dist);
						theta = curPoint.pointPQ.getPriority();
					}
				}
			}
		}
	}

	public void findOverLapTreeNodesAndSearch_del(prQuadInternal prRoot, float[] expectedRange, MetricObject curPoint,
			int K,HashMap<Long, MetricObject> KnnPoints) {

		ArrayList<prQuadLeaf> prQuadLeafList = new ArrayList<prQuadLeaf>();
		Stack<prQuadInternal> prQuadTreeInternal = new Stack<prQuadInternal>();
		prQuadTreeInternal.push(prRoot);
		
		//if(((Record)curPoint.getObj()).getRId()==4001){
			while (!prQuadTreeInternal.empty()) {
				prQuadInternal curPRNode = prQuadTreeInternal.pop();
				
				// traverse 4 childs and save to the stack if inside the expecting
				// range

				for (prQuadNode tempNode : curPRNode.getChildNodes()) {
					
					// check range
					if (!checkRange(expectedRange, tempNode.getCoordinates())){
						//System.out.println("cheack");

						continue;
					}
						
					if (tempNode.getClass().getName().endsWith("prQuadInternal")) {
						prQuadTreeInternal.push((prQuadInternal) tempNode);
					} else { // leaf
						
						prQuadLeafList.add((prQuadLeaf) tempNode);
					}
				}
			}
		
//			if(((Record)curPoint.getObj()).getRId()==4001){
//				System.out.println("leaf   "+prQuadLeafList.size());
//			}
		// traverse each leaf nodes to find kNN
		float dist = 0.0f;
		float theta;
		if (curPoint.pointPQ.size() == K){
			theta = curPoint.pointPQ.getPriority();
		}	
		else
			theta = Float.POSITIVE_INFINITY;
		for (prQuadLeaf currentLeaf : prQuadLeafList) {
			for (MetricObject o_S : currentLeaf.getListOfPoints()) {
				if (((Record) o_S.getObj()).getRId() == ((Record) curPoint.getObj()).getRId()) {
					continue;
				} else if (o_S.getType() == 'C')
					continue;
				else {
					try {
						dist = metric.dist(curPoint.getObj(), o_S.getObj());
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					if (curPoint.pointPQ.size() < K) {
						//System.out.println("metricSpace.getID(o_S.getObj())");
						o_S.insert(curPoint);
						curPoint.insertKNN(metricSpace.getID(o_S.getObj()), dist);
						curPoint.pointPQ.insert(metricSpace.getID(o_S.getObj()), dist);
						theta = curPoint.pointPQ.getPriority();
					} else if (dist < theta) {
					
						KnnPoints.get(curPoint.pointPQ.getValue()).delete(curPoint);
						curPoint.deleteKNN(curPoint.pointPQ.getValue());
						curPoint.pointPQ.pop();
						curPoint.insertKNN(metricSpace.getID(o_S.getObj()), dist);
						curPoint.pointPQ.insert(metricSpace.getID(o_S.getObj()), dist);
						theta = curPoint.pointPQ.getPriority();
						o_S.insert(curPoint);
					}
				}
			}
		}
	}
	
	public prQuadLeaf findLeafWithSmallCellIndex(prQuadInternal prRoot, int[] indexForSmallCells,
			int[] independentDims) {
		Stack<prQuadInternal> prQuadTreeInternal = new Stack<prQuadInternal>();
		prQuadTreeInternal.push(prRoot);
		// System.out.println(xx + "," + yy);
		while (!prQuadTreeInternal.empty()) {
			prQuadInternal curPRNode = prQuadTreeInternal.pop();
			// traverse childs and save to the stack if inside the expecting
			// range
			for (prQuadNode tempNode : curPRNode.getChildNodes()) {
				boolean tempFlag = true;
				int[] indexRangeOfSmallCell = tempNode.getIndexInSmallCell();
				for (int i = 0; i < independentDims.length; i++) {
					if (indexRangeOfSmallCell[2 * i] > indexForSmallCells[i]
							|| indexForSmallCells[i] > indexRangeOfSmallCell[2 * i + 1]) {
						tempFlag = false;
						break;
					}
				}
				if (tempFlag) {
					if (tempNode.getClass().getName().endsWith("prQuadInternal")) {
						prQuadTreeInternal.push((prQuadInternal) tempNode);
					} else { // leaf
						return (prQuadLeaf) tempNode;
					}
				} else
					continue;
			}
		}
		return null;
	}

	public prQuadLeaf RangeQuery(prQuadInternal prRoot, float[] expectedRange) {
		Stack<prQuadInternal> prQuadTreeInternal = new Stack<prQuadInternal>();
		prQuadTreeInternal.push(prRoot);
		while (!prQuadTreeInternal.empty()) {
			prQuadInternal curPRNode = prQuadTreeInternal.pop();
			// traverse 4 childs and save to the stack if inside the expecting
			// range
			for (prQuadNode tempNode : curPRNode.getChildNodes()) {
				// check range
				if (!checkRange(expectedRange, tempNode.getCoordinates()))
					continue;
				if (tempNode.getClass().getName().endsWith("prQuadInternal")) {
					prQuadTreeInternal.push((prQuadInternal) tempNode);
				} else { // leaf
					return (prQuadLeaf) tempNode;
				}
			}
		}
		return null;
	}
	/**
	 * 2019.01.03 LF 
	 * @param large_cell_store
	 * @param curPoint
	 * @param kdist
	 * @param independentDims
	 * @return
	 */
	public ArrayList<LargeCellStore> calExtendDistanceSli(LargeCellStore curentSlide,LinkQueue<LargeCellStore> queueWin, MetricObject curPoint, float kdist,
			int[] independentDims) {
		float extendDist = 0.0f;		
		ArrayList<LargeCellStore> needCalSlide=new ArrayList<>();
		float[] curPointCoor = ((Record) curPoint.getObj()).getValue();
		
		ListNode aa1=queueWin.getFront();
		while(aa1!=null){
			boolean flag=true;
			if(((LargeCellStore)aa1.getE())==curentSlide){//Judging coincidence with the deleted slide
				
			}else{
				
					float[] largeCellCoor = ((LargeCellStore)aa1.getE()).coordinates;
					for (int i = 0; i < independentDims.length; i++) {

						if((curPointCoor[independentDims[i]] - kdist)>largeCellCoor[2 * i+1]||
								largeCellCoor[2 * i]>(curPointCoor[independentDims[i]] + kdist)){
							flag=false;
							break;	
						}
					}
				
					if(flag){
						needCalSlide.add(((LargeCellStore)aa1.getE()));
					}
				}
				
			
			aa1=aa1.getNext();
		}//end while
		
		return needCalSlide;
	}

	//Determine whether the new slide creates rtree----Yes
	
	public void findKnnsForRtreeSlide(LargeCellStore slide,LinkQueue<LargeCellStore> queueWin,int[] independentDims,
			int K,float []partition_store,int num_dims){//int[] independentDims=SQConfig.indexOfIndependentDims
		for (prQuadLeaf curLeaf : slide.prLeaves) {
			for (MetricObject curPoint : curLeaf.getListOfPoints()) {
				
				findKnnsForOnePointRtreeInsideSlide(curPoint,curLeaf, this, SQConfig.K,independentDims);
				float kdist=curPoint.getKdist();
			
				findKnnsForOnePointOutsideSlide(curPoint, slide, queueWin,partition_store, SQConfig.K,num_dims,independentDims);
				
			}
		}
	}
	
		//Determine whether the new slide creates rtree----No
		public void findKnnsForLargeCellInsideSlide(LargeCellStore slide,LinkQueue<LargeCellStore> queueWin,
				int K, int[] independentDims,float []partition_store,int num_dims) {
			// find knns for each point in the large cell
			for (MetricObject curPoint : slide.getListOfPoints()) {
				
				//findKnnsForOnePointInLargeCellInsideBucket(TrueKnnPoints, curPoint, large_cell_store, K, independentDims);
				findKnnsForOnePointInLargeCellInsideSlide(curPoint,slide,  K, independentDims);
				//float kdist=curPoint.getKdist();
				float kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
					findKnnsForOnePointOutsideSlide(curPoint, slide, queueWin,partition_store, SQConfig.K,num_dims,independentDims);
				
			}
		}
			
		//Determine whether the new slide establishes rtree-----No, each new entry rNN finds it when it asks for kNN
		public void findKnnsForLargeCellInsideSlide_add(LargeCellStore slide,LinkQueue<LargeCellStore> queueWin,
				int K, int[] independentDims,float []partition_store,int num_dims,HashMap<Long, MetricObject> KnnPoints) {
			// find knns for each point in the large cell
			for (MetricObject curPoint : slide.getListOfPoints()) {
				
				
				//findKnnsForOnePointInLargeCellInsideBucket(TrueKnnPoints, curPoint, large_cell_store, K, independentDims);
				findKnnsForOnePointInLargeCellInsideSlide_add(curPoint,slide,  K, independentDims,KnnPoints);
				//float kdist=curPoint.getKdist();
				float kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
					findKnnsForOnePointOutsideSlide_add(curPoint, slide, queueWin,partition_store, SQConfig.K,
							num_dims,independentDims,KnnPoints);
				
				
				
			}
		}
		
		/**
		 * 2019.01.04 LF
		 * bulid rtree in slide
		 * @param curPoint
		 * @param curLeaf
		 * @param currentLeafNode
		 * @param K
		 * @param independentDims
		 */
		public void findKnnsForOnePointRtreeInsideSlide(MetricObject curPoint,
				prQuadLeaf curLeaf, LargeCellStore currentLeafNode, int K, int[] independentDims) {
			curPoint.setInsideKNNfind(true);
			float kdist = Float.POSITIVE_INFINITY;
			// first find kNNs within the large cell and bound a partition area for
			// largeCell
			// first find kNNs within the leaf
			findKnns(curPoint, curLeaf, K);
			kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;

			// then find KNNs within the large cell
			findKnnsWithinOneCell(curPoint, curLeaf, currentLeafNode, K, independentDims);
			kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
			// System.out.println("old kdistance: " + kdist);
			// check if kNNs exceeds the large cell
			//curPoint.setLargeCellExpand(calExtendDistance(currentLeafNode, curPoint, kdist, independentDims));
			curPoint.setKdist(kdist);
			// if not exceed the large cell, don't need to traverse other large
			// cells
//			if (curPoint.getLargeCellExpand() <= 1e-9) {
//				savePriorityQueueToKNN(curPoint, K);
//			}

			// if not exceed the large cell, don't need to traverse other large
			// cells
			
		}
		
		/**
		 * 2019.01.04 LF
		 * not bulid rtree in slide
		 * @param curPoint
		 * @param large_cell_store
		 * @param K
		 * @param independentDims
		 */
		public void findKnnsForOnePointInLargeCellInsideSlide(MetricObject curPoint, LargeCellStore large_cell_store, int K, int[] independentDims) {
			curPoint.setInsideKNNfind(true);
			float kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
			// first traverse the large cell
			traverseLargeCell(curPoint, large_cell_store, K);
			kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
			curPoint.setKdist(kdist);
		}
		
		public void findKnnsForOnePointInLargeCellInsideSlide_add(MetricObject curPoint, LargeCellStore large_cell_store, int K, int[] independentDims,HashMap<Long, MetricObject> KnnPoints) {
			curPoint.setInsideKNNfind(true);
			float kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
			// first traverse the large cell
			traverseLargeCell_del(curPoint, large_cell_store, K,KnnPoints);
			kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
			curPoint.setKdist(kdist);

		}
	
	public float calExtendDistance(LargeCellStore large_cell_store, MetricObject curPoint, float kdist,
			int[] independentDims) {
		float extendDist = 0.0f;
		float[] curPointCoor = ((Record) curPoint.getObj()).getValue();
		float[] largeCellCoor = large_cell_store.coordinates;
		for (int i = 0; i < independentDims.length; i++) {
			extendDist = (float) Math.max(extendDist,
					largeCellCoor[2 * i] - (curPointCoor[independentDims[i]] - kdist));
			extendDist = (float) Math.max(extendDist,
					(curPointCoor[independentDims[i]] + kdist) - largeCellCoor[2 * i + 1]);
		}
		return extendDist;
	}
	
	public void findKnnsWithinPRTreeInsideBucket(HashMap<Long, MetricObject> TempTrueKnnPoints,
			LargeCellStore large_cell_store, int K, int[] independentDims) {
		// find kNNs for each point in the prLeaves
		for (prQuadLeaf curLeaf : this.prLeaves) {
			for (MetricObject curPoint : curLeaf.getListOfPoints()) {
				findKnnsForOnePointInsideBucket(TempTrueKnnPoints, curPoint, curLeaf, large_cell_store, K,
						independentDims);
			}
		}
	}

	public void findKnnsWithinPRTreeOutsideBucket(HashMap<Long, MetricObject> TrueKnnPoints,
			ArrayList<LargeCellStore> large_cell_store, int xx, partitionTreeNode ptn, float[] partition_store, int K,
			int num_dims, int[] independentDims) {
		// find kNNs for each point in the prLeaves
		for (prQuadLeaf curLeaf : this.prLeaves) {
			for (MetricObject curPoint : curLeaf.getListOfPoints()) {
				if (curPoint.getType() == 'F' && !curPoint.isCanPrune())
					findKnnsForOnePointOutsideBucket(TrueKnnPoints, curPoint, large_cell_store,
							large_cell_store.get(xx), ptn, partition_store, K, num_dims, independentDims);
			}
		}
	}

	public void findKnnsForLargeCellInsideBucket(HashMap<Long, MetricObject> TrueKnnPoints,
			LargeCellStore large_cell_store, int K, int[] independentDims) {
		// find knns for each point in the large cell
		for (MetricObject curPoint : large_cell_store.getListOfPoints()) {
			findKnnsForOnePointInLargeCellInsideBucket(TrueKnnPoints, curPoint, large_cell_store, K, independentDims);
		}
	}

	public void findKnnsForLargeCellOutsideBucket(HashMap<Long, MetricObject> TrueKnnPoints,
			ArrayList<LargeCellStore> large_cell_store, int xx, partitionTreeNode ptn, float[] partition_store, int K,
			int num_dims, int[] independentDims) {
		// find knns for each point in the large cell
		for (MetricObject curPoint : large_cell_store.get(xx).getListOfPoints()) {
			if (curPoint.getType() == 'F')
				findKnnsForOnePointOutsideBucket(TrueKnnPoints, curPoint, large_cell_store, large_cell_store.get(xx),
						ptn, partition_store, K, num_dims, independentDims);
		}
		
	}

	public void findKnnsForOnePointInLargeCellInsideBucket(HashMap<Long, MetricObject> TrueKnnPoints,
			MetricObject curPoint, LargeCellStore large_cell_store, int K, int[] independentDims) {
		curPoint.setInsideKNNfind(true);
		float kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
		// first traverse the large cell
		traverseLargeCell(curPoint, large_cell_store, K);
		kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
		// check if kNNs exceeds the large cell
		curPoint.setLargeCellExpand(calExtendDistance(large_cell_store, curPoint, kdist, independentDims));
		if (curPoint.getLargeCellExpand() <= 1e-9) {
			savePriorityQueueToKNN(curPoint, K);
			TrueKnnPoints.put(((Record) curPoint.getObj()).getRId(), curPoint);
		}
	}

	
	public boolean isBreakIntoSmallCells() {
		return breakIntoSmallCells;
	}

	public void setBreakIntoSmallCells(boolean breakIntoSmallCells) {
		this.breakIntoSmallCells = breakIntoSmallCells;
	}

	/** print Large cell store */
	public String printCellStoreWithSupport() {
		String str = "";
		for (float x : coordinates) {
			str += x + ",";
		}
		str += "\n";
		str += "number of points:" + numOfPoints + "\n" + "closest pair distance: " + cpDist + "\n"
				+ "Points in detail: ";
		for (Iterator<MetricObject> itr = listOfPoints.iterator(); itr.hasNext();) {
			str = str + itr.next().getObj().toString() + "\n";

		}
		return str.substring(0, str.length());
	}

	public String printCellStoreDetailedInfo() {
		String str = "";
		for (float x : coordinates) {
			str += x + ",";
		}
		str += "\n";
		str += "number of points:" + numOfPoints + "\n" + "closest pair distance: " + cpDist + "\n";
		str += "isbreakup?" + breakIntoSmallCells + "\n";
		str += "small cell size = " + smallCellSize;
		return str.substring(0, str.length());
	}

	public float getCpDist() {
		return cpDist;
	}

	public void setCpDist(float cpDist) {
		this.cpDist = cpDist;
	}

	public float getSmallCellSize() {
		return smallCellSize;
	}

	public void setSmallCellSize(float smallCellSize) {
		this.smallCellSize = smallCellSize;
	}

	public IMetric getMetric() {
		return metric;
	}

	public void setMetric(IMetric metric) {
		this.metric = metric;
	}

	public int getNumOfPoints() {
		return numOfPoints;
	}

	public void setNumOfPoints(int numOfPoints) {
		this.numOfPoints = numOfPoints;
	}

	public float[] getCoordinates() {
		return this.coordinates;
	}

	public prQuadInternal getRootForPRTree() {
		return rootForPRTree;
	}

	public void setRootForPRTree(prQuadInternal rootForPRTree) {
		this.rootForPRTree = rootForPRTree;
	}

	public ArrayList<prQuadLeaf> getPrLeaves() {
		return prLeaves;
	}

	public void setPrLeaves(ArrayList<prQuadLeaf> prLeaves) {
		this.prLeaves = prLeaves;
	}

	public double getBucketPriority() {
		return bucketPriority;
	}

	public void setBucketPriority(double bucketPriority) {
		this.bucketPriority = bucketPriority;
	}

}
