package metricspace;

import java.util.HashMap;
import java.util.Map;

import util.PriorityQueue;

@SuppressWarnings("rawtypes")
public class MetricObject {

	private Object obj;
	// private Map<Long,Float> knnInDetail = new HashMap<Long,Float>();
	public PriorityQueue pointPQ = new PriorityQueue(PriorityQueue.SORT_ORDER_DESCENDING);
	public PriorityQueue pointPQR = new PriorityQueue(PriorityQueue.SORT_ORDER_DESCENDING);
	// private Map<Long, coreInfoKNNs> knnMoreDetail = new HashMap<Long,
	// coreInfoKNNs>();
	private float kdist = -1;
	private float lrdValue = -1;
	private float lofValue = -1;
	private char type = 'F';
    private float kde=-1;
    private float kof=-1;
    private float rKOF=-1;
    private float kmax=-1;
    private float kmin=-1;
    private float ub=-1;
    private HashMap<Long, MetricObject> RnnList=new HashMap<>();
    private HashMap<Long, Float> KnnList=new HashMap<>();
    
    
	private float nearestNeighborDist = Float.MAX_VALUE;
	private boolean canPrune = false;
	private int[] indexForSmallCell;
	private int indexOfCPCellInList = -1; // index of which cell it is in
	private float largeCellExpand = 0.0f;
	private boolean insideKNNfind = false;
	private int cellIndexPoint=-1;
	public MetricObject(Object obj) {
		this.obj = obj;
	}

	
	
	public float getrKOF() {
		return rKOF;
	}
	public void setrKOF(float rKOF) {
		this.rKOF = rKOF;
	}



	public void insert(MetricObject o_S){
		this.RnnList.put(((Record)o_S.getObj()).getRId(), o_S);
	}
	public void delete(MetricObject o_S){
		this.RnnList.remove(((Record)o_S.getObj()).getRId());
	}
	
	public void insertKNN(long id,float dist){
		this.KnnList.put(id, dist);
	}
	public void deleteKNN(long id){
		
		this.KnnList.remove(id);
	}
	public float getKde() {
		return kde;
	}

	public void setKde(float kde) {
		this.kde = kde;
	}

	public float getKof() {
		return kof;
	}

	public void setKof(float kof) {
		this.kof = kof;
	}

	public float getKmax() {
		return kmax;
	}

	public void setKmax(float kmax) {
		this.kmax = kmax;
	}



	public HashMap<Long, Float> getKnnList() {
		return KnnList;
	}

	public void setKnnList(HashMap<Long, Float> knnList) {
		KnnList = knnList;
	}

	public float getKmin() {
		return kmin;
	}

	public void setKmin(float kmin) {
		this.kmin = kmin;
	}

	public float getUb() {
		return ub;
	}

	public void setUb(float ub) {
		this.ub = ub;
	}

	public void setRnnList(HashMap<Long, MetricObject> rnnList) {
		RnnList = rnnList;
	}

	public HashMap<Long, MetricObject> getRnnList() {
		return RnnList;
	}

	public float getNearestNeighborDist() {
		return nearestNeighborDist;
	}

	public void setNearestNeighborDist(float nearestNeighborDist) {
		this.nearestNeighborDist = nearestNeighborDist;
	}

	public char getType() {
		return type;
	}

	public void setType(char type) {
		this.type = type;
	}

	public float getLrdValue() {
		return lrdValue;
	}

	public void setLrdValue(float lrdValue) {
		this.lrdValue = lrdValue;
	}

	public float getLofValue() {
		return lofValue;
	}

	public void setLofValue(float lofValue) {
		this.lofValue = lofValue;
	}

	public int[] getIndexForSmallCell() {
		return indexForSmallCell;
	}

	public void setIndexForSmallCell(int[] indexForSmallCell) {
		this.indexForSmallCell = indexForSmallCell;
	}

	public boolean isCanPrune() {
		return canPrune;
	}

	public void setCanPrune(boolean canPrune) {
		this.canPrune = canPrune;
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		Record r = (Record) obj;
		sb.append(", Knn in detail: ");
		return sb.toString();
	}

	public Object getObj() {
		return obj;
	}

	public void setObj(Object obj) {
		this.obj = obj;
	}

	public float getKdist() {
		return kdist;
	}

	public void setKdist(float kdist) {
		this.kdist = kdist;
	}

	@SuppressWarnings("unchecked")
	public static void main(String[] args) {

	}

	public int getIndexOfCPCellInList() {
		return indexOfCPCellInList;
	}

	public void setIndexOfCPCellInList(int indexOfCPCellInList) {
		this.indexOfCPCellInList = indexOfCPCellInList;
	}

	public float getLargeCellExpand() {
		return largeCellExpand;
	}

	public void setLargeCellExpand(float largeCellExpand) {
		this.largeCellExpand = largeCellExpand;
	}

	public PriorityQueue getPointPQ() {
		return pointPQ;
	}

	public void setPointPQ(PriorityQueue pointPQ) {
		this.pointPQ = pointPQ;
	}

	public boolean isInsideKNNfind() {
		return insideKNNfind;
	}

	public void setInsideKNNfind(boolean insideKNNfind) {
		this.insideKNNfind = insideKNNfind;
	}

	public PriorityQueue getPointPQR() {
		return pointPQR;
	}

	public void setPointPQR(PriorityQueue pointPQR) {
		this.pointPQR = pointPQR;
	}



	public int getCellIndexPoint() {
		return cellIndexPoint;
	}



	public void setCellIndexPoint(int cellIndexPoint) {
		this.cellIndexPoint = cellIndexPoint;
	}
	
	
}