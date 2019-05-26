package cellpruning.lof.pruning;

import java.util.HashMap;

import metricspace.MetricObject;
import metricspace.Record;

public class ComputerUb {
	
	/**
	 * create 2019.01.08LF
	 * */
	
	public static float PI=3.14f;
	
	public static void calKDEMax(MetricObject o_S,HashMap<Long, MetricObject> KnnPoints,int dim,float fixh){
		
		long[] KNN_moObjectsID = o_S.getPointPQ().getValueSet();
		float[] moDistToKNN = o_S.getPointPQ().getPrioritySet();
		float minDist=Float.MAX_VALUE;
		long id=-1;
		float kdeMax=0f;
		float h=0f;
		for(int i=0;i<moDistToKNN.length;i++){
			if(minDist>moDistToKNN[i]){
				minDist=moDistToKNN[i];
				id=KNN_moObjectsID[i];
			}
		}
		//System.out.println("ID  "+id+" juli "+minDist);
		MetricObject q=KnnPoints.get(id);
		h=fixh*q.getKdist();
		float a1=(float) (1/(Math.pow(h, dim)*Math.pow(2*Math.PI, dim/2)));
		double a2=-(Math.pow(minDist, 2)/(2*Math.pow(h, 2)));
		 //System.out.println("temp_kNNKey "+temp_kNNKey+"   a1*Math.pow(Math.E, a2)   "+a1*Math.pow(Math.E, a2)+" daikuan  "+h+" 1/hd a1  "+a1+" e�Ĳ��� "+a2+"  Math.pow(Math.E, -a2)  "+Math.pow(Math.E, a2)+"  "+Math.E);
		kdeMax=(float) (a1*Math.pow(Math.E, a2));
		if (Float.isNaN(kdeMax) || Float.isInfinite(kdeMax))
			kdeMax=0;
		o_S.setKmax(kdeMax);
	}
	
	
public static void calKDEMin(MetricObject o_S,HashMap<Long, MetricObject> KnnPoints,int dim,float fixh){
		
		long[] KNN_moObjectsID = o_S.getPointPQ().getValueSet();
		float[] moDistToKNN = o_S.getPointPQ().getPrioritySet();
		long id=KNN_moObjectsID[0];
		float kdeMin=0f;
		float h=0f;	
	//	System.out.println("ID  "+id+" juli "+minDist);
		MetricObject q=KnnPoints.get(id);
		
		h=fixh*q.getKdist();
		float a1=(float) (1/(Math.pow(h, dim)*Math.pow(2*Math.PI, dim/2)));
		double a2=-(Math.pow(o_S.getKdist(), 2)/(2*Math.pow(h, 2)));
		 //System.out.println("temp_kNNKey "+temp_kNNKey+"   a1*Math.pow(Math.E, a2)   "+a1*Math.pow(Math.E, a2)+" daikuan  "+h+" 1/hd a1  "+a1+" e�Ĳ��� "+a2+"  Math.pow(Math.E, -a2)  "+Math.pow(Math.E, a2)+"  "+Math.E);
		kdeMin=(float) (a1*Math.pow(Math.E, a2));
		if (Float.isNaN(kdeMin) || Float.isInfinite(kdeMin))
			kdeMin=0;
		o_S.setKmin(kdeMin);
		
		
			
	}
	
public static void calUB(MetricObject o_S,HashMap<Long, MetricObject> KnnPoints,int dim,float fixh,int k){
	float ub=0f;
	float  temp_os=0.0f;
	if(o_S.getKde()!=-1){
		temp_os=o_S.getKde();	
	}
	else if(o_S.getKmin()==-1){
		calKDEMin(o_S,KnnPoints,dim,fixh);
		temp_os=o_S.getKmin();
	}else{
		temp_os=o_S.getKmin();
	}
		
	
	if(temp_os==0){
		ComputerKOF.calKOF(o_S, KnnPoints, dim, fixh, k);
		ub=o_S.getKof();
	}else{
		long[] KNN_moObjectsID = o_S.getPointPQ().getValueSet();
		float[] moDistToKNN = o_S.getPointPQ().getPrioritySet();
		for(int i=0;i<KNN_moObjectsID.length;i++){
			float temp_kMax=0f;
			long temp_kNNKey = KNN_moObjectsID[i];
			//float temp_dist = moDistToKNN[i];
			//System.out.println("hjhj  "+temp_kNNKey+"  "+((Record)o_S.getObj()).getRId());
			if(KnnPoints.get(temp_kNNKey).getKde()!=-1){
				temp_kMax=KnnPoints.get(temp_kNNKey).getKde();
				if (temp_kMax== 0)
					continue;
				else
				ub+=temp_kMax/temp_os;
			}else{
				float kMax=KnnPoints.get(temp_kNNKey).getKmax();
				if(kMax==-1){
					calKDEMax(KnnPoints.get(temp_kNNKey), KnnPoints,dim,fixh);
					temp_kMax=KnnPoints.get(temp_kNNKey).getKmax();
					if (temp_kMax== 0)
						continue;
					else{
						
						ub+=temp_kMax/temp_os;
					}
				}else{
					temp_kMax=KnnPoints.get(temp_kNNKey).getKmax();
					if (temp_kMax== 0)
						continue;
					else{
						
						ub+=temp_kMax/temp_os;
					}
				}
			}
		
		}
		
	}
	ub=ub/k;
	if (Float.isNaN(ub) || Float.isInfinite(ub)){
		ComputerKOF.calKOF(o_S, KnnPoints, dim, fixh, k);
		ub=o_S.getKof();
	}
		
	o_S.setUb(ub);
}



}
