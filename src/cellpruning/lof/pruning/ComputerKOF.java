package cellpruning.lof.pruning;

import java.util.HashMap;

import metricspace.MetricObject;

public class ComputerKOF {
	public static float PI=3.14f;
	
	
	/**
	 * 2019.01.08LF
	 * @param o_S
	 * @param KnnPoints
	 * @param dim
	 * @param fixh
	 * @param k
	 */
	public static void calKDE(MetricObject o_S,HashMap<Long, MetricObject> KnnPoints,int dim,float fixh,int k){
		float kde=0.0f;
		long[] KNN_moObjectsID = o_S.getPointPQ().getValueSet();
		float[] moDistToKNN = o_S.getPointPQ().getPrioritySet();
		for(int i=0;i<KNN_moObjectsID.length;i++){
			long temp_kNNKey = KNN_moObjectsID[i];
			float temp_dist = moDistToKNN[i];
			float h=fixh*KnnPoints.get(temp_kNNKey).getKdist();
			float a1=(float) (1/Math.pow(h, dim));
			double a2=-(Math.pow(temp_dist, 2)/(2*Math.pow(h, 2)));
			 //System.out.println("temp_kNNKey "+temp_kNNKey+"   a1*Math.pow(Math.E, a2)   "+a1*Math.pow(Math.E, a2)+" daikuan  "+h+" 1/hd a1  "+a1+" eµÄ²ÎÊý "+a2+"  Math.pow(Math.E, -a2)  "+Math.pow(Math.E, a2)+"  "+Math.E);
			 kde+=(float) (a1*Math.pow(Math.E, a2));
//			if(temp_kNNKey==7){
//				System.out.println("ceshi "+(float) (a1*Math.pow(Math.E, a2))/Math.pow(2*PI, dim/2));
//			}
		}
		
		float a3=(float) (k*Math.pow(2*Math.PI, dim/2));
		kde=kde/a3;
//		if (Float.isNaN(kde) || Float.isInfinite(kde))
//			kde = 0;
		o_S.setKde(kde);
		
	}
	
	/**
	 * 2019.01.08LF
	 * @param o_S
	 * @param KnnPoints
	 * @param k
	 */
	public static void calKOF(MetricObject o_S,HashMap<Long, MetricObject> KnnPoints,int dim,float fixh,int k){
		float kof=0.0f;
		long[] KNN_moObjectsID = o_S.getPointPQ().getValueSet();
		float[] moDistToKNN = o_S.getPointPQ().getPrioritySet();
		if(o_S.getKde()==-1){
			calKDE(o_S, KnnPoints,dim,fixh,k);
		
		}
		if(o_S.getKde()==0){
			kof=0;
		}else{
			for(int i=0;i<KNN_moObjectsID.length;i++){
				long temp_kNNKey = KNN_moObjectsID[i];
				//float temp_dist = moDistToKNN[i];
				float kde=KnnPoints.get(temp_kNNKey).getKde();
				if(kde==-1){
					calKDE(KnnPoints.get(temp_kNNKey), KnnPoints,dim,fixh,k);
					float temp_kde=KnnPoints.get(temp_kNNKey).getKde();
					if (temp_kde== 0 || o_S.getLrdValue() == 0)
						continue;
					else
					kof+=temp_kde / o_S.getKde()*1.0f;
					//System.out.println("j123 "+temp_kNNKey+" kde "+KnnPoints.get(temp_kNNKey).getKde());
				}else{
					float temp_kde=KnnPoints.get(temp_kNNKey).getKde();
					if (temp_kde== 0 || o_S.getLrdValue() == 0)
						continue;
					else
					kof+=temp_kde / o_S.getKde()*1.0f;
				}
				
			}
		}
		
		//System.out.println("vgv KOF "+kof);
		if (Float.isNaN(kof) || Float.isInfinite(kof))
			kof = 0;
		kof=kof/k;
		o_S.setKof(kof);
	}
	
	
}
