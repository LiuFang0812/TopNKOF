package cellpruning.lof.pruning;

import java.util.HashMap;
import java.util.Map;

public class publicFunction {
	public static long maxValueForTopn(HashMap<Long,Float> list){
		//System.out.println("list³¤¶È  "+list.size());
		long id=0;
		float value=-Float.MIN_VALUE;
		for(Map.Entry<Long, Float> entry:list.entrySet()){
			//System.out.println("value "+entry.getKey()+"  "+entry.getValue()+"   "+value);
			if(entry.getValue()>value){
				value=entry.getValue();
				id=entry.getKey();
			//	System.out.println("bhb  id   "+id+"   "+entry.getKey());
			}
		}
	//	System.out.println("value   "+value);
		return id;
	}
	
	public static long minValueForTopn(HashMap<Long,Float> list){
		long id=0;
		float value=0.0f;
		for(Map.Entry<Long, Float> entry:list.entrySet()){
			if(entry.getValue()<value){
				value=entry.getValue();
				id=entry.getKey();
			}
		}
		return id;
		
	}
	
	
	
	
	
}
