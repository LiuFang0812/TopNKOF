Experimental Code for KBS: Scalable KDE-based Top-n Local Outlier Detection over Large-Scale DataStreams

(1)Main Methods
The proposed method consists two versions: UKOF and LUKOF method. 

 Main Class for UKOF method: cellpruning.lof.pruning.ComputeTopNKDE

 Main Class for LUKOF method: cellpruning.lof.pruning.ComputeTopNKDE_LazyUpdate
 
(2)Build and Use the Software Artifact
1.Open Eclipse

2.Import the code named "TopNKOF"

3.Set parameters in "util.SQConfig", such as the number of nearest neighbors k, top outliers n,
window size w and slide size s.

4.Run the corresponding main methods for UKOF and LUKOF, namely "cellpruning.lof.pruning.ComputeTopNKDE" and "cellpruning.lof.pruning.ComputeTopNKDE_LazyUpdate"
