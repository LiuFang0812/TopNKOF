Experimental Code for KBS: Scalable KDE-based Top-n Local Outlier Detection over Large-Scale DataStreams
==
Overview
--
The detection of local outliers over high-volume data streams is critical for diverse real-time applications in the real world, where the distributions in different subsets of the data tend to be skewed. However, existing methods are not scalable to large-scale high-volume data streams owing to the high complexity of the re-detection of data updates. In this work, we propose a top-$n$ local outlier detection method based on Kernel Density Estimation (KDE) over large-scale high-volume data streams. 


Main Methods
--
The proposed method consists two versions: UKOF and LUKOF method. 

 Main Class for UKOF method: cellpruning.lof.pruning.ComputeTopNKDE

 Main Class for LUKOF method: cellpruning.lof.pruning.ComputeTopNKDE_LazyUpdate
 
 Environment
 --
* Eclipse
 
* Build and Use the Software Artifact

>>>1.Open Eclipse

>>>2.Import the code named "TopNKOF"

>>>3.Set parameters in "util.SQConfig", such as the number of nearest neighbors k, top outliers n,
window size w and slide size s.

>>>4.Run the corresponding main methods for UKOF and LUKOF, namely "cellpruning.lof.pruning.ComputeTopNKDE" and "cellpruning.lof.pruning.ComputeTopNKDE_LazyUpdate"

Dataset
--
 extensive experiments are conducted on ten real-world and synthetic datasets. The real-world datasets are extracted from  UCI Machine Learning Repository.
 
 synthetic dataset:Interchanging RBF,Moving Squares,Mixture RBF
 
 real-word datasets: Vowels,KDDCup,Subhttp,Smtp,ForestCover,Mobike,GeoLife
