### Code to reproduce figures/tables in our manuscript
Please find below pointers to code in this repository that can be used to reproduce the figures/tables in the main text and supplement of our manuscript: _"Robust discovery of causal gene networks via measurement error estimation and correction" by Rahul Biswas, Brintha V P, Amol Dumrewal and Manikandan Narayanan._

#### Figure/Table to Code mappings for Main text:

----Figures

Fig. 1 -> shotnoise.R - Modeling of shot noise using Poisson distribution in the Overview figure

Fig. 2 -> 2a, 2b, 2c -> simulation_yeast.manuscript.R - causal.model.plot.analysis()

Fig. 3 -> 3a findr.analysis.manuscript.R - pAUPR.all()

Fig. 3 -> 3b findr.analysis.manuscript.R - DESeq_alternative_using_package()

Fig. 3 -> 3c findr.analysis.manuscript.R - pAUPR.noisy()

Fig. 3 -> 3d,e findr.analysis.manuscript.R - precision.recall.noisy.fdr()

Fig. 4 -> 4a The orientations obtained by running allworkflow.gtex.withintissue.R script are saved in cytoscape.muscle.txt and then visualized using the cytoscape tool. 

Fig. 4 -> 4b - allworkflow.gtex.withintissue.R - boxplot.analysis() inside meqtl.cit.ecit function 


#### Figure/Table to Code mappings for Supplement:

----Suppl Figures

S1 Fig: simulation_yeast.manuscript.R - independent.model.plot.analysis()

S2 Fig: findr.analysis.manuscript.R - pAUPR.all() - obtained by setting the error.ratio variable to 0.3 and 0.5.

S3 Fig: allworkflow.gtex.withintissue.R - boxplot.analysis() inside meqtl.cit.ecit function 

S4 Fig: The entire network obtained by running allworkflow.gtex.withintissue.R script is saved in cytoscape.muscle.txt and then visualized using the cytoscape tool.


