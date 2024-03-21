# ReverseTimeInference
This repository demonstrates inference in target state aligned ensembles for the example of cytokinetic ring constriction. Please cite our corresponding publication 

"Reverse Time Analysis of Target-State Directed Processes in Living Systems" - https://arxiv.org/abs/2304.03226

Contents of this repository:
- cytokinesis inference:  TSAInferenceExampleForCytokinesis.ipynb
- inference subsampling: TSAEnsembleVsPathInference.ipynb
- code to generate an ensemble of SDEs with power law force
- code to generate the cytokinesis ensemble


## Getting started

prerequisites:
- numpy
- scipy
- matplotlib
- dynesty

optional:
- seaborn


## Acknowledgements
For inference we used dynesty: https://github.com/joshspeagle/dynesty
