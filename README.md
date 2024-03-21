# ReverseTimeInference
This repository demonstrates inference in target state aligned ensembles for the example of cytokinetic ring constriction. Please cite our publication:

"Reverse Time Analysis of Target-State Directed Processes in Living Systems" - https://arxiv.org/abs/2304.03226


## Getting started

python prerequisites:
- numpy
- scipy
- matplotlib
- dynesty
optional:
- seaborn

Inference explained in:
-  TSAInferenceExampleForCytokinesis.ipynb

Path Inference Vs Ensemble Inference:
- TSAEnsembleVsPathInference.ipynb

Langevin simulation and corresponding data:
- /generatePowerlawEnsemble/

Cytokinesis simulation and corresponding data:
- /generateCytokinEnsemble/

To generate trajectories compile the C++ code on your own machine (make clean, make) and then execute via python script.


## Acknowledgements
For inference we used dynesty: https://github.com/joshspeagle/dynesty
