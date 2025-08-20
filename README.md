# ReverseTimeInference
This repository contains code for the paper "The Landscape of Terminating Processes in Living Systems" by Nicolas Lenner, Matthias Häring, Stephan Eule, Jörg Großhans, and Fred Wolf. It demonstrates reverse-time inference using target state aligned ensembles for the example of cytokinetic ring constriction. 


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

To generate trajectories execute the already compiled c code via python script. In case the program does not run on your machine, you have to compile yourself (make clean, make) and then execute via the python script.

## Contributors

Matthias Häring - Göttingen Campus Institute for Dynamics of Biological Networks

Nicolas Lenner - Institute for Advanced Study Princeton

## Acknowledgements
For nested sampling inference we used dynesty: https://github.com/joshspeagle/dynesty
