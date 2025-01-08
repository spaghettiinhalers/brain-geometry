CONTENTS

 - ```data/```: data needed for analysis
 - ```code/```: scripts
 - ```results/```: stores output
 - ```demo/```: common functions used in the scripts in code/, since most eigenmode calculation and analysis processes are very similar
 - ```tuxerva.png```: required file

DEPENDENCIES

 - ```numpy v1.22.4```
 - ```brainspace v0.1.1```
 - ```vtk v9.1.0```
 - ```lapy v1.0.0```

i do not recommend using updated versions of these libraries due to function deprecation issues i initially encountered. there may also be other dependencies not listed here (especially within matlab), but it should otherwise run properly based on what software and packages have been pre-installed on the mpi cbs system.

REPLICATION

1. download ```code/``` and ```data/```, and set up an empty ```results/``` folder with the same structure as the one provided in this repo
2. set up environments, dependencines and missing ```data/``` files as appropriate
3. run every script in ```code/``` to produce an output in ```results/``` identical to the folder in this repo
4. adjust scripts as needed (using ```demo/``` will be helpful)

MORE

please see the repository for the original paper on which this project was based ([https://www.nature.com/articles/s41586-023-06098-1](https://www.nature.com/articles/s41586-023-06098-1)) for further detail on data sources, script sources, hpc dependencies, eigenmode analysis, etc: [https://github.com/NSBLab/BrainEigenmodes](https://github.com/NSBLab/BrainEigenmodes)
