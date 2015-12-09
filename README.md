![Gene2farm project] (https://github.com/bioinformatics-ptp/Zanardi/wiki/Images/logog2f.png "Gene2farm project logo")

GenoSim: _an open-source multi-population simulator_
=========
ref: E.L. Nicolazzi (Fondazione Parco Tecnologico Padano) - Via Einstein, Loc. Cascina Codazza (26900) Lodi (Italy).


This is a simulator of genomic/phenotypic data, heavily relying on the work from Hossein Jorjani (Interbull Centre, 2009). 
From its original version, GenoSim maintains general functionallity (and routines) but it has been significantly modified such that the current version has the following features: 
  - i) multiple populations can be generated from a single common (base) population that has reached a mutation-drift equilibrium; 
  - ii) the simulation accommodates multiple traits, with three different heritabilities (low, medium and high), positive and negative genetic correlations between traits and traits with two genetic structures (high and low number of QTL, dependent on the number of independent chromosome segments – Me -); 
  - iii) real pedigrees can be included to study the effect of population structure on multiple/across-breed evaluations; 
  - iv) selection can be applied at different intensities (e.g. for male/female populations) and for different traits, and; 
  - v) the number of markers can reach densities approaching full sequence.
  - vi) **THIS VERSION allows admixture between _2_ populations (do not use more than 2!)**

**Full documentation (including how to install and run the program) can be found [here](Manual/GenoSim_software.pdf)**
*Important note: the manual is referred to the base GenoSim program. The only difference here is that this version contains one more line in the parameter file (user_options.txt) in order to account for a user-defined degree of admixture between populations.*

For bug report, feedback and questions contact ezequiel [dot] nicolazzi [at] ptp [dot] it.
