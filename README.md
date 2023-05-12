# CONE
Code to run the Clustered Onion Network Ensemble from network data
Hébert-Dufresne et al. (2023)

### Codes

The integration of the ODEs is done in C++ using the Gnu Scientific Library with the Boost library for random number generation and data structures.

File run.sh is a simple bash script that takes the path to a network edgelist and runs everything. It uses the following pipeline:

File decomposition.py produces the CONE parametrization based on clique detection and onion decomposition of a network edgelist.

The parametrization are then stored in the matrices folder.

File tevol_source_sis_cone.cpp contains the main code and takes the transmission rate as its single argument.

File dyn_sis_cone.hpp contains the differential equations to be integrated.

Results are then stored in the results folder.

File plot_figure_localization.py produces a plot of the bifurcation diagram of the CONE.

The plot is stored in different format in the figures folder.

### Data

The file edges_friendship.csv contains the network data used as a case study in the paper.