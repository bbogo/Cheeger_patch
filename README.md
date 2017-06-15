# Cheeger_patch
Cheeger sets and optimal Cheeger patches in 2D

Example of implementation of the search for an optimal alpha-Cheeger cluster using finite elements in 2D. The theoretical aspects can be found in the article "Phase Field Approach to Optimal Packing Problems and Related Cheeger Clusters" by Beiamin Bogosel, Dorin Bucur and Ilaria Fragala (2017).

Contents: 
1. A simple implementation of the Kawohl, Lachand-Robert algorithm for finding Cheeger sets of convex domains in 2D. This uses the toolbox Clipper together with a Matlab interface.
2. An optimization algorithm based on a Gamma-convergence relaxation for finding the Cheeger sets for general domains.
3. A Gamma-convergence approach for finding optimal Cheeger cluster and optimal packings
4. A  gradient-free local optimizing algorithm for a circle/sphere packing.

Requirements:
1. Optimization algorithms (necessary in order to perform the optimization):
- LBFGS wrapper http://www.cs.toronto.edu/~liam/software.shtml 
- minConf optimizer https://www.cs.ubc.ca/~schmidtm/Software/minConf.html 
2. Polygon offset algorithms (necessary for the implementation of the Kawohl, Lachand-Robert algorithm)
- Clipper toolbox http://www.angusj.com/delphi/clipper.php
- Matlab interface for Clipper https://fr.mathworks.com/matlabcentral/fileexchange/61329-new-polygon-clipping-and-offsetting

Install: 
Just copy all the files into one folder and make sure this folder is in the Matlab path. 

Usage: main functions are listed below. Type 'help function_name' to find detailed instructions and examples
- FEMcheeger2.m    : the optimization algorithm
- cheeger_poly.m   : Kawohl, Lachand-Robert algorithm
- CH_testing.m     : testing of the KLR algorithm
