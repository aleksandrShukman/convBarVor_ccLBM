# convBarVor_ccLBM
This C++ code implements the Lattice Boltzmann Method based on an athermal D3Q19 velocity set and served as a foundation in a study of spurious aero-acoustic emissions created by the interaction of a convected barotropic vortex and cell-centered grid refinement interface, located in the middle of a double-periodic quasi-2D domain, cf. https://doi.org/10.3390/fluids10020031.

<img width="550" height="323" alt="grafik" src="https://github.com/user-attachments/assets/a94b241d-e0ad-47cc-9a5c-9017ab0926de" />
<br />
<br />

Included are several established collision operators:
+ Classical Bhatnagar-Gross-Krook model, cf. https://doi.org/10.1103/PhysRev.94.511
+ Multiple Relaxation Time model (Gram-Schmidt basis and raw-moments basis) , cf. https://doi.org/10.1103/PhysRevE.100.033305
+ Recursive-regularization model, cf. https://doi.org/10.48550/arXiv.1505.06900
+ Hybrid recursive-regularization model, cf. https://doi.org/10.1080/14685248.2018.1540879

... as well as a cell-centered grid refinement scheme, based on:
+ Uniform explosion, cf. https://doi.org/10.1002/fld.1140
+ Linear interpolation during explosion of populations parallel to the refinement interface, cf. https://doi.org/10.1016/j.physa.2005.09.036

Case parameters and model types have to be selected by the user within BoundaryConditions.txt by commenting/uncommenting the desired feature. The absolute path to this file needs to be set once within main (look for inputFile = "/absolute-path-to-BC-file/BoundaryConditions.txt";).
