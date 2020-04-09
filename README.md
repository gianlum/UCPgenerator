# UCPgenerator

Python program to derive urban canopy parameters for running COSMO-BEP-Tree.

Further details can be found at:

Mussetti, G., Brunner, D., Henne, S., Allegrini, J., Krayenhoff, E. S., Schubert, S., Feigenwinter, C., Vogt, R., Wicki, A., and Carmeliet, J.: COSMO-BEP-Tree v1.0: a coupled urban climate model with explicit representation of street trees, Geosci. Model Dev., 13, 1685–1710, https://doi.org/10.5194/gmd-13-1685-2020, 2020.


**Inputs**

- building geometry dataset (.shp file, LoD1 preferred)
- soil sealing dataset (.geotiff)
- tree dataset (.geotiff)
- COSMO static input file (laf file)

**Output**
- COSMO static input file updated with the field necessary to run COSMO-BEP-Tree
