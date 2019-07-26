# UCPgenerator

Python program to derive urban canopy parameters for running COSMO-BEP-Tree.

Further details can be found at:

Gianluca Mussetti et al. *Street trees in urban climate modelling – coupling the multi-layerurban canopy model BEP-Tree with the regional climate modelCOSMO-CLM (v5.0_clm2.1).* Submitted to Geoscientific Model Development. 2019


**Inputs**

- building geometry dataset (.shp file, LoD1 preferred)
- soil sealing dataset (.geotiff)
- tree dataset (.geotiff)
- COSMO static input file (laf file)

**Output**
- COSMO static input file updated with the field necessary to run COSMO-BEP-Tree
