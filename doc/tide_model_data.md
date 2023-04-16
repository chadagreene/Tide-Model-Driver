[&larr; Back to TMD3.0 Main Page](../README.md)

# Getting Tide Models
After you [add TMD3.0 to MATLAB](installing_tmd.md), you'll need to get one or more tide model files described below. For TMD3.0, these will always be NetCDF ( \*.nc) files, but when you download them, they may be compressed into files with a .zip extension, in which case, unzip the folder and you'll be ready to start using TMD3.0.
 
## Barotropic Ocean Models
Barotropic ocean models contain complex tidal coefficients for water height and depth-integrated transport, for a set of tidal constituents. Depth-averaged tidal current velocity coefficients are calculated as the predicted transport values divided by water column thickness. Ocean models are either *global* or *regional*. 

### Global Ocean Models
* **`TPXO9_atlas_v5.nc`** by Gary Egbert and Svetlana Erofeeva.
	* **Variables:** Complex height & transport coefficients.
	* **Resolution:** 1/30 degree.
	* **Constituents:** 15: 2n2 k1 k2 m2 m4 mf mm mn4 ms4 n2 o1 p1 q1 s1 s2.
	* **Restrictions:** For noncommercial use only (contact Gary Egbert or Svetlana Erofeeva otherwise). 
	* **Data Access:** [Download here](https://www.chadagreene.com/tide_data/TPXO9_atlas_v5.nc.zip) [1.6 GB].
	* **Website:** [https://www.tpxo.net/global/tpxo9-atlas](https://www.tpxo.net/global/tpxo9-atlas)
	* **Data Citation:** Egbert, G. D., and Erofeeva, S. Y. (2002). Efficient inverse modeling of barotropic ocean tides. *Journal of Atmospheric and Oceanic Technology* 19(2), 183-204.
	
* **`EOT20_ocean.nc`** by Hart-Davis et al. EOT20 is the latest in a series of global ocean tide models developed by *Hart-Davis et al*. at DGFI-TUM. EOT20 is created using residual tidal analysis of multi-mission altimetry data from 1992-2019. Eleven satellite altimetry missions are used and the FES2014b tide model used as the reference model for the residual tidal analysis. The model extends throughout the global ocean, with EOT20 ranging from 66°S to 66°N with FES2014b being used to fill in the higher latitudes.
	* **Variables:** Complex height only.
	* **Resolution:** 1/8 degree.
	* **Constituents:** 17: 2n2 j1 k1 k2 m2 m4 mf mm n2 o1 p1 q1 s1 s2 sa ssa t2.
	* **Restrictions:** Users must cite Hart-Davis et al., 2021. 
	* **Data Access:** [Download here](https://www.chadagreene.com/tide_data/EOT20_ocean.nc.zip) [26 MB].
	* **Website:** [SEANOE data page](https://doi.org/10.17882/79489).
	* **Data Citation:** Hart-Davis, M., Piccioni, G., Dettmering, D., Schwatke, C., Passaro, M., and Seitz, F. (2021). EOT20 - A global Empirical Ocean Tide model from multi-mission satellite altimetry. SEANOE.  [doi:10.17882/79489](https://doi.org/10.17882/79489).
	* **Literature Reference:** Hart-Davis, M., Piccioni, G., Dettmering, D., Schwatke, C., Passaro,M., and Seitz, F.  (2021). EOT20: a global ocean tide model from multi-mission satellite altimetry. *Earth System Science Data*, 13(8), 3869-3884. [doi:10.5194/essd-13-3869-2021](https://doi.org/10.5194/essd-13-3869-2021).

### Regional Ocean Models 
* **`CATS2008_v2023.nc`** Circum-Antarctic Tidal Simulation by Lana Erofeeva, Chad Greene, Susan Howard, and Laurie Padman. Updates from the original CATS2008 are described [here](cats2008_updates.md). 
	* **Variables:** Complex height & transport coefficients.
	* **Resolution:** 2 km.
	* **Constituents:** 10: m2 s2 n2 k2 k1 o1 p1 q1 mf mm.
	* **Restrictions:** This work is licensed under the Creative Commons Attribution 4.0 International License. 
	* **Data Access:** [Download here](https://www.chadagreene.com/tide_data/CATS2008_v2023.nc.zip) [519 MB].
	* **Website:** .
	* **Data Citation:** .
	
* **`Gr1kmTM_v1.nc`** by Susan Howard and Laurie Padman is a regional model of Greenland, developed using the Regional Ocean Modeling System (ROMS) on a standard north polar stereographic projection. 
	* **Variables:** Complex height & transport coefficients.
	* **Resolution:** 1 km.
	* **Constituents:** 8: m2 s2 k1 o1 n2 p1 k2 q1.
	* **Restrictions:** This work is licensed under the Creative Commons Attribution 4.0 International License.. 
	* **Data Access:** [Download here](https://www.chadagreene.com/tide_data/Gr1kmTM_v1.nc.zip) [286 MB].
	* **Website:**  [https://arcticdata.io/catalog/portals/ArcticTides]( https://arcticdata.io/catalog/portals/ArcticTides).
	* **Data Citation:** Howard, S. L., and Padman, L. (2021). Gr1kmTM: Greenland 1 kilometer Tide Model, 2021. Arctic Data Center. [doi:10.18739/A2B853K18](https://doi.org/10.18739/A2B853K18)
	
* **`Arc2kmTM_v1.nc`** by Susan Howard and Laurie Padman is a regional model of the Arctic Ocean, developed using the Regional Ocean Modeling System (ROMS) on a standard north polar stereographic projection.  
	* **Variables:** Complex height & transport coefficients.
	* **Resolution:** 2 km.
	* **Constituents:** 8: m2 s2 k1 o1 n2 p1 k2 q1.
	* **Restrictions:** This work is licensed under the Creative Commons Attribution 4.0 International License.
	* **Data Access:** [Download here](https://www.chadagreene.com/tide_data/Arc2kmTM_v1.nc.zip) [293 MB].
	* **Website:**  [https://arcticdata.io/catalog/portals/ArcticTides]( https://arcticdata.io/catalog/portals/ArcticTides).
	* **Data Citation:** Howard, S. L., and Padman, L. (2021). Arc2kmTM: Arctic 2 kilometer Tide Model, 2021. Arctic Data Center. [doi:10.18739/A2D21RK6K](https://doi.org/10.18739/A2D21RK6K).
	
* **`Arc5km2018.nc`** by Svetlana Erofeeva and Gary Egbert is an update to the original Arctic Ocean Tidal Inverse Model (AOTIM5) model developed in 2004, described by [Padman and Erofeeva (2004)](https://doi.org/10.1029/2003GL019003). The consolidated NetCDF data file for TMD3.0 has also been interpolated from the original custom grid to a polar stereographic grid with a standard longitude at 0°E.
	* **Variables:** Complex height & transport coefficients.
	* **Resolution:** 5 km.
	* **Constituents:** 12: m2 s2 n2 k2 k1 o1 p1 q1 m4 mn4 ms4 2n2.
	* **Restrictions:** This work is licensed under the Creative Commons Attribution 4.0 International License. 
	* **Data Access:** [Download here](https://www.chadagreene.com/tide_data/Arc5km2018.nc.zip) [47 MB].
	* **Website:**  [https://arcticdata.io/catalog/portals/ArcticTides]( https://arcticdata.io/catalog/portals/ArcticTides).
	* **Data Citation:** Erofeeva, S. and Egbert, G. 2020. Arc5km2018: Arctic Ocean Inverse Tide Model on a 5 kilometer grid, 2018. Arctic Data Center. [doi:10.18739/A21R6N14K](https://doi.org/10.18739/A21R6N14K).

* **`AOTIM5.nc`** Arctic Ocean Tidal Inverse Model by Padman, Erofeeva, & Howard is a barotropic model created with the OSU Tidal Inversion Software (OTIS) package. AOTIM5 used AODTM5 as a “prior” model, then assimilated coastal and benthic tide gauges, and TOPEX/Poseidon and ERS satellite radar altimetry, to improve the 4 largest-amplitude constituents, M2, S2, K1 and O1. The consolidated NetCDF data file for TMD3.0 has also been interpolated from the original custom grid to a polar stereographic grid with a standard longitude at 0°E.
	* **Variables:** Complex height & transport coefficients.
	* **Resolution:** 5 km.
	* **Constituents:** 8: m2 s2 n2 k2 k1 o1 p1 q1.
	* **Restrictions:** This work is licensed under the Creative Commons Attribution 4.0 International License.
	* **Data Access:** [Download here](https://www.chadagreene.com/tide_data/AOTIM5.nc.zip) [38 MB].
	* **Website:**  [https://arcticdata.io/catalog/portals/ArcticTides]( https://arcticdata.io/catalog/portals/ArcticTides).
	* **Data Citation:** Padman, L., Erofeeva, S., and Howard, S. (2020). AOTIM5: Arctic Ocean Inverse Tide Model, on 5 kilometer grid, developed in 2004. [doi:10.18739/A2S17SS80](https://doi.org/10.18739/A2S17SS80).
	* **Literature Reference:** Padman, L., and S. Erofeeva (2004). A barotropic inverse tidal model for the Arctic Ocean. *Geophysical Research Letters* 31(2). [doi:10.1029/2003GL019003](https://doi.org/10.1029/2003GL019003).
	
* **`AODTM5.nc`** Arctic Ocean Dynamics-based Tide Model by Padman, Erofeeva, & Howard is a barotropic model model that solves the depth-integrated shallow water equations, with forcing at open ocean boundaries from the TOPEX/Poseidon global barotropic tidal solution version 6.2 (TPXO.6.2) and local astronomical forcing (“potential tides”). Each constituent in AODTM5 was tuned, separately, to Arctic tide height data by optimizing the linear drag coefficient. The consolidated NetCDF data file for TMD3.0 has also been interpolated from the original custom grid to a polar stereographic grid with a standard longitude at 0°E.
	* **Variables:** Complex height & transport coefficients.
	* **Resolution:** 5 km.
	* **Constituents:** 8: m2 s2 n2 k2 k1 o1 p1 q1.
	* **Restrictions:** This work is licensed under the Creative Commons Attribution 4.0 International License.. 
	* **Data Access:** [Download here](https://www.chadagreene.com/tide_data/AODTM5.nc.zip) [28 MB].
	* **Website:**  [https://arcticdata.io/catalog/portals/ArcticTides]( https://arcticdata.io/catalog/portals/ArcticTides).
	* **Data Citation:** : Padman, L., Erofeeva, S., and Howard, S. 2020. AODTM5: Arctic Ocean Dynamics-based Tide Model, on 5 kilometer grid, developed in 2004. Arctic Data Center. [doi:10.18739/A2901ZG3N](https://doi.org/10.18739/A2901ZG3N).
	* **Literature Reference:** Padman, L., and S. Erofeeva (2004). A barotropic inverse tidal model for the Arctic Ocean. *Geophysical Research Letters*, 31(2). [doi:10.1029/2003GL019003](https://doi.org/10.1029/2003GL019003).


## Load Tides
Load tides refer to the deformation of the Earth's crust that occurs in response to the weight of ocean water when tides are high or low. A tide gauge such as a bottom-pressure recorder rests on the seafloor, so it moves up and down with the load tide, and therefore only measures the ocean tide. A satellite altimeter, on the other hand, measures the ocean's surface height, which is the sum of the load tide and the ocean tide. 

These load tide model data files are compatible with TMD3.0: 

* **`EOT20_load.nc`** by Hart-Davis et al. EOT20 is the latest in a series of global ocean tide models developed by *Hart-Davis et al*. at DGFI-TUM. EOT20 is created using residual tidal analysis of multi-mission altimetry data from the period of 1992 to 2019. Eleven satellite altimetry missions are used and the FES2014b tide model used as the reference model for the residual tidal analysis. The model extends throughout the global ocean, with EOT20 ranging from 66°S to 66°N with FES2014b being used to fill in the higher latitudes.
	* **Variables:** Complex height only.
	* **Resolution:** 1/8 degree.
	* **Constituents:** 17: 2n2 j1 k1 k2 m2 m4 mf mm n2 o1 p1 q1 s1 s2 sa ssa t2.
	* **Restrictions:** Users must cite Hart-Davis et al., 2021. 
	* **Data Access:** [Download here](https://www.chadagreene.com/tide_data/EOT20_load.nc.zip) [23 MB].
	* **Website:** [SEANOE data page](https://doi.org/10.17882/79489).
	* **Data Citation:** Hart-Davis, M., Piccioni, G., Dettmering, D., Schwatke, C., Passaro, M., and Seitz, F. (2021). EOT20 - A global Empirical Ocean Tide model from multi-mission satellite altimetry. SEANOE.  [doi:10.17882/79489](https://doi.org/10.17882/79489).
	* **Literature Reference:** Hart-Davis, M., Piccioni, G., Dettmering, D., Schwatke, C., Passaro,M., and Seitz, F. (2021). EOT20: a global ocean tide model from multi-mission satellite altimetry. *Earth System Science Data*, 13(8), 3869-3884. [doi:10.5194/essd-13-3869-2021](https://doi.org/10.5194/essd-13-3869-2021).
	
# Tide Model Intercomparison
For a quick intercomparison of tide model performance, check out [this intercomparison](tide_model_intercomparison.md) page. Only a few models are compared on that page, but it provides a template for how you may explore model performance on your own. 

# Author Info
This page was written by [Chad A. Greene](https://www.chadagreene.com), June 2022. 
