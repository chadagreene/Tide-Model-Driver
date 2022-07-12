[&larr; Back to TMD3.0 Main Page](../README.md)

# Getting Tide Model Data
After you [add TMD3.0 to MATLAB](installing_tmd.md), you'll need to get some tide model data. Click on any of the links below to download tide model data files.

## Barotropic Ocean Models
Barotropic ocean models contain tidal coefficients for water height and depth-integrated transport. Depth-averaged tidal current velocities are calculated as the predicted transport values divided by water column thickness. Ocean models are either *global* or *regional*. 

### Global Ocean Models
* [**`TPXO9_atlas_v5.nc`**](https://www.chadagreene.com/tide_data/TPXO9_atlas_v5.nc) [963 MB] by Gary Egbert and Svetlana Erofeeva.
	* **Variables:** height & transport.
	* **Resolution:** 1/30 degree.
	* **Constituents:** 15: 2n2 k1 k2 m2 m4 mf mm mn4 ms4 n2 o1 p1 q1 s1 s2.
	* **Restrictions:** For noncommercial use only (contact Gary Egbert or Svetlana Erofeeva otherwise). 
	* **Data Access:** [Download here](https://www.chadagreene.com/tide_data/TPXO9_atlas_v5.nc) [963 MB].
	* **Website:** [https://www.tpxo.net/global/tpxo9-atlas](https://www.tpxo.net/global/tpxo9-atlas)
	* **Data Citation:** Egbert, Gary D., and Svetlana Y. Erofeeva. "Efficient inverse modeling of barotropic ocean tides." Journal of Atmospheric and Oceanic Technology 19.2 (2002): 183-204.
* [**`EOT20_ocean.nc`**](https://www.chadagreene.com/tide_data/EOT20_ocean.nc) [23 MB] by Hart-Davis et al. EOT20 is the latest in a series of global ocean tide models developed by *Hart-Davis et al*. at DGFI-TUM. EOT20 is created using residual tidal analysis of multi-mission altimetry data from the period of 1992 to 2019. Eleven satellite altimetry missions are used and the FES2014b tide model used as the reference model for the residual tidal analysis. The model extends throughout the global ocean, with EOT20 ranging from 66°S to 66°N with FES2014b being used to fill in the higher latitudes.
	* **Variables:** height only.
	* **Resolution:** 1/8 degree.
	* **Constituents:** 17: 2n2 j1 k1 k2 m2 m4 mf mm n2 o1 p1 q1 s1 s2 sa ssa t2.
	* **Restrictions:** Users must cite Hart-Davis et al., 2021. 
	* **Data Access:** [Download here](https://www.chadagreene.com/tide_data/EOT20_ocean.nc) [23 MB].
	* **Website:** [SEANOE data page](https://doi.org/10.17882/79489).
	* **Literature Reference:** Hart-Davis et al., 2021 ESSD [doi:10.5194/essd-13-3869-2021](https://doi.org/10.5194/essd-13-3869-2021).
	* **Data Citation:** Hart-Davis Michael, Piccioni Gaia, Dettmering Denise, Schwatke Christian, Passaro Marcello, Seitz Florian (2021). EOT20 - A global Empirical Ocean Tide model from multi-mission satellite altimetry. SEANOE. [doi:10.17882/79489](https://doi.org/10.17882/79489).

### Regional Ocean Models 
* [**`CATS2008_update_2022-06-11.nc`**](https://www.chadagreene.com/tide_data/CATS2008_update_2022-06-11.nc) [355 MB] Circum-Antarctic Tidal Simulation by Susan Howard and Laurie Padman. Updates from the original CATS2008 are described [here](cats2008_updates.md). 
	* **Variables:** height & transport.
	* **Resolution:** 2 km.
	* **Constituents:** 10: m2 s2 n2 k2 k1 o1 p1 q1 mf mm.
	* **Restrictions:** . 
	* **Data Access:** [Download here](https://www.chadagreene.com/tide_data/CATS2008_update_2022-06-11.nc) [355 MB].
	* **Website:** .
	* **Literature Reference:** .
	* **Data Citation:** .
* [**`Gr1kmTM_v1.nc`**](https://www.chadagreene.com/tide_data/Gr1kmTM_v1.nc) [193 MB] by Susan Howard and Laurie Padman is a regional model of Greenland, developed using the Regional Ocean Modeling System (ROMS) on a standard north polar stereographic projection. 
	* **Variables:** height & transport.
	* **Resolution:** 1 km.
	* **Constituents:** 8: m2 s2 k1 o1 n2 p1 k2 q1.
	* **Restrictions:** . 
	* **Data Access:** [Download here](https://www.chadagreene.com/tide_data/Gr1kmTM_v1.nc) [193 MB].
	* **Website:** .
	* **Literature Reference:** .
	* **Data Citation:** Susan L. Howard and Laurie Padman. 2021. Gr1kmTM: Greenland 1 kilometer Tide Model, 2021. [doi:10.18739/A2251FM3S](https://www.doi.org/doi:10.18739/A2251FM3S).
* [**`Arc2kmTM_v1.nc`**](https://www.chadagreene.com/tide_data/Arc2kmTM_v1.nc) [206 MB] by Susan Howard and Laurie Padman is a regional model of the Arctic Ocean, developed using the Regional Ocean Modeling System (ROMS) on a standard north polar stereographic projection.  
	* **Variables:** height & transport.
	* **Resolution:** 2 km.
	* **Constituents:** 8: m2 s2 k1 o1 n2 p1 k2 q1.
	* **Restrictions:** . 
	* **Data Access:** [Download here](https://www.chadagreene.com/tide_data/Arc2kmTM_v1.nc) [206 MB].
	* **Website:** .
	* **Literature Reference:** .
	* **Data Citation:** Susan L. Howard and Laurie Padman. 2021. Arc2kmTM: Arctic 2 kilometer Tide Model, 2021. Arctic Data Center. [doi:10.18739/A2PV6B79W](https://doi.org/10.18739/A2PV6B79W).
* [**`Arc5km2018.nc`**](https://www.chadagreene.com/tide_data/Arc5km2018.nc) [40 MB] by Susan Howard and Laurie Padman is an update to the original Arctic Ocean Tidal Inverse Model (AOTIM5) model developed in 2004, described by [Padman and Erofeeva (2004)](https://doi.org/10.1029/2003GL019003). The consolidated NetCDF data file for TMD3.0 has also been interpolated from the original custom grid to a polar stereographic grid with a standard longitude at 0°E.
	* **Variables:** height & transport.
	* **Resolution:** 5 km.
	* **Constituents:** 12: m2 s2 n2 k2 k1 o1 p1 q1 m4 mn4 ms4 2n2.
	* **Restrictions:** . 
	* **Data Access:** [Download here](https://www.chadagreene.com/tide_data/Arc5km2018.nc) [40 MB].
	* **Website:** .
	* **Literature Reference:** .
	* **Data Citation:** Erofeeva S and Egbert G. 2020. **<--FIX?** Arc5km2018: Arctic Ocean Inverse Tide Model on a 5 kilometer grid, 2018 [doi:10.18739/A21R6N14K](https://doi.org/10.18739/A21R6N14K).
* [**`AOTIM5.nc`**](https://www.chadagreene.com/tide_data/AOTIM5.nc) [29 MB] Arctic Ocean Tidal Inverse Model by Padman, Erofeeva, & Howard is a barotropic model forced at open ocean boundaries by the TOPEX/Poseidon global barotropic tidal solution version 6.2 (TPXO.6.2), and by local astronomical forcing (“potential tides”). Each constituent in AODTM5 was tuned, separately, to Arctic tide height data by optimizing the linear drag coefficient. AOTIM5 used AODTM5 as a “prior” model, then assimilated coastal and benthic tide gauges, and TOPEX/Poseidon and ERS satellite radar altimetry, to improve the 4 largest-amplitude constituents, M2, S2, K1 and O1. The consolidated NetCDF data file for TMD3.0 has also been interpolated from the original custom grid to a polar stereographic grid with a standard longitude at 0°E.
	* **Variables:** height & transport.
	* **Resolution:** 5 km.
	* **Constituents:** 8: m2 s2 n2 k2 k1 o1 p1 q1.
	* **Restrictions:** . 
	* **Data Access:** [Download here](https://www.chadagreene.com/tide_data/AOTIM5.nc) [29 MB].
	* **Website:** .
	* **Literature Reference:** Padman, L., and S. Erofeeva. "A barotropic inverse tidal model for the Arctic Ocean." *Geophysical Research Letters* 31.2 (2004) [doi:10.1029/2003GL019003](https://doi.org/10.1029/2003GL019003).
	* **Data Citation:** Padman L , Erofeeva S , and Howard S. 2020. AOTIM5: Arctic Ocean Inverse Tide Model, on 5 kilometer grid, developed in 2004. [doi:10.18739/A2S17SS80](https://doi.org/10.18739/A2S17SS80).
* [**`AODTM5.nc`**](https://www.chadagreene.com/tide_data/AODTM5.nc) [23 MB] Arctic Ocean Dynamics-based Tide Model by Padman, Erofeeva, & Howard is a barotropic model model that solves the depth-integrated shallow water equations, with forcing at open ocean boundaries from the TOPEX/Poseidon global barotropic tidal solution version 6.2 (TPXO.6.2) and local astronomical forcing (“potential tides”). Each constituent in AODTM5 was tuned, separately, to Arctic tide height data by optimizing the linear drag coefficient. The consolidated NetCDF data file for TMD3.0 has also been interpolated from the original custom grid to a polar stereographic grid with a standard longitude at 0°E.
	* **Variables:** height & transport.
	* **Resolution:** 5 km.
	* **Constituents:** 8: m2 s2 n2 k2 k1 o1 p1 q1.
	* **Restrictions:** . 
	* **Data Access:** [Download here](https://www.chadagreene.com/tide_data/AODTM5.nc) [23 MB].
	* **Website:** .
	* **Literature Reference:** Padman, L., and S. Erofeeva. "A barotropic inverse tidal model for the Arctic Ocean." *Geophysical Research Letters* 31.2 (2004) [doi:10.1029/2003GL019003](https://doi.org/10.1029/2003GL019003).
	* **Data Citation:** Padman L , Erofeeva S , and Howard S. 2020. AODTM5: Arctic Ocean Dynamics-based Tide Model, on 5 kilometer grid, developed in 2004. [doi:10.18739/A2901ZG3N](https://doi.org/10.18739/A2901ZG3N).
* **`placeholder.nc`** . 
	* **Variables:** height & transport.
	* **Resolution:** .
	* **Constituents:** .
	* **Restrictions:** . 
	* **Data Access:** [Download here](https://www.chadagreene.com/tide_data/.nc) [MB].
	* **Website:** .
	* **Literature Reference:** .
	* **Data Citation:** .



## Load Tides
Load tides refer to the deformation of the Earth's crust that occurs in response to the weight of ocean water when tides are high or low. A tide gauge such as a bottom-pressure recorder rests on the seafloor, so it moves up and down with the load tide, and therefore only measures the ocean tide. A satellite altimeter, on the other hand, measures the ocean's surface height, which is the sum of the load tide and the ocean tide. 

These load tide model data files are compatible with TMD3.0: 

* [**`EOT20_load.nc`**](https://www.chadagreene.com/tide_data/EOT20_load.nc) [20 MB] by Hart-Davis et al. EOT20 is the latest in a series of global ocean tide models developed by *Hart-Davis et al*. at DGFI-TUM. EOT20 is created using residual tidal analysis of multi-mission altimetry data from the period of 1992 to 2019. Eleven satellite altimetry missions are used and the FES2014b tide model used as the reference model for the residual tidal analysis. The model extends throughout the global ocean, with EOT20 ranging from 66°S to 66°N with FES2014b being used to fill in the higher latitudes.
	* **Variables:** height only.
	* **Resolution:** 1/8 degree.
	* **Constituents:** 17: 2n2 j1 k1 k2 m2 m4 mf mm n2 o1 p1 q1 s1 s2 sa ssa t2.
	* **Restrictions:** Users must cite Hart-Davis et al., 2021. 
	* **Data Access:** [Download here](https://www.chadagreene.com/tide_data/EOT20_load.nc) [20 MB].
	* **Website:** [SEANOE data page](https://doi.org/10.17882/79489).
	* **Literature Reference:** Hart-Davis et al., 2021 ESSD [doi:10.5194/essd-13-3869-2021](https://doi.org/10.5194/essd-13-3869-2021).
	* **Data Citation:** Hart-Davis Michael, Piccioni Gaia, Dettmering Denise, Schwatke Christian, Passaro Marcello, Seitz Florian (2021). EOT20 - A global Empirical Ocean Tide model from multi-mission satellite altimetry. SEANOE. [doi:10.17882/79489](https://doi.org/10.17882/79489).

# Tide Model Intercomparison
For a quick intercomparison of tide model performance, check out [this intercomparison](tide_model_intercomparison.md) page. Only a few models are compared on that page, but it provides a template for how you may explore model performance on your own. 

# Author Info
This page was written by [Chad A. Greene](https://www.chadagreene.com), June 2022. 