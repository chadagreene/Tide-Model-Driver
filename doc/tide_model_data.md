[&larr; Back to TMD3.0 Main Page](../README.md)

# Getting Tide Model Data
After you [add TMD3.0 to MATLAB](installing_tmd.md), you'll need to get some tide model data. Click on any of the links below to download tide model data files.

# Barotropic Ocean Models
Barotropic ocean models contain tidal coefficients for water height and depth-integrated transport. Depth-averaged tidal current velocities are calculated as the predicted transport values divided by water column thickness. Ocean models are either *global* or *regional*. 

### Global Ocean Models 
* [`TPXO9_atlas_v5.nc`](https://www.chadagreene.com/tide_data/TPXO9_atlas_v5.nc) [963 MB] 1/30 degree global model with 15 tidal consituents. Includes tide height and transport variables. Developed by Gary Egbert and Svetlana Erofeeva. More information on TPXO9-atlas models can be found [here](https://www.tpxo.net/global/tpxo9-atlas).
* [`EOT20_ocean.nc`](https://www.chadagreene.com/tide_data/EOT20_ocean.nc) [23 MB] 1/8 degree resolution global model with 17 tidal constituents (height only). EOT20 is the latest in a series of global ocean tide models developed by *Hart-Davis et al*. at DGFI-TUM. EOT20 is created using residual tidal analysis of multi-mission altimetry data from the period of 1992 to 2019. Eleven satellite altimetry missions are used and the FES2014b tide model used as the reference model for the residual tidal analysis. The model extends throughout the global ocean, with EOT20 ranging from 66°S to 66°N with FES2014b being used to fill in the higher latitudes. Tidal constituents are: 2N2, J1, K1, K2, M2, M4, MF, MM, N2, O1, P1, Q1, S1, S2, SA, SSA and T2. More information on EOT20 can be found at the [SEANOE data page](https://doi.org/10.17882/79489) and in the [*Hart-Davis et al., 2021* ESSD paper](https://doi.org/10.5194/essd-13-3869-2021).

### Regional Ocean Models 
* [`CATS2008_update_2022-06-11.nc`](https://www.chadagreene.com/tide_data/CATS2008_update_2022-06-11.nc) [355 MB] 2 km Antarctic regional model with 10 tidal constituents, from Susan Howard and Laurie Padman. Updates from the original CATS2008 are described [here](cats2008_updates.md). 
* [`Gr1kmTM_v1.nc`](https://www.chadagreene.com/tide_data/Gr1kmTm_v1.nc) [193 MB] 1 km resolution, 8 constituent regional model surrounding Greenland, from Susan Howard and Laurie Padman. The Greenland 1 kilometer Tide Model (Gr1kmTM) is a barotropic ocean tide model on a 1 km x 1 km polar stereographic grid, developed using the Regional Ocean Modeling System (ROMS). Gr1kmTM consists of spatial grids of complex amplitude coefficients for sea surface height and depth-integrated currents (“volume transports”) for 8 principal tidal constituents: 4 semidiurnal (M2, S2, K2, N2) and 4 diurnal (K1, O1, P1, Q1). More information on Gr1kmTM_v1 can be found [here](https://doi.org/10.18739/A2B853K18). 
* [`Arc2kmTM_v1.nc`](https://www.chadagreene.com/tide_data/Arc2kmTm_v1.nc) [206 MB] 2 km resolution, 8 constituent regional model of the Arctic Ocean by Susan Howard and Laurie Padman. The Arctic 2 km (kilometer) Tide Model (Arc2kmTM) is a barotropic ocean tide model on a 2x2 km polar stereographic grid, developed using the Regional Ocean Modeling System (ROMS). Arc2kmTM consists of spatial grids of complex amplitude coefficients for sea surface height and depth-integrated currents (“volume transports”) for 8 principal tidal constituents: 4 semidiurnal (M2, S2, K2, N2) and 4 diurnal (K1, O1, P1, Q1). More information on Arc2kmTM_v1 can be found [here](https://doi.org/10.18739/A2PV6B79W). 

# Load Tides
Load tides refer to the deformation of the Earth's crust that occurs in response to the weight of ocean water when tides are high or low. A tide gauge such as a bottom-pressure recorder rests on the seafloor, so it moves up and down with the load tide, and therefore only measures the ocean tide. A satellite altimeter, on the other hand, measures the ocean's surface height, which is the sum of the load tide and the ocean tide. 

These load tide model data files are compatible with TMD3.0: 

* [`EOT20_load.nc`](https://www.chadagreene.com/tide_data/EOT20_load.nc) [20 MB] 1/8 degree resolution global model with 17 tidal constituents. EOT20 is the latest in a series of global ocean tide models developed by *Hart-Davis et al*. at DGFI-TUM. EOT20 is created using residual tidal analysis of multi-mission altimetry data from the period of 1992 to 2019. Eleven satellite altimetry missions are used and the FES2014b tide model used as the reference model for the residual tidal analysis. The model extends throughout the global ocean, with EOT20 ranging from 66°S to 66°N with FES2014b being used to fill in the higher latitudes. Tidal constituents are: 2N2, J1, K1, K2, M2, M4, MF, MM, N2, O1, P1, Q1, S1, S2, SA, SSA and T2. More information on EOT20 can be found at the [SEANOE data page](https://doi.org/10.17882/79489) and in the [*Hart-Davis et al., 2021* ESSD paper](https://doi.org/10.5194/essd-13-3869-2021).

# Tide Model Intercomparison
For a quick intercomparison of tide model performance, check out [this intercomparison](tide_model_intercomparison.md) page.


# Author Info
This page was written by [Chad A. Greene](https://www.chadagreene.com), June 2022. 