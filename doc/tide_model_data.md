[&larr; Back to TMD3.0 Main Page](../README.md)

# Getting Tide Model Data
After you [add TMD3.0 to MATLAB](installing_tmd.md), you'll need to get some tide model data. Click on any of the links below to download tide model data files.

# Ocean Models
Ocean models contain water height predictions and the barotropic transport variables used to calculate tidal current velocities. Ocean models are either *global* or *regional*. 

### Global Ocean Models 
* [`TPXO9_atlas_v5.nc`](https://www.chadagreene.com/tide_data/TPXO9_atlas_v5.nc) [963 MB] 1/30 degree global model with 15 tidal consituents, from Gary Egbert and Svetlana Erofeeva. More information on TPXO9-atlas models can be found [here](https://www.tpxo.net/global/tpxo9-atlas).

### Regional Ocean Models 
* [`CATS2008_update_2022-06-11.nc`](https://www.chadagreene.com/tide_data/CATS2008_update_2022-06-11.nc) [355 MB] 2 km Antarctic regional model with 10 tidal constituents, from Susan Howard and Laurie Padman. Updates from the original CATS2008 are described [here](cats2008_updates.md). 
* [`Gr1kmTM_v1.nc`](https://www.chadagreene.com/tide_data/Gr1kmTm_v1.nc) [193 MB] 1 km resolution, 8 constituent regional model surrounding Greenland, from Susan Howard and Laurie Padman. The Greenland 1 kilometer Tide Model (Gr1kmTM) is a barotropic ocean tide model on a 1 km x 1 km polar stereographic grid, developed using the Regional Ocean Modeling System (ROMS). Gr1kmTM consists of spatial grids of complex amplitude coefficients for sea surface height and depth-integrated currents (“volume transports”) for 8 principal tidal constituents: 4 semidiurnal (M2, S2, K2, N2) and 4 diurnal (K1, O1, P1, Q1). More information on Gr1kmTM_v1 can be found [here](https://doi.org/10.18739/A2B853K18). 
* [`Arc2kmTM_v1.nc`](https://www.chadagreene.com/tide_data/Arc2kmTm_v1.nc) [206 MB] 2 km resolution, 8 constituent regional model of the Arctic Ocean by Susan Howard and Laurie Padman. The Arctic 2 km (kilometer) Tide Model (Arc2kmTM) is a barotropic ocean tide model on a 2x2 km polar stereographic grid, developed using the Regional Ocean Modeling System (ROMS). Arc2kmTM consists of spatial grids of complex amplitude coefficients for sea surface height and depth-integrated currents (“volume transports”) for 8 principal tidal constituents: 4 semidiurnal (M2, S2, K2, N2) and 4 diurnal (K1, O1, P1, Q1). More information on Arc2kmTM_v1 can be found [here](https://doi.org/10.18739/A2PV6B79W). 

# Load Tides
🚧 Check back soon for load tides.🚧

# Tide Model Intercomparison
For a quick intercomparison of tide model performance, check out [this intercomparison](doc/tide_model_intercomparison.md) page.

# Author Info
This page was written by [Chad A. Greene](https://www.chadagreene.com), June 2022. 