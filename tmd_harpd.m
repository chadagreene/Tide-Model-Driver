function hhat = tmd_harpd(t, hc, constituents, tt_ut1=0.0007)
% tmd_harpd predicts tidal time series using harmonic constants. 
% Nodal corrections are included. 
% 
%% Syntax
% 
%  hhat = tmd_harpd(t,hc,constituents)
% 
%% Description 
% 
% hhat = tmd_harpd(t,hc,constituents,p,N) returns the time series hhat, which is
% the sum of tidal constituents given by hc for times t and has the units of hc.
% The frequencies of the equilibrium tides are computed using Doodson coefficients.
% 
% Inputs:
%    t:  time in Matlab's datenum format (size Tx1, where T is # of timesteps),
%    hc: complex harmonic constants (size Cx1, where C is # of constituents), 
%    constituents: cell array, e.g., {'m2','k1'} (size 1xC), 
%    tt_ut1: delta time (in days) between Terrestrial Time (TT) and Universal Time (UT1).
% 
% Outputs: 
%    hhat: reconstructed tide (size Tx1).
%
%% References 
%
% Doodson, A. T. (1941). Admiralty manual of tides. His Majesty's Stationery Office.
%
% Schureman, P. (1958). Manual of harmonic analysis and prediction of tides.
% Special Edition No. 98. US Coast and Geodetic Survey, United States Government
% Printing Office, Washington, DC.
%
% Foreman, M. G. G., & Henry, R. F. (1989). The harmonic analysis of tidal 
% model time series. Advances in water resources, 12(3), 109-120.
% https://doi.org/10.1016/0309-1708(89)90017-1
%
%% Author Info
% TMD was written by Lana Erofeeva (serofeev@coas.oregonstate.edu) and is maintained by
% Laurie Padman (padman@esr.org). 
%
% See also tmd_predict, tmd_astrol, and tmd_harp. 

[pu,pf,G] = tmd_arguments(t, constituents, tt_ut1=tt_ut1);

th = G*pi/180.0 + pu;
hhat = sum(pf.*(real(hc)'.*cos(th) + imag(hc)'.*sin(th)),2);

end
 
