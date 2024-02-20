function hhat = tmd_harp(t,hc,constituents,p,N,ph,omega)
% tmd_harp predicts tidal time series using harmonic constants. 
% Nodal corrections are included. 
% 
%% Syntax
% 
%  hhat = tmd_harp(t,hc,constituents,p,N,ph,omega)
% 
%% Description 
% 
% hhat = tmd_harp(t,hc,constituents,p,N) returns the time series hhat, which is
% the sum of tidal constituents given by hc for times t and has the units of hc.
% 
% Inputs:
%    t:  time in Matlab's datenum format (size Tx1, where T is # of timesteps),
%    hc: complex harmonic constants (size Cx1, where C is # of constituents), 
%    constituents: cell array, e.g., {'m2','k1'} (size 1xC), 
%    p: lunar perigee p given by the tmd_astrol function. (size Tx1), 
%    N: ascending lunar node N given by the tmd_astrol function (size Tx1), 
%    ph: astronomical phase given by the tmd_constit function (size Cx1), 
%    omega: frequency given by the tmd_constit function (size Cx1) units s^-1. 
% 
% Outpus: 
%    hhat: reconstructed tide (size Tx1).
% 
%% Author Info
% TMD was written by Lana Erofeeva (serofeev@coas.oregonstate.edu) and is maintained by
% Laurie Padman (padman@esr.org). 
%
% In February 2016, Chad A. Greene modified harp1 to increase speed and 
% readability. Changes made by Chad include: 
%    * Removed unused variables. 
%    * Replaced function end with function return. 
%    * Replaced constit loop with a single call to consit, which is newly re-written to allow multiple constituents. 
%    * Rewrote hhat loop to be a little less cryptic. 
%
% In April 2022, Chad Greene 
%    * rewrote the function to remove loops, 
%    * changed to datenum input format, 
%    * combined harp and harp1 into tmd_harp. 
%    * added p, N, ph, and omega as inputs (so we don't have to call tmd_astrol and tmd_constit every time tmd_harp is called). 
% 
% See also tmd_predict and tmd_astrol. 

time_s = (t - datenum(1992,1,1))*86400; % time (days since Jan 1, 1992) is converted to seconds by multiplying by 86400 s/yr. 

[pu,pf] = tmd_nodal(constituents,p,N);

tmp = omega'.*time_s + ph' + pu; 
hhat = sum(pf.*(real(hc)'.*cos(tmp) + imag(hc)'.*sin(tmp)),2);

% % Previous implementation ignored bad constituents, but in TMD3.0 there should be no bad constituents in the TMD3.0 compatible models: 
% igood=ispec~=-1;
% tmp = omega(igood)'.*time_s + ph(igood)' + pu; 
% hhat = sum(pf.*(real(hc(igood,:))'.*cos(tmp) + imag(hc(igood,:))'.*sin(tmp)),2);

end
 
