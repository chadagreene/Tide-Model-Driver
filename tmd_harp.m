function hhat = tmd_harp(t,hc,constituents,astrol_p,astrol_N)
% harp1 predicts tidal time series using harmonic constants. This is the time series version of harp. 
% Nodal corrections are included. 
% 
% harp1 is called by tmd_tide_pred. 
% 
%% Syntax
% 
% hhat = harp1(time,hc,con)
% 
%% Description 
% 
% hhat = harp1(time,hc,con) returns the time series hhat reconstructed using HC. Inputs are: 
% 
%     t datetime or datenum
%     con(nc,4) - char*4 tidal constituent IDs 
%     hc(nc) - harmonic constant vector  (complex)
% 
%% Author Info
% TMD was written by Lana Erofeeva (serofeev@coas.oregonstate.edu) and is maintained by
% Laurie Padman (padman@esr.org). In February 2016, Chad A. Greene modified harp1 to
% increase speed and readability. Changes made by Chad include: 
% * Removed unused variables. 
% * Replaced function end with function return. 
% * Replaced constit loop with a single call to consit, which is newly re-written to allow multiple constituents. 
% * Rewrote hhat loop to be a little less cryptic. 
% 
% See also harp, nodal, and tmd_tide_pred. 

% if isdatetime(t)
%    t = datenum(t); 
% end

time_s = (t - datenum(1992,1,1))*86400; % time (days since Jan 1, 1992) is converted to seconds by multiplying by 86400 s/yr. 

[ispec,~,ph,omega,~,~] = tmd_constit(constituents);

igood=ispec~=-1;

[pu,pf] = tmd_nodal(constituents(igood),astrol_p,astrol_N); % time (days since Jan 1, 1992) is converted to MJD by adding 48622

tmp = omega(igood)'.*time_s + ph(igood)' + pu; 
hhat = sum(pf.*(real(hc(igood,:))'.*cos(tmp) - imag(hc(igood,:))'.*sin(tmp)),2);


end
 
