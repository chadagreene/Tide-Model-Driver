function hhat = tmd_harp(t,hc,constituents,astrol_p,astrol_N,ispec,ph,omega)
% tmd_harp predicts tidal time series using harmonic constants. 
% Nodal corrections are included. 
% 
%% Syntax
% 
% hhat = tmd_harp(t,hc,constituents,astrol_p,astrol_N)
% 
%% Description 
% 
% hhat = harp(time,hc,con) returns the time series hhat reconstructed using HC.
% Inputs are: 
% 
%     t  time in Matlab's datenum format
%     hc complex harmonic constants  
%     constituents cell array, e.g., {'m2','k1'} 
% 
%% Author Info
% TMD was written by Lana Erofeeva (serofeev@coas.oregonstate.edu) and is maintained by
% Laurie Padman (padman@esr.org). In February 2016, Chad A. Greene modified harp1 to
% increase speed and readability. Changes made by Chad include: 
% * Removed unused variables. 
% * Replaced function end with function return. 
% * Replaced constit loop with a single call to consit, which is newly re-written to allow multiple constituents. 
% * Rewrote hhat loop to be a little less cryptic. 
% In April 2022, Chad Greene rewrote the function to remove loops and
% accept input times in datenum format. 
% 
% See also tmd_predict and tmd_astrol. 

% if isdatetime(t)
%    t = datenum(t); 
% end

time_s = (t - datenum(1992,1,1))*86400; % time (days since Jan 1, 1992) is converted to seconds by multiplying by 86400 s/yr. 

igood=ispec~=-1;

[pu,pf] = tmd_nodal(constituents(igood),astrol_p,astrol_N); % time (days since Jan 1, 1992) is converted to MJD by adding 48622

tmp = omega(igood)'.*time_s + ph(igood)' + pu; 
hhat = sum(pf.*(real(hc(igood,:))'.*cos(tmp) + imag(hc(igood,:))'.*sin(tmp)),2);


end
 
