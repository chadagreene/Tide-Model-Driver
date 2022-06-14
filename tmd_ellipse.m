function [umajor,uminor,uphase,uincl] = tmd_ellipse(filename,constituent,lati,loni)
% tmd_ellipse gives tidal ellipse parameters at specified location(s).  
% 
%% Syntax
% 
%  [umajor,uminor,uphase,uincl] = tmd_ellipse(filename,constituent,lati,loni)
% 
%% Description
% 
% [umajor,uminor,uphase,uincl] = tmd_ellipse(filename,constituent,lati,loni)
% gives the major and minor axes of tidal current velocities. 
%
% Inputs: 
%    filename: TMD3.0 compatible tide model ending in .nc. 
%    constituent: string constituent (e.g., 'm2') 
%    lati,loni: geographic coordinates (can be any size)
% 
% Outputs: 
%    umajor: Major axis, the largest current for the constituent (m/s). Always positive. 
%    uminior: Minor axis, the smallest current (m/s). Negative uminor indicates clockwise flow.  
%    uphase: Reference point on the ellipse (degrees)
%    uincl: Inclination
%
%% Example 
% % 
% fn = 'CATS2008_update_2022-06-12.nc'; 
% lat = -71.9102; 
% lon =  172.6590;
% 
% [umaji,umini,uphasei,uincli] = tmd_ellipse(fn,'o1',lat,lon)
% umaji =
%     0.2060
% umini =
%    -0.0886
% uphasei =
%    36.9514
% uincli =
%    56.3513
% 
% % The above indicates clockwise flow with a max velocity of about 20
% cm/s.
%
%% References 
% 
% Foreman, M. G. G., & Henry, R. F. (1989). The harmonic analysis of tidal 
% model time series. Advances in water resources, 12(3), 109-120.
% https://doi.org/10.1016/0309-1708(89)90017-1
% 
%% Version History 
% TMD release 2.02: 21 July 2010
% 
% June 2022: Chad Greene 
%   * Changed name from tmd_get_ellipse to tmd_ellipse.
%   * Changed behavior to calculate ellipses at user-specified location(s).
% 
% See also: tmd_predict, tmd_interp, and tmd_data. 

%% Input checks

narginchk(4,4)
assert(isequal(size(lati),size(loni)),'Dimensions of lati,loni must agree.')

%% Load data 

u = tmd_interp(filename,'u',lati,loni,'constituents',constituent);
v = tmd_interp(filename,'v',lati,loni,'constituents',constituent); 

%% Calculate ellipses 

% change to polar coordinates 
% Chad Greene swapped +/- to adapt to TMD3.0's complex number convention.  
t1p = (real(u) + imag(v));
t2p = (real(v) - imag(u));
t1m = (real(u) - imag(v));
t2m = (real(v) + imag(u));

% ap, am - amplitudes of positevely and negatively
% rotated vectors
ap = hypot(t1p,t2p)/2;
am = hypot(t1m,t2m)/2;

% ep, em - phases of positively and negatively rotating vectors
ep = atan2( t2p, t1p);
ep = ep + 2 * pi * (ep < 0.0);
ep = 180. * ep / pi;
em = atan2( t2m, t1m);
em = em + 2 * pi * (em < 0.0);
em = 180. * em / pi;

% determine the major and minor axes, phase and inclination using Foreman's formula 
umajor = (ap + am); 
uminor = (ap - am);
uincl = 0.5*(em + ep);
uincl = uincl - 180. *  (uincl > 180);
uphase = 0.5*(em - ep) ;
uphase = uphase + 360. * (uphase < 0);
uphase = uphase - 360. * (uphase >= 360);

end
