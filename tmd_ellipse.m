function [umajor,uminor,uphase,uincl] = tmd_ellipse(filename,constituent,lati,loni)
% tmd_get_ellipse extracts tidal ellipse grids from a model.  
% 
% usage:
% [x,y,umaj,umin,uphase,uincl]=tmd_get_ellipse(Model,cons);
% 
% Model - control file name for a tidal model, consisting of lines
%         <elevation file name>
%         <transport file name>
%         <grid file name>
%         <function to convert lat,lon to x,y>
% 4th line is given only for models on cartesian grid (in km)
% All model files should be provided in OTIS format
% cons - tidal constituent given as char* 
%
% output:
% umaj,umin - major and minor ellipse axis (cm/s)
% uphase, uincl - ellipse phase and inclination degrees GMT
% x,y - grid coordinates
%
% sample call:
% [x,y,umaj,umin,uphase,uincl]=tmd_get_ellipse('DATA/Model_Ross_prior','k1');
%
% TMD release 2.02: 21 July 2010
% TMD release 3.00: August 2018, changes by Chad Greene include: 
%   - Turned TideEl.m into a subfunction of tmd_get_ellipse.
% 
% 
%

%% Input checks

narginchk(4,4)
assert(isequal(size(lati),size(loni)),'Dimensions of lati,loni must agree.')

%% Load data 

u = tmd_interp(filename,'u',lati,loni,'constituents',constituent); 
v = tmd_interp(filename,'v',lati,loni,'constituents',constituent); 


% TideEl calculates tidal ellipse parameters for the arrays of
% u and v - COMPLEX amplitudes of EW and NS currents of
% a given tidal constituent
% land should be set to 0 or NaN in u,v prior to calling tideEl
% usage: [umajor,uminor,uincl,uphase]=TideEl(u,v);
% Version History 
% TMD2.04: ??
% TMD3.00: Chad Greene, August 2018.
%   - Renamed (case change) function from tideEl to TideEl. 
%   - Changed function return to function end. 
%   - Turned this into a subfunction of tmd_get_ellipse. 
%   

% change to polar coordinates 
% in Robin's was - + + -, this is as in Foreman's
t1p = (real(u) - imag(v));
t2p = (real(v) + imag(u));
t1m = (real(u) + imag(v));
t2m = (real(v) - imag(u));

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

%  determine the major and minor axes, phase and inclination using Foreman's formula 
umajor = (ap + am); 
uminor = (ap - am);
uincl = 0.5 * (em + ep);
uincl = uincl - 180. *  (uincl > 180);
uphase = - 0.5*(ep-em) ;
uphase = uphase + 360. * (uphase < 0);
uphase = uphase - 360. * (uphase >= 360);

end
