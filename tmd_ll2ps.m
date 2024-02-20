function [x,y] = tmd_ll2ps(lat,lon,SLAT,SLON,HEMI)
% tmd_ll2ps converts geographic coordinates to polar stereographic km.
% This is a simple renaming of the old TMD function mapll.m. 
% 
% This function is equivalent to the AMT function ll2ps, but returns
% kilometers instead of meters, and is written for single precision (even
% overwrites the pi variable to single precision!) so you may find that
% tmd_ll2ps and AMT's ll2ps produce results that differ by up to a few
% millimeters. 
% 
% The original description is as follows: 
%                                                                         
%    DESCRIPTION:
%                                                                         
%    This function converts from geodetic latitude and longitude to Polar 
%    Stereographic (X,Y) coordinates for the polar regions.  The equations 
%    are from Snyder, J. P., 1982,  Map Projections Used by the U.S.      
%    Geological Survey, Geological Survey Bulletin 1532, U.S. Government  
%    Printing Office.  See JPL Technical Memorandum 3349-85-101 for further
%    details.                                                             
%                                                                       
%    ARGUMENTS:                                                         
%                                                                       
%    Variable     I/O    Description                          
%                                                                        
%    lat           I     Geodetic Latitude (degrees, +90 to -90) Any input dimensions are allowed.  
%    lon           I     Geodetic Longitude (degrees, 0 to 360) Any input dimensions are allowed.  
%    SLAT          I     Standard latitude (typ. 71, or 70)
%    SLON          I  
%    HEMI          I     Hemisphere (char*1: 'N' or 'S' (not
%                                    case-sensitive)
%    x             O     Polar Stereographic X Coordinate (km)
%    y             O     Polar Stereographic Y Coordinate (km)
%                                                                      
%
% FORTRAN CODE HISTORY
%    Written by C. S. Morris - April 29, 1985             
%    Revised by C. S. Morris - December 11, 1985          
%    Revised by V. J. Troisi - January 1990               
%       SGN - provides hemisphere dependency (+/- 1)         
%    Revised by Xiaoming Li - October 1996                     
%		Corrected equation for RHO 
%
%    Converted from FORTRAN to Matlab by L. Padman - 25-Oct-2006
%    Updated for SLON                 by L. Padman - 21-Nov-2006
% 
%    Renamed from mapll.m to tmd_ll2ps.m by C. Greene - April 2022. 
%
% Sample call:
%             [x,y] = tmd_ll2ps(lat,lon,SLAT,SLON,HEMI);
%
%*************************************************************************
%                                                                        
%    DEFINITION OF CONSTANTS:                                            
%                                                                        
%    Conversion constant from degrees to radians = 57.29577951.          
      CDR=57.29577951;
      E2 = 6.694379852e-3;           % Eccentricity squared
      E  = sqrt(E2);
      pi = 3.141592654;
      %RE = 6378.273;                % Original value
      RE=6378.1370;                  % Updated 2/11/08 (see email from
      %                                  Shad O'Neel
%                                                                        
%*************************************************************************

if(abs(SLAT)==90)
    RHO=2*RE/((1+E).^(1+E)*(1-E).^(1-E)).^(E/2);
else
    SL  = abs(SLAT)/CDR;
    TC  = tan(pi/4-SL/2)/((1-E*sin(SL))./(1+E*sin(SL))).^(E/2);
    MC  = cos(SL)/sqrt(1-E2*(sin(SL).^2));
    RHO = RE.*MC./TC;
end
lat = abs(lat)/CDR;
T   = tan(pi/4-lat/2)./((1-E*sin(lat))./(1+E*sin(lat))).^(E/2);
lon =-(lon-SLON)/CDR;
x   =-RHO.*T.*sin(lon);
y   = RHO.*T.*cos(lon);

if(upper(HEMI)=='N')
   y=-y; 
end

end
