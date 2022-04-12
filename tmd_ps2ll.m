function [lat,lon] = tmd_ps2ll(X,Y,SLAT,SLON,HEMI)
% tmd_ps2ll converts geographic coordinates to polar stereographic km.
% This is a simple renaming of the old TMD function mapxy.m. 
% 
% This function is equivalent to the AMT function ps2ll, but accepts
% kilometers instead of meters, and is written for single precision (even
% overwrites the pi variable to single precision!) so you may find that
% tmd_ps2ll and AMT's ps2ll produce differ by up to a few millimeters. 
% 
% The original description is as follows: 
%                                                                         
%                                                                         
%    DESCRIPTION:                                                         
%                                                                         
%    This subroutine converts from Polar Stereographic (X,Y) coordinates  
%    to geodetic latitude and longitude for the polar regions. The
%    equations are from Snyder, J. P., 1982,  Map Projections Used by the U.S.       
%    Geological Survey, Geological Survey Bulletin 1532, U.S. Government   
%    Printing Office.  See JPL Technical Memorandum 3349-85-101 for further
%    details.                                                              
%                                                                          
%                                                                          
%    ARGUMENTS:                                                           
%                                                                         
%    Variable       I/O    Description                          
%                                                                        
%    X               I     Polar Stereographic X Coordinate (km) 
%    Y               I     Polar Stereographic Y Coordinate (km)
%    SLAT            I     Standard latitude
%    SLON            I     Standard longitude
%    HEMI            I     Hemisphere (char*1, 'S' or 'N', 
%                                      not case-sensitive)
%    lat             O     Geodetic Latitude (degrees, +90 to -90)
%    lon             O     Geodetic Longitude (degrees, 0 to 360) 
%                                                                          
% FORTRAN HISTORY                                                                        
%    Written by C. S. Morris - April 29, 1985               
%    Revised by C. S. Morris - December 11, 1985            
%    Revised by V. J. Troisi - January 1990
%      SGN - provide hemisphere dependency (+/- 1)
% MATLAB HISTORY
%    Converted from FORTRAN to Matlab     by L. Padman - 25-Oct-2006
%    Updated for SLON                     by L. Padman - 21-Nov-2006
%    Updated to keep lon between +/-180   by L. Padman - 01-Nov-2007
%    Updated w/ logical indexing          by C. Greene - 31-Jul-2018
%    Renamed tmd_ps2ll                    by C. Greene - 11-Apr-2022
%
%   Sample call:
%                [lat,lon]=tmd_ps2ll(X,Y,SLAT,SLON,HEMI);
%
%*************************************************************************

%                                                                         
%    DEFINITION OF CONSTANTS:                                             
%                                                                         
%    Conversion constant from degrees to radians = 57.29577951.           
      CDR= 57.29577951;
      E2 = 6.694379852e-3;           % Eccentricity squared
      E  = sqrt(E2);
      pi = 3.141592654;
      %RE = 6378.273;                % Original value
      RE=6378.1370;                  % Updated 2/11/08 (see email from
      %                                  Shad O'Neel
%                                                                         
%*************************************************************************
   
if(upper(HEMI)=='S')
    SGN=-1;
else
    SGN=+1;
end
if(upper(HEMI)=='N')
   Y=-Y;
end
SLAT=abs(SLAT);
SL  = SLAT/CDR;
RHO = sqrt(X.^2+Y.^2);
if(RHO<0.1)           % Don't calculate if on the equator
    lat=90.*SGN;
    lon=0.0;
    return
else
    CM=cos(SL)./sqrt(1.0-E2.*(sin(SL).^2));
    T=tan((pi/4.0)-(SL/2))./((1.0-E.*sin(SL))./(1.0+E.*sin(SL))).^(E/2.0);
    if(abs(SLAT-90.)<1.e-5)
        T=RHO*sqrt((1.+E).^(1+E)*(1-E).^(1-E))/2/RE;
    else
        T=RHO.*T./(RE.*CM);
    end
    a1 =  5*E2^2 / 24;
    a2 =    E2^3 / 12;
    a3 =  7*E2^2 / 48;
    a4 = 29*E2^3 /240;
    a5 =  7*E2^3 /120;

    CHI= (pi/2)-2*atan(T);
    lat= CHI+((E2/2) + a1 + a2).*sin(2*CHI)+(a3 + a4).*sin(4*CHI)+ ...
                              a5*sin(6*CHI);
    lat= SGN*lat*CDR;
    %lon= SGN*(atan2(SGN*X,-SGN*Y)*CDR)+SLON;  %Original
    lon= -(atan2(-X,Y)*CDR)+SLON;
    lon(lon<-180)=lon(lon<-180)+360;
    lon(lon>+180)=lon(lon>+180)-360;
end

end
