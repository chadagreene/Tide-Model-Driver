function [pu,pf] = tmd_nodal(constituents,p,N)
% This is a Matlab remake of ARGUMENTS by Lana Erofeeva, Jan 2003.
% ARGUMENTS and ASTROL FORTRAN subroutines SUPPLIED by RICHARD RAY, March 1999.
% NOTE - "no1" in constit.h corresponds to "M1" in arguments.
% 
% This function is called by tmd_harp. 
% 
%% Syntax 
% 
%  [pu,pf] = tmd_nodal(constituents,p,N) 
%
%% Description 
% 
% [pu,pf] = tmd_nodal(constituents,p,N) takes input constituents 
% as cell array (1xN constituents). Inputs p and N are the lunar perigee p 
% and ascending lunar node N given by the tmd_astrol function, and are 
% dimesnions Mx1, with one row per timestep. Nodal correction outputs pu,pf 
% have units of radians and are MxN dimensions, corresponding to M 
% timesteps and N constituents. 
% 
%% Author Info
% This function is part of the Tide Model Driver (TMD), which was written by Lana Erofeeva
% and is maintained by Laurie Padman. This function is a Matlab remake of ARGUMENTS by Lana 
% Erofeeva, Jan 2003. ARGUMENTS and ASTROL FORTRAN subroutines SUPPLIED by RICHARD RAY, March 1999.
% 
% In February 2016, Chad A. Greene modified nodal for computational efficiency and readability. 
%    * Function return replaced by function end. 
%    * Removed several dozen lines of unused variables. 
%    * Removed several loops which were used for indexing; replaced with logicals and ismember.  
% 
% In April 2022, Chad Greene rewrote this as tmd_nodal, which replaces both the
% old nodal.m and nodal1.m, but note that the inputs to tmd_nodal are different 
% from the original nodal.m function. Now requires p and N as inputs, to
% prevent calling tmd_astrol every time tmd_nodal is called. 
% 
% See also tmd_harp. 
     
cid0 = {'sa';'ssa';'mm';'msf';'mf';'mt';'alpha1';'2q1';'sigma1';'q1';
        'rho1';'o1';'tau1';'m1';'chi1';'pi1';'p1';'s1';'k1';'psi1';'phi1';
        'theta1';'j1';'oo1';'2n2';'mu2';'n2';'nu2';'m2a';'m2';'m2b';'lambda2';
        'l2';'t2';'s2';'r2';'k2';'eta2';'mns2';'2sm2';'m3';'mk3';'s3';'mn4';
        'm4';'ms4';'mk4';'s4';'s5';'m6';'s6';'s7';'s8'};
 
% Determine equilibrium arguments
rad=pi/180;

nT = length(p);

% Determine nodal corrections f and u 
sinn = sin(N*rad);
cosn = cos(N*rad);
sin2n = sin(2*N*rad);
cos2n = cos(2*N*rad);
sin3n = sin(3*N*rad);

%% Define constants

f=ones(nT,53);
% f(:,1) = 1;                                     % Sa
% f(:,2) = 1;                                     % Ssa
f(:,3) = 1 - 0.130*cosn;                        % Mm
% f(:,4) = 1;                                     % MSf
f(:,5) = 1.043 + 0.414*cosn;                    % Mf
f(:,6) = sqrt((1+.203*cosn+.040*cos2n).^2 + ...
              (.203*sinn+.040*sin2n).^2);        % Mt

% f(:,7) = 1;                                     % alpha1
f(:,8) = sqrt((1.+.188*cosn).^2+(.188*sinn).^2);% 2Q1
f(:,9) = f(:,8);                                % sigma1
f(:,10) = f(:,8);                               % q1
f(:,11) = f(:,8);                               % rho1
f(:,12) = sqrt((1.0+0.189*cosn-0.0058*cos2n).^2 + ...
                 (0.189*sinn-0.0058*sin2n).^2);% O1
f(:, 13) = 1;                                   % tau1
% tmp1  = 2.*cos(astrol_p*rad)+.4*cos((astrol_p-astrol_N)*rad);
% tmp2  = sin(astrol_p*rad)+.2*sin((astrol_p-astrol_N)*rad);% Doodson's
tmp1  = 1.36*cos(p*rad)+.267*cos((p-N)*rad);% Ray's
tmp2  = 0.64*sin(p*rad)+.135*sin((p-N)*rad);
f(:,14) = hypot(tmp1,tmp2);                % M1
f(:,15) = sqrt((1.+.221*cosn).^2+(.221*sinn).^2);% chi1
% f(:,16) = 1;                                    % pi1
% f(:,17) = 1;                                    % P1
% f(:,18) = 1;                                    % S1
f(:,19) = sqrt((1.+.1158*cosn-.0029*cos2n).^2 + ...
                (.1554*sinn-.0029*sin2n).^2);  % K1
% f(:,20) = 1;                                    % psi1
% f(:,21) = 1;                                    % phi1
% f(:,22) = 1;                                    % theta1
f(:,23) = sqrt((1.+.169*cosn).^2+(.227*sinn).^2); % J1
f(:,24) = sqrt((1.0+0.640*cosn+0.134*cos2n).^2 + ...
                (0.640*sinn+0.134*sin2n).^2 ); % OO1
f(:,25) = sqrt((1.-.03731*cosn+.00052*cos2n).^2 + ...
                (.03731*sinn-.00052*sin2n).^2);% 2N2
f(:,26) = f(:,25);                                % mu2
f(:,27) = f(:,25);                                % N2
f(:,28) = f(:,25);                                % nu2
% f(:,29) = 1;                                    % M2a
f(:,30) = f(:,25);                                % M2
% f(:,31) = 1;                                    % M2b
% f(:,32) = 1;                                    % lambda2
temp1 = 1.-0.25*cos(2*p*rad)-0.11*cos((2*p-N)*rad)-0.04*cosn;
temp2 = 0.25*sin(2*p*rad)+0.11*sin((2*p-N)*rad)+ 0.04*sinn;
f(:,33) = sqrt(temp1.^2 + temp2.^2);              % L2
% f(:,34) = 1;                                    % T2
% f(:,35) = 1;                                    % S2
% f(:,36) = 1;                                    % R2
f(:,37) = sqrt((1.+.2852*cosn+.0324*cos2n).^2 + ...
                (.3108*sinn+.0324*sin2n).^2);  % K2
f(:,38) = sqrt((1.+.436*cosn).^2+(.436*sinn).^2); % eta2
f(:,39) = f(:,30).^2;                            % MNS2
f(:,40) = f(:,30);                              % 2SM2
% f(:,41) = 1;   % wrong                          % M3
f(:,42) = f(:,19).*f(:,30);                     % MK3
% f(:,43) = 1;                                    % S3
f(:,44) = f(:,30).^2;                           % MN4
f(:,45) = f(:,44);                              % M4
f(:,46) = f(:,30);                              % MS4
f(:,47) = f(:,30).*f(:,37);                     % MK4
% f(:,48) = 1;                                    % S4
% f(:,49) = 1;                                    % S5
f(:,50) = f(:,30).^3;                           % M6
% f(:,51) = 1;                                    % S6
% f(:,52) = 1;                                    % S7
% f(:,53) = 1;                                    % S8

u=zeros(nT,53);
% u(:, 1) = 0;                                       % Sa
% u(:, 2) = 0;                                       % Ssa
% u(:, 3) = 0;                                       % Mm
% u(:, 4) = 0;                                       % MSf
u(:, 5) = (-23.7*sinn + 2.7*sin2n - 0.4*sin3n)*rad;      % Mf
u(:, 6) = atan(-(.203*sinn+.040*sin2n)./...
             (1+.203*cosn+.040*cos2n));        % Mt
% u(:, 7) = 0;                                       % alpha1
u(:, 8) = atan(.189*sinn./(1.+.189*cosn));     % 2Q1
u(:, 9) = u(:,8);                                  % sigma1
u(:,10) = u(:,8);                                  % q1
u(:,11) = u(:,8);                                  % rho1
u(:,12) = (10.8*sinn - 1.3*sin2n + 0.2*sin3n)*rad;       % O1
% u(:,13) = 0;                                       % tau1
u(:,14) = atan2(tmp2,tmp1);                    % M1
u(:,15) = atan(-.221*sinn./(1.+.221*cosn));    % chi1
% u(:,16) = 0;                                       % pi1
% u(:,17) = 0;                                       % P1
% u(:,18) = 0;                                       % S1
u(:,19) = atan((-.1554*sinn+.0029*sin2n)./...
           (1.+.1158*cosn-.0029*cos2n));       % K1
% u(:,20) = 0;                                       % psi1
% u(:,21) = 0;                                       % phi1
% u(:,22) = 0;                                       % theta1
u(:,23) = atan(-.227*sinn./(1.+.169*cosn));    % J1
u(:,24) = atan(-(.640*sinn+.134*sin2n)./...
           (1.+.640*cosn+.134*cos2n));         % OO1
u(:,25) = atan((-.03731*sinn+.00052*sin2n)./ ...
           (1.-.03731*cosn+.00052*cos2n));     % 2N2
u(:,26) = u(:,25);                                 % mu2
u(:,27) = u(:,25);                                 % N2
u(:,28) = u(:,25);                                 % nu2
% u(:,29) = 0;                                       % M2a
u(:,30) = u(:,25);                                   % M2
% u(:,31) = 0;                                       % M2b
% u(:,32) = 0;                                       % lambda2
u(:,33) = atan(-temp2./temp1) ;                % L2
% u(:,34) = 0;                                       % T2
% u(:,35) = 0;                                       % S2
% u(:,36) = 0;                                       % R2
u(:,37) = atan(-(.3108*sinn+.0324*sin2n)./ ...
             (1.+.2852*cosn+.0324*cos2n));     % K2
u(:,38) = atan(-.436*sinn./(1.+.436*cosn));    % eta2
u(:,39) = u(:,30)*2;                               % MNS2
u(:,40) = u(:,30);                                 % 2SM2
u(:,41) = 1.5d0*u(:,30);                           % M3
u(:,42) = u(:,30) + u(:,19);                       % MK3
% u(:,43) = 0;                                       % S3
u(:,44) = u(:,30)*2;                               % MN4
u(:,45) = u(:,44);                                 % M4
u(:,46) = u(:,30);                                 % MS4
u(:,47) = u(:,30)+u(:,37);                         % MK4
% u(:,48) = 0;                                       % S4
% u(:,49) = 0;                                       % S5
u(:,50) = u(:,30)*3;                               % M6
% u(:,51) = 0;                                       % S6
% u(:,52) = 0;                                       % S7
% u(:,53) = 0;                                       % S8

%% Trim to only the user-requested constituents

[~,Lib] = ismember(constituents,cid0); 

pu = u(:,Lib); 
pf = f(:,Lib); 

end

