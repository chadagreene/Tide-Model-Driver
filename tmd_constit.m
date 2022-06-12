function [ispec,amp,ph,omega,alpha] = tmd_constit(constituents)
% tmd_constit returns amplitude, phase, frequency, alpha, species for tidal 
% constituents. 
% 
%% Syntax
% 
%  [ispec,amp,ph,omega,alpha] = constit(constituents)
% 
%% Description 
% 
% [ispec,amp,ph,omega,alpha] = constit(constituents) returns constants for tidal 
% constituents, specified by a string (e.g., 'm2') or a cell array of
% strings (e.g., {'m2','s2'}). 
% 
%    ispec: species type (spherical harmonic dependence of quadropole
%    potential). Species are 0 (long period and overtides such as M4), 1 
%    (diurnal), 2(semidiurnal)
% 
%    amp: amplitude of equilibrium tide in m for tidal constituent.
%
%    phase: (radians) Astronomical arguments (relative to t0 = 1 Jan 0:00 1992].
%
%    omega: frequency (s^-1).
%
%    alpha: loading love number ... frequncy dependance here is suspect.
% 
%% Example
% Get constants for the k1 constituent: 
% 
% [ispec,amp,ph,omega,alpha] = tmd_constit('k1')
% ispec =
%      1
% amp =
%     0.1416
% ph =
%     0.1730
% omega =
%    7.2921e-05
% alpha =
%     0.7360
% 
%% References 
% 
% Egbert and Erofeeva, "Efficient Inverse Modeling of Barotropic Ocean Tides", 
% Journal of Atmospheric and Oceanic Technology, (2002).
%
%% Author Info
% TMD was written by Lana Erofeeva (serofeev@coas.oregonstate.edu) and is maintained by
% Laurie Padman (padman@esr.org). 
% 
% In February 2016, Chad A. Greene modified constit to:
% increase speed and readability. Changes made by Chad include: 
%    * Changed c_data from char array to cell array. 
%    * Previously a loop was used to match c with its corresponding index in c_data --the loop has been replaced by strcmpi.
%
% In June 2022, Chad Greene: 
%    * changed name to tmd_constit
%    * removed constitNum output (because no one knew what it was.)
%    * added some constituents to match constit.h. 
%

%% Define constants: 

c_data = { 'm2';  's2';  'k1';  'o1';
           'n2';  'p1';  'k2';  'q1';
          '2n2'; 'mu2'; 'nu2';  'l2';
           't2';  'j1';  'm1'; 'oo1';
         'rho1';  'mf';  'mm'; 'ssa';
           'm4'; 'ms4'; 'mn4';  'm6';
           'm8'; 'mk3';  's6';'2sm2';
         '2mk3';  's1'; '2q1';  'm3'};
 
% ispec: species type (spherical harmonic dependence of quadropole potential)
ispec_data = [2;2;1;1; 
              2;1;2;1; 
              2;2;2;2; 
              2;1;1;1; 
              1;0;0;0; 
              0;0;0;0; 
              0;0;0;0;
              0;1;1;3]; 
           
% amp: amplitudes
amp_data = [0.2441;0.112743;0.141565;0.100661;
          0.046397;0.046848;0.030684;0.019273;
          0.006141;0.007408;0.008811;0.006931;
          0.006608;0.007915;0.007915;0.004338;
          0.003661;0.042041;0.022191;0.019567;
          0;0;0;0;
          0;0;0;0;
          0;7.6464e-04;0.002565;0.003192];

% phase: Astronomical arguments (relative to t0 = 1 Jan 0:00 1992]
% Richard Ray subs: "arguments" and "astrol"
phase_data = [1.731557546;0.000000000;0.173003674;1.558553872;
              6.050721243;6.110181633;3.487600001;5.877717569;
              4.086699633;3.463115091;5.427136701;0.553986502;
              0.052841931;2.137025284;2.436575100;1.929046130;
              5.254133027;1.756042456;1.964021610;3.487600001;
              3.463115091;1.731557546;1.499093481;5.194672637;
              6.926230184;1.904561220;0.000000000;4.551627762;
              3.809122439;          0;   3.913707;5.738991];

% omega: frequencies
omega_data = [1.405189e-04;1.454441e-04;7.292117e-05;6.759774e-05;
              1.378797e-04;7.252295e-05;1.458423e-04;6.495854e-05;
              1.352405e-04;1.355937e-04;1.382329e-04;1.431581e-04;
              1.452450e-04;7.556036e-05;7.028195e-05;7.824458e-05;
              6.531174e-05;0.053234e-04;0.026392e-04;0.003982e-04;
              2.810377e-04;2.859630e-04;2.783984e-04;4.215566e-04;
              5.620755e-04;2.134402e-04;4.363323e-04;1.503693e-04;
              2.081166e-04; 7.2722e-05;0.6231934E-04;2.107783523e-04];

% alpha: loading love number ... frequncy dependance here is suspect
alpha_data = [0.693;0.693;0.736;0.695;
              0.693;0.706;0.693;0.695;
              0.693;0.693;0.693;0.693;
              0.693;0.695;0.695;0.695;
              0.695;0.693;0.693;0.693;
              0.693;0.693;0.693;0.693;
              0.693;0.693;0.693;0.693;
              0.693;0.693;0.693;0.802];

% % constitNum: why are constituents 9-23 zero?            
% constitNum_data = [1,2,5,6,...
%                    3,7,4,8,...
%                    0,0,0,0,...
%                    0,0,0,0,...
%                    0,0,0,0,...
%                    0,0,0,0,...
%                    0,0,0,0,...
%                    0,0,0,0];
%                  
     
%% Return only constants associated with constituent c: 

[~,kk] = ismember(constituents,c_data);

% Define outputs: 
if any(kk)
   ispec = ispec_data(kk);
   amp = amp_data(kk);
   ph = phase_data(kk);
   omega = omega_data(kk);
   alpha = alpha_data(kk);
else
   ispec = -1;
   amp=0;
   ph=0;
   omega=0;
   alpha=0;
end

end

