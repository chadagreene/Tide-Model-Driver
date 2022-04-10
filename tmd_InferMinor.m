function dh = tmd_InferMinor(hc,constituents,t,aS,aH,aP,aN)
% InferMinor returns correction for 16 minor tidal constiuents.  Zeros 
% are returned if not enough input for inference. This is based on Richard
% Ray's code perth2. 
% 
%% Syntax
% 
%  dh = tmd_InferMinor(hc,constituents,t,aS,aH,aP,aN)
% 
%% Description 
% 
% dh = tmd_InferMinor(hc,constituents,t,aS,aH,aP,aN) returns the correction dh given 
% complex hc for GIVEN constituents/points, constituents are in cell format, 
% t is serial date in Matlab's datenum or datetime format. The variables aS,aH,aP,aN 
% are the outputs of tmd_astrol, and must have the same dimensions as t. 
% 
%% Function history
% This function was written in Matlab by Lana Erofeeva, Oct 2004, based on 
% Richard Ray's perth2 code.  In January 2016, Chad A. Greene made the following 
% changes: 
% * Added documentation in function header. 
% * Removed variables which go unused. 
% * Replaced summation loop with a simple call to sum. 
% * Replaced function return with function end. 
% * Replaced if-statement error checks with assertions. 
% * Added input checks. 
% * Replaced nested loops previously used to create z8 with ismember. 
% In April 2022, Chad Greene changed the inputs for speed and readability.
% 
% See also tmd_predict and tmd_astrol.

%% Initial error checks: 

narginchk(7,7) 
assert(isnumeric(hc)==1,'InferMinor input error: hc must be numeric.') 
assert(iscell(constituents),'InferMinor input error: cid must be a cell array of constituents.') 
assert(isequal(size(t),size(aS),size(aH),size(aP),size(aN)),'Time and astrol inputs must all have the same dimensions.')

%% Check sizes: 

M = size(hc,2); 
N = numel(t); 
if N>1
   assert(M==1 | M==N,'Error: If neither input times nor location are scalar, they must be vectors of the same length.')
   M=N; 
end

%% Create data based on input dimensions: 

if isdatetime(t)
   t = datenum(t); 
end

t_1992 = t(:)-datenum(1992,1,1); % time (days) relatively Jan 1 1992 GMT corresponds to 48622mjd
t_hr = (t_1992 - floor(t_1992))*24;
t1 = 15.*t_hr;
t2 = 30.*t_hr;

cid8 = {'q1';'o1';'p1';'k1';'n2';'m2';'s2';'k2'};

% Reorder zmaj to correspond to cid8: 
[Lia,Lib] = ismember(cid8,constituents); 
z8 = hc(Lib,:); 

assert(sum(Lia)>5,'Not enough constituents for inference.')

%% Solve constants: 

rad=pi/180;
PP=282.8;

zmin=zeros(18,M);   
zmin(1,:)  = 0.263 *z8(1,:) - 0.0252*z8(2,:);   %2Q1
zmin(2,:)  = 0.297 *z8(1,:) - 0.0264*z8(2,:);   %sigma1
zmin(3,:)  = 0.164 *z8(1,:) + 0.0048*z8(2,:);   %rho1 +
zmin(4,:)  = 0.0140*z8(2,:) + 0.0101*z8(4,:);   %M1
zmin(5,:)  = 0.0389*z8(2,:) + 0.0282*z8(4,:);   %M1
zmin(6,:)  = 0.0064*z8(2,:) + 0.0060*z8(4,:);   %chi1
zmin(7,:)  = 0.0030*z8(2,:) + 0.0171*z8(4,:);   %pi1
zmin(8,:)  =-0.0015*z8(2,:) + 0.0152*z8(4,:);   %phi1
zmin(9,:)  =-0.0065*z8(2,:) + 0.0155*z8(4,:);   %theta1
zmin(10,:) =-0.0389*z8(2,:) + 0.0836*z8(4,:);   %J1 +
zmin(11,:) =-0.0431*z8(2,:) + 0.0613*z8(4,:);   %OO1 +
zmin(12,:) = 0.264 *z8(5,:) - 0.0253*z8(6,:);   %2N2 +
zmin(13,:) = 0.298 *z8(5,:) - 0.0264*z8(6,:);   %mu2 +
zmin(14,:) = 0.165 *z8(5,:) + 0.00487*z8(6,:);  %nu2 +
zmin(15,:) = 0.0040*z8(6,:) + 0.0074*z8(7,:);   %lambda2
zmin(16,:) = 0.0131*z8(6,:) + 0.0326*z8(7,:);   %L2 +
zmin(17,:) = 0.0033*z8(6,:) + 0.0082*z8(7,:);   %L2 +
zmin(18,:) = 0.0585*z8(7,:);                      %t2 + 

arg=zeros(18,M);
arg(1,:) = t1 - 4.*aS + aH + 2.*aP - 90;     % 2Q1
arg(2,:) = t1 - 4.*aS + 3.*aH - 90;         % sigma1
arg(3,:) = t1 - 3.*aS + 3.*aH - aP - 90;     % rho1
arg(4,:) = t1 - aS + aH - aP + 90;           % M1
arg(5,:) = t1 - aS + aH + aP + 90;           % M1
arg(6,:) = t1 - aS + 3.*aH - aP + 90;        % chi1
arg(7,:) = t1 - 2.*aH + PP - 90;           % pi1
arg(8,:) = t1 + 3.*aH + 90;                % phi1
arg(9,:) = t1 + aS - aH + aP + 90;           % theta1
arg(10,:) = t1 + aS + aH - aP + 90;          % J1
arg(11,:) = t1 + 2.*aS + aH + 90;           % OO1
arg(12,:) = t2 - 4.*aS + 2.*aH + 2.*aP;       % 2N2
arg(13,:) = t2 - 4.*aS + 4.*aH;              % mu2
arg(14,:) = t2 - 3.*aS + 4.*aH - aP;          % nu2
arg(15,:) = t2 - aS + aP + 180;           % lambda2
arg(16,:) = t2 - aS + 2.*aH - aP + 180;    % L2
arg(17,:) = t2 - aS + 2.*aH + aP;             % L2
arg(18,:) = t2 - aH + PP;                   % t2

% determine nodal corrections f and u:
sinn = sin(aN*rad);
cosn = cos(aN*rad);
sin2n = sin(2.*aN*rad);
cos2n = cos(2.*aN*rad);

f = ones(18,M);
f(1,:) = hypot(1.0 + 0.189*cosn - 0.0058*cos2n, 0.189*sinn - 0.0058*sin2n);
f(2,:) = f(1,:);
f(3,:) = f(1,:); % was mistakenly f(3,:) = f(1); in 2.04 
f(4,:) = hypot(1.0 + 0.185*cosn, 0.185*sinn);
f(5,:) = hypot(1.0 + 0.201*cosn, 0.201*sinn);
f(6,:) = hypot(1.0 + 0.221*cosn, 0.221*sinn);
f(10,:) = hypot(1.0 + 0.198*cosn, 0.198*sinn);
f(11,:) = hypot(1.0 + 0.640*cosn + 0.134*cos2n,0.640*sinn + 0.134*sin2n);
f(12,:) = hypot(1.0 - 0.0373*cosn,0.0373*sinn);
f(13,:) = f(12,:);
f(14,:) = f(12,:);
f(16,:) = f(12,:);
f(17,:) = hypot(1.0 + 0.441*cosn,0.441*sinn);

u = zeros(18,M);
u(1,:) = atan2(0.189*sinn - 0.0058*sin2n,...
             1.0 + 0.189*cosn - 0.0058*sin2n)/rad;
u(2,:) = u(1,:);
u(3,:) = u(1,:);
u(4,:) = atan2( 0.185*sinn, 1.0 + 0.185*cosn)/rad;
u(5,:) = atan2(-0.201*sinn, 1.0 + 0.201*cosn)/rad;
u(6,:) = atan2(-0.221*sinn, 1.0 + 0.221*cosn)/rad;
u(10,:) = atan2(-0.198*sinn, 1.0 + 0.198*cosn)/rad;
u(11,:) = atan2(-0.640*sinn - 0.134*sin2n,...
              1.0 + 0.640*cosn + 0.134*cos2n)/rad;
u(12,:) = atan2(-0.0373*sinn, 1.0 - 0.0373*cosn)/rad;
u(13,:) = u(12,:);
u(14,:) = u(12,:);
u(16,:) = u(12,:);
u(17,:) = atan2(-0.441*sinn, 1.0 + 0.441*cosn)/rad;

%% Sum over all tides: 

tmp = (arg+u)*rad; % This tmp variable prevents performing the same operation twice on the next line.
dh = permute(sum(real(zmin).*f.*cos(tmp) + imag(zmin).*f.*sin(tmp)),[2 3 1]); 

end

