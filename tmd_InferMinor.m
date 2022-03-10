function dh = tmd_InferMinor(zmaj,cid,SDtime)
% InferMinor returns correction for 16 minor tidal constiuents.  Zeros 
% are returned if not enough input for inference. This is based on Richard
% Ray's code perth2. 
% 
% Modes:
%      Time series: zmaj(ncon,1),  SDtime(nt,1)  ->Output: dh(nt,1)
%      Drift Track: zmaj(ncon,nt), SDtime(nt,1)  ->Output: dh(nt,1)
%      Map:         zmaj(ncon,N,M),SDtime(1,1)   ->Output: dh(N,M)
% 
% This function is called by tmd_tide_pred. 
% 
%% Syntax
% 
% dh = InferMinor(zmaj,cid,SDtime)
% 
%% Description 
% 
% dh = InferMinor(zmaj,cid,SDtime) returns the correction dh given complex HC zmaj
% for GIVEN constituents/points, cid are the given constituents in Nx4 char format,
% and SDtime is serial date in Matlab's datenum format. 
% 
%% Example 
% In tmd_tide_pred, InferMinor is called by 
% 
% dh = InferMinor(hc,conList,t); 
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
% 
% See also tmd_tide_pred.

%% Initial error checks: 

narginchk(3,3) 
assert(isnumeric(zmaj)==1,'InferMinor input error: zmaj must be numeric.') 
assert(isnumeric(cid)==0,'InferMinor input error: cid must be in Nx4 char format.') 
assert(isnumeric(SDtime)==1,'InferMinor input error: SDtime must be in Matlab''s numeric datenum format. Type doc datenum for help.') 

%% Check sizes: 

assert(size(cid,2)==4,'InferMinor call: Wrong constituents') 

ncon = size(cid,1);
[n1,n2,n3]=size(zmaj);
nt=length(SDtime);

% Try reordering: 
if n1~=ncon
   zmaj=conj(zmaj');
   [n1,n2,n3]=size(zmaj);
   
   % If reordering did not work, throw an error: 
   assert(n1==ncon,'InferMinor call: Wrong zmaj size.');
end

if n2==1 & n3==1 % time series 
   N=nt;
   M=1;
elseif n3==1 & n2>1 & nt==n2 % drift track
   N=nt;
   M=1;
elseif n2>1 & nt==1 % map solution
   N=n2;
   M=n3;
else
   error('InferMinor call: Wrong zmaj/SDtime size\n');
end

%% Create data based on input dimensions: 

time=SDtime-datenum(1992,1,1); % time (days) relatively Jan 1 1992 GMT corresponds to 48622mjd
%time_mjd=48622+time;
cid8 = ['q1  ';'o1  ';'p1  ';'k1  ';'n2  ';'m2  ';'s2  ';'k2  '];

% Reorder zmaj to correspond to cid8: 
[Lia,Lib] = ismember(cid8,cid,'rows'); 
z8 = zmaj(Lib,:,:); 

assert(sum(Lia)>5,'Not enough constituents for inference.')

%% Solve constants: 

rad=pi/180;
PP=282.8;

zmin=zeros(18,N,M);   
zmin(1,:,:)  = 0.263 *z8(1,:,:) - 0.0252*z8(2,:,:);   %2Q1
zmin(2,:,:)  = 0.297 *z8(1,:,:) - 0.0264*z8(2,:,:);   %sigma1
zmin(3,:,:)  = 0.164 *z8(1,:,:) + 0.0048*z8(2,:,:);   %rho1 +
zmin(4,:,:)  = 0.0140*z8(2,:,:) + 0.0101*z8(4,:,:);   %M1
zmin(5,:,:)  = 0.0389*z8(2,:,:) + 0.0282*z8(4,:,:);   %M1
zmin(6,:,:)  = 0.0064*z8(2,:,:) + 0.0060*z8(4,:,:);   %chi1
zmin(7,:,:)  = 0.0030*z8(2,:,:) + 0.0171*z8(4,:,:);   %pi1
zmin(8,:,:)  =-0.0015*z8(2,:,:) + 0.0152*z8(4,:,:);   %phi1
zmin(9,:,:)  =-0.0065*z8(2,:,:) + 0.0155*z8(4,:,:);   %theta1
zmin(10,:,:) =-0.0389*z8(2,:,:) + 0.0836*z8(4,:,:);   %J1 +
zmin(11,:,:) =-0.0431*z8(2,:,:) + 0.0613*z8(4,:,:);   %OO1 +
zmin(12,:,:) = 0.264 *z8(5,:,:) - 0.0253*z8(6,:,:);   %2N2 +
zmin(13,:,:) = 0.298 *z8(5,:,:) - 0.0264*z8(6,:,:);   %mu2 +
zmin(14,:,:) = 0.165 *z8(5,:,:) + 0.00487*z8(6,:,:);  %nu2 +
zmin(15,:,:) = 0.0040*z8(6,:,:) + 0.0074*z8(7,:,:);   %lambda2
zmin(16,:,:) = 0.0131*z8(6,:,:) + 0.0326*z8(7,:,:);   %L2 +
zmin(17,:,:) = 0.0033*z8(6,:,:) + 0.0082*z8(7,:,:);   %L2 +
zmin(18,:,:) = 0.0585*z8(7,:,:);                      %t2 + 

hour = (time - floor(time))*24.D0;
t1 = 15.*hour;
t2 = 30.*hour;
[S,H,P,omega]=tmd_astrol(SDtime); % nsites x nt

arg=zeros(18,N,M);
arg(1,:,:) = t1 - 4.*S + H + 2.*P - 90.;     % 2Q1
arg(2,:,:) = t1 - 4.*S + 3.*H - 90.;         % sigma1
arg(3,:,:) = t1 - 3.*S + 3.*H - P - 90.;     % rho1
arg(4,:,:) = t1 - S + H - P + 90.;           % M1
arg(5,:,:) = t1 - S + H + P + 90.;           % M1
arg(6,:,:) = t1 - S + 3.*H - P + 90.;        % chi1
arg(7,:,:) = t1 - 2.*H + PP - 90.;           % pi1
arg(8,:,:) = t1 + 3.*H + 90.;                % phi1
arg(9,:,:) = t1 + S - H + P + 90.;           % theta1
arg(10,:,:) = t1 + S + H - P + 90.;          % J1
arg(11,:,:) = t1 + 2.*S + H + 90.;           % OO1
arg(12,:,:) = t2 - 4.*S + 2.*H + 2.*P;       % 2N2
arg(13,:,:) = t2 - 4.*S + 4.*H;              % mu2
arg(14,:,:) = t2 - 3.*S + 4.*H - P;          % nu2
arg(15,:,:) = t2 - S + P + 180.D0;           % lambda2
arg(16,:,:) = t2 - S + 2.*H - P + 180.D0;    % L2
arg(17,:,:) = t2 - S + 2.*H + P;             % L2
arg(18,:,:) = t2 - H + PP;                   % t2

% determine nodal corrections f and u:
sinn = sin(omega*rad);
cosn = cos(omega*rad);
sin2n = sin(2.*omega*rad);
cos2n = cos(2.*omega*rad);
f = ones(18,N,M);
f(1,:,:) = sqrt((1.0 + 0.189*cosn - 0.0058*cos2n).^2 +...
                (0.189*sinn - 0.0058*sin2n).^2);
f(2,:,:) = f(1,:,:);
f(3,:,:) = f(1,:,:); % was mistakenly f(3,:,:) = f(1); in 2.04 
f(4,:,:) = sqrt((1.0 + 0.185*cosn).^2 + (0.185*sinn).^2);
f(5,:,:) = sqrt((1.0 + 0.201*cosn).^2 + (0.201*sinn).^2);
f(6,:,:) = sqrt((1.0 + 0.221*cosn).^2 + (0.221*sinn).^2);
f(10,:,:) = sqrt((1.0 + 0.198*cosn).^2 + (0.198*sinn).^2);
f(11,:,:) = sqrt((1.0 + 0.640*cosn + 0.134*cos2n).^2 +...
             (0.640*sinn + 0.134*sin2n).^2 );
f(12,:,:) = sqrt((1.0 - 0.0373*cosn).^2 + (0.0373*sinn).^2);
f(13,:,:) = f(12,:,:);
f(14,:,:) = f(12,:,:);
f(16,:,:) = f(12,:,:);
f(17,:,:) = sqrt((1.0 + 0.441*cosn).^2 + (0.441*sinn).^2);

u = zeros(18,N,M);
u(1,:,:) = atan2(0.189*sinn - 0.0058*sin2n,...
             1.0 + 0.189*cosn - 0.0058*sin2n)/rad;
u(2,:,:) = u(1,:,:);
u(3,:,:) = u(1,:,:);
u(4,:,:) = atan2( 0.185*sinn, 1.0 + 0.185*cosn)/rad;
u(5,:,:) = atan2(-0.201*sinn, 1.0 + 0.201*cosn)/rad;
u(6,:,:) = atan2(-0.221*sinn, 1.0 + 0.221*cosn)/rad;
u(10,:,:) = atan2(-0.198*sinn, 1.0 + 0.198*cosn)/rad;
u(11,:,:) = atan2(-0.640*sinn - 0.134*sin2n,...
              1.0 + 0.640*cosn + 0.134*cos2n)/rad;
u(12,:,:) = atan2(-0.0373*sinn, 1.0 - 0.0373*cosn)/rad;
u(13,:,:) = u(12,:,:);
u(14,:,:) = u(12,:,:);
u(16,:,:) = u(12,:,:);
u(17,:,:) = atan2(-0.441*sinn, 1.0 + 0.441*cosn)/rad;

%% Sum over all tides: 

% permute gets rid of the singleton dimension and is faster than squeeze. 
tmp = (arg+u)*rad; % This tmp variable prevents performing the same operation twice on the next line.
dh = permute(sum(real(zmin).*f.*cos(tmp) - imag(zmin).*f.*sin(tmp)),[2 3 1]); 

%% Columnate for consistent behavior: 
 
if isrow(dh)
   dh = dh'; 
end
      
end

