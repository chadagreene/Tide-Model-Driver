function hhat = tmd_harp1(t,hc,constituents)
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


if isrow(t)
   t = t';
end

if isdatetime(t)
   t = datenum(t); 
end

t_1992 = t - datenum(1992,1,1);
time_s = t_1992*86400; % time (days since Jan 1, 1992) is converted to seconds by multiplying by 86400 s/yr. 

nc = length(constituents); % number of constituents. 

[ispec,~,ph,omega,~,~] = constit(con);

igood=ispec~=-1;
con1=con(igood,:);

[pu1,pf1] = tmd_nodal(constituents,astrol_p,astrol_N); % time (days since Jan 1, 1992) is converted to MJD by adding 48622

% phase variables: 
L = length(t_1992);
pu=zeros(L,nc);
pf=ones(L,nc);
pu(:,igood)=pu1; 
pf(:,igood)=pf1;

hhat=zeros(size(time_s));

hc_real = real(hc); 
hc_imag = imag(hc); 

for k=1:nc
   
   % Calculate just the frequency and phase so we won't have to calculate it twice: 
   tmp = omega(k)*time_s + ph(k)+pu(:,k); 
   
   arg = pf(:,k).*(hc_real(k).*cos(tmp) - hc_imag(k).*sin(tmp));
   hhat=hhat+arg;
end



end
 
