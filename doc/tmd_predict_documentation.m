

addpath(genpath('/Users/cgreene/Downloads/TMD3.00_alpha'))
addpath(genpath('/Users/cgreene/Downloads/CATS2008'))


fn = '/Users/cgreene/Documents/GitHub/Tide-Model-Driver/doc/h189.nc'; 
lat = ncread(fn,'lat');
lon = ncread(fn,'lon');
t = ncdateread(fn,'time');
sl = ncread(fn,'sea_level')/1000; 
sl = sl-mean(sl,'omitnan'); 

plot(t,sl)

%%

profile clear 
profile on
sl_old = tmd_tide_pred('Model_CATS2008',datenum(t),lat*ones(size(t)),lon*ones(size(t)),'z'); 
profile viewer 

hold on
plot(t,sl_old)

isf = isfinite(sl) & isfinite(sl_old);
%corr(sl(isf),sl_old(isf)).^2


%%

sl_new = tmd_predict('CATS2022_02-21.nc',lat,lon,t); 


%%

filename = 'CATS2022_02-21.nc'; 
ptype='h'; 

hc = squeeze(tmd_interp(filename,ptype,lat,lon));



%tmp1 = tmd_harp1(t,hc,conList) ;
%tmp2 = tmd_InferMinor(hc,conList,t);

%%

%[pu,pf] = nodal(datenum(t([1 end])) - datenum(1992,1,1) + 48622,cid)

[astrol_s,astrol_h,astrol_p,astrol_N] = tmd_astrol(t);

[pu,pf] = tmd_nodal(conList,astrol_p,astrol_N); 

%%

dh = InferMinor(hc,cid,datenum(t)); 

dh2 = tmd_InferMinor(hc,conList,t,astrol_s,astrol_h,astrol_p,astrol_N);

%%

dh3 = tmd_InferMinor(hc,conList,t,astrol_s,astrol_h,astrol_p,astrol_N);

