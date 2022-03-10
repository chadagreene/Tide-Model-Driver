

fn = '/Users/cgreene/Documents/GitHub/Tide-Model-Driver/doc/h189.nc'; 
lat = ncread(fn,'lat');
lon = ncread(fn,'lon');
t = ncdateread(fn,'time');
sl = ncread(fn,'sea_level')/1000; 
sl = sl-mean(sl,'omitnan'); 

plot(t,sl)

%%

sl_old = tmd_tide_pred('Model_CATS2008',datenum(t),lat,lon,'z'); 

hold on
plot(t,sl_old)

isf = isfinite(sl) & isfinite(sl_old);
corr(sl(isf),sl_old(isf)).^2


%%

sl_new = tmd_predict('CATS2022_02-21.nc',lat,lon,t); 


%%

filename = 'CATS2022_02-21.nc'; 
ptype='h'; 


