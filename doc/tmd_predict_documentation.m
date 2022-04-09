

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


ptype = 'h'; 
filename = 'CATS2022_02-21.nc';

%%


%sl_old = tmd_tide_pred('Model_CATS2008',datenum(t),lat*ones(size(t)),lon*ones(size(t)),'z'); 
tic
sl_old = tmd_tide_pred('Model_CATS2008',datenum(t),lat,lon,'z'); 
toc 

tic
profile clear 
profile on
sl_new = tmd_predict('CATS2022_02-21.nc',lat,lon,t);
profile viewer
toc

hold on
plot(t,sl_old)
plot(t,sl_new)

isf = isfinite(sl) & isfinite(sl_old);
%corr(sl(isf),sl_old(isf)).^2

return
%%
tic
sl_old = tmd_tide_pred('Model_CATS2008',datenum(t),lat*ones(size(t)),lon*ones(size(t)),'z'); 
toc 

tic
sl_new2 = tmd_predict('CATS2022_02-21.nc',lat*ones(size(t)),lon*ones(size(t)),t); 
toc

%%

profile clear 
profile on
sl_new2 = tmd_predict('CATS2022_02-21.nc',lat,lon,t); 
profile viewer 

plot(t,sl_new2)
%%



ptype = 'h'; 
filename = 'CATS2022_02-21.nc';

[lat,lon] = psgrid('south pole',5000,2); 

t = datenum(1984,5,27,2,17,0); 

ptype = 'h'; 


tic
sl_old = tmd_tide_pred('Model_CATS2008',t,lat,lon,'z'); 
toc % 21 sec

%
tic
z = tmd_predict('CATS2022_02-21.nc',lat,lon,t); 
toc

%%

figure
subplot(1,2,1)
pcolorps(lat,lon,sl_old)
colorbar
bedmachine

subplot(1,2,2)
pcolorps(lat,lon,z)
colorbar
bedmachine

%%


addpath(genpath('/Users/cgreene/Downloads/TMD3.00_alpha'))
addpath(genpath('/Users/cgreene/Downloads/CATS2008'))


ptype = 'h'; 

%%

filename = 'CATS2022_02-21.nc';

[lat,lon] = psgrid('amundsen sea',800,1); 

%t = datenum(1901,2,6,12,0,0); 

t = datenum(1984,5,27,2,17,0);


t = [t t+.2]; 


wct = tmd_interp(filename,'wct',lat,lon); 

profile clear 
profile on 
z = tmd_predict(filename,lat,lon,t);
profile viewer 


figure
surfps(lat,lon,-wct,z)
view(2)
axis off
shadem(-3)

%%


ptype = 'h'; 
conList = strsplit(ncreadatt(filename,'cons','long_name')); 

if isdatetime(t)
   t = datenum(t); 
end

% Columnate input coordinates for consistent behavior: 
InputGridSize = size(lat); 
lat = reshape(lat,[],1); 
lon = reshape(lon,[],1); 

InputTimeSize = size(t); 
t = reshape(t,[],1); 


hc = tmd_interp(filename,ptype,lat,lon,'constituents',conList);

hc = permute(hc,[3 1 2])/1000;

hc = complex(real(hc),-imag(hc)); 

[astrol_s,astrol_h,astrol_p,astrol_N] = tmd_astrol(t);

%%
