

hc = tmd_interp('CATS2008_update_2022-04-22.nc','h',lat,lon,'constituents',{'m2','mf'});
hc = permute(hc,[3 1 2]);


hg = tmd_interp('/Users/cgreene/Downloads/TPXO9_atlas_v5/TPXO9_atlas30_update_2022-05-06.nc','h',lat,lon,'constituents',{'m2','mf'});
hg = permute(hg,[3 1 2]);


filename = '/Users/cgreene/Downloads/TPXO9_atlas_v5/TPXO9_atlas30_update_2022-05-06.nc';
ph = ncread(filename,'phase');
[angle(hc) angle(conj(hg)) -ph([4 6])]
%[abs(hg) angle(hg)]

%%

addpath(genpath('/Users/cgreene/Downloads/TMD3.00_alpha'))
addpath(genpath('/Users/cgreene/Downloads/CATS2008'))
addpath(genpath('/Users/cgreene/Downloads/Gr1kmTM'))


 fn = '/Users/cgreene/Documents/GitHub/Tide-Model-Driver/doc/h127.nc';% bad
fn = '/Users/cgreene/Documents/GitHub/Tide-Model-Driver/doc/h189.nc'; % good
%fn = '/Users/cgreene/Documents/GitHub/Tide-Model-Driver/doc/h700.nc'; 
lat = ncread(fn,'lat');
lon = ncread(fn,'lon');
t = ncdateread(fn,'time');
sl = ncread(fn,'sea_level')/1000; 

plot(t,sl-mean(sl,'omitnan'))

sl_old = tmd_tide_pred('Model_CATS2008',datenum(t),lat,lon,'z'); 
hold on
plot(t,sl_old-mean(sl_old,'omitnan'))

%sl_new = tmd_predict('CATS2022_02-21.nc',lat,lon,t);

sl_cats2 = tmd_predict('CATS2008_update_2022-04-22.nc',lat,lon,t,'h','coast','unmask');
plot(t,sl_cats2-mean(sl_cats2,'omitnan'),'--')

%%
clf
plot(t,sl_cats2-mean(sl_cats2,'omitnan'))
hold on
sl_tpxo = tmd_predict('/Users/cgreene/Downloads/TPXO9_atlas_v5/TPXO9_atlas30_update_2022-05-06.nc',lat,lon,t,'h','coast','unmask'); 

plot(t,sl_tpxo-mean(sl_tpxo,'omitnan'))

xlim(datetime([datenum('Jun 20, 2009, 12:00:56') datenum('Jul 02, 2009, 09:44:26')],'convertfrom','datenum'))


%%

sl_cats2 = tmd_predict('CATS2008_update_2022-04-22.nc',lat,lon,t,'h','coast','unmask','constituents',{'m2'});
sl_tpx = tmd_predict('/Users/cgreene/Downloads/TPXO9_atlas_v5/TPXO9_atlas30_update_2022-05-03.nc',lat,lon,t,'h','coast','unmask','constituents',{'m2'});
plot(t,sl_cats2)
hold on
plot(t,sl_tpx)

%%
figure
plot(t,sl)
hold on
plot(t,sl-sl_new)

%% Example: NaNs near the coast 
% In some cases, you may be interested in a tide gauge close to the coast
% that gets NaN'd out because it's closest to a land mask pixel by the
% default nearest-neighbor interpolation. Here's an example:

fn = '/Users/cgreene/Documents/GitHub/Tide-Model-Driver/doc/h127.nc';% bad

lat = ncread(fn,'lat');
lon = ncread(fn,'lon');
t = ncdateread(fn,'time');
sl = ncread(fn,'sea_level')/1000; 

% Predict tides at the location of interest: 
sl_predict = tmd_predict('CATS2008_update_2022-04-22.nc',lat,lon,t,'h');

% Look how many tide predictions are NaNs: 
sum(isnan(sl_predict))

%%
figure
mapzoomps(lat,lon,'mapwidth',100,'inset','ne') % centers a map on the region of interest
plotps(lat,lon,'ko','markerfacecolor','y') % plots the tide gauge location
modismoaps('contrast','low') % background image

%%

figure
plot(t,sl-mean(sl,'omitnan'))

sl_predict = tmd_predict('CATS2008_update_2022-04-22.nc',lat,lon,t,'h','coast','unmask');
hold on
plot(t,sl_predict)
plot(t,sl - mean(sl,'omitnan') - sl_predict)
legend('observations','prediction','detided observation','location','northwest')


xlim(datetime([datenum('Sep 07, 2001, 17:03:49') datenum('Sep 21, 2001, 22:29:18')],'convertfrom','datenum'))


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

[lat,lon] = psgrid('amundsen sea',900,1); 

t = datenum(1984,5,27,2,17,0); 

ptype = 'h'; 


tic
sl_old = tmd_tide_pred('Model_CATS2008',t,lat,lon,'z'); 
toc % 21 sec

%
tic
z = tmd_predict('CATS2008_update_2022-04-11.nc',lat,lon,t,'h','flexure'); 
toc


%%

m2 = tmd_interp('CATS2008_update_2022-04-11.nc','h',lat,lon);

%%

mx = max(abs([sl_old(:);z(:)])); 

figure
subplot(1,2,1)
pcolorps(lat,lon,sl_old)
colorbar
axis image off
bedmachine
ax(1) = gca; 
caxis([-1 1]*mx)
cmocean bal
shading interp

subplot(1,2,2)
pcolorps(lat,lon,z)
colorbar
axis image off
bedmachine
ax(2) = gca; 
linkaxes(ax,'xy')
caxis([-1 1]*mx)
cmocean bal
shading interp

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

%% Greenland 

fn = 'h820_nuuk.nc';
lat = ncread(fn,'lat');
lon = ncread(fn,'lon');
t = ncdateread(fn,'time');
sl = ncread(fn,'sea_level')/1000; 

plot(t,sl-mean(sl,'omitnan'))


sl_old = tmd_tide_pred('Model_Gr1kmTM_v1',datenum(t),lat,lon,'z'); 
sl_predict2 = tmd_predict('/Users/cgreene/Downloads/Gr1kmTM/Gr1kmTM_update_2022-05-15.nc',lat,lon,t,'h','coast','unmask');

hold on

%plot(t,sl)
plot(t,sl_old,'k','linewidth',2)
plot(t,sl_predict2,'y--','linewidth',2)

%%
filename = '/Users/cgreene/Downloads/Gr1kmTM/Gr1kmTM_update_2022-05-06.nc'; 
Model='Model_Gr1kmTM_v1'


%%
figure
plot(t,sl)
hold on
plot(t,sl-sl_predict)


