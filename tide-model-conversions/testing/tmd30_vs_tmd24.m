% This script provides some visual checks to ensure converted tide model
% files with TMD3.0 are providing the same results as the old TMD2.4 (within
% some tolerance for numerical rounding and slight differences in
% interpolation). 
% 

addpath(genpath('/Users/cgreene/Downloads/TMD3.00_alpha'))

%% Model_Gr1kmTM_v1

% Load reference data: 
fn = 'h820_nuuk.nc';
lat = ncread(fn,'lat');
lon = ncread(fn,'lon');
t = ncdateread(fn,'time');
sl = ncread(fn,'sea_level')/1000; 

% Predict tides: 
tic
sl_old = tmd_tide_pred('Model_Gr1kmTM_v1',datenum(t),lat,lon,'z'); 
old_time = toc 
tic
sl_new = tmd_predict('Gr1kmTM_update_2022-05-15.nc',lat,lon,t,'h','coast','unmask');
new_time = toc

% Plot
figure
subplot(2,1,1)
plot(t,sl-mean(sl,'omitnan'),'color',.6*[1 1 1])
hold on
plot(t,sl_old,'k','linewidth',2)
plot(t,sl_new,'y--','linewidth',2)
legend('observations','TMD 2.4','TMD 3.0')
axis tight
xlim(datetime([2015 2015],[6 7],[28 5]))

% Do it again for "drift track" solution: 
lat = repmat(lat,size(t)); 
lon = repmat(lon,size(t)); 

% Predict tides: 
tic
sl_old = tmd_tide_pred('Model_Gr1kmTM_v1',datenum(t),lat,lon,'z'); 
old_time = toc 
tic
sl_new = tmd_predict('Gr1kmTM_update_2022-05-15.nc',lat,lon,t,'h','coast','unmask');
new_time = toc


% Plot
subplot(2,1,2)
plot(t,sl-mean(sl,'omitnan'),'color',.6*[1 1 1])
hold on
plot(t,sl_old,'k','linewidth',2)
plot(t,sl_new,'y--','linewidth',2)
legend('observations','TMD 2.4','TMD 3.0')
axis tight
xlim(datetime([2015 2015],[6 7],[28 5]))