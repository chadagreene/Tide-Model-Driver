% This script provides some visual checks to ensure converted tide model
% files with TMD3.0 are providing the same results as the old TMD2.4 (within
% some tolerance for numerical rounding and slight differences in
% interpolation). 
% 
%% Three model names for intercomparison: 

model1 = 'TPXO9_atlas_v5.nc';
model2 = 'AOTIM5.nc';
model3 = 'AODTM5.nc';

%% Compare wct of three models: 

[wct,x,y] = bedmachine_data('wct','greenland');

skip = 10; 
wct = wct(1:skip:end,1:skip:end); 
x = x(1:skip:end); 
y = y(1:skip:end); 
[X,Y] = meshgrid(x,y); 
[Lat,Lon] = psn2ll(X,Y); 

wct1 = tmd_interp(model1,'wct',Lat,Lon); 
wct2 = tmd_interp(model2,'wct',Lat,Lon); 
wct3 = tmd_interp(model3,'wct',Lat,Lon); 

figure
subplot(3,3,1) 
imagescn(x,y,wct)
axis image off
title bedmachine

subplot(3,3,2) 
imagescn(x,y,wct1)
axis image off
title(model1,'interpreter','none')

subplot(3,3,3) 
imagescn(x,y,wct1-wct)
axis image off
title diff
caxis([-1 1]*1000)
cmocean balance 


subplot(3,3,4) 
imagescn(x,y,wct)
axis image off
title bedmachine
subplot(3,3,5) 
imagescn(x,y,wct2)
axis image off
title(model2,'interpreter','none')
subplot(3,3,6) 
imagescn(x,y,wct2-wct)
axis image off
title diff
caxis([-1 1]*1000)
cmocean balance 


subplot(3,3,7) 
imagescn(x,y,wct)
axis image off
title bedmachine
subplot(3,3,8) 
imagescn(x,y,wct3)
axis image off
title(model3,'interpreter','none')
subplot(3,3,9) 
imagescn(x,y,wct3-wct)
axis image off
title diff
caxis([-1 1]*1000)
cmocean balance 

%% Compare masks of three models: 

figure
subplot(1,3,1) 
imagescn(x,y,tmd_interp(model1,'mask',Lat,Lon))
title(model1,'interpreter','none')
axis image off

subplot(1,3,2) 
imagescn(x,y,tmd_interp(model2,'mask',Lat,Lon))
title(model2,'interpreter','none')
axis image off

subplot(1,3,3) 
imagescn(x,y,tmd_interp(model3,'mask',Lat,Lon))
title(model3,'interpreter','none')
axis image off

%% Tide height observations near Nuuk, Greenland 
% Note, the Nuuk tide gauge is too close to shore for the Arc2km model, 
% so we have to "unmask" the solution there. 

% Load example tide gauge data (found in the doc/example_data folder.)
fn = 'h820_nuuk.nc'; 
lat = ncread(fn,'lat');
lon = ncread(fn,'lon');
t = ncread(fn,'time')+datenum(1800,1,1,0,0,0); % units = 'days since 1800-01-01 00:00:00' 
sl = ncread(fn,'sea_level')/1000; 

% Predict tides at the tide gauge location: 
sl_model1 = tmd_predict(model1,lat,lon,t,'h','coasts','unmask');
sl_model2 = tmd_predict(model2,lat,lon,t,'h','coasts','unmask');
sl_model3 = tmd_predict(model3,lat,lon,t,'h','coasts','unmask');
u_model1 = tmd_predict(model1,lat,lon,t,'U');
u_model2 = tmd_predict(model2,lat,lon,t,'U');
u_model3 = tmd_predict(model3,lat,lon,t,'U');
v_model1 = tmd_predict(model1,lat,lon,t,'V');
v_model2 = tmd_predict(model2,lat,lon,t,'V');
v_model3 = tmd_predict(model3,lat,lon,t,'V');

% Plot observed and predicted tides: 
figure
subsubplot(3,1,1)
p(1)=plot(t,sl-mean(sl,'omitnan'),'k','linewidth',2);
hold on
p(2)=plot(t,sl_model1,'linewidth',1);
p(3)=plot(t,sl_model2,'linewidth',1);
p(4)=plot(t,sl_model3,'linewidth',1);
ylabel 'tide height (m)'
legend(p(1:4),'observations',model1,model2,model3,...
   'interpreter','none','location','best')
legend boxoff
box off
axis tight
xlim([736976.01     736981.12])
datetick('x','keeplimits')

% Plot tidal currents:
subsubplot(3,1,2)
hold on
plot(t,u_model1,'linewidth',1,'color',p(2).Color);
plot(t,u_model2,'linewidth',1,'color',p(3).Color);
plot(t,u_model3,'linewidth',1,'color',p(4).Color);
ylabel 'zonal transport (m^2/s)'
box off
axis tight
xlim([736976.01     736981.12])
datetick('x','keeplimits')

subsubplot(3,1,3)
hold on
plot(t,v_model1,'linewidth',1,'color',p(2).Color);
plot(t,v_model2,'linewidth',1,'color',p(3).Color);
plot(t,v_model3,'linewidth',1,'color',p(4).Color);
ylabel 'meridional transport (m^2/s)'
box off
axis tight
xlim([736976.01     736981.12])
datetick('x','keeplimits')

%% Compare a model to TMD4 version of the same model and observations
% At a random location 

t = datetime('Jan 1, 2020'):hours(1):datetime('Mar 1, 2020'); 
lat = 79.5568; 
lon = -69.3935; 

addpath(genpath('/Users/cgreene/Downloads/TMD3.00_alpha'))
%cd('/Users/cgreene/Documents/data/tides/Arc5km2018')

% Predict tides: 
% tic
% sl_old = tmd_tide_pred('Model_Arc2kmTM_v1',datenum(t),lat,lon,'z'); 
% old_time = toc 
% tic
% sl_new = tmd_predict('Arc2kmTM_v1.nc',lat,lon,t,'h','coast','unmask');
% new_time = toc
tic
sl_old = tmd_tide_pred('Model_Arc5km2018',datenum(t),lat,lon,'z'); 
old_time = toc 
tic
sl_new = tmd_predict('Arc5km2018.nc',lat,lon,t,'h','coast','unmask');
new_time = toc

% Plot
figure
subplot(2,1,1)
hold on
plot(t,sl_old,'k','linewidth',2)
plot(t,sl_new,'y--','linewidth',2)
legend('TMD 2.4','TMD 3.0')
axis tight

%% AOTIM5

cd('/Users/cgreene/Documents/data/tides/AOTIM5')

% Load reference data: 
fn = 'h820_nuuk.nc';
lat = ncread(fn,'lat');
lon = ncread(fn,'lon');
t = ncdateread(fn,'time');
sl = ncread(fn,'sea_level')/1000; 

% Predict tides: 
tic
sl_old = tmd_tide_pred('Model_AOTIM5',datenum(t),lat,lon,'z'); 
old_time = toc 
tic
sl_new = tmd_predict('AOTIM5.nc',lat,lon,t,'h','coast','unmask');
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

%% AODTM5

%cd('/Users/cgreene/Documents/data/tides/AODTM5')

% Load reference data: 
fn = 'h820_nuuk.nc';
lat = ncread(fn,'lat');
lon = ncread(fn,'lon');
t = ncdateread(fn,'time');
sl = ncread(fn,'sea_level')/1000; 

% Predict tides: 
tic
sl_old = tmd_tide_pred('/Users/cgreene/Documents/data/tides/AODTM5/Model_AODTM5',datenum(t),lat,lon,'z'); 
old_time = toc 
tic
sl_new = tmd_predict('AODTM5.nc',lat,lon,t,'h','coast','unmask');
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

%% Drift track 

% Do it again for "drift track" solution: 
lat = repmat(lat,size(t)); 
lon = repmat(lon,size(t)); 

% Predict tides: 
tic
sl_old = tmd_tide_pred('Model_Arc5km2018',datenum(t),lat,lon,'z'); 
old_time = toc 
tic
sl_new = tmd_predict(model3,lat,lon,t,'h','coast','unmask');
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


%% Intercomparison 


%% Tide height observations near Nuuk, Greenland AODTM5
% Note, the Nuuk tide gauge is too close to shore for the Arc2km model, 
% so we have to "unmask" the solution there. 

% Load example tide gauge data (found in the doc/example_data folder.)
fn = 'h820_nuuk.nc'; 
lat = ncread(fn,'lat');
lon = ncread(fn,'lon');
t = ncread(fn,'time')+datenum(1800,1,1,0,0,0); % units = 'days since 1800-01-01 00:00:00' 
sl = ncread(fn,'sea_level')/1000; 

% Predict tides at the tide gauge location: 
sl_model1 = tmd_predict(model1,lat,lon,t,'h');
sl_Gr1km = tmd_predict('Gr1kmTM_v1.nc',lat,lon,t,'h');
sl_model2 = tmd_predict('AODTM5.nc',lat,lon,t,'h','coasts','unmask');
u_model1 = tmd_predict(model1,lat,lon,t,'U');
u_model2 = tmd_predict('Gr1kmTM_v1.nc',lat,lon,t,'U');
u_model2 = tmd_predict('AODTM5.nc',lat,lon,t,'U');
v_model1 = tmd_predict(model1,lat,lon,t,'V');
v_model2 = tmd_predict('Gr1kmTM_v1.nc',lat,lon,t,'V');
v_model2 = tmd_predict('AODTM5.nc',lat,lon,t,'V');

% Plot observed and predicted tides: 
figure
subsubplot(3,1,1)
p(1)=plot(t,sl-mean(sl,'omitnan'),'k','linewidth',2);
hold on
p(2)=plot(t,sl_model1,'linewidth',1);
p(3)=plot(t,sl_Gr1km,'linewidth',1);
p(4)=plot(t,sl_model2,'linewidth',1);
ylabel 'tide height (m)'
legend('observations','TPXO9_atlas_v5','Gr1kmTm','AODTM5',...
   'interpreter','none','location','best')
legend boxoff
box off
axis tight
xlim([736976.01     736981.12])
datetick('x','keeplimits')

% Plot tidal currents:
subsubplot(3,1,2)
hold on
plot(t,u_model1,'linewidth',1,'color',p(2).Color);
plot(t,u_model2,'linewidth',1,'color',p(3).Color);
plot(t,u_model2,'linewidth',1,'color',p(4).Color);
ylabel 'zonal transport (m^2/s)'
box off
axis tight
xlim([736976.01     736981.12])
datetick('x','keeplimits')

subsubplot(3,1,3)
hold on
plot(t,v_model1,'linewidth',1,'color',p(2).Color);
plot(t,v_model2,'linewidth',1,'color',p(3).Color);
plot(t,v_model2,'linewidth',1,'color',p(4).Color);
ylabel 'meridional transport (m^2/s)'
box off
axis tight
xlim([736976.01     736981.12])
datetick('x','keeplimits')

%% Tide height observations near Nuuk, Greenland AOTIM5
% Note, the Nuuk tide gauge is too close to shore for the Arc2km model, 
% so we have to "unmask" the solution there. 

% Load example tide gauge data (found in the doc/example_data folder.)
fn = 'h820_nuuk.nc'; 
lat = ncread(fn,'lat');
lon = ncread(fn,'lon');
t = ncread(fn,'time')+datenum(1800,1,1,0,0,0); % units = 'days since 1800-01-01 00:00:00' 
sl = ncread(fn,'sea_level')/1000; 

% Predict tides at the tide gauge location: 
sl_model1 = tmd_predict(model1,lat,lon,t,'h');
sl_Gr1km = tmd_predict('Gr1kmTM_v1.nc',lat,lon,t,'h');
sl_model2 = tmd_predict('AOTIM5.nc',lat,lon,t,'h','coasts','unmask');
u_model1 = tmd_predict(model1,lat,lon,t,'U');
u_model2 = tmd_predict('Gr1kmTM_v1.nc',lat,lon,t,'U');
u_model2 = tmd_predict('AOTIM5.nc',lat,lon,t,'U');
v_model1 = tmd_predict(model1,lat,lon,t,'V');
v_model2 = tmd_predict('Gr1kmTM_v1.nc',lat,lon,t,'V');
v_model2 = tmd_predict('AOTIM5.nc',lat,lon,t,'V');

% Plot observed and predicted tides: 
figure
subsubplot(3,1,1)
p(1)=plot(t,sl-mean(sl,'omitnan'),'k','linewidth',2);
hold on
p(2)=plot(t,sl_model1,'linewidth',1);
p(3)=plot(t,sl_Gr1km,'linewidth',1);
p(4)=plot(t,sl_model2,'linewidth',1);
ylabel 'tide height (m)'
legend('observations','TPXO9_atlas_v5','Gr1kmTm','AOTIM5',...
   'interpreter','none','location','best')
legend boxoff
box off
axis tight
xlim([736976.01     736981.12])
datetick('x','keeplimits')

% Plot tidal currents:
subsubplot(3,1,2)
hold on
plot(t,u_model1,'linewidth',1,'color',p(2).Color);
plot(t,u_model2,'linewidth',1,'color',p(3).Color);
plot(t,u_model2,'linewidth',1,'color',p(4).Color);
ylabel 'zonal transport (m^2/s)'
box off
axis tight
xlim([736976.01     736981.12])
datetick('x','keeplimits')

subsubplot(3,1,3)
hold on
plot(t,v_model1,'linewidth',1,'color',p(2).Color);
plot(t,v_model2,'linewidth',1,'color',p(3).Color);
plot(t,v_model2,'linewidth',1,'color',p(4).Color);
ylabel 'meridional transport (m^2/s)'
box off
axis tight
xlim([736976.01     736981.12])
datetick('x','keeplimits')

%% Check where currents are known to be high
% The location of the Nuuk tide gauge may be bad for testing, 
% so check where transport is known better. 

[UAm,lon,lat] = tmd_data(model1,'UAm'); 
mask = tmd_data(model1,'mask');
tmp = sum(abs(UAm),3); 

figure
h = imagescn(lon,lat,tmp); 
h.AlphaData = mask; 
axis image  

loni = 294.44; 
lati = 60.96; 
hold on
plot(loni,lati,'ro')

%%


lon = 277.6393; 
lat = 74.1257; 
t = datetime('jan 1, 2020'):minutes(1):datetime('feb 28, 2020'); 

% Predict tides at the tide gauge location: 
sl_model1 = tmd_predict(model1,lat,lon,t,'h');
sl_model2 = tmd_predict(model2,lat,lon,t,'h');
sl_model3 = tmd_predict(model3,lat,lon,t,'h');
u_model1 = tmd_predict(model1,lat,lon,t,'U');
u_model2 = tmd_predict(model2,lat,lon,t,'U');
u_model3 = tmd_predict(model3,lat,lon,t,'U');
v_model1 = tmd_predict(model1,lat,lon,t,'V');
v_model2 = tmd_predict(model2,lat,lon,t,'V');
v_model3 = tmd_predict(model3,lat,lon,t,'V');

% Plot observed and predicted tides: 
figure
subsubplot(3,1,1)
hold on
p(2)=plot(t,sl_model1,'linewidth',1);
p(3)=plot(t,sl_model2,'linewidth',1);
p(4)=plot(t,sl_model3,'linewidth',1);
ylabel 'tide height (m)'
legend(p(2:4),model1,model2,model3,...
   'interpreter','none','location','best')
legend boxoff
box off
axis tight
ax(1) = gca; 

% Plot tidal currents:
subsubplot(3,1,2)
hold on
plot(t,u_model1,'linewidth',1,'color',p(2).Color);
plot(t,u_model2,'linewidth',1,'color',p(3).Color);
plot(t,u_model3,'linewidth',1,'color',p(4).Color);
ylabel 'zonal transport (m^2/s)'
box off
axis tight
ax(2) = gca; 

subsubplot(3,1,3)
hold on
plot(t,v_model1,'linewidth',1,'color',p(2).Color);
plot(t,v_model2,'linewidth',1,'color',p(3).Color);
plot(t,v_model3,'linewidth',1,'color',p(4).Color);
ylabel 'meridional transport (m^2/s)'
box off
axis tight
ax(3) = gca; 
linkaxes(ax,'x')
