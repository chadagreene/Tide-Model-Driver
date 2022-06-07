
%% Load measured data 
% For this example, we use the same data that's described in the <tutorial_currents.html Tidal Current Tutorial>.
% Begin by loading the data: 

t = ncread('ADCP_S112.nc','time') + datenum(1950,1,1,0,0,0); 
lat = ncread('ADCP_S112.nc','lat'); 
lon = ncread('ADCP_S112.nc','lon'); 
z = ncread('ADCP_S112.nc','z'); 
u = ncread('ADCP_S112.nc','u'); 
v = ncread('ADCP_S112.nc','v'); 

% Calculate depth-averaged time series 
u_bar = mean(u,'omitnan'); 
v_bar = mean(v,'omitnan'); 

% Indices of two weeks in May 2012:
ind = t>datenum('may 16, 2012') & t<datenum('may 30, 2012'); 
ind = t>datenum('may 22, 2012') & t<datenum('may 26, 2012'); 

figure
plot(t(ind),u_bar(ind))
hold on
plot(t(ind),v_bar(ind))
axis tight 
box off
datetick('x','keeplimits')
ylabel 'depth-averaged current (m/s)'
legend('u','v')

%% 
% Here's a different way to plot the same data, with u and v plotted on
% their respective spatial dimensions, and time represented as color: 

figure
scatter(u_bar(ind),v_bar(ind),25,t(ind),'filled')

%%

filename = '/Users/cgreene/Documents/data/tides/CATS2008_update_2022-06-05.nc';
[umajor,uminor,uphase,uincl] = tmd_ellipse(filename,'o1',lat,lon);

%%

hold on
h=plot_ellipse(umajor,uminor,uincl,uphase); 

