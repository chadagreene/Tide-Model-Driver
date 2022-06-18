[&larr; Back to TMD3.0 Main Page](../README.md)

# Animating tidal time series 
This script creates the animated TMD logo gif. It uses the `earthimage`, `cmocean`, `shadem`, and `gif` functions, which are all part of the [Climate Data Toolbox for MATLAB](https://github.com/chadagreene/CDT). 

As an alternative to the gif function, you may prefer MATLAB's built-in `exportgraphics` function (after R2022a) for gifs, or the `VideoWriter` function for other video formats. 
 
Another example of creating a gif like this can be found [in the documentation](tmd_predict_documentation.md#example-time-series-of-maps) for the `tmd_predict` function. 

*Note*: If the grid you're solving is particularly large and/or you have lots of timesteps, you may need to predict tides one map at a time. That will be slower than solving the entire cube at once, because loading and interpolating the tidal constituent data is takes a fair bit of time, but memory might become an issue for very large data cubes. 

This script takes about 2 minutes to run on my laptop from 2019. 

```matlab
% Model filename: 
fn = 'TPXO9_atlas_v5.nc';

% Create 0.2 degree resolution grid: 
[Lat,Lon] = cdtgrid(0.2); 

% Create an hourly time array for 25 hours: 
t = datetime('mar 29, 2017'):hours(1):datetime('mar 30, 2017'); 

% Get water column thickness for the grid: 
wct = tmd_interp(fn,'wct',Lat,Lon); 

% Predict tides: 
Z = tmd_predict(fn,Lat,Lon,t,'h'); 

figure('position',[10 10 560 280],'color','k')
hs = surf(Lon,Lat,-wct,Z(:,:,1)); % the first "slice" of Z data 
hold on
he = earthimage('bottom');
shading interp
view(2)
axis image off
cmocean balance
caxis([-1 1]*2)
shadem(-12) % tiny bit of hillshade
hs.ZData = hs.ZData-min(hs.ZData(:))+1; % raises above earthimage 
set(gca,'position',[0 0 1 1])

txt = text(0.5,0.5,{'Tide Model Driver';'for MATLAB'},'horiz','center',...
   'vert','middle','units','normalized','fontweight','bold',...
   'fontangle','italic','fontsize',36); 

% Write the first frame: 
gif('tmd_logo.gif','delaytime',1/8,'resolution',200)

% Loop through the remaining frames
for k = 2:length(t) 
   hs.CData = Z(:,:,k); 
   gif % writes this frame
end
```
![TMD logo](markdown_figures/tmd_logo.gif)

## Author Info
This script was written by Chad A. Greene, June 2022. 