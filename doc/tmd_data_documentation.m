

%% Example 1 

[h,x,y] = tmd_data('CATS2022_02-21.nc','h'); 

whos h x y

%% Example 2 

[wct,lon,lat] = tmd_data('CATS2022_02-21.nc','wct','geo'); 

figure
pcolorps(lat,lon,wct) 
bedmachine

%% Example 3

[hAm,x,y] = tmd_data('CATS2022_02-21.nc','hAm','constituents','s2');  
hPh = tmd_data('CATS2022_02-21.nc','hPh','constituents','s2'); 

figure
imagesc(x,y,hAm); 
axis xy image
hold on
contour(x,y,rad2deg(hPh),-180:30:180,'k')
xlabel 'easting (km)' 
ylabel 'northing (km)' 

%% 

[u1,x1,y1] = tmd_data('CATS2022_02-21.nc','U','bounds',[-2498 3049]);

%%

