


[Lat,Lon] = psgrid('scar inlet',1500,1); 

wct = tmd_interp('CATS2022_02-21.nc','wct',Lat,Lon); 

pcolorps(Lat,Lon,wct)
axis tight off
bedmachine
cmocean deep 
shadem(4)

return

%%


%mask = tmd_interp('CATS2022_02-21.nc','mask',Lat,Lon); 


h = tmd_interp('CATS2022_02-21.nc','hAm',Lat,Lon,'constituent','k2'); 

hold on
contourps(Lat,Lon,h,'k')

%%
U = tmd_interp('CATS2022_02-21.nc','uAm',Lat,Lon,'constituent','k2'); 
V = tmd_interp('CATS2022_02-21.nc','vAm',Lat,Lon,'constituent','k2'); 


[Vx,Vy] = uv2vxvy(Lat,Lon,U,V); 


[X,Y] = ll2ps(Lat,Lon); 

sc = 0.02; 
Xr = imresize(X,sc); 
Yr = imresize(Y,sc); 
Vxr = imresize(Vx,sc); 
Vyr = imresize(Vy,sc); 


hold on
q = quiver(Xr,Yr,Vxr,Vyr,'r'); 

%%

