%% Description 
% This script creates the animated TMD logo gif. It uses the |earthimage|, 
% |cmocean|, |shadem|, and |gif| functions, which are all part of the Climate Data 
% Toolbox for MATLAB. 
% 
% As an alternative to the gif function, you may prefer MATLAB's built-in
% |exportgraphics| function (after R2022a) for gifs, or the |VideoWriter|
% function for other video formats. 
% 
% Another example of creating a gif like this can be found in the
% documentation for the tmd_predict function. One distinct difference
% between this example and the <tmd_predict_documentation.html |tmd_predict|> example, is that here we only
% predict one map at a time, whereas the tmd_predict example creates a
% cube, where the third dimension corresponds to time. Generally for small
% regions, for a small number of timesteps, and for small model files (with
% few grid points and/or few constituents), the cube approach will be
% faster, because it only loads and interpolates the data once. In this
% case, we solve one frame at a time (which is slower, overall) because
% memory can become an issue for very large cubes. 
% 

%% Solve tides 
% Solving for the whole cube takes about 2 minutes to run on my laptop from 2019. 
% 
%  % Model filename: 
%  fn = 'TPXO9_atlas_v5.nc';
%  
%  % Create 0.2 degree resolution grid: 
%  [Lat,Lon] = cdtgrid(0.2); 
%  
%  % Create an hourly time array for 25 hours: 
%  t = datetime('mar 29, 2017'):hours(1):datetime('mar 30, 2017'); 
%  
%  % Get water column thickness for the grid: 
%  wct = tmd_interp(fn,'wct',Lat,Lon); 
%  
%  % Predict tides: 
%  Z = tmd_predict(fn,Lat,Lon,t,'h'); 

%% Animate
% 
%  figure('position',[10 10 560 280],'color','k')
%  hs = surf(Lon,Lat,-wct,Z(:,:,1)); % the first "slice" of Z data
%  hold on
%  he = earthimage('bottom');
%  shading interp
%  view(2)
%  axis image off
%  cmocean balance
%  caxis([-1 1]*2)
%  shadem(-12) % tiny bit of hillshade
%  hs.ZData = hs.ZData-min(hs.ZData(:))+1; % raises above earthimage 
%  set(gca,'position',[0 0 1 1])
%  
%  txt = text(0.5,0.5,{'Tide Model Driver';'for MATLAB'},'horiz','center',...
%     'vert','middle','units','normalized','fontweight','bold',...
%     'fontangle','italic','fontsize',36); 
%  
%  % Write the first frame: 
%  gif('tmd_logo.gif','delaytime',1/8,'resolution',200)
%  
%  % Loop through the remaining frames
%  for k = 2:length(t) 
%     hs.CData = Z(:,:,k); 
%     gif % writes this frame
%  end

%% Author Info
% This script was written by Chad A. Greene, June 2022. 