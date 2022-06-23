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
%
%% 
% <<../markdown_figures/tmd_logo.gif>>
% 
%% Not seasick yet?
% If the animation above doesn't make you want to lose your lunch, perhaps
% this one will. This example closely follows the one above, but in addition 
% to tidal height, we'll also predict tidal transports. 
% 
% Some notes: 
% 
% # We'll solve transports on a 0.2 degree grid, but we'll only plot
% vectors on a 3 degree grid. This requires anti-aliasing. 
% # Technically, we should solve tides on the native 1/30 degree grid, then
% lowpass filter to 6 degrees or a little more, then interpolate to a 3
% degree grid. That would produce anti-aliased transports that meet the
% Nyquist criteria. Instead, we'll do this sort of hybrid approach, where
% we'll solve transports on a 0.2 degree grid, then lowpass filter to 6
% degrees (using the |filt2| function in the Climate Data Toolbox)
% and interpolate. So we might end up with a little bit of noise,
% but let's not be pedants.
% # Matlab's |text| function doesn't allow bordered text, so below I create
% eight white text objects, each offset by some combination of 0.3 degrees
% in longitude or 0.15 degrees in latitude. The end result produces the
% appearance of outlined text. 
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
%  U = tmd_predict(fn,Lat,Lon,t,'U'); 
%  V = tmd_predict(fn,Lat,Lon,t,'V'); 
% 
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
%  ax = axis; % gets axis limits 
%  
%  % Downsampled grid: 
%  skip = 15; 
%  Lonf = Lon(floor(skip/2):skip:end,floor(skip/2):skip:end); 
%  Latf = Lat(floor(skip/2):skip:end,floor(skip/2):skip:end); 
% 
%  % Land mask for the downsampled grid: 
%  land = ~tmd_interp(fn,'mask',Latf,Lonf); 
% 
%  % Anti-aliased (by lowpass filtering) U and V: 
%  tmpu = interp2(Lon,Lat,filt2(U(:,:,1),1,skip*2,'lp'),Lonf,Latf); 
%  tmpv = interp2(Lon,Lat,filt2(V(:,:,1),1,skip*2,'lp'),Lonf,Latf); 
% 
%  % Ensure land pixels are NaN: 
%  tmpu(land) = nan; 
% 
%  % Plot arrows (with a Z component to get above surface topography): 
%  q = quiver3(Lonf,Latf,9*ones(size(Lonf)),tmpu,tmpv,zeros(size(Lonf)),'k'); 
%  q.LineWidth = 0.1; 
% 
%  % Plot background text: 
%  dl = 0.3; 
%  txtw(1) = text(dl,0,10,{'Tide Model Driver';'for MATLAB'},'horiz','center',...
%     'vert','middle','fontweight','bold',...
%     'fontangle','italic','fontsize',36,'color','w'); 
%   txtw(2) = text(-dl,0,10,{'Tide Model Driver';'for MATLAB'},'horiz','center',...
%     'vert','middle','fontweight','bold',...
%     'fontangle','italic','fontsize',36,'color','w'); 
%   txtw(3) = text(dl,dl/2,10,{'Tide Model Driver';'for MATLAB'},'horiz','center',...
%     'vert','middle','fontweight','bold',...
%     'fontangle','italic','fontsize',36,'color','w'); 
%   txtw(4) = text(dl,-dl/2,10,{'Tide Model Driver';'for MATLAB'},'horiz','center',...
%     'vert','middle','fontweight','bold',...
%     'fontangle','italic','fontsize',36,'color','w'); 
%   txtw(5) = text(-dl,dl/2,10,{'Tide Model Driver';'for MATLAB'},'horiz','center',...
%     'vert','middle','fontweight','bold',...
%     'fontangle','italic','fontsize',36,'color','w'); 
%   txtw(6) = text(-dl,-dl/2,10,{'Tide Model Driver';'for MATLAB'},'horiz','center',...
%     'vert','middle','fontweight','bold',...
%     'fontangle','italic','fontsize',36,'color','w'); 
%   txtw(7) = text(0,dl/2,10,{'Tide Model Driver';'for MATLAB'},'horiz','center',...
%     'vert','middle','fontweight','bold',...
%     'fontangle','italic','fontsize',36,'color','w'); 
%   txtw(8) = text(0,-dl/2,10,{'Tide Model Driver';'for MATLAB'},'horiz','center',...
%     'vert','middle','fontweight','bold',...
%     'fontangle','italic','fontsize',36,'color','w'); 
%  
%  % Plot
%  txt = text(0,0,10,{'Tide Model Driver';'for MATLAB'},'horiz','center',...
%     'vert','middle','fontweight','bold',...
%     'fontangle','italic','fontsize',36); 
%  
%  col = cmocean('phase',length(t)); 
%  txt.Color=col(1,:); 
% 
%  axis(ax); % ensures axis limits match in every frame 
% 
%  % Write the first frame: 
%  gif('tmd_logo_v2.gif','delaytime',1/8,'resolution',200)
%  
%  % Loop through the remaining frames
%  for k = 2:length(t) 
%     
%     % Update ocean color (height data): 
%     hs.CData = Z(:,:,k); 
%     
%     % Get anti-aliased U and V: 
%     tmpu = interp2(Lon,Lat,filt2(U(:,:,k),1,skip*2,'lp'),Lonf,Latf); 
%     tmpv = interp2(Lon,Lat,filt2(V(:,:,k),1,skip*2,'lp'),Lonf,Latf); 
%     tmpu(land) = nan; 
% 
%     % Update quiver arrows: 
%     q.UData = tmpu; 
%     q.VData = tmpv; 
%     
%     txt.Color=col(k,:); 
%     axis(ax) % ensures axis limits match in every frame
%     gif % writes this frame
%  end
% 
%% 
% 
% <<../markdown_figures/tmd_logo_v2.gif>>
%
%% Author Info
% This script was written by Chad A. Greene, June 2022. 