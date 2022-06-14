%% |tmd_conlist| documentation
% |tmd_conlist| returns a list of tidal constituents in a TMD3.0 compatible
% consolidated NetCDF tide model file. 
% 
% <TMD_Contents.html Back to Tide Model Driver Contents>.
%% Syntax
% 
%  conList = tmd_conlist(filename)
% 
%% Description 
% 
% |conList = tmd_conlist(filename)| returns a cell array of tidal
% constituents in the specified model file. 
% 
%% Example 
% Get a list of constituents in the updated CATS2008 model: 

conList = tmd_conlist('CATS2008_update_2022-06-11.nc')

%% Author Info 
% The |tmd_conlist| function and its documentation were written by Chad A.
% Greene, June 2022. 