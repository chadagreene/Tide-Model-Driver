function conList = tmd_conlist(filename)
% tmd_conlist returns a list of tidal constituents in a TMD3.0 compatible
% consolidated NetCDF tide model file. 
% 
%% Syntax
% 
%  conList = tmd_conlist(filename)
% 
%% Description 
% 
% conList = tmd_conlist(filename) returns a cell array of tidal
% constituents in the specified model file. 
% 
%% Example 
% 
% 
%% Author Info
% This function was written by Chad A. Greene, June 2022. 
% 
% See also: tmd_data and tmd_interp. 

%% Error checks 

assert(contains(filename,'.nc'),'Input filename must end in .nc.')
assert(exist(filename,'file'),['Cannot find ',filename,'. Check the path and try again.'])

% Ensure the model file is TMD3.0 compatible: 
try 
   tmd_version = ncreadatt(filename,'/','tmd_version'); 
   assert(tmd_version>=3.0,[filename,' is not compatible with TMD3.0+.'])
catch 
   error([filename,' is not compatible with TMD3.0+.'])
end

%% Get constituents

conList = strsplit(ncreadatt(filename,'constituents','constituent_order')); 

end

