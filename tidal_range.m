function R = tidal_range(filename_or_hc,conList,mask)
% tidal_range calculates the peak-to-peak tidal height range. (May take
% several hours.) 
% 
%% Syntax
% 
%  R = tidal_range(filename)
%  R = tidal_range(hc,conList,mask)
% 
%% Description
% 
% R = tidal_range(filename) returns the peak-to-peak tidal height range for
% a given TMD3.0 tide model file. Tidal range is calculated as the maximum
% minus minimum predicted values over a 14 day period, calculated at 30
% minute temporal resolution, including major and minor constituents. 
% 
% R = tidal_range(hc,conList,mask) computes the tidal range for the
% complex coefficients hc corresponding to the constituents conList, and a
% binary mask that is true for any pixels to be solved. hc should be MxNxP,
% where M and N are spatial dimensions and P is the number of constituents.
% conList is a Px1 cell array of consitutents, and mask is MxN. 
%
%% Example 
% 
% R = tidal_range('/Users/cgreene/data/CATS2008_update_2022-04-22.nc');
% 
%% Author Info 
% Written by Chad A. Greene, June 2022. 

%% 

% A month of 30 minute timesteps: 
t = (datenum(2000,1,1):1/48:datenum(2000,1,31))';

%% Load data

if isnumeric(filename_or_hc)
   assert(nargin==3,'If the first input is numeric, the inputs must be hRe, hIm, conList, mask.')
   hc = filename_or_hc; 
else
   conList = strsplit(ncreadatt(filename_or_hc,'cons','long_name')); 
   mask = tmd_data(filename_or_hc,'mask');
   hc = tmd_data(filename_or_hc,'h');
end

hc = cube2rect(hc,mask); % cube2rect is in Climate Data Toolbox for Matlab, included as a subfunction here. 

[astrol_s,astrol_h,astrol_p,astrol_N] = tmd_astrol(t);
[ispec,~,ph,omega,~] = tmd_constit(conList);

%% Solve

% Start with an assumption of zero tides: 
z_min = zeros(size(hc,2),1);
z_max = zeros(size(hc,2),1);

% Solve for each timestep:
w = waitbar(0,'Calculating tidal range (may take a few hours).');
for k = 1:length(t)

   % Major constiuents: 
   hhat = tmd_harp(t(k),hc,conList,astrol_p(k),astrol_N(k),ispec,ph,omega);

   % Minor constiuents:
   d_minor = tmd_InferMinor(hc,conList,t(k),astrol_s(k),astrol_h(k),astrol_p(k),astrol_N(k)); 

   % Just calculate min and max values (thus far) bc all we need is the total range:  
   z_max = max(z_max, d_minor + hhat); 
   z_min = min(z_min, d_minor + hhat); 
   
   waitbar(k/length(t),w,'Calculating tidal range (may take a few hours).')
end

close(w)

R = rect2cube((z_max-z_min)',mask); % rect2cube is in Climate Data Toolbox for Matlab, included as a subfunction here.
R(isnan(R)) = 0;

end


function A2 = cube2rect(A3,mask)
% cube2rect reshapes a 3D matrix for use with standard Matlab functions. 
% 
% This function enables easy, efficient, vectorized analysis instead of using 
% nested loops to operate on each row and column of a matrix. 
% 
%% Syntax
% 
%  A2 = cube2rect(A3) 
%  A2 = cube2rect(A3,mask)
% 
%% Description 
% 
% A2 = cube2rect(A3) for a 3D matrix A3 whose dimension correspond to space x space x time, 
% cube2rect reshapes the A3 into a 2D matrix A2 whose dimensions correspond to time x space. 
% 
% A2 = cube2rect(A3,mask) uses only the grid cells corresponding to true values in a 2D mask. 
% This option can reduce memory requirements for large datasets where some regions (perhaps all land
% or all ocean grid cells) can be neglected in processing. 
% 
%% Example 
% For examples type 
% 
%   cdt cube2rect
% 
%% Author Info
% This function was written by Chad A. Greene of the University of Texas 
% Institute for Geophysics (UTIG). 
% 
% See also: rect2cube, permute, and reshape.

%% Error checks: 

narginchk(1,2)
if nargin>1
   assert(isequal(size(mask),[size(A3,1),size(A3,2)]),'Error: The dimensions of the input mask must match the first two dimensions of A3') 
   assert(islogical(mask),'Error: mask must be logical.') 
end

% Reshape
A2 = permute(A3,[3 1 2]); % brings the "time" dimenion to dimension 1, which makes time go down the rows. 
A2 = reshape(A2,size(A2,1),size(A2,2)*size(A2,3)); % combines both spatial dimensions into one.  

% Remove masked-out values: 
if nargin>1
   A2 = A2(:,mask(:)); 
end

end

function A3 = rect2cube(A2,gridsize_or_mask) 
% rect2cube is the complement of cube2rect. It reshapes and permutes
% a 2D matrix into a 3D cube. 
% 
% This function enables easy, efficient, vectorized analysis instead of using 
% nested loops to operate on each row and column of a matrix. 
% 
%% Syntax
% 
%  A3 = rect2cube(A2,gridsize)
%  A3 = rect2cube(A2,mask)
% 
%% Description
% 
% A3 = rect2cube(A2,gridsize) reshapes 2D matrix A2 into a 3D matrix whose first
% two dimensions are spatial (e.g., lat x lon or lon x lat) and whose third 
% dimension is time or perhaps ocean depth or some variable along which operations
% are performed. The final dimensions of A3 are specified by gridsize, which 
% may be a complete 3 element array describing A3, or gridsize may be just a 2 
% element array containing the first two dimensions of A3. 
% 
% A3 = rect2cube(A2,mask) reshapes the elements of A2 into the true grid
% cells in a 2D matrix mask. 
% 
%% Examples
% For examples, type 
% 
%  cdt rect2cube
% 
%% Author Info
% This function was written by Chad A. Greene of the University 
% of Texas at Austin. 
% 
% See also: cube2rect, reshape, and permute. 

%% Error checks: 

assert(nargin==2,'Error: rect2cube requires 2 inputs.') 

%% Input parsing: 

% maskMethod means the user has input a 2D mask.  
if numel(gridsize_or_mask)<=3
   maskMethod = false; % the user has input only the sizes
   gridsize = gridsize_or_mask; 
   assert(gridsize(1)*gridsize(2)==size(A2,2),'Error: Dimensions of the gridsize must match the number of columns in A2.')
else
   maskMethod = true;  % the user has input a full 2D mask
   mask = gridsize_or_mask; 
   assert(islogical(mask),'Error: gridsize_or_mask must be either a 1D array describing grid size or a 2D logical mask.')
   gridsize = size(mask); 
end

if maskMethod
   
   % Create a full-size 2D NaN array: 
   if islogical(A2)
      A2full = false(size(A2,1),numel(mask)); 
   else
      A2full = NaN(size(A2,1),numel(mask)); 
   end
   
   % Fill in the "true" elements from the mask with A2: 
   A2full(:,mask) = A2; 
   
   % Unreshape back to original size of A: 
   A3 = ipermute(reshape(A2full,[],gridsize(1),gridsize(2)),[3 1 2]); 
   
else
   % Unreshape back to original size of A: 
   A3 = ipermute(reshape(A2,[],gridsize(1),gridsize(2)),[3 1 2]); 
end


end