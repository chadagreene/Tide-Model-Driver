function tmd(TMDFunctionName)
% tmd displays help files for the Tide Model Driver (TMD) software. 
% 
%% Syntax
% 
%  tmd
%  tmd(TMDFunctionName) 
% 
%% Description 
% 
% tmd displays the TMD Getting Started page. 
% 
% tmd(TMDFunctionName) displays help for any of the TMD function names. 
% 
%% Example
% To display help for the tmd_tide_pred function simply type this into 
% your Command Window: 
% 
%  tmd tmd_predict
% 
%% Author Info
% TMD version 3.00; August 2018; Chad A. Greene (chad@chadagreene.com). 
% Modified from TMD version 2.04 & 2.05 by Lana Erofeeva (Oregon State 
% University) and Laurie Padman (Earth & Space Research). 

if nargin==0
   TMDFunctionName = 'tmd'; 
end

switch TMDFunctionName
   case 'tmd'
      showdemo TMD_function_list
      
   case {'tide_pred','tmd_tide_pred'}
      showdemo tmd_tide_pred_documentation
      
   case {'get_coeff','tmd_get_coeff'}
      showdemo tmd_get_coeff_documentation
      
   case {'get_ellipse','tmd_get_ellipse'}
      showdemo tmd_get_ellipse_documentation
      
   case {'extract_HC','tmd_extract_HC'}
      showdemo tmd_extract_HC_documentation
      
   case {'get_bathy','tmd_get_bathy'}
      showdemo tmd_get_bathy_documentation
      
   otherwise 
      warning(['Unrecognized TMD function ',TMDFunctionName,'.'])
      showdemo TMD_function_list

end

