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
      
   case {'predict','tmd_predict'}
      showdemo tmd_predict_documentation
      
   case {'interp','tmd_interp'}
      showdemo tmd_interp_documentation
      
   case {'data','tmd_data'}
      showdemo tmd_data_documentation
      
   case {'conlist','tmd_conlist'}
      showdemo tmd_extract_HC_documentation
      
   case {'ellipse','tmd_ellipse'}
      showdemo tmd_ellipse_documentation
      
   otherwise 
      warning(['Unrecognized TMD function ',TMDFunctionName,'.'])
      showdemo TMD_function_list

end

