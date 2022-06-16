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
% This function was written by Chad A. Greene.

if nargin==0
   TMDFunctionName = 'tmd'; 
end

switch TMDFunctionName
   case 'tmd'
      showdemo TMD_main_page
      
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
      warning(['Sorry, there''s no fancy documentation for ',TMDFunctionName,'. Here''s what we''ve got, or you can try typing ',TMDFunctionName,'.'])
      showdemo TMD_main_page

end

