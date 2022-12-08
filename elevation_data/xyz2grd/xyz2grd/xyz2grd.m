%% file definition
clc;
clear;
xyz2grd_func('USGS_one_meter_x58y384_AR_R6_JeffersonCO_2015.xyz')

%% 

function xyz2grd_func(file)
% Function to convert  xyz file (obtained from Global mapper) into a 
% GRD file (Golden Software Surfer, ASCII format)
%
% xyz2grd(file)
%
% Input:
%      file = name of the file to be written (include ".xyz" extension)
% Output:
%      grd file in current directory
%
% Coded by :    Siva Srinivas Kolukula, PhD      
%               Indian Tsunami Early Warning Centre (ITEWC)
%               Advisory Services and Satellite Oceanography Group (ASG)
%               Indian National Centre for Ocean Information Services (INCOIS)
%               Hyderabad, INDIA
% E-mail   :    allwayzitzme@gmail.com                                        
% web-link :    https://sites.google.com/site/kolukulasivasrinivas/   

% Check for the input file format
[pathstr,name,ext] = fileparts(file) ;
if ~strcmpi(ext,'.xyz')
    error('Input file is not xyz')
end
data = load(file) ;                         % Load the file 
% Get longitude and latitude vectors 
x = unique(data(:,1)) ;
y = unique(data(:,2)) ;
z = unique(data(:,3)) ;
% dimensions of the Bathymetry 
nx = length(x) ; 
ny = length(y) ;
% Frame matrix of grid 
matrix = reshape(data(:,3),[ny,nx]) ;
% flip bathymetry matrix to adjust for plot
matrix = flipud(matrix) ; % Required if .xyz is from global mapper 
% Write the matrix into grd file 
grdfile = strcat(name,'.grd') ;
grd_write(matrix,min(x),max(x),min(y),max(y),grdfile)
end