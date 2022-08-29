%% This script estimates the scaling of Dmax with Lmax for the surface ruptures in Southern California
%% contained in the FDHI database 
% this script requires:
% the ll2utm function (https://www.mathworks.com/matlabcentral/fileexchange/45699-ll2utm-and-utm2ll), 
% and the Matlab Mapping Toolbox

% this script also requires:
% shapefiles of the primary fractures for CA events contained in the raw_data folder
% the spreadsheet 'data_FDHI.xlsx' containing info from the FDHI database, also contained 
% in the raw_data folder

% this script must be run in the raw_data shapefiles_CA directory

clear all; close all 

%% select CA events from FDHI database 
data = readtable('data_FDHI.xlsx'); 
type = data.fps_meas_type;
idxb =  find(strcmp(type,'field')); % select field data only
data = data(idxb,:); 

ID = data.EQ_ID; 
event = data.eq_name; 
slip = data.recommended_net_preferred_for_analysis_meters;
mag = data.magnitude;
region = data.region;

idx_CA = find(strcmp(region,'California'));
slip = slip(idx_CA);
slip = slip(~isnan(slip));
mag = mag(idx_CA);
mag = mag(~isnan(mag));
region = region(idx_CA); 
event = event(idx_CA); 
ID = ID(idx_CA);
ID = ID(~isnan(ID)); 

IDn = unique(ID); 

%% calculate maximum length and displacement per event 
shapefiles = dir('*.shp'); % access all shapefile names in the folder

L_max = [];
nametrack = {}; 
max_slip = [];

for i=1:length(shapefiles)
% read shapefile name
shapename = shapefiles(i).name;
fullname = strsplit(shapename,{'_','.'}); % string containing shapefile name
name = fullname{1};
lines = shaperead(shapename);

% calculate length of all fractures in shapefile
L = [];

for n = 1:length(lines) 
  L(n) =  measure_length(lines(n).X,lines(n).Y); 
end

% save length of longest fracture
L_max(i) = max(L); 
nametrack{i} = {name}; 

% find maximum slip for event from the FDHI database
idx = find(strcmp(event,name));
max_slip(i) = max(slip(idx));
end

%% compile data into table
maxdata = table(nametrack', L_max', max_slip') ;
maxdata.Properties.VariableNames = {'Name',...
    'L_max',...
    'Max_slip'};

%% note EMC must be added later 
data = readtable('data_FDHI.xlsx');
name = data.eq_name; 

EMC = shaperead('EMC_secondary_fractures.shp'); 
for n = 1:length(EMC) 
  L_EMC(n) =  measure_length(EMC(n).X,EMC(n).Y); 
end 
L_max_EMC = max(L_EMC); 

idx = find(strcmp(name,'EMC'));
data = data(idx,:);

type = data.fps_meas_type;
idxb =  find(strcmp(type,'field'));
data = data(idxb,:); 

slip = data.recommended_net_preferred_for_analysis_meters;
Dmax_EMC = max(slip); 

%% add EMC to max data 
EMC_T = table({'EMC'},L_max_EMC, Dmax_EMC); 
EMC_T.Properties.VariableNames = {'Name',...
    'L_max',...
    'Max_slip'};
Tdata_all = [maxdata;EMC_T];


%% get scaling ratio
Scaling_ratio = log10(Tdata_all.Max_slip./Tdata_all.L_max);
mean_Scaling_ratio = mean(Scaling_ratio);
std_Scaling_ratio = std(Scaling_ratio);
var_Scaling_ratio = var(Scaling_ratio); 
Scaling_ratio = table(Scaling_ratio); 
Scaling_ratio.Properties.VariableNames = {'Scaling_ratio'};

% add scaling ratio to table
Tdata_all = [Tdata_all,Scaling_ratio]; 
Tdata_all.Properties.VariableNames = {'Name',...
    'L_max',...
    'Max_slip','Scaling_ratio'};

%% write table as csv file 
%writetable(Tdata_all,'max_slip_length_scaling.csv') 

%% function dumpster
function [L] = measure_length(fault_x,fault_y)
fault_x = fault_x(~isnan(fault_x)); % removes NaN artifact at end of each fault in shapefile
fault_y =fault_y(~isnan(fault_y));
if fault_y<90
[fault_x,fault_y]=ll2utm(fault_y,fault_x);
else 
end
% calculate length
x_1 = fault_x(1:end-1);
x_2 = fault_x(2:end);
y_1 = fault_y(1:end-1);
y_2 = fault_y(2:end);
segment_length = sqrt((x_1-x_2).^2+(y_1-y_2).^2); % note transformation to local coordinate system 
L = sum(segment_length);
end 

