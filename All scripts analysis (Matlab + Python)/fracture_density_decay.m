%% This script generates a fracture density decay with fault-perpendicular distance using a simplified rupture/fault trace and a fracture map

%% Inputs required
% shapefile of the main rupture trace or fault
% shapefile of the fracture map 

% External functions required and Mathworks download links
% ll2utm
% (https://www.mathworks.com/matlabcentral/fileexchange/45699-ll2utm-and-utm2ll)
% interparc
% (https://www.mathworks.com/matlabcentral/fileexchange/34874-interparc)

% This script requires the Matlab Mapping Toolbox

% This script outputs a two column text file with fault-perpendicular
% distance on column 1 and fracture density on column 2. This text file
% serves as the input for the Jupyter Notebook
% MCMC_density_decay_Poisson.ipynb where the decay is fit following an MCMC
% approach. 

%% load simplified main rupture trace

main_rupture = shaperead('mainrupture_mainshock.shp'); % load main trace shapefile 
pt_x_main = [];
pt_y_main = [];

% break-down the main rupture into evenly spaced points at 1m spacing
for n=1:numel(main_rupture)
    [pt_x_i,pt_y_i] = subdivide_points(main_rupture(n).X,main_rupture(n).Y) ;
    pt_x_main = [pt_x_main; pt_x_i]; % save discretized points to array
    pt_y_main = [pt_y_main; pt_y_i];
end 

%% load high-resolution fracture map

lines_secondary = shaperead('mainshock_fractures_Ridgecrest.shp'); % load fracture shapefile

D = [];

% break-down fractures into evenly spaced points at 1m spacing and measure
% distance to nearest point on discretized main rupture

for n=1:numel(lines_secondary)
    % generate spline of fracture and resample at 1 m spaced points
    [pt_x_i,pt_y_i] = subdivide_points(lines_secondary(n).X,lines_secondary(n).Y);
    % calculate Euclidean distance between each point in discretized line
    % and each point in main rupture
    mindist = pdist2([pt_x_i pt_y_i],[pt_x_main pt_y_main],'euclidean');
    % minimum of each row is the minimum distance between each point in
    % discretized line and the nearest point in the main rupture
    Di = min(mindist,[],2);
    D = [D; Di];
    disp(n) % keeps track of iteration #
end 

% locations = CFMRidgecrestdistancesHauksson(:,29:44);
% depthi = CFMRidgecrestdistancesHauksson(:,10);
% depth_idx = find(depthi<=5); 
% error_loc = CFMRidgecrestdistancesHauksson(:,23);
% min_loc = min(locations');
% min_loc = min_loc(depth_idx);
% error_loc = error_loc(depth_idx);
% idx = find(abs(error_loc)<0.1);
% min_loc = min_loc(idx);
% distance_EQ = min_loc;

%% generate and plot fault-perpendicular decay

% subdivide datasets into log-spaced bins
nbin = 100;
edges = logspace(0,log10(max([D])), nbin);
% note that if comparing several datasets, the maximum distance for binning
% should be determined by the dataset reching the furthest out.
% For the case of Ridgecrest, this is the aftershocks. Edges should be the
% same for all datasets to ensure consistent binning. 

hist_fr = histcounts(D, 'Binedges',edges);

% measure length of main rupture for density normalization 
L_main_rupture = zeros(1,length(main_rupture));
for n=1:numel(main_rupture)
L_main_rupture(n) = measure_length(main_rupture(n).X,main_rupture(n).Y);
end 

% normalize histogram by length of the main rupture
normalize_main = sum(L_main_rupture); % sum lengths of all segments on main rupture to estimate total rupture length
normalized_fr = (hist_fr./diff(edges))/normalize_main; % normalize each density by bin size and total rupture length(yields units fr/m^2)
normalized_fr(isnan(normalized_fr)) = 0; 
normalized_fr(isinf(normalized_fr)) = 0;

% normalization factors for sediment and bedrock Ridgecrest (0.7044 and
% 0.2956 respectively) from portion of rupture in material

% plot histogram of fracture density with fault-perpendicular distance
figure
h_fr = histogram('Binedges',edges,'BinCounts',normalized_fr,'FaceColor',[0.8000    0.8000    0.8000],'FaceAlpha',0.8)
set(gca,'YScale','log','XScale','log')
ylabel('Frequency')
xlabel('Distance away from fault (m)')
ax = gca;
ax.FontSize = 12; 
hold on 

% fit middle of each box in histogram to generate decay line
xvals_fr = (h_fr.BinEdges(2:end)+h_fr.BinEdges(1:end-1))/2; % find midpoint of each box in histogram
frval_fr = h_fr.Values; % retrieve value of each box in the histogram
figure
plot(xvals_fr, h_fr.Values, 'Color',[0.8510    0.3255    0.0980],'linewidth',1.5);
hold on 
set(gca,'YScale','log','XScale','log')
ylabel('Fractures/m^{2}')
xlabel('Distance away from fault (m)')
hold on
set(gca,'FontSize',14)

% Note that if the footprint over which the data was mapped is not even on
% the edges, the last bit of the decay that spans the distance with
% incomplete coverage may need to be cut off

%% output data in text file for MCMC in Jupyter Notebook (MCMC_density_decay_Poisson.ipynb)
% column 1: fault-perpendicular distance
% column 2: fracture density
yvals = h_fr.Values;
results = [xvals_fr' yvals'];
results = results(all(results,2),:); % removes empty rows

%writematrix(results, '.txt') 


%% function dumpster
% functions called in this script live here 
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
euclidean_matrix = sqrt((x_1-x_2).^2+(y_1-y_2).^2); % note transformation to local coordinate system 
L = sum(euclidean_matrix);
end 
function [pt_x,pt_y] = subdivide_points(fault_x,fault_y)
if length(fault_x)>0 % if statement to deal with empty lines in shapefile
 fault_x = fault_x(~isnan(fault_x));
 fault_y = fault_y(~isnan(fault_y));
if fault_y<90
[fault_x,fault_y]=ll2utm(fault_y,fault_x);
 else
end
 %% measure length
x_1 = fault_x(1:end-1);
x_2 = fault_x(2:end);
y_1 = fault_y(1:end-1);
y_2 = fault_y(2:end);
euclidean_matrix = sqrt((x_1-x_2).^2+(y_1-y_2).^2); % note transformation to local coordinate system
L = sum(euclidean_matrix);
%% create spline of fault 
spacing = 1; % meters
if length(fault_x) > 2
pt = interparc(0:(spacing/L):1,fault_x,fault_y,'linear'); 
pt_x = pt(:,1);
pt_y = pt(:,2);
%spacetest = sqrt((pt_x(1)-pt_x(2)).^2 + (pt_y(1)-pt_y(2)).^2) % check
%that spacing is 1 m 
else
    pt_x = fault_x(1);
    pt_y = fault_y(1);
end
end
end
