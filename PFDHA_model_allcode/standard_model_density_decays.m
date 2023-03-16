%% standard model, using average of all events 

% External functions required: 
% wgs2utm
% interparc
% distance2curve


% all the inputs here are derived in the other scripts
clear all; close all;

%% distances for each event
D_Landers = []; 
D_HectorMine = []; 
D_EMC = []; 
D_Ridgecrest1 = []; 
D_Ridgecrest2 = []; 
Dtotal = [];

shapefiles_secondary = {'Landers_secondary_fractures.shp', 'HectorMine_secondary_fractures.shp',...
    'EMC_secondary_fractures.shp', 'Ridgecrest1_secondary_fractures.shp',...
    'Ridgecrest2_secondary_fractures.shp'};

shapefiles_main = {'Landers_main_rupture.shp', 'HectorMine_main_rupture.shp',...
    'EMC_main_rupture.shp', 'Ridgecrest1_main_rupture.shp',...
    'Ridgecrest2_main_rupture.shp'};

for n = 1:length(shapefiles_secondary)
    shapefiles_secondaryi = shapefiles_secondary{n};
    shapefiles_maini = shapefiles_main{n};
    secondarylines = shaperead(shapefiles_secondaryi);
    mainlines = shaperead(shapefiles_maini);
    Dit = measure_distances(secondarylines,mainlines);
    Dtotal = [Dtotal; Dit];
    disp(n)
end 

%% rupture length of each event
L_Ridgecrest1 = measure_lengthi(shaperead('Ridgecrest1_main_rupture.shp')); 
L_Ridgecrest2 = measure_lengthi(shaperead('Ridgecrest2_main_rupture.shp')); 
L_Landers = measure_lengthi(shaperead('Landers_main_rupture.shp')); 
L_HectorMine = measure_lengthi(shaperead('HectorMine_main_rupture.shp')); 
L_EMC =  measure_lengthi(shaperead('EMC_main_rupture.shp')); 


L_FDHI = L_Ridgecrest1+ L_Ridgecrest2 + L_Landers + L_HectorMine + L_EMC; 


%% calculate decay 

% subdivide datasets into log-spaced bins
nbin = 100;
edges = logspace(0,log10(max(Dtotal)), nbin);
% note that if comparing several datasets, the maximum distance for binning
% should be determined by the dataset reching the furthest out.
% For the case of Ridgecrest, this is the aftershocks. Edges should be the
% same for all datasets to ensure consistent binning. 

hist_fr = histcounts(Dtotal, 'Binedges',edges);

% normalize histogram by length of the main rupture
normalize_main = L_FDHI; % sum lengths of all segments on main rupture to estimate total rupture length
normalized_fr = (hist_fr./diff(edges))/(L_FDHI); % normalize each density by bin size and total rupture length(yields units fr/m^2)
normalized_fr(isnan(normalized_fr)) = 0; 
normalized_fr(isinf(normalized_fr)) = 0;


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

% output data in text file for MCMC in Jupyter Notebook (MCMC_density_decay_Poisson.ipynb)
% column 1: fault-perpendicular distance
% column 2: fracture density
yvals = h_fr.Values;
results = [xvals_fr' yvals'];
results = results(all(results,2),:); % removes empty rows

writematrix(results, 'general_density_decay.txt')  ;


%% function dumpster

function D = measure_distances(lines_secondary,main_rupture)
pt_x = [];
pt_y = [];

for n=1:numel(lines_secondary)
%     % generate spline of fracture and resample at 1 m spaced points
    [pt_x_i,pt_y_i] = subdivide_points(lines_secondary(n).X,lines_secondary(n).Y);
    pt_x = [pt_x; pt_x_i];
    pt_y = [pt_y; pt_y_i];
end


    for n=1:length(main_rupture)
        [coords_refx, coords_refy] =  wgs2utm(main_rupture(n).Y,main_rupture(n).X,11,'N');
        curvexy = [coords_refx' coords_refy'];
        curvexy = rmmissing(curvexy);
        [xy,distance(:,n),t_a] = distance2curve(curvexy,[pt_x pt_y],'linear');
    end
    
 D = min(distance,[],2);
 
end

function L_main_rupturei = measure_lengthi(main_rupture)

L_main_rupture = zeros(1,length(main_rupture));

for n=1:numel(main_rupture)
L_main_rupture(n) = measure_length(main_rupture(n).X,main_rupture(n).Y);
end 

L_main_rupturei = sum(L_main_rupture);

end 

function [L] = measure_length(fault_x,fault_y)
fault_x = fault_x(~isnan(fault_x)); % removes NaN artifact at end of each fault in shapefile
fault_y =fault_y(~isnan(fault_y));
[fault_x,fault_y]=wgs2utm(fault_y,fault_x,11,'N');
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
[fault_x,fault_y]=wgs2utm(fault_y,fault_x,11,'N');
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

