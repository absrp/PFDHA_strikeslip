% this script requires:

% External functions required 
% wgs2utm
% interparc
% distance2curve

% this script requires the Matlab Mapping Toolbox

% this script requires five shapefiles from the FDHI rupture database
% appendix (https://www.risksciences.ucla.edu/girs-reports/2021/08) and the
% main rupture shapefiles from the appendix of Rodriguez Padilla and Oskin
% (202X). We also use the file 'data_FDHI.xlsx', which contains the
% data from the database, including displacement. 

%% generate displacement distribution
close all; clear;

data = readtable('data_FDHI.xlsx');

events = {'Landers','EMC', 'HectorMine','Ridgecrest1','Ridgecrest2'}; 


  %  figure
for i=1:length(events)
    event = events{i};
    
    % extract pre-saved color and file name choices for event
    info = event_info(event);
    c = info{1}; % event color
    str1 = info{2}; % event name for titles
    str2 = info{3}; % event-based name for pdf saving figure 2
    str3 = info{4}; % event-based name for pdf saving figure 3

    % subset spreadsheet to event data 
    name = data.eq_name; 
    idx = find(strcmp(name,event));
    subset_data = data(idx,:);
    type = subset_data.fps_meas_type;
    field = find(strcmp(type,'field'));
    subset_data = subset_data(field,:);
    slip = subset_data.recommended_net_preferred_for_analysis_meters;
    slipidx = find(slip>0); 
    slip = slip(slipidx);
    coordsx = subset_data.longitude_degrees(slipidx,:);
    coordsy = subset_data.latitude_degrees(slipidx,:); 
    EQ_style = subset_data.style;
    measurement_style = subset_data.fps_style;
    [coords_refx, coords_refy] = wgs2utm(coordsy,coordsx,11,'N');
    coords_ref = [coords_refx coords_refy];
    
    % load reference primary fault trace
    strname = '_main_rupture.shp';
    combined_str_main = append(event,strname);
    main_rupture = shaperead(combined_str_main); 
    
    distance = zeros(length(coords_ref),length(main_rupture)); 
 

    if i == 5
        p == 4;
    else 
        p = i;
    end
   % 
   %  subplot(2,2,p)
   %  scatter(coordsx,coordsy,20,'MarkerFaceColor',c,'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none') 
   %  hold on
   %  for n=1:length(main_rupture)
   %      plot(main_rupture(n).X,main_rupture(n).Y,'k','linewidth',1.5)
   %      hold on
   %      %set(gca,'FontSize',14)
   %  end
   %  axis equal
   %  % box on
   %  xlabel('Lon')
   %  ylabel('Lat')
   % 
   %  if p == 1
   %      text(-116.95, 34.68, 'Landers', 'VerticalAlignment', 'top')
   %  elseif p == 2
   %     text(-115.5, 32.715, 'El Mayor-Cucapah', 'VerticalAlignment', 'top')
   %  elseif p == 3
   %     text(-116.2, 34.76, 'Hector Mine', 'VerticalAlignment', 'top')
   %  elseif p == 4
   %     text(-117.59, 35.8, 'Mainshock', 'VerticalAlignment', 'top','Color', [0.4941    0.1843    0.5569])
   %     text(-117.7, 35.659, 'Foreshock', 'VerticalAlignment', 'top','Color', [0.8706    0.4902         0])
   %     text(-117.46, 35.92, 'Ridgecrest', 'VerticalAlignment', 'top')
   %  end
   % 
   % namestrfig = ('displacement_map_ref.pdf');
   % saveas(gcf,namestrfig);

    for n=1:length(main_rupture)
        [curvexyx, curvexyy] = wgs2utm(main_rupture(n).Y,main_rupture(n).X,11,'N');
        curvexy = [curvexyx' curvexyy'];
        curvexy = rmmissing(curvexy); 
        [xy,distance(:,n),t_a] = distance2curve(curvexy,coords_ref,'linear');
    end
    
    
    dist = min(distance,[],2); 
    
    
    % bin displacement-decay data and get mean and std
    nbin = 40;
    edges = logspace(0,log10(max(dist)),nbin);
    % find indices in each bin 
    
    [N,distvals,bins] = histcounts(dist,edges); % bins contains the bin index for each dist value
    centers = (edges(1:end-1) + edges(2:end)) / 2;
    
    % for each value of bin, calculate average and standard deviation of
    % slip magnitude 
    binsID = unique(bins); 
    
    meanslip = [];
    stdslip = [];
    medianslip = [];
    xval = [];
    xval_all = [];
    slip_all = [];
    
    for s = 1:length(binsID)-1
       idxbins =  find(bins == binsID(s)); 
       slipbins = slip(idxbins); 
%        if s<7
%        figure
%        histogram(slipbins)
%        title(distvals(s))
%        end
       meanslip(s) = mean(slipbins); 
       stdslip(s) = std(slipbins); 
       medianslip(s) = median(slipbins); 
       xval(s) = centers(s);
       xval_alli = repelem(xval(s),length(slipbins));
       xval_all = [xval_all; xval_alli'];
       slip_all = [slip_all; slipbins];
    end 
    
       
      strbase = 'all_displacements_binned.csv';
      strexport = strcat(event,'_',strbase);
      tablei = table(xval_all, slip_all,'VariableNames', {'x', 'disp'}); 
      writetable(tablei,strexport);

    
    figure(1)
        if i<5
        subplot(3,2,i)
    elseif i==5
        subplot(3,2,5.5)  
    else
        end

    scatter(xval,meanslip,'Filled')
    hold on
    scatter(xval,stdslip,'filled','MarkerFaceAlpha',0.3)
    legend('Mean','Standard deviation')
    title(event)
    set(gca,'XScale','log','YScale','log')
    ylabel('Slip (m)') 
    xlabel('Fault-perpendicular distance (m)')
    if i==2
        ylim([5*10^-3 3])
    else
    end

%     figure(2)
%     if i<5
%         subplot(3,2,i)
%     elseif i==5
%         subplot(3,2,5.5)  
%     else
%     end
%         
%     scatter(xval,meanslip,'MarkerFaceColor',c,'MarkerEdgeColor','none')
%     hold on
%     title(str1)
%     set(gca,'XScale','log','YScale','log')
%     ylabel('Slip (m)') 
%     xlabel('Fault-perpendicular distance (m)')
%  
%     % fit pwl through xval, meanslip
%     
%      % Define Start points, fit-function and fit curve
% 
%    
      strbase = 'mean_displacement_decay.txt';
      strexport = strcat(event,'_',strbase);
      writematrix([xval',meanslip'],strexport);
%     
    
end 

saveas(gcf,'dispmap.pdf');


function event_info = event_info(event) 
if strcmp(event,'Landers')
    c = [0.6353    0.0784    0.1843];
    str1 = 'Landers';
    str2 = 'Landers_Scracks.txt';
    str3 = 'Landers_Dpoints.txt';
elseif strcmp(event,'HectorMine')
    c = [0.1647    0.3843    0.2745];
    str1 = 'Hector Mine';
    str2 = 'HectorMine_Scracks.txt';
    str3 = 'HectorMine_Dpoints.txt';
elseif strcmp(event,'Ridgecrest1')
    str1 = 'Ridgecrest foreshock';
    c = [0.8706    0.4902         0];
    str2 = 'Ridgecrest1_Scracks.txt';
    str3 = 'Ridgecrest1_Dpoints.txt';
elseif strcmp(event,'Ridgecrest2')
    str1 = 'Ridgecrest mainshock';
    c = [0.4941    0.1843    0.5569];
    str2 = 'Ridgecrest2_Scracks.txt';
    str3 = 'Ridgecrest2_Dpoints.txt';
elseif strcmp(event,'EMC')
    c = [0    0.6000    0.6000];
    str1 = 'El Mayor Cucapah';
    str2 = 'EMC_Scracks.txt';
    str3 = 'EMC_Dpoints.txt';
else
    slip('Event name does not match database name')
end
    event_info = {c str1 str2 str3};
end 

