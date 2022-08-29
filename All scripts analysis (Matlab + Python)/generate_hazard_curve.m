% this script requires:
% the ll2utm function (https://www.mathworks.com/matlabcentral/fileexchange/45699-ll2utm-and-utm2ll), 
% the logfit function (https://www.mathworks.com/matlabcentral/fileexchange/29545-power-law-exponential-and-logarithmic-fit),
% the polyparci function (https://www.mathworks.com/matlabcentral/fileexchange/39126-polyparci?s_tid=FX_rc2_behav)
% and the Matlab Mapping Toolbox

% this script also requires four shapefiles to run:
% the shapefile for the Ridgecrest fracture maps is available at https://sandbox.zenodo.org/record/902426#.YSAyW9NKiDU
% rupture maps for the remaining events are available from the FDHI rupture database appendix (https://www.risksciences.ucla.edu/girs-reports/2021/08)

% the referenced equations are in Rodriguez Padilla and Oskin (2022)
%% load fracture maps 
close all; clear all;

EMC = shaperead('EMC_noN_noS_radar.shp'); 
Landers = shaperead('Landers_secondary_fractures.shp');
HM = shaperead('HM_secondary.shp'); 
Ridgecrest = shaperead('Ridgecrest_highres_maps.shp'); 

%% measure length of secondary fractures in each population
L_EMC = []; 
L_Landers = [];
L_HM = [];
L_Ridgecrest = []; 

for n = 1:length(EMC) 
  L_EMC(n) =  measure_length(EMC(n).X,EMC(n).Y); 
end 
for n = 1:length(Landers) 
  L_Landers(n) =  measure_length(Landers(n).X,Landers(n).Y); 
end 
for n = 1:length(HM) 
  L_HM(n) =  measure_length(HM(n).X,HM(n).Y); 
end 
for n = 1:length(Ridgecrest) 
  L_Ridgecrest(n) =  measure_length(Ridgecrest(n).X,Ridgecrest(n).Y); 
end 

%% generate fracture distribution 
% equations 1-2

% minimum and maximum length of completeness estimated visually from
% histograms
Lo_EMC = 20;
Lo_Ridgecrest = 2;
Lo_Landers = 60;
Lo_HM = 8;

Lmax_EMC = 1100;
Lmax_Ridgecrest = 93;
Lmax_Landers = 1500;
Lmax_HM = 150;

figure
[struct_length_distribution_EMC] = length_distribution(L_EMC,Lo_EMC,Lmax_EMC,[0    0.6000    0.6000],1);
[struct_length_distribution_Landers] = length_distribution(L_Landers,Lo_Landers,Lmax_Landers,[0.6353    0.0784    0.1843],2);
[struct_length_distribution_HM] = length_distribution(L_HM,Lo_HM,Lmax_HM,[0.1647    0.3843    0.2745],3);
[struct_length_distribution_Ridgecrest] = length_distribution(L_Ridgecrest,Lo_Ridgecrest,Lmax_Ridgecrest,[0.8706    0.4902         0],4); 

hold on 
subplot(2,2,2)
text(Lo_Landers*1.1,3*10^3,'L_{o}','Color',[0.6353    0.0784    0.1843],'FontSize',10)
text(Lmax_Landers*1.1,3*10^3,'L_{max}','Color',[0.6353    0.0784    0.1843],'FontSize',10)

%saveas(gcf,'Lpopulations.pdf')

% fit equation 1 using least-squares method

[struct_length_fit_EMC] = population_fit(struct_length_distribution_EMC{1,4},struct_length_distribution_EMC{1,3});
[struct_length_fit_Landers] = population_fit(struct_length_distribution_Landers{1,4},struct_length_distribution_Landers{1,3});
[struct_length_fit_HM] = population_fit(struct_length_distribution_HM{1,4},struct_length_distribution_HM{1,3});
[struct_length_fit_Ridgecrest] = population_fit(struct_length_distribution_Ridgecrest{1,4},struct_length_distribution_Ridgecrest{1,3});

%% generate displacement density decay
% parameters from MCMC fit in Python following method in Rodriguez Padilla
% et al. (2022) and available in fracture_density_decay folder

% Dmin fixed for minimum fracture length of 1 meter
Dmin = 10^-3;

% equation 10

Do = [0.1,0.5,1];
symbol = {':','-','--'};

x = 1:0.05:20000; 
m = 2;

% parameters for equation 8 (second term in equation 10) obtained using MCMC code in Python folder

bestfit_EMC = load('EMC_noradar_best_fit_parameters.txt'); 
bestfit_HM = load('HM_best_fit_parameters.txt'); 
bestfit_Landers = load('Landers_best_fit_parameters.txt'); 
bestfit_Ridgecrest = load('Ridgecrest_imagery_best_fit_parameters.txt'); 

v_o_EMC = bestfit_EMC(:,1); d_EMC = bestfit_EMC(:,2); gamma_EMC = bestfit_EMC(:,3); 
v_o_HM = bestfit_HM(:,1); d_HM = bestfit_HM(:,2); gamma_HM = bestfit_HM(:,3); 
v_o_Landers = bestfit_Landers(:,1); d_Landers = bestfit_Landers(:,2); gamma_Landers = bestfit_Landers(:,3); 
v_o_Ridgecrest = bestfit_Ridgecrest(:,1); d_Ridgecrest = bestfit_Ridgecrest(:,2); gamma_Ridgecrest = bestfit_Ridgecrest(:,3);

figure 

for n=1:length(Do)
estimate_PD_decay(Do(n),Dmin,struct_length_fit_EMC{1,3},x,m,v_o_EMC,d_EMC,gamma_EMC,[0    0.6000    0.6000],symbol{n});
estimate_PD_decay(Do(n),Dmin,struct_length_fit_Landers{1,3},x,m,v_o_Landers,d_Landers,gamma_Landers,[0.6353    0.0784    0.1843],symbol{n});
estimate_PD_decay(Do(n),Dmin,struct_length_fit_HM{1,3},x,m,v_o_HM,d_HM,gamma_HM,[0.1647    0.3843    0.2745],symbol{n});
estimate_PD_decay(Do(n),Dmin,struct_length_fit_Ridgecrest{1,3},x,m,v_o_Ridgecrest,d_Ridgecrest,gamma_Ridgecrest,[0.8706    0.4902         0],symbol{n});
end 

hold on 
f1 = plot(nan,nan,'LineStyle',':','Color',[0.2 0.2 0.2])
f2 = plot(nan,nan,'LineStyle','-','Color',[0.2 0.2 0.2])
f3 = plot(nan,nan,'LineStyle','--','Color',[0.2 0.2 0.2])

legend([f1,f2,f3],'D_{o} = 0.1 m','D_{o} = 0.5 m','D_{o} = 1 m')

ylabel('P(D>D_{o})')
xlabel('Fault perpendicular distance (m)')
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'FontSize',14)
xlim([min(x) max(x)])

%saveas(gcf,'P_D_Do_allsample.pdf')

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

function [struct_fracture_distribution] = length_distribution(L,Lo,Lmax,d,nsub)
% bin fracture lengths
nbins = 30;
bins = logspace(0,log10(max(L)),nbins); 
% make histogram of fracture lengths
subplot(2,2,nsub)
h = histogram(L,bins,'FaceColor',d,'EdgeColor','none','EdgeAlpha',0.2,'FaceAlpha',0.6)
hold on 
xline(Lo,'Color',d,'LineStyle',':','LineWidth', 1.5) 
hold on 
xline(Lmax,'Color',d,'LineStyle',':','LineWidth', 1.5) 
set(gca,'XScale','log','YScale','log')
ylabel('Frequency')
xlabel('Fracture length (m)')
set(gca,'FontSize',12)
ylim([1 5*10^3])
xlim([0 5*10^3])
% account for fracture being discretized into 1 m spaced segments in density decays
freq = h.Values;
% generate displacement PMF from fracture length histogram
midpoint = (h.BinEdges(2:end)+h.BinEdges(1:end-1))/2;
idx = find(midpoint > Lo & midpoint < Lmax); 
midpoint = midpoint(idx);
freq = freq(idx);% 10^-4 is the constant relating fracture length and maximum displacement/aperture per LEFM
struct_fracture_distribution = {bins, h, freq, midpoint};
end 

function [struct_fracture_fit] = population_fit(lfrac,freq)
figure
x = lfrac; 
y = freq;
% fit power-law through truncated displacement 
[slope, intercept,MSE, R2, S]  = logfit(x,y,'powerlaw');
[yfit,std_error] = polyval([slope, intercept],x,S);
ci = polyparci([slope, intercept],S,0.95);
error = ci(1)-slope
struct_fracture_fit = {x,y,slope,intercept,yfit,error};
set(gca,'YScale','log','XScale','log')
end 

function [] = estimate_PD_decay(Do,Dmin,n,x,m,v_o,d,gamma,c,style)
n = abs(n); % hack to make exponent + since negative sign was accounted for in derivation
PD = (Do/Dmin)^(1-n);
lambda = (v_o * (d.^m./(x.^m+d^m)).^(gamma/m)); 
PD_Do= lambda*PD;
plot(x,PD_Do,'linewidth',1.5,'LineStyle',style,'Color',c)
set(gca,'YScale','log','XScale','log')
hold on 
end