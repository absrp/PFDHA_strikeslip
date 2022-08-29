%clear all; close all;
%% load event specific data
% from MonteCarlo density decay sampling
Markov_chain = load('Landers_chain.txt');
bestfit = load('Landers_best_fit_parameters.txt');
str = 'P_D_Do_error_Landers.pdf'; % figure name for pdf to be exported

% color associated with that event
c = [0.6353    0.0784    0.1843];

% from least-squares fit to inverse power-law fracture population
n = 1.48;
sigma = 0.17; % 0.17 0.12 0.08 0.06  L H E R


% user input Do (displacement threshold)
Do = 0.1;
% Error estimates 
% extract 1000 iterations of one walker 
% columns 1,2,3 for one walker 
onewalker = Markov_chain(:,1:3);

% draw 1000 random samples (rows) without replacement from MCMC Markov_chain
rng('default')
index = randsample(1:length(onewalker), 5000);
selection = onewalker(index,:);

vo = selection(:,1); 
d = selection(:,2); 
gamma = selection(:,3); 

% draw 1000 random samples from normally distributed n based on linear fit
R = normrnd(n,sigma,5000);

% displacements in meters (min equivalent to 1 m long fracture
Dmin = 10^-3;

% other parameters from MCMC fracture density decay
m = 2; % following Powers and Jordan (2010) 

% plot Montecarlo 
figure
subplot(2,4,[1 2])
x = linspace(1,20000,30000);

PD_Do = [];
hold on 

for i = 1:length(gamma) 
lambda = (vo(i) * (d(i).^m./(x.^m+d(i)^m)).^(gamma(i)/m)); 
PD = (Do^(1-R(i)))/(Dmin^(1-R(i)));
PD_Do= lambda*PD;
lh = plot(x,PD_Do,'Color',c);
lh.Color = [lh.Color 0.01];
disp(i)
end 

set(gca,'XScale','log','YScale','log')
%%
ylabel('P(D>D_{o})')
set(gca,'XScale','log','YScale','log','Xticklabel',[])
xlim([min(x) max(x)])

voR = bestfit(:,1); dR = bestfit(:,2); gammaR = bestfit(:,3);
lambda = (voR * (dR.^m./(x.^m+dR^m)).^(gammaR/m)); 
PD = (Do^(1-n))/(Dmin^(1-n));
PD_Dopref = lambda*PD;

plot(x,PD_Dopref,'Color',c,'linewidth',2)

subplot(2,4,3)
histogram(R,25,'FaceColor',[0.8000    0.8000    0.8000],'EdgeColor','none')
xlabel('n')
subplot(2,4,4)
histogram(vo,25,'FaceColor',[0.8000    0.8000    0.8000],'EdgeColor','none')
xlabel('v_{o}')
subplot(2,4,7)
histogram(gamma,25,'FaceColor',[0.8000    0.8000    0.8000],'EdgeColor','none')
xlabel('\gamma')
subplot(2,4,8)
histogram(d,25,'FaceColor',[0.8000    0.8000    0.8000],'EdgeColor','none')
xlabel('d')

% plot residuals 
subplot(2,4,[5 6])

res_struc = {};

hold on
for i = 1:length(gamma) 
lambda = (vo(i) * (d(i).^m./(x.^m+d(i)^m)).^(gamma(i)/m)); 
PD = (Do^(1-R(i)))/(Dmin^(1-R(i)));
PD_Do= lambda*PD;
lh = plot(x,log10(PD_Dopref)-log10(PD_Do),'Color',c);
res_struc{i} = log10(PD_Dopref)-log10(PD_Do);
lh.Color = [lh.Color 0.01];
end 

res = vertcat(res_struc{:});
res_sorted = sort(res,'ascend');

% pick 2.5 and 97.5 percentiles
prc_975 = prctile(res_sorted,97.5); 
prc_025 = prctile(res_sorted,2.5);
prc_034 = prctile(res_sorted,34); 
prc_068 = prctile(res_sorted,68);

hold on 
plot(x,log10(PD_Dopref)-log10(PD_Dopref),'Color',c,'linewidth',2)
xlabel('Fault-perpendicular distance (m)')
ylabel('\Delta log(P(D>D_{o}))')
plot(x,prc_975,'Color',c,'linewidth',1.5) 
plot(x,prc_025,'Color',c,'linewidth',1.5) 

% plot(x,prc_034,'Color',c,'linewidth',1.5) 
% hold on 
% plot(x,prc_068,'Color',c,'linewidth',1.5) 
set(gca,'XScale','log')
xlim([min(x) max(x)])

%% save uncertainty plots as pdf
print(gcf, '-dpdf', '-r300', str)
