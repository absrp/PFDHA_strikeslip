% this script requires the file 'FDHI_data.xlsx'
% this script fits a linear regression through the mean displacement data
% for strike-slip events in the FDHI database. The mean slip is calculated
% for each event, and then the regresion is done using package curvefit
clear all; close all;

%% load data 
dataslip = readtable('data_FDHI.xlsx');
dataSRL = readtable()


% create mean slip and magnitude variables for regression
meanslip_SS = []; 
mag_SS = [];
event_names = [];

%% calculate mean displacement per event

% select style strike-slip only
style = data.style;
style_select = find(strcmp(style,'Strike-Slip'));
strike_slip = data(style_select,:);

% group per event
eventID = strike_slip.EQ_ID; 
uniqueIDs = unique(eventID);

for n=1:length(uniqueIDs)
    % select all rows associated with the event 
    eventIDi = uniqueIDs(n);
    idx_event = find(strike_slip.EQ_ID == eventIDi);
    subset_event = strike_slip(idx_event,:);
    event_name = subset_event.eq_name;
    event_name = event_name(1);

    % select all rows associated with field and preferred displacements
    % only
    % field = find(strcmp(subset_event.fps_meas_type,'field'));
    % subset_data = subset_event(field,:);
    slip = subset_event.recommended_net_preferred_for_analysis_meters;
    slipidx = find(slip>0); 
    slip = slip(slipidx);

    % calculate mean slip 
    mean_slip_event = mean(slip);

    % get event magnitude
    mag = subset_event.magnitude;
    mag = mag(1); 

    mag_SS = [mag_SS;mag];
    meanslip_SS = [meanslip_SS;mean_slip_event];
    event_names = [event_names; event_name];

end 

%% fit regression through data 

% Take the base 10 logarithm of 'meanslip_SS'
logged_meanslip_SS = log10(meanslip_SS);

% Fit a linear regression model
[coefficients,S] = polyfit(mag_SS, logged_meanslip_SS, 1);

% Compute the covariance using the coefficients
covariance = coefficients(2) * var(mag_SS);

% Extract the slope (m) and intercept (b) of the linear regression line
m = coefficients(1);
b = coefficients(2);

% Create a function for the fitted line
fitted_line = @(x) m*x+b;

% Generate points for the fitted line
x_fit = linspace(min(mag_SS), max(mag_SS), 100);
y_fit = fitted_line(x_fit); % expected value of y given regression

% standard error of the slope and intercept
residuals = logged_meanslip_SS - polyval(coefficients, mag_SS); % residuals of each expected value and preferred y
residuals_sq = residuals.^2;
sum_residuals_sq = sum(residuals_sq);
N = length(mag_SS);
mse = sum_residuals_sq/(N-2);
rmse = sqrt(mse);
 
mean_mag = mean(mag_SS);
diff_sqrs = (mag_SS-mean_mag).^2;
sum_diff_sqrs = sum(diff_sqrs);
other = sqrt((1/N)+(mean_mag^2)/sum_diff_sqrs^2);

se_intercept = rmse*other;
se_slope = sqrt(mse/sum_diff_sqrs);

% plot standard error as a function of magnitude
mwcheck = linspace(5,8,100);

semw = [];

for n=1:length(mwcheck)
other_mw = sqrt((1/N)+((mwcheck(n)-mean_mag)^2)/sum_diff_sqrs^2);
se_mw = rmse*other_mw;
semw = [semw; se_mw];
end

figure
plot(mwcheck,semw,'linewidth',4)
ylabel('Standard error')
xlabel('Magni')


figure;
subplot(3,1,[1 2])
scatter(mag_SS, logged_meanslip_SS, 'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','none');
hold on;
plot(x_fit, y_fit, 'Color',[0.6353    0.0784    0.1843], 'LineWidth', 2);
set(gca,'FontSize',14)
hold off;

xlabel('Magnitude');
ylabel(' Mean Slip (log_{10})');
% Create the equation string with the regression coefficients
equation_str = sprintf('log_{10}(MS) = %.2f x %.2f', m, b);

% Add the equation to the legend
legend({'FDHI data', [equation_str]}, 'Location', 'best','fontsize',14);

subplot(3,1,3)
scatter(mag_SS, residuals,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','none');
hold on;
plot(mag_SS, zeros(size(mag_SS)), 'Color',[0.6353    0.0784    0.1843], 'LineWidth', 2); % Plot a horizontal line at y = 0
xlabel('Magnitude');
ylabel('Residuals');
set(gca,'FontSize',14)
