% Load HCP age data
load('HCP_individual_age.mat')

% Define group indices: D = Developmental, YA = Young Adult, A = Adult
idx_D = 1:633;         % Developmental group indices
idx_YA = 634:1028;     % Young Adult group indices  
idx_A = 1029:1321;     % Adult group indices


% Create group labels: 1 = Developmental/Adult, 2 = Young Adult
Group = zeros(size(Age));
Group(idx_D) = 1;Group(idx_YA) = 2;Group(idx_A) = 1;

% Load indicator data (e.g., brain connectivity metrics)
load('your_indicator_path.mat');
data = table(Age, Indicator, Group);
data.Group = categorical(data.Group);     % Convert group to categorical variable

% Perform linear regression: Indicator ~ Age + Group
mdl = fitlm(data, 'Indicator ~ Age + Group');

% Calculate age-adjusted indicator by removing systematic bias
% Add residuals back to age effect to get adjusted values
adjusted_Indicator = mdl.Residuals.Raw + mdl.Coefficients{"Age","Estimate"}*Age;

mdl.Coefficients("Age","pValue")
[~,idx] = sort(Age);
y_fit = mdl.Coefficients{"Age","Estimate"}*Age;

% Create visualization plot
% Scatter plot of age vs adjusted indicator
figure;
plot(Age,adjusted_Indicator,'.','color',[ 115 165 162]/255,'MarkerSize',6);
hold on;
plot(Age(idx),y_fit(idx),'color',[37   125    139]/255, 'LineWidth',0.8);
hold off;
axis square
set(gca, 'Box', 'off');
set(gcf,'unit','centimeters','position',[0.3 5 3.2 3.2])
set(gca,'FontUnits','points','FontName','Arial',"FontSize",7);