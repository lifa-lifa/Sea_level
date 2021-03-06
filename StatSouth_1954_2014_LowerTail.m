%% Oeresund A102102
% * LIFA 2017-09-06
% * Statistics for southern storms, lower tail
close all
clear all
clc

%% Get data
ReadFolder = strcat(pwd,'\input\');  % data folder
dirList = dir(strcat(ReadFolder,'*.txt'));  %list all txt files
fReadName = dirList(1).name; 
full_name = [ReadFolder fReadName];
fid = fopen(full_name, 'r');
data_temp = textscan(fid, '%f', 'Delimiter', ',');
fclose(fid);
data_temp = cell2mat(data_temp);  % convert cell array to matrix


%% Inputs
rng default;  % fixes the rng random seed to Mersenne Twister 5489
threshold = 111;
alpha = 0.10;  % significance level
nSim = 10000;  % number of simulations, typ. 10000

% return periods and years
yearT = [10, 20, 50, 100, 250, 500, ...
         1000, 2000, 3000, 5000, 10000, 100000, 1000000]';  % return periods
yearT_plot = [50, 100, 250, 1000, 10000]';  % return periods for plotting
obsYearsRange = 2016 - 1956 + 1; % range of years
msYears = [1990, 2017, 2025, 2050, 2080, 2100, 2115]'; % milestone years
num_msYears = length(msYears);  % number of milestone years
num_yearT = length(yearT); % number of return period years

% quantile values
quantileValues = [0.025, 0.05, 0.16, 0.5, 0.84, 0.95, 0.975]';

% adjustment rates
rIA = 0.00145; % rate of isostatic adjustment (in meter/year)
rSC = 0; % rate of storm contribution (in meter/year)
rOB = 0.0039; % rate of observed SLR

%% SLR - simulated sea level rise for milestone years
slr_sim_msYears = f_SLR_Sim(nSim, msYears, rIA, rSC, rOB);

%% Water level calculations
%Data preparation
obs = data_temp(data_temp > threshold); % data above threshold
obsSorted = sort(obs, 'descend');  % sort, largest on top (descending order)
obsOverThresh = obsSorted - threshold; % minus threshold
nobs = length(obsOverThresh);

% Poisson intensity = number of obs / number of years in obs
lambdaPoisson = nobs / obsYearsRange;

% Plotting position, select option in plotting position function
% calculate pp for each number of observations, transpose results to column vector
pp = arrayfun(@(i) f_PP('haz',nobs,i),1:nobs)'; 
ariPP=1./(lambdaPoisson*pp); % element-wise division %%LIFA moved lambdaPoisson multiplication from pp to ariPP

% Maximum likelihood estimator
% estimate the Weibull distribution parameters for scale and shape
% Extract scale and shape variables from parameter estimate output
mleParam = mle(obsOverThresh, 'distribution', 'Weibull');
wblScale = mleParam(1);  
wblShape = mleParam(2);
% Compute the Weibull values for our return periods
qSim = 1-1./(lambdaPoisson*yearT);  % get quantiles
% get wbl values for calculated quantiles, and add threshold back in
wbl_yearT = wblinv(qSim, wblScale, wblShape) + threshold; 
% Compute the Weibull estimates for observations
qSimObs = 1-pp; % get quantiles for ariPP %%LIFA much clearer than (1-1./(lambdaPoisson*ariPP))
wbl_obs = wblinv(qSimObs, wblScale, wblShape) + threshold; % get wbl values, add threshold


%% Monte Carlo simulations for WL
% generate random values Poisson distr. with lambda parameter (avg) of nobs,
% to array of size nSim rows and 1 column. This generated avg number of events
randPoisson = poissrnd(nobs, nSim, 1);
% preallocate cell array (needed cell because different result lengths
% corresponding to randPoisson values)
wbl_MC = cell(nSim, 1);
% for each avg num of events, generate rand values from Weibull distribution
for i = 1:length(randPoisson);
     MC_temp = wblrnd(wblScale, wblShape, [1 randPoisson(i)])'; % sim results given avg num of events
     wbl_MC(i) = {MC_temp};  % store all results in cell array
end

% Determine variation in Weibull parameters based on MC simulations
re3 = zeros(nSim, 2); % preallocate to store wlb scale and shape param
wbl_est = zeros(nSim, num_yearT); % preallocate
for i = 1:nSim;
    % Extract Weibull scale and shape variables for each set of MC generated values
    re3(i,:) = mle(wbl_MC{i}, 'distribution', 'Weibull'); 
    % evaluate the Weibull values corresponding to quantiles, given MC based scale and shape parameters
    wbl_est(i,:) = wblinv(qSim, re3(i,1), re3(i,2)) + threshold;
end

% Get quantile values from the MC generated values at the desired quantile values
% each column of wbl_est is for one return period
quantiles_wbl = zeros(num_yearT, length(quantileValues)); % preallocate
for i = 1:num_yearT;
    quantiles_wbl(i,:) = quantile(wbl_est(:,i), quantileValues);
end

%% Combined SLR and WL
sim_wl_slr = zeros(nSim, num_yearT, num_msYears); % preallocate, 3D array
% convert sim slr from m to cm
slr_sim_msYears_cm = 100*slr_sim_msYears;
for i = 1:num_msYears;
    % use broadcast function to element-wise add simulated SLR to Weibull
    % estimated values. Loop through all milestone years
    sim_wl_slr(:,:,i) = bsxfun(@plus, wbl_est, slr_sim_msYears_cm(:,i));
end

% Get quantile values from the combined WL and SLR data
quantile_wl_slr = zeros(length(quantileValues), num_yearT, num_msYears); % preallocate
for i = 1:num_msYears;
    for j = 1:num_yearT;
        quantile_wl_slr(:,j,i) = quantile(sim_wl_slr(:,j,i), quantileValues);
    end
end

% %% Save outputs to file
% writeFolder = strcat(pwd,'\output\');  % save to subfolder called output
% % save simulated SLR
%     fWriteName = [writeFolder 'SLR.txt'];
%     save(fWriteName, 'slr_sim','-ascii');


%% Plots
% %% Return period vs water level
% hfig = figure(1);
%     % set figure appearances
%     set(gca, 'XScale', 'log',...
%              'XLim',[0 max(yearT)]);
%     % plot things
%     hold on
%     plot(ariPP, obsSorted, '+'); % plot observations
%     plot(yearT, quantile_wl_slr(4,:,1)); % plot median
%    
%     
%     hold off
    
%% QQ plots
% for msYear_idx = 1:num_msYears;
%     hfig = figure(msYear_idx);
%     title_generated = ['QQ-Plot Storms South Copenhagen 1954-2014 CRES yr ' num2str(msYears(msYear_idx))];
%     % set figure appearances
%         set(gca, 'XLim',[110 max(quantile_wl_slr(7,:,msYear_idx))],...
%                  'YLim',[110 max(quantile_wl_slr(7,:,msYear_idx))]);
%         box on
%         ylabel('Observed quantile [cm]');
%         xlabel('Theoretical quantile [cm]');
%         title(title_generated);
%         % plot things
%             hold on
%             plot(wbl_obs, obsSorted, '+'); % observations
%             plot([0,400],[0,400], 'k', 'HandleVisibility','off'); % plot diagonal
%             plot(quantile_wl_slr(4,:,msYear_idx), quantile_wl_slr(4,:,msYear_idx), 'k'); % median
%             % 68% confidence
%             plot(quantile_wl_slr(4,:,msYear_idx), quantile_wl_slr(3,:,msYear_idx), 'b');
%             plot(quantile_wl_slr(4,:,msYear_idx), quantile_wl_slr(5,:,msYear_idx), 'b', 'HandleVisibility','off');
%             % 90% confidence
%             plot(quantile_wl_slr(4,:,msYear_idx), quantile_wl_slr(2,:,msYear_idx), 'g');
%             plot(quantile_wl_slr(4,:,msYear_idx), quantile_wl_slr(6,:,msYear_idx), 'g', 'HandleVisibility','off');
%             % 95% confidence
%             plot(quantile_wl_slr(4,:,msYear_idx), quantile_wl_slr(1,:,msYear_idx), 'r');
%             plot(quantile_wl_slr(4,:,msYear_idx), quantile_wl_slr(7,:,msYear_idx), 'r', 'HandleVisibility','off');
%         legend('obs', 'median', '68% conf.', '90% conf.', '95% conf.');
%         hold off
%         
%         % save plots
%         writeFolder = strcat(pwd,'\output\');  % save to subfolder called output
%         fWriteName = [writeFolder 'QQ-Plot Storms South CPH 1954-2014 CRES yr ' num2str(msYears(msYear_idx))];
%         print(fWriteName,'-dpng')
% end

%% Print for Excel
printForExcelTrue = 1.0;
if printForExcelTrue == 1.0;
    for i = 1:num_msYears;
        msYears(i)
        printForExcel = [quantile_wl_slr(4,:,i)' ... % median
                         quantile_wl_slr(5,:,i)' quantile_wl_slr(3,:,i)' ... % 68% upper/lower
                         quantile_wl_slr(6,:,i)' quantile_wl_slr(2,:,i)']    % 90% upper/lower
    end;
end;























