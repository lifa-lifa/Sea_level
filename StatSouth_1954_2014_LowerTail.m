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
threshold = 111;
alpha = 0.10;  % significance level
nSim = 1000;  % number of simulations, typ. 10000

% return periods and years
yearT = [10, 20, 50, 100, 250, 500, ...
         1000, 2000, 3000, 5000, 10000, 100000, 1000000]';  % return periods
yearT_plot = [50, 100, 250, 1000, 10000]';  % return periods for plotting
obsYearsRange = 2016 - 1956 + 1; % range of years
msYears = [1990, 2017, 2050, 2080, 2100]'; % milestone years
num_msYears = length(msYears);  % number of milestone years
num_yearT = length(yearT); % number of return period years

% quantile values
quantileValues = [0.025, 0.05, 0.16, 0.5, 0.84, 0.95, 0.975]';

% adjustment rates
rIA = 0.0011; % rate of isostatic adjustment (in meter/year)
rSC = 0.0012; % rate of storm contribution (in meter/year)

% SLR values from CRES 
%(for reference see Fig 3.7 of Storebælt Østtunnel, Klimavurdering og Sikring, Ramperne 2015)
slr_values=[0:0.125:3]';  % stepped SLR 0 to 3 in steps of 0.125m
% Probabilities for this scenario
slr_prob =[0, 0.007488233, 0.025673941, 0.067394095, 0.132648695,...
    0.177578092, 0.174368849, 0.133718442, 0.084510056, 0.054557125,...
    0.036371416, 0.023534446, 0.017115961, 0.013906718, 0.011232349,...
    0.009092854, 0.008023107, 0.006953359, 0.005348738, 0.003637142,...
    0.002674369, 0.00203252, 0.001390672, 0.000641849, 0.000106975]';

%% Sea level rise calculations
slr_prob_cum = cumsum(slr_prob);
% generate random values between 0 and 1 from uniform distribution, to array of size nSim rows and 1 column
rand0_1 = rand(nSim,1);
% anonymous function to find the index of the cum prob that just exceeds (x). 
% Example, if cum_prob =[0, 0.19, 0.22, 0.3...], then fHandle(0.2)= 3.  'first' is stated explicitly
fHandle = @(x) find(slr_prob_cum > x, 1, 'first');
% find indices of all simulations
indices = arrayfun(fHandle, rand0_1);
% get the slr values corresponding to the indices
slr_sim = slr_values(indices);

% evaluate climate factor for milestone years
climateFactor = f_SLR_Norm(msYears, rIA, rSC);
% simulated slr including climate change for milestone years
% array columns are msYears, rows are each slr_sim. 
% transposed climateFactor for multiplication to match dimensions
slr_sim_msYears = slr_sim * climateFactor';

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


%% Monte Carlo simulations 
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
wbl_quantiles = zeros(num_yearT, length(quantileValues)); % preallocate
for i = 1:num_yearT;
    wbl_quantiles(i,:) = quantile(wbl_est(:,i), quantileValues);
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



















