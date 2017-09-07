%% Oeresund A102102
% * LIFA 2017-09-06
% * Test to replicate Northern storm matlab file, to check that Southern
% storm is replicated correctly in matlab
close all
clear all
clc

%% Get data
ReadFolder = strcat(pwd,'\input\');  % data folder
dirList = dir(strcat(ReadFolder,'*.txt'));  %list all txt files
fReadName = dirList(2).name; 
full_name = [ReadFolder fReadName];
fid = fopen(full_name, 'r');
data_temp = textscan(fid, '%f', 'Delimiter', ',');
fclose(fid);
data_temp = cell2mat(data_temp);  % convert cell array to matrix


%% Inputs
threshold = 119;
alpha = 0.10;  % significance level
nSim = 1000;  % number of simulations, typ. 10000

% return periods and years
yearT = [2, 5, 10, 20, 50, 100, 250, 500, 1000, 10000]';  % return periods
yearT_plot = [50, 100, 250, 1000, 10000]';  % return periods for plotting
obsYearsRange = 185; % range of years
msYears = [1990, 2100]'; % milestone years
num_msYears = length(msYears);  % number of milestone years
num_yearT = length(yearT); % number of return period years

% adjustment rates
rIA = 0.0011; % rate of isostatic adjustment (in meter/year)
rSC = 0.0012; % rate of storm contribution (in meter/year)

% SLR values from CRES 
%(for reference see Fig 3.7 of Storeb�lt �sttunnel, Klimavurdering og Sikring, Ramperne 2015)
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

% Plotting position, select Weibull option in plotting position function
% calculate pp for each number of observations, transpose results to column vector
pp = (lambdaPoisson)*arrayfun(@(i) f_PP('wbl',nobs,i),1:nobs)'; 
ariPP=1./pp; % element-wise division

% Maximum likelihood estimator
% estimate the Weibull distribution parameters for scale and shape
% Extract scale and shape variables from parameter estimate output
mleParam = mle(obsOverThresh, 'distribution', 'Weibull');
wblScale = mleParam(1);  
wblShape = mleParam(2);
% Compute the Weibull quantiles for our return periods
qSim = 1-1./(lambdaPoisson*yearT);  % get quantiles
wbl_yearT = wblinv(qSim, wblScale, wblShape); % get return periods for given quantiles
whatsit = wbl_yearT + threshold;



%% Combined SLR and WL















