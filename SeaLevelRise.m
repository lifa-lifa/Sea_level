%% Oeresund A102102
% * LIFA 2017-09-19
% * script to generate Sea Level Rise
close all
clc

%% Inputs
nSim = 10000;  % number of simulations, typ. 10000

% return periods and years
msYears = [1990, 2017, 2050, 2080, 2100]'; % milestone years
num_msYears = length(msYears);  % number of milestone years

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



% %% Save outputs to file
% writeFolder = strcat(pwd,'\output\');  % save to subfolder called output
% % save simulated SLR
%     fWriteName = [writeFolder 'SLR.txt'];
%     save(fWriteName, 'slr_sim','-ascii');

























