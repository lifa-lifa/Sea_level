%% Oeresund A102102
% * LIFA 2017-09-06
% * Statistics for southern storms, lower tail
%  

close all
clear all
clc

%% Get data
ReadFolder = strcat(pwd,'\input\');
dirList = dir(strcat(ReadFolder,'*.txt'));  %list all txt files
fReadName = dirList(1).name;
full_name = [ReadFolder fReadName];
fid = fopen(full_name, 'r');
data_temp = textscan(fid, '%f', 'Delimiter', ',');
fclose(fid);

data_temp = cell2mat(data_temp);


%% Inputs
threshold = 111;
alpha = 0.10;  % significance level
NSim = 1000;  % number of simulations, typ. 10000

% return periods and years
TYear = [10, 20, 50, 100, 250, 500, ...
         1000, 2000, 3000, 5000, 10000, 100000, 1000000];  % return periods
TYearPlot = [50, 100, 250, 1000, 10000];  % return periods for plotting
TRange = 2016 - 1956 + 1; % range of years
MileStoneYears = [1990, 2017, 2050, 2080, 2100];

% adjustment rates
IsoAdjust = 0.0011; % rate of isostatic adjustment (in meter/year)
StormContrib = 0.0012; % rate of storm contribution (in meter/year)

% SLR values from CRES
seaLevelRiseValues=[0:0.125:3];  % stepped SLR 0 to 3 in steps of 0.125m
% Probabilities for this scenario (LIFA -> NEJO  Where is this from?)
probabilities=[0, 0.007488233, 0.025673941, 0.067394095, 0.132648695,...
    0.177578092, 0.174368849, 0.133718442, 0.084510056, 0.054557125,...
    0.036371416, 0.023534446, 0.017115961, 0.013906718, 0.011232349,...
    0.009092854, 0.008023107, 0.006953359, 0.005348738, 0.003637142,...
    0.002674369, 0.00203252, 0.001390672, 0.000641849, 0.000106975];

%% Calcs
obs = data_temp(data_temp > threshold); % data above threshold
obsSorted = sort(obs, 'descend');  % sort, largest on top (descending order)
obsOverThresh = obsSorted - threshold; % minus threshold
xLength = length(obsOverThresh);












