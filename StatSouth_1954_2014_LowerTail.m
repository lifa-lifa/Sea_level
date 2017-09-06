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
TYear = [10, 20, 50, 100, 250, 500, ...
         1000, 2000, 3000, 5000, 10000, 100000, 1000000];  % return periods
TYearPlot = [50, 100, 250, 1000, 10000];  % return periods for plotting
TRange = 2016 - 1956 + 1; % range of years
MileStoneYears = [1990, 2017, 2050, 2080, 2100];

%% Calcs
obs = data_temp(data_temp > threshold); % data above threshold
obsSorted = sort(obs, 'descend');  % sort, largest on top (descending order)
obsOverThresh = obsSorted - threshold; % minus threshold
xLength = length(obsOverThresh);












