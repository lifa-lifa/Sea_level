%% Oeresund A102102
% * LIFA 2017-09-11
% * Statistics for southern storms, middle segment
% * Based on NEJO's Wakeby-Simulation-1044-2014-ThreePeriodModel - Køge-Gedser-Travemunde-RepetitionFactor-Mixture.nb
close all
clear all
clc

%% Get data
ReadFolder = strcat(pwd,'\input\forWakeby\');  % data folder
dirList = dir(strcat(ReadFolder,'*.txt'));  %list all txt files
% loop through all files, save to cell array
for i = 1:length(dirList);
    fReadName = dirList(i).name; 
    full_name = [ReadFolder fReadName];
    fid = fopen(full_name, 'r');
    data_temp(i) = textscan(fid, '%f', 'Delimiter', ',');
    fclose(fid);
%     data_temp(i) = cell2mat(data_temp);  % convert cell array to matrix
end

% move cell array data to matrix, and sort
obsP1 = sort(cell2mat(data_temp(1)), 'descend');
obsP2 = sort(cell2mat(data_temp(2)), 'descend');
obsP3 = sort(cell2mat(data_temp(3)), 'descend');

%% Inputs
thresholdP1 = 270;
thresholdP2 = 240;
repetitionFactor = 100;
nSim = 10000;  % number of simulations, typ. 10000
lambdaPoisson = 1;

% return periods and years
yearT = [10, 20, 50, 100, 250, 500, ...
         1000, 2000, 3000, 5000, 10000, 100000]';  % return periods

%% Prepare data
% repeat data by repetition factor
xP1 = repelem(obsP1, repetitionFactor);
xP2 = repelem(obsP2, repetitionFactor);
xP3 = repelem(obsP3, repetitionFactor);

% quantiles
q = 1-1./(lambdaPoisson * yearT);

%% calculate L-moments
% results from lmom function stored in array
L_P3 = f_lmom(xP3, 5);  
% calculate L-CV = t = l_2/l_1
t_P3 = L_P3(2)/L_P3(1);
% calculate L-moment ratios (t_r = l_r / l_2)
t3_P3 = L_P3(3)/L_P3(2);
t4_P3 = L_P3(4)/L_P3(2);
t5_P3 = L_P3(5)/L_P3(2);

%% Wakeby
% compute Wakeby Constants
% takes L-moments as input, outputs constats as [D1 D2 D3]
computed_WC = f_WakebyConst(L_P3); 

l2normP3 = t_P3;
l3normP3 = L_P3(3)/L_P3(1);

% compute Wakeby parameters
[alphaW, betaW, gammaW, deltaW, xiW] = f_Wakeby(computed_WC(1), computed_WC(2), computed_WC(3),...
    l2normP3, l3normP3);

%% Generate values
nYrGen = 1824-1500+1-length(obsP2);
    





