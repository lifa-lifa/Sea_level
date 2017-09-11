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
obs1044_1499 = sort(cell2mat(data_temp(1)), 'descend');
obs1500_1824 = sort(cell2mat(data_temp(2)), 'descend');
obs1825_2014 = sort(cell2mat(data_temp(3)), 'descend');

%% Inputs
threshold1044_1499 = 270;
threshold1500_1824 = 240;
repetitionFactor = 100;
nSim = 10000;  % number of simulations, typ. 10000
lambdaPoisson = 1;

% return periods and years
yearT = [10, 20, 50, 100, 250, 500, ...
         1000, 2000, 3000, 5000, 10000, 100000]';  % return periods

%% Prepare data
% repeat data by repetition factor
x1044_1499 = repelem(obs1044_1499, repetitionFactor);
x1500_1824 = repelem(obs1500_1824, repetitionFactor);
x1825_2014 = repelem(obs1825_2014, repetitionFactor);

% quantiles
q = 1-1./(lambdaPoisson * yearT);

%% calculate L-moments
% results from lmom function stored in array
L_1825_2014 = f_lmom(x1825_2014, 5);  
% calculate L-CV = t = l_2/l_1
t_1825_2014 = L_1825_2014(2)/L_1825_2014(1);
% calculate L-moment ratios (t_r = l_r / l_2)
t3_1825_2014 = L_1825_2014(3)/L_1825_2014(2);
t4_1825_2014 = L_1825_2014(4)/L_1825_2014(2);
t5_1825_2014 = L_1825_2014(5)/L_1825_2014(2);

%% Wakeby
% compute Wakeby Constants
% takes L-moments as input, outputs constats as [D1 D2 D3]
computed_WC = f_WakebyConst(L_1825_2014); 

l2norm1825_2014 = t_1825_2014;
l3norm1825_2014 = L_1825_2014(3)/L_1825_2014(1);

% compute Wakeby parameters
[alphaW, betaW, gammaW, deltaW, xiW] = f_Wakeby(computed_WC(1), computed_WC(2), computed_WC(3),...
    l2norm1825_2014, l3norm1825_2014);


    





