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
obsP1 = sort(cell2mat(data_temp(1)), 'descend');  % P1 is 1044-1499
obsP2 = sort(cell2mat(data_temp(2)), 'descend');  % P2 is 1500-1824
obsP3 = sort(cell2mat(data_temp(3)), 'descend');  % P3 is 1825-2014

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
% output named using capital L, because confusion between l and 1 
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
% output contains 5 parameters for Wakeby distribution with shape 
% parameters B and D, scale parameters A and G, and location parameter M (or X).
% output is in the order of [A, B, G, D, X]
[alphaW, betaW, gammaW, deltaW, xiW] = f_Wakeby(computed_WC(1), computed_WC(2), computed_WC(3),...
    l2normP3, l3normP3);

%% Generate values
% find num years that have no data, and therefore need generated values
nYrsGenP2 = (1824-1500+1-length(obsP2))*repetitionFactor;
% Generate values for P2 from Wakeby, multiply by mean (l1) of P3
wakebyEstP2 = L_P3(1)*f_wkbrnd(alphaW, betaW, gammaW, deltaW, xiW, [nYrsGenP2,1]);
wakebyEstP2 = sort(wakebyEstP2, 'descend');    

    % print out max and min to check
    str = sprintf('Before replacement max and min: %.1f and %.1f', [max(wakebyEstP2), min(wakebyEstP2)]);
    disp(str);

% num of random values that exceed threshold and need to be re-generated
nRandReplace = sum(wakebyEstP2>thresholdP2);
nTemp = 1;  % temp variable for while loop
while nTemp > 0; % do until we have a list with only values below threshold
    listToReplace = L_P3(1)*f_wkbrnd(alphaW, betaW, gammaW, deltaW, xiW, [nRandReplace,1]);
    listToReplace = sort(listToReplace, 'descend');
    nTemp = sum(listToReplace>thresholdP2);
end

% replace all values above threshold, by the re-generated list
wakebyEstP2(wakebyEstP2>thresholdP2) = listToReplace;
wakebyEstP2 = sort(wakebyEstP2, 'descend');  % sort data again

    % print out max and min to check
    str = sprintf('After  replacement max and min: %.1f and %.1f', [max(wakebyEstP2), min(wakebyEstP2)]);
    disp(str);
    
 
    
    

