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
obsP1 = sort(cell2mat(data_temp(3)), 'descend');  % P1 is 1825-2014
obsP2 = sort(cell2mat(data_temp(2)), 'descend');  % P2 is 1500-1824
obsP3 = sort(cell2mat(data_temp(1)), 'descend');  % P3 is 1044-1499


%% Inputs
rng default;  % sets random seed. If uncommented, randoms are repeatable
thresholdP3 = 270;
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

% quantiles (Probability of non-exceedance, F = 1-1/T)
q = 1-1./(lambdaPoisson * yearT);

%% Calculate L-moments
% results from lmom function stored in array
% output named using capital L, because confusion between l and 1 
L_P1 = f_lmom(xP1, 5);  
% calculate L-CV = t = l_2/l_1
t_P1 = L_P1(2)/L_P1(1);
% calculate L-moment ratios (t_r = l_r / l_2)
t3_P1 = L_P1(3)/L_P1(2);
t4_P1 = L_P1(4)/L_P1(2);
t5_P1 = L_P1(5)/L_P1(2);

%% Wakeby from P1
% compute Wakeby Constants
% takes L-moments as input, outputs constats as [D1 D2 D3]
wc_P1 = f_WakebyConst(L_P1); 

l2normP1 = t_P1;
l3normP1 = L_P1(3)/L_P1(1);

% compute Wakeby parameters
% output contains 5 parameters for Wakeby distribution with shape 
% parameters B and D, scale parameters A and G, and location parameter M (or X).
% output is in the order of [A, B, G, D, X]
[A_P1, B_P1, G_P1, D_P1, X_P1] = f_Wakeby(wc_P1(1), wc_P1(2), wc_P1(3),...
    l2normP1, l3normP1);

%% Generate values for P2
% find num years that have no data, and therefore need generated values
nYrsGenP2 = (1824-1500+1-length(obsP2))*repetitionFactor;
% Generate values for P2 from Wakeby, multiply by mean (l1) of P3
wakebyEstP2 = L_P1(1)*f_wkbrnd(A_P1, B_P1, G_P1, D_P1, X_P1, [nYrsGenP2,1]);
wakebyEstP2 = sort(wakebyEstP2, 'descend');    

    % print out max and min to check
    str = sprintf('P2: Before replacement max and min: %.1f and %.1f', [max(wakebyEstP2), min(wakebyEstP2)]);
    disp(str);

% num of random values that exceed threshold and need to be re-generated
nReplaceP2 = sum(wakebyEstP2>thresholdP2);
nTemp = 1;  % temp variable for while loop
while nTemp > 0; % do until we have a list with only values below threshold
    listReplaceP2 = L_P1(1)*f_wkbrnd(A_P1, B_P1, G_P1, D_P1, X_P1, [nReplaceP2,1]);
    nTemp = sum(listReplaceP2>thresholdP2);
end

% replace all values above threshold, by the re-generated list
wakebyEstP2(wakebyEstP2>thresholdP2) = listReplaceP2;
wakebyEstP2 = sort(wakebyEstP2, 'descend');  % sort data again

    % print out max and min to check
    str = sprintf('P2: After  replacement max and min: %.1f and %.1f', [max(wakebyEstP2), min(wakebyEstP2)]);
    disp(str);
    
%% Combined P2, P1 data 
% join together observed data and generated data
xP12 = sort(cat(1, xP2, xP1, wakebyEstP2), 'descend');  % concatenate along vertical dim
% calculate L-moments and ratios
L_P12 = f_lmom(xP12, 5);  
% calculate L-CV = t = l_2/l_1
t_P12 = L_P12(2)/L_P12(1);
% calculate L-moment ratios (t_r = l_r / l_2)
t3_P12 = L_P12(3)/L_P12(2);
t4_P12 = L_P12(4)/L_P12(2);
t5_P12 = L_P12(5)/L_P12(2);

% compute Wakeby Constants
wc_P12 = f_WakebyConst(L_P12); 
l2normP12 = t_P12;
l3normP12 = L_P12(3)/L_P12(1); 

% compute Wakeby parameters [A, B, G, D, X]
[A_P12, B_P12, G_P12, D_P12, X_P12] = f_Wakeby(wc_P12(1), wc_P12(2), wc_P12(3),...
    l2normP12, l3normP12);

%% Generate values for P3
% find num years that have no data, and therefore need generated values
nYrsGenP3 = (1499-1044+1-length(obsP3))*repetitionFactor;
% Generate values for P1 from Wakeby, multiply by mean (l1) of combined P23
wakebyEstP3 = L_P12(1)*f_wkbrnd(A_P12, B_P12, G_P12, D_P12, X_P12, [nYrsGenP3,1]);
wakebyEstP3 = sort(wakebyEstP3, 'descend'); 

    % print out max and min to check
    str = sprintf('P3: Before replacement max and min: %.1f and %.1f', [max(wakebyEstP3), min(wakebyEstP3)]);
    disp(str);

% num of random values that exceed threshold and need to be re-generated
nReplaceP3 = sum(wakebyEstP3>thresholdP3);
nTemp = 1;  % temp variable for while loop
while nTemp > 0; % do until we have a list with only values below threshold
    listReplaceP3 = L_P12(1)*f_wkbrnd(A_P12, B_P12, G_P12, D_P12, X_P12, [nReplaceP3,1]);
    nTemp = sum(listReplaceP3>thresholdP3);
end    

% replace all values above threshold, by the re-generated list
wakebyEstP3(wakebyEstP3>thresholdP3) = listReplaceP3;
wakebyEstP3 = sort(wakebyEstP3, 'descend');  % sort data again

    % print out max and min to check
    str = sprintf('P3: After  replacement max and min: %.1f and %.1f', [max(wakebyEstP3), min(wakebyEstP3)]);
    disp(str);

%% Combined P1, P2, P3
% join together observed data and generated data
xP123 = sort(cat(1, xP1, xP2, xP3, wakebyEstP3, wakebyEstP2), 'descend');  % concatenate along vertical dim
% calculate L-moments and ratios
L_P123 = f_lmom(xP123, 5);  
% calculate L-CV = t = l_2/l_1
t_P123 = L_P123(2)/L_P123(1);
% calculate L-moment ratios (t_r = l_r / l_2)
t3_P123 = L_P123(3)/L_P123(2);
t4_P123 = L_P123(4)/L_P123(2);
t5_P123 = L_P123(5)/L_P123(2);

% compute Wakeby Constants
wc_P123 = f_WakebyConst(L_P123); 
l2normP123 = t_P123;
l3normP123 = L_P123(3)/L_P123(1); 

% compute Wakeby parameters [A, B, G, D, X]
[A_P123, B_P123, G_P123, D_P123, X_P123] = f_Wakeby(wc_P123(1), wc_P123(2), wc_P123(3),...
    l2normP123, l3normP123);

%% Wakeby distribution
% plot check
% x_plot = (0:0.01:4);
% plot(x_plot, f_wkbpdf(x_plot , A_P123, B_P123, G_P123, D_P123, X_P123));

% quantiles
quantEst = f_wkbinv(q, A_P1, B_P1, G_P1, D_P1, X_P1);
quantEst2 = f_wkbinv(q, A_P12, B_P12, G_P12, D_P12, X_P12);
quantEst3 = f_wkbinv(q, A_P123, B_P123, G_P123, D_P123, X_P123);

QE1 = L_P1(1)*quantEst;
QE2 = L_P12(1)*quantEst2;
QE3 = L_P123(1)*quantEst3;

%% Print results for comparison

L_P1
t_val_P1 = [t_P1 t3_P1 t4_P1 t5_P1]
L_P12
t_val_P12 = [t_P12 t3_P12 t4_P12 t5_P12]
L_P123
t_val_P123 = [t_P123 t3_P123 t4_P123 t5_P123]


