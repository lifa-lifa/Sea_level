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
computed_WC1 = f_WakebyConst(L_1825_2014); % takes L-moments as input

% compute Wakeby Constants again
l2norm1825_2014 = t_1825_2014;
l3norm1825_2014 = L_1825_2014(3)/L_1825_2014(1);
computed_WC2 = f_WakebyConst([computed_WC1 l2norm1825_2014 l3norm1825_2014]);

% solve polynomial of D1z^2 + D2z + D3 == 0 for z
z = roots(computed_WC1);
l1norm = 1.0;

% beta and negative delta are the roots of the quadratic eq. 
% beta being the larger of the two roots
betaE = max(z); 
deltaE = -min(z);
if (betaE + deltaE) < 0; deltaE = 0.00001-betaE; else deltaE = deltaE; end
if deltaE < 0; deltaE = 0.00001; else deltaE = deltaE; end
% the other parameter, alpha, gamma, xi
alphaE = (1+betaE)*(2+betaE)*(3+betaE)*((1+deltaE)*l2norm1825_2014...
    -(3-deltaE)*l3norm1825_2014)/(4*(betaE+deltaE));
gammaE = -(1-deltaE)*(2-deltaE)*(3-deltaE)*((1-betaE)*l2norm1825_2014...
    -(3+betaE)*l3norm1825_2014)/(4*(betaE+deltaE));
if gammaE < 0; gammaE = 0.00001; else gammaE = gammaE; end
if gammaE < 0.000011; deltaE = 0.00001; else deltaE = deltaE; end
xiE = l1norm - alphaE/(1+betaE) - gammaE/(1-deltaE);






