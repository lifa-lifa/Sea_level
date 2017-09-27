%% Oeresund A102102
% * LIFA 2017-09-22
% * Statistics for southern storms, upper segment (Exponential)
% * Based on NEJO's Statistics South-1044-2013-UpperTail.nb
% * refer to Sec 6.3 of CR-CSJV-GEN=GSh-GC-DES-REP-212004  
close all
clear all
clc

%% Inputs
% rng default;  % sets random seed. Matlab default is Mersenne Twister seed 5498
nSim = 10;
% return periods and years
msYears = [1990, 2025, 2050, 2080, 2100]'; % milestone years
yearT = [250, 500, 1000, 2000, 5000, 10000, 100000]';  % return periods
% truncation
gammaWE = 240;
% take the larger of gammaWE or historial truncation level
gammaP1 = max(gammaWE, 200); % assumed historical truncation for P1
gammaP2 = max(gammaWE, 240); % ditto P2
gammaP3 = max(gammaWE, 270); % ditto P3
sigmaP1 = 15; % std deviation in P1, 15 cm
sigmaP2 = 45;
sigmaP3 = 60;
% duration of data periods
EndP1 = 2013; % NEJOs mathematica code inconsistently uses 2013 and 2015?
EndP2 = 1824;
EndP3 = 1499;
StartP3 = 1044;
% TE = EndP1-StartP3+1; % duration of all periods
TE = 972;
TP1 = EndP1 - EndP2;
TP2 = EndP2 - EndP3;
TP3 = EndP3 - StartP3 + 1;


%% Data
sampleObs = [243.513, 257.182, 290.481, 250.69, 283.059, 311.376, 331.04, 304.774, 358.54];
sampleObsP12 = [243.513, 257.182, 290.481, 250.69, 283.059, 311.376, 331.04];
truncation = [0, 0, 0, 0, 0, 0, 0, 30, 30];


%% Distribution estimate
% Exponential (Matlab returns mu, Mathematica returns lambda = 1/mu)
estExpMle = mle(sampleObs-gammaWE, 'distribution', 'Exponential')
estExpMle_lambda = 1/estExpMle;  % inverse, for comparison with Mathematica results


%% Conditional MLE estimate    
% prepare the input data in cell array, necessary at MLE only takes one input arg. 
% i.e. cannot input samples and truncation as two arguments
data = cell(2,1);
data{1} = (sampleObs-gammaWE);
data{2} = truncation;

% custom negative log likelihood function, we want to miminize this function
% recall mu = 1/lambda. Matlab uses mu to define exponential distribution
custnloglf = @(mu, data, cens, freq) -nansum(log((pdf('Exponential', data{1}, mu))./...
                                                     (1-cdf('Exponential', data{2}, mu))));
estExpMleTrunc = mle(data, 'nloglf', custnloglf, 'start', estExpMle);
1/estExpMleTrunc % Mathematica output check
LLExpMle = custnloglf(estExpMleTrunc, data)

NPoissonSimulation = round(TE/(TP1+TP2)*length(sampleObsP12))
lambdaPoisson = NPoissonSimulation/TE;
lambdaPoissonMin = ceil(TE/min(yearT))/TE
qSim = 1-1./(lambdaPoisson*yearT);  % get quantiles

% random number from Poisson distribution
k = poissrnd(NPoissonSimulation, 1)
% cannot have too low number
kmin = ceil(TE/min(yearT));
kval = max(k, kmin);
kval = 13; % test!!!!!

% new sample from Exponential
lambdaPoissonSimNewSample = kval/TE
qSimNewSample = 1-1./(lambdaPoissonSimNewSample*yearT);  % get quantiles
% generate k numbers from the estimated Exponential distribution
newSampleEx = exprnd(estExpMleTrunc, kval, 1);
%!!!!! for test, overwrite random generated new sample
newSampleEx = [46.6563, 76.9496, 9.64744, 44.7102, 62.7202, 125.892, 23.0628, 60.1462, 27.1351, 4.99122,...
41.8752, 90.3281, 27.3091];

% determine Exponential - MLE
estExpMleNewSample = mle(newSampleEx, 'distribution', 'Exponential')
estExpMleNewSample_lambda = 1/estExpMleNewSample  % inverse, for comparison with Mathematica results

% values at quantiles
qEstExpMle = expinv(qSimNewSample, estExpMleNewSample) + gammaWE

%% KS test
ExpCFD = makedist('Exponential', estExpMleTrunc);
[h,p] = kstest(data{1}, ExpCFD)






