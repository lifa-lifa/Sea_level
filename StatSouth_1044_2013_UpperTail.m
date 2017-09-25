%% Oeresund A102102
% * LIFA 2017-09-22
% * Statistics for southern storms, upper segment (Exponential)
% * Based on NEJO's Statistics South-1044-2013-UpperTail.nb
% * refer to Sec 6.3 of CR-CSJV-GEN=GSh-GC-DES-REP-212004  
close all
clear all
clc

%% Inputs
rng default;  % sets random seed. Matlab default is Mersenne Twister seed 5498
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
TE = EndP1-StartP3+1; % duration of all periods
% TE = 972;
TP1 = EndP1 - EndP2;
TP2 = EndP2 - EndP3;
TP3 = EndP3 - StartP3 + 1;


%% Data
% Data preparation
obsP1 = sort([193, 216, 220, 230, 235, 286], 'descend'); % 1825-2013? or 2015?
obsP2 = sort([258, 290, 301, 366], 'descend'); % 1500-1824
obsP3 = sort([276, 286, 343], 'descend'); % 1044-1500
% truncate data
obsP1Trunc = obsP1(obsP1 > gammaP1);
obsP2Trunc = obsP2(obsP2 > gammaP2);
obsP3Trunc = obsP3(obsP3 > gammaP3);
% join
obsAll = sort([obsP1 obsP2 obsP3], 'descend');
obsAllTrunc = sort([obsP1Trunc obsP2Trunc obsP3Trunc], 'descend');
% Poisson intensity, lambda Poisson
nobs = length(obsAllTrunc); % num of points above gammaWE per 1000 yrs
lambdaPoisson = nobs/TE;

% Compute the Weibull values for our return periods
qSim = 1-1./(lambdaPoisson*yearT);  % get quantiles

%% Simulate uncertainty
% generate rand values between 0 and 1 from uniform distribution, to array of size nSim rows and 1 column
rand0_1 = normrnd(10, 2, nSim, 1);



%% Distribution estimate
% Weibull
wbMleParam = mle(obsP3Trunc, 'distribution', 'Weibull');
% Exponential
exMleParam = mle(obsP3Trunc, 'distribution', 'Exponential');





%% direct comparison to NEJO mathematica MLS KS test
% obsminusThresh = [258, 286.01, 290, 301, 366, 246, 256, 313];
% estExMLE = mle((obsminusThresh-gammaWE), 'distribution', 'Exponential');




