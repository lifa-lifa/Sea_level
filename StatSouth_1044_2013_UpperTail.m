%% Oeresund A102102
% * LIFA 2017-09-22
% * Statistics for southern storms, upper segment (Exponential)
% * Based on NEJO's Statistics South-1044-2013-UpperTail.nb
% * refer to Sec 6.3 of CR-CSJV-GEN=GSh-GC-DES-REP-212004  
close all
clear all
clc
tic % start timer
%% Options
test = 0; % if set to 1, random generated values are replaced by test values
rng default;  % sets random seed. Matlab default is Mersenne Twister seed 5498

%% Inputs
nSim = 10000;
% return periods and years
msYears = [1990, 2017, 2025, 2050, 2080, 2100, 2115]'; % milestone years
yearT = [500, 1000, 2000, 5000, 10000]';  % return periods
num_yearT = length(yearT);
num_msYears = length(msYears);
% truncation
gammaWE = 240; % take the larger of gammaWE or historial truncation level
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
if test ==1; TE = 972; end;
TP1 = EndP1 - EndP2;
TP2 = EndP2 - EndP3;
TP3 = EndP3 - StartP3 + 1;

% adjustment rates
rIA = 0.00145; % rate of isostatic adjustment (in meter/year)
rSC = 0; % rate of storm contribution (in meter/year)
rOB = 0.0039; % rate of observed SLR

% quantile values
quantileValues = [0.025, 0.05, 0.16, 0.5, 0.84, 0.95, 0.975]';

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

%% preallocate arrays for Monte Carlo sim
sampleObs = cell(nSim, 1);
sampleObsP12 = cell(nSim, 1);
truncation = cell(nSim,1);
qEstExpMle = zeros(nSim, length(yearT));
LLExpMle = zeros(nSim, 1);
ksExpMLE = zeros(nSim, 1);
lambda = cell(nSim, 1);

%% Simulate uncertainty
% generate values from normal distribution, truncate using gamma
for n = 1:nSim
    for i= 1:3
       % eval using uncertainty simulation, and do for each period using that
       % period's observations, sigma, and gamma. Generates variables uncObsPx, x=1,2,3
       % e.g.  uncObsP1 = f_UncertaintySim(obsP1, sigmaP1, gammaP1);
       eval(['uncObsP' num2str(i) '=f_UncertaintySim(obsP' num2str(i)...
                                                  ', sigmaP' num2str(i)...
                                                  ', gammaP' num2str(i) ');']); 
    end
    % join outputs
    sampleObs{n} = [uncObsP1 uncObsP2 uncObsP3];
    sampleObsP12{n} = [uncObsP1 uncObsP2];
    if test == 1; sampleObs{1} = [243.513, 257.182, 290.481, 250.69, 283.059, 311.376, 331.04, 304.774, 358.54]; end;

    % create truncation array, of same length as sampleObs
    for i= 1:3
       eval(['TruncP' num2str(i) '= gammaP' num2str(i) '- gammaWE;']); 
       eval(['TruncP' num2str(i) '= repmat(TruncP' num2str(i) ',1, length(uncObsP' num2str(i) '));'])
    end
    truncation{n} = [TruncP1 TruncP2 TruncP3];
    if test == 1; truncation{1} = [0, 0, 0, 0, 0, 0, 0, 30, 30]; end;


    %% Distribution estimate
    % Weibull (parambers in order of scale (beta), shape (alpha))
    wbMleParam = mle((sampleObs{n}-gammaWE), 'distribution', 'Weibull');
    % Exponential (Matlab returns mu, Mathematica returns lambda = 1/mu)
    estExpMle = mle(sampleObs{n}-gammaWE, 'distribution', 'Exponential');
    estExpMle_lambda = 1/estExpMle;  % inverse, for comparison with Mathematica results

    %% Conditional MLE estimate    
    % prepare the input data in cell array, necessary at MLE only takes one input arg. 
    % i.e. cannot input samples and truncation as two arguments
    data = cell(2,1);
    data{1} = (sampleObs{n}-gammaWE);
    data{2} = truncation{n};

    % custom negative log likelihood function, we want to miminize this function
    % recall mu = 1/lambda. Matlab uses mu to define exponential distribution
    custnloglf = @(mu, data, cens, freq) -nansum(log((pdf('Exponential', data{1}, mu))./...
                                                         (1-cdf('Exponential', data{2}, mu))));
    estExpMleTrunc = mle(data, 'nloglf', custnloglf, 'start', estExpMle);
    estExpMleTrunc_lambda = 1/estExpMleTrunc; % Mathematica output check
    LLExpMle(n) = custnloglf(estExpMleTrunc, data);

    NPoissonSimulation = round(TE/(TP1+TP2)*length(sampleObsP12{n}));
    lambdaPoisson = NPoissonSimulation/TE;
    lambdaPoissonMin = ceil(TE/min(yearT))/TE;
    qSim = 1-1./(lambdaPoisson*yearT);  % get quantiles

    % random number from Poisson distribution
    k = poissrnd(NPoissonSimulation, 1);
    % cannot have too low number
    kmin = ceil(TE/min(yearT));
    kval = max(k, kmin);
    if test == 1; kval = 13; end;

    % new sample from Exponential
    lambdaPoissonSimNewSample = kval/TE;
    qSimNewSample = 1-1./(lambdaPoissonSimNewSample*yearT);  % get quantiles
    % generate k numbers from the estimated Exponential distribution
    newSampleEx = exprnd(estExpMleTrunc, kval, 1);
    %!!!!! for test, overwrite random generated new sample
    if test == 1; newSampleEx = [46.6563, 76.9496, 9.64744, 44.7102, 62.7202, 125.892, 23.0628, 60.1462, 27.1351, 4.99122,...
    41.8752, 90.3281, 27.3091]; end;

    % determine Exponential - MLE
    estExpMleNewSample = mle(newSampleEx, 'distribution', 'Exponential');
    estExpMleNewSample_lambda = 1/estExpMleNewSample;  % inverse, for comparison with Mathematica results

    % values at quantiles
    qEstExpMle(n,:) = expinv(qSimNewSample, estExpMleNewSample) + gammaWE;

    %% KS test
    ExpCFD = makedist('Exponential', estExpMleTrunc);
    [~,p] = kstest(data{1}, ExpCFD); % only save p value
    

    %% print output option
    print = 1.0;
    if print == 1
        n
        sampleObs{n}
%         truncation{n}
%         estExpMleTrunc_lambda
%         kmin
%         kval
%         estExpMleNewSample_lambda
%         qEstExpMle
%         p
    end

end % end of nSim loop

%% SLR - simulated sea level rise for milestone years
slr_sim_msYears = f_SLR_Sim(nSim, msYears, rIA, rSC, rOB);
slr_sim_msYears_cm = 100*slr_sim_msYears;

%% add SLR and WL together
qEstExpMLESLR = zeros(nSim, num_yearT, num_msYears); % preallocate
for i = 1:num_msYears;
   % use broadfast function to element-wise add SLR to Est Exp
   % loop through all milestone years
   qEstExpMLESLR(:,:,i) = bsxfun(@plus, qEstExpMle, slr_sim_msYears_cm(:,i));
end

% get quantile values from the combined SLR and WL data
% this is quantile values from joint distribution of SLR and Exponential estimates together
quantile_wl_slr = zeros(length(quantileValues), num_yearT, num_msYears); % preallocate
for i = 1:num_msYears;
    for j = 1:num_yearT;
        quantile_wl_slr(:,j,i) = quantile(qEstExpMLESLR(:,j,i), quantileValues);
    end
end


%% Print for Excel
printForExcelTrue = 1.0;
if printForExcelTrue == 1.0;
    format longG
    for i = 1:num_msYears;
        msYears(i)
        printForExcel = [quantile_wl_slr(4,:,i)' ... % median
                         quantile_wl_slr(5,:,i)' quantile_wl_slr(3,:,i)' ... % 68% upper/lower
                         quantile_wl_slr(6,:,i)' quantile_wl_slr(2,:,i)']    % 90% upper/lower
    end;
end;


toc % stop timer, print duration


