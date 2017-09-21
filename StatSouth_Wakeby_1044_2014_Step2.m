%% Oeresund A102102
% * LIFA 2017-09-21
% * Statistics for southern storms, middle segment
% * Based on NEJO's Statistics South-1044-2014-Køge-Wakeby-Mixture.nb
% * Takes Wakeby parameters and performs Monte Carlo, adds SLR
close all
clear all
clc

%% Get data
% load Wakeby parameters (no need to assign to variable, loads req parameters from Step 1)
load(strcat(pwd,'\output\WakebyParam\WakebyParam_P123.mat'));

% get observations (only using 1984-2014 data)
ReadFolder = strcat(pwd,'\input\forWakeby\');  % data folder
full_name = [ReadFolder 'Obs_1825_2014.txt'];
fid = fopen(full_name, 'r');
obsP1 = textscan(fid, '%f', 'Delimiter', ',');
fclose(fid);
% move cell array data to matrix, and sort
obsP1 = sort(cell2mat(obsP1), 'descend');  % P1 is 1825-2014

%% Inputs
rng default;  % sets random seed. Matlab default is Mersenne Twister seed 5498
nSim = 100;
rIA = 0.0011;
rSC = 0.0012;
msYears = [1990 2100];
num_yearT = length(yearT);
xMax = 10; % used to clip the Wakeby estimates to 10 times the avg value

% quantiles (Probability of non-exceedance, F = 1-1/T)
q = 1-1./(1 * yearT); 

%% SLR - simulated sea level rise for milestone years
slr_sim_msYears = f_SLR_Sim(nSim, msYears, rIA, rSC);

%% Monte Carlo sim for WL - based on Wakeby distribution
nobs = length(obsP1); % based on AMS, therefore use nobs, instead of Poisson dist.
wkb_MC = cell(nSim, 1); % preallocate

% Generate rand values from Wakeby distribution for nSim times
for i = 1:nSim;
    % use Wakeby parameters, each simulation contains nobs number of events
    MC_temp = f_wkbrnd(wp_P123{:}, [nobs,1]); 
    wkb_MC(i) = {MC_temp}; % store results in cell array
end

% For each MC simulation, fit the Wakeby distribution again
wp_MC = cell(nSim, 1); % preallocate
xval = zeros(nSim, num_yearT); % preallocate
for i = 1:nSim;
   % L-moments
   L_MC = f_lmom(wkb_MC{i}, 5);
   % t values
   t_MC = f_tvalues(L_MC);
   % Wakeby constants
   wc_MC = f_WakebyConst(L_MC);
   l2norm = t_MC(1);
   l3norm = L_MC(3)/L_MC(1);
   % Wakeby parameters
   wp_MC{i} = f_Wakeby(wc_MC(1), wc_MC(2), wc_MC(3), l2norm, l3norm);
   
   % given quantiles, return value from Wakeby distribution
   xval(i,:) = f_wkbinv(q, wp_MC{i}{:}); % acceess all 5 param from row i of wp_MC
end

% clip the data so wkb_est does not exceed xMax
xval = min(xval, xMax);

% multiply l1 of observed data back in for storm surge
wkb_est = mean(obsP1)*xval;





