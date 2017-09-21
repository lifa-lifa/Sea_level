%% Oeresund A102102
% * LIFA 2017-09-21
% * Statistics for southern storms, middle segment
% * Based on NEJO's Statistics South-1044-2014-Køge-Wakeby-Mixture.nb
% * Takes Wakeby parameters and performs Monte Carlo, adds SLR
close all
clear all
clc

%% Get data
ReadFolder = strcat(pwd,'\output\WakebyParam\');  % data folder
% load Wakeby parameters (no need to assign to variable, loads as wp_P123)
load([ReadFolder 'WakebyParam_P123.mat']);

%% Inputs
rIA = 0.0011;
rSC = 0.0012;
msYears = [1990 2100];
nSim = 100;

%% SLR - simulated sea level rise for milestone years
slr_sim_msYears = f_SLR_Sim(nSim, msYears, rIA, rSC);

