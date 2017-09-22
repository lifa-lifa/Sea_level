%% Oeresund A102102
% * LIFA 2017-09-22
% * Statistics for southern storms, upper segment (Exponential)
% * Based on NEJO's Statistics South-1044-2013-UpperTail.nb
close all
clear all
clc

%% Inputs
msYears = [1990, 2025, 2050, 2080, 2100]'; % milestone years
% truncation
gammaWE = 240;
% take the larger of gammaWE or historial truncation level
gammaP1 = max(gammaWE, 200); % assumed historical truncation for P1
gammaP2 = max(gammaWE, 240); % ditto P2
gammaP3 = max(gammaWE, 270); % ditto P3

sigmaP1 = 15; % std deviation in P1, 15 cm
sigmaP2 = 45;
sigmaP3 = 60;

%% Data
obsP1 = sort([193, 216, 220, 230, 235, 286], 'descend'); % 1825-2015
obsP2 = sort([258, 290, 301, 366], 'descend'); % 1500-1824
obsP3 = sort([276, 286, 343], 'descend'); % 1044-1500

% truncate data
obsP1Trunc = obsP1(obsP1 > gammaP1);
obsP2Trunc = obsP2(obsP2 > gammaP2);
obsP3Trunc = obsP3(obsP3 > gammaP3);

% join
obsAll = sort([obsP1Trunc obsP2Trunc obsP3Trunc], 'descend');






