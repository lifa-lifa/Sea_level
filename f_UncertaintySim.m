%% Oeresund A102102
% * LIFA 2017-09-25
% * script to add noise to observations
% INPUTS:
% mu - array of obs considered average values
% sigma - the considered standard deviation
% trunc - truncation threshold, simulated values below trunc will be
% deleted


function [returnValues] = f_UncertaintySim(mu, sigma, trunc)

    % generate values from normal distribution
    normalValues = normrnd(mu, sigma);
    % keep values greater than truncation 
    normalValues = normalValues(normalValues>trunc);

    returnValues = normalValues;    

end
























