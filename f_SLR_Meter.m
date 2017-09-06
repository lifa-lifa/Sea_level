% LIFA copied on 2017-09-06 from:
% NEJO matlab function created for Nothern storms, 2017-04-01
% This function return the sea level rise in meter
% With factor1=0.0071 and factor2=0.1054, the sea level rise
% end up to be 0.975 m from 1990 to 2100. This is the average value
% determined by CRES in the report from the environmental ministry 
% published year 2014
% It deduct the isostatic movement of the coast line (input to program).
% It adds the storm contribution (if any).
% 
function y=f_SLR_Meter(year,rateOfIsostaticAdjustmentMeter,...
                            rateOfStormContributionMeter)
numYears = length(year);
% factors in parabolic equation
factor1=0.0071;
factor2=0.1054;
y = zeros(numYears, 1); % initialize array

for i = 1:numYears;                        
    if (year(i) < 1990)
        disp('Error: No climate change effect considered before year 1990');
        % Year is then set to 1990, equal to no climate change effect
        year(i)=1990;
    end

    y(i)=(factor1*(year(i) - 1990)^2 + factor2*(year(i) - 1990))/100-...
        rateOfIsostaticAdjustmentMeter*(year(i)-1990)+...
        rateOfStormContributionMeter*(year(i)-1990);
end
end