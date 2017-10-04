% LIFA modified to select maximum of parabolic or linear functions
% LIFA copied on 2017-09-06 from:
% NEJO matlab function created for Nothern storms, 2017-04-01
% This function return the normalized sea level rise
% It deduct the isostatic movement of the coast line
% and add the storm contribution (if any).
% The sea level rise in year 2100 is set to factor 1.0
% The sea level rise in year 1990 is set to factor 0.0
% 
function y=f_SLR_Norm_ParaLinear(year,rateOfIsostaticAdjustmentMeter,...
                            rateOfStormContributionMeter, rateOfObservedSLR)
                        
numYears = length(year);
y = zeros(numYears, 1); % initialize array                        

for i = 1:numYears;
    % Parabolic function
    parabVal=f_SLR_Meter(year(i),rateOfIsostaticAdjustmentMeter,...
                            rateOfStormContributionMeter);
                        
    parab2100val=f_SLR_Meter(2100,rateOfIsostaticAdjustmentMeter,...
                            rateOfStormContributionMeter);
       
    parabNorm = parabVal/parab2100val;                    
                        
    % Linear function, normalized to 2100 value from Parabolic function                    
    linearVal = rateOfObservedSLR*(year(i) - 1990);
    linearNorm = linearVal/parab2100val;                  
              
    % Chose the max of either function
    y(i)=max(parabNorm, linearNorm);
end 

end