% This function return the normalized sea level rise
% It deduct the isostatic movement of the coast line
% and add the storm contribution (if any).
% The sea level rise in year 2100 is set to factor 1.0
% The sea level rise in year 1990 is set to factor 0.0
% 
function y=f_SLR_Norm(year,rateOfIsostaticAdjustmentMeter,...
                            rateOfStormContributionMeter)
    yyval=seaLevelRiseMeter(year,rateOfIsostaticAdjustmentMeter,...
                            rateOfStormContributionMeter);
    y2100val=seaLevelRiseMeter(2100,rateOfIsostaticAdjustmentMeter,...
                            rateOfStormContributionMeter);
    y=yyval/y2100val;
end