% LIFA 2017-09-11 
% estimating Wakeby parameters
% based on 1997 Hosking & Wallis, Regional Frequency Analysis, see A108
% INPUTS: Wakeby constants D1, D2, D3
% OUTPUTS: Wakeby parameters


function [alphaE, betaE, gammaE, deltaE, xiE] = f_Wakeby(D1, D2, D3, l2, l3)
    
    % solve polynomial of D1z^2 + D2z + D3 == 0 for z
    z = roots([D1 D2 D3]);
    l1norm = 1.0;

    % beta and negative delta are the roots of the quadratic eq. 
    % beta being the larger of the two roots
    betaE = max(z); 
    deltaE = -min(z);
    if (betaE + deltaE) < 0; deltaE = 0.00001-betaE; else deltaE = deltaE; end
    if deltaE < 0; deltaE = 0.00001; else deltaE = deltaE; end
    
    % the other parameters: alpha, gamma, xi
    alphaE = (1+betaE)*(2+betaE)*(3+betaE)*((1+deltaE)*l2...
        -(3-deltaE)*l3)/(4*(betaE+deltaE));
   
    gammaE = -(1-deltaE)*(2-deltaE)*(3-deltaE)*((1-betaE)*l2...
        -(3+betaE)*l3)/(4*(betaE+deltaE));

    if gammaE < 0; gammaE = 0.00001; else gammaE = gammaE; end
    if gammaE < 0.000011; deltaE = 0.00001; else deltaE = deltaE; end
    
    xiE = l1norm - alphaE/(1+betaE) - gammaE/(1-deltaE);
    
% P = [alphaE, betaE, gammaE, deltaE, xiE];    
end