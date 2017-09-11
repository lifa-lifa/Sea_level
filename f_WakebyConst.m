% LIFA 2017-09-11 
% INPUTS: array of L-moments, l1 to l5
% OUTPUTS: array of D (Wakeby constants)


function [D] = f_WakebyConst(L)
    
    N1 = 3*L(2) - 25*L(3) + 32*L(4);
    N2 = -3*L(2) + 5*L(3) + 8*L(4);
    N3 = 3*L(2) + 5*L(3) + 2*L(4);
    
    C1 = 7*L(2) - 85*L(3) + 203*L(4) - 125*L(5);
    C2 = -7*L(2) + 25*L(3) + 7*L(4) - 25*L(5);
    C3 = 7*L(2) + 5*L(3) - 7*L(4) - 5*L(5);
    
    D1 = N2*C3 - N3*C2;
    D2 = N1*C3 - N3*C1;
    D3 = N1*C2 - N2*C1;
    
    D = [D1, D2, D3];    
end