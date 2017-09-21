% LIFA 2017-09-21 
% calculates t values from L-moments

function [T] = f_tvalues(L)
    
    % calculate L-CV = t = l_2/l_1
    t = L(2)/L(1);
    % calculate L-moment ratios (t_r = l_r / l_2)
    t3 = L(3)/L(2);
    t4 = L(4)/L(2);
    t5 = L(5)/L(2);

T = [t t3 t4 t5];
    
end



