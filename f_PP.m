%This function computes the plotting position
%The type can be Weibull, Hazen, Gringoroton and Beard
%
function y=plottingPosition(type,nobs,i)
    switch (type)
        case 'wbl'
            %Weibull plotting positions
            val=i/(nobs+1);
        case 'haz'
            %Hazen plotting positions
            val=(i-0.5)/(nobs+1);
        case 'gri'
            %Gringoroton plotting positions
            val=(i-0.44)/(nobs+0.12);
        case 'bea'
            %Beard plotting positions
            val=(i-0.31)/(nobs+0.38);
        otherwise
            %Keep Weibull default for now, but display warning
            val=i/(nobs+1);
            sprintf('Warning: the plotting position type is not recognized');
    end        
    y=val;
end