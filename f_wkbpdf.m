function y = f_wkbpdf(x,a,b,g,d,m)
%WBKPDF Wakeby probability distribution function
%   y = WKBPDF(X,A,B,G,D,M) returns the PDF of the Wakeby distribution
%   with shape parameters B and D, scale parameters A and G, and location
%   parameter M.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 28-Mar-2011


if nargin < 6
    error('wkbpdf:TooFewInputs',...
          'Requires at least six input arguments.'); 
end

[errorcode, x, a, b, g, d, m] = distchck(6,x,a,b,g,d,m);

if errorcode > 0
    error('wkbpdf:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end


%Using Eq. (12.87) from "Continuous Univariate Distributions" Kotz
%The pdf p(x) of X is related to the cdf F(x) by the formula
%p(x)=[a(1-F)^(b-1) + g(1-F)^(-d-1)]^(-1)
F = f_wkbcdf(x,a,b,g,d,m);
S = 1-F;
t = S.^(b+d);
f = (S.^(d+1))./(a.*t+g); %pdf
y = f;


end

