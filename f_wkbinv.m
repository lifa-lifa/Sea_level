function x = f_wkbinv(p,a,b,g,d,m)
%WBKINV Inverse of the Wakeby cumulative distribution function
%   X = WKBPDF(P,A,B,G,D,M) returns the inverse of the Wakeby probability
%   density function with shape parameters B and D, scale parameters A
%   and G, and location parameter M.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 27-Mar-2011


if nargin < 6
    error('wkbinv:TooFewInputs',...
          'Requires at least six input arguments.');
end

[errorcode, p, a, b, g, d, m] = distchck(6,p,a,b,g,d,m);

if errorcode > 0
    error('wkbinv:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end

% Initialize x to NaN
if isa(p,'single') || isa(a,'single') || isa(b,'single') || isa(g,'single') || isa(d,'single') || isa(m,'single')
   x = zeros(size(p),'single');
else
   x = zeros(size(p));
end

%Edge Cases
k1=(p==0 & a>0 & b>0 & g>0 & d>0);
if any(k1)
    x(k1)=m(k1);  %The CDF is 0 at mu
end
k2=(p==1 & a>0 & b>0 & g>0 & d>0);
if any(k2)
    x(k2)=Inf;  %The CDF is 1 at Inf
end

%Other valid cases
k=(p>0 & p<1 & a>0 & b>0 & g>0 & d>0);
if any(k)
    pk=p(k);
    ak=a(k);
    bk=b(k);
    gk=g(k);
    dk=d(k);
    mk=m(k);
    x(k)=mk+((ak./bk).*(1-(1-pk).^bk))-((gk./dk).*(1-(1-pk).^(-dk)));
end


end