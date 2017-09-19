function p = f_wkbcdf(x,a,b,g,d,m)
%WBKCDF Wakeby cumulative distribution function
%   P = WKBCDF(X,A,B,G,D,M) returns the CDF of the Wakeby distribution
%   with shape parameters B and D, scale parameters A and G, and location
%   parameter M.
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   WKBCDF uses Newton's method to converge to the solution, as there is no
%   closed form expression other than for WKBINV.

%   Mike Sheppard
%   Last Modified 27-Mar-2011


%   Algorithm based on 'betainv' revision 2.11.2.5 (2004/12/06)
%   Which uses the Newton's Method

if nargin < 6
    error('wkbcdf:TooFewInputs',...
          'Requires at least six input arguments.'); 
end

[errorcode, x, a, b, g, d, m] = distchck(6,x,a,b,g,d,m);

if errorcode > 0
    error('wkbstat:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end

% Initialize p to NaN
if isa(x,'single') || isa(a,'single') || isa(b,'single') || isa(g,'single') || isa(d,'single') || isa(m,'single')
   p = zeros(size(x),'single');
   seps=sqrt(eps('single'));
else
   p = zeros(size(x));
   seps=sqrt(eps);
end


%The cdf of Inf is 1, and the CDF of m is 0
k0=find(x==Inf & a>0 & b>0 & g>0 & d>0);
if any(k0)
    p(k0)=1;
end
k1=find(x==m & a>0 & b>0 & g>0 & d>0);
if any(k1);
    p(k1)=0;
end


%Newton's Method
%Permit no more than count_limit iterations
count_limit=100;
count=0;

k=find(x>m & isfinite(x) & a>0 & b>0 & g>0 & d>0);
if isempty(k)
    return;
end
xk=x(k);


% Use 0.5 a starting guess
pk=0.5.*ones(size(xk));
if isa(x,'single')
    pk=single(pk);
end



h=ones(size(xk));
crit=seps;

%Break out of the iteration loop for the following:
% 1) The last update is very small (compared to p).
% 2) The last update is very small (compared to 100*eps)
% 3) There are more than 100 iterations. This should NEVER happen.

while(any(abs(h)>crit*abs(pk)) && max(abs(j))>crit && count<count_limit),
    count=count+1;
    %The derivative of the inverse is explicitly known, subfunction below
    h=(f_wkbinv(pk,a(k),b(k),g(k),d(k),m(k)) - xk) ./ wbkinvderv(pk,a(k),b(k),g(k),d(k));
    pnew=pk-h;
    
    %Make sure that the values stay inside the bounds (0<=p<=1)
    %Initially, Newton's Method may take big steps
    psmall = find(pnew<=0);
    plarge = find(pnew>=1);
    if any(psmall) || any(plarge)
        pnew(psmall) = pk(psmall) / 10;
        pnew(plarge) = 1 - (1-pk(plarge))/10;
    end
    pk=pnew;
end

%Return the converged value(s).
p(k)=pk;

if count==count_limit
    fprintf('\nWarning: WKBCDF did not converge.\n');
    str='The last step was: ';
    outstr=sprintf([str,'%13.8f\n'],max(h(:)));
    fprintf(outstr);
end


%Reshape to original input
p=reshape(p,size(x));

end


function dinv=wbkinvderv(p,a,b,g,d)
%Parameters have passed error checking
%Note: m is not used in the derivative
dinv=(a.*((1-p).^(-1+b)))+(g.*((1-p).^(-1-d)));
end