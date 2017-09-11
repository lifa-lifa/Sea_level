function r = f_wkbrnd(a,b,g,d,m,varargin)
%WKBRND Random arrays from the Wakeby distribution
%   R = WKBRND(A,B,G,D,M) returns an array of random numbers from the
%   Wakeby distribution with shape parameters B and D, scale parameters 
%   A and G, and location parameter M.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%   R = WKBRND(A,B,G,D,M,P,Q,R,...) or R = WKBRND(A,B,G,D,M,[P,Q,R,...])
%   returns a P-by-Q-by-R-by-... arrray.
%

%   Mike Sheppard
%   Last Modified 27-Mar-2011

if nargin < 5
    error('wkbrnd:TooFewInputs','Requires at least five input arguments.');
end

if isempty(varargin), varargin={1}; end %Scalar

%Use wkbinv
r = f_wkbinv(rand(varargin{:}),a,b,g,d,m);

end