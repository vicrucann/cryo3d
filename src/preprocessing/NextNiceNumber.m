function [n factors]=NextNiceNumber(x,f,step)
% function n=NextNiceNumber(x,f,step)
% Find the next integer >= x that is a multiple of step and whose largest
% prime factor is f (default=5).  This is used for finding a nice
% dimension of vector for use with mixed-radix FFTs and other algorithms.
% The default value for step is 2, so an even result is returned.
% If step is negative, the returned value
% is the next lower nice number that is a multiple of abs(step).

if (nargin<3) || step==0
    step=2;
end;

if nargin<2
    f=5;
end;

% start with n being a multiple of mult.
as=abs(step);
factors=factor(as);
if max(factors)>f
    error(['NextNiceNumber: step size =' num2str(step) ...
          ' has a factor larger than f =' num2str(f)]);
end;

if step>0
    n=as*ceil(x/as);
else
    n=as*floor(x/as);
end;

% Increment n until the largest factor is <= f.
factors=factor(n);
while max(factors)>f
    n=n+step;
    factors=factor(n);
end;
