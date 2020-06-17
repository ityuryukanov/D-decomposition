function y = mp_polyval(p,x)
%POLYVAL Evaluate polynomial.
%   Y = POLYVAL(P,X) returns the value of a polynomial P evaluated at X. P
%   is a vector of length N+1 whose elements are the coefficients of the
%   polynomial in descending powers.
%
%       Y = P(1)*X^N + P(2)*X^(N-1) + ... + P(N)*X + P(N+1)
%
%   If X is a matrix or vector, the polynomial is evaluated at all
%   points in X.  See POLYVALM for evaluation in a matrix sense.
%
%   Example:
%      Evaluate the polynomial p(x) = 3x^2+2x+1 at x = 5,7, and 9:
%
%      p = [3 2 1];
%      polyval(p,[5 7 9])%
%
%   Class support for inputs P, X: mp, float: double

%   Copyright 1984-2010 The MathWorks, Inc.
%   $Revision: 5.16.4.11 $  $Date: 2010/11/17 11:29:33 $
%   2014/04/22: Modified for mp input arguments

% Check p
if ~isvector(p) || ~isa(p, 'mp')
    error('ERROR(mp_polyval.m): p must be an mp vector');
end
% Check x
if ~isvector(x) || ~(isa(x, 'mp') || isa(x, 'double'))
    error('ERROR(mp_polyval.m): x must be a vector');
end

% Use Horner's method for general case where X is an array.
siz_x = size(x);
nc = length(p);
if nc==1
    y = p(1)*ones(siz_x); 
    return; 
else
    y = mp(zeros(siz_x));
end 
if nc>0, y = p(1); end
for i=2:nc
    y = x .* y + p(i);
end