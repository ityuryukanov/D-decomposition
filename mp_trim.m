function r = mp_trim(p)
% POLYTRIM Trim polynomial by stripping off leading zeros.
% A modified version of the Peter John Acklam's function

% Check number of input arguments.
error(nargchk(1, 1, nargin));

% Check array class.
if ~isa(p,'mp')
  error('ERROR(mp_trim.m): P must be an mp array.');
end

s=double(p);  % double for speed =) 1e-50 is 1e-50 also in double, so zeros will be zeros only if they are really zeros!

% Check array size.
if (ndims(s) ~= 2) || (size(s, 2) ~= 1)
  error('ERROR(mp_trim.m): P must be column vector.');
end

k = find(s);         % Find non-zero coefficients.
if isempty(k)        % If non were found...
  r = mp(0);         % ...return zero polynomial... %%??? may be in mp() format?
else                 % or else...
  k = min(k);       % ... get index of first ...
  n = length(s);    % ... get length of vector ...
  r = p(k:n);       % ... and assign output polynomial
end
