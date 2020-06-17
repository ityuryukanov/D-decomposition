function out_vec=mp_trunc(in_vec, tol)
% this function truncates very small polynomial coefficients in the in_vec
% input vector that might appear due to round-off errors. The "smallness"
% is determined by the value of tol.

% tol is an estimate of the largest "round-off-error" coefficient (i.e.
% everything larger than tol could be a "normal coefficient"), or the
% largest realistic round-off error (i.e. everything smaller than tol
% can't be a "normal coefficient").

% tol should be positive! For symbolic, vpa or mp in_vec the truncation
% tolerance can be much smaller (e.g. 1e-25 and lower).

% Check array class
if ~isa(in_vec, 'mp')
  error('ERROR(mp_trunc.m): in_vec must be an mp array.');
end

% Check number of input arguments.
error(nargchk(1, 1, nargin));

in_vec_dbl=double(in_vec);  % use double for speed if possible

% Check array size.
if (ndims(in_vec_dbl) ~= 2) || (size(in_vec_dbl, 2) ~= 1)
  error('ERROR(mp_trunc.m): in_vec must be column vector.');
end

if nargin<2
  tol=max(abs(double(in_vec(1:end-1))))*10^(ceil(-0.5*mp.Digits));  % Tolerance for round-off errors. tol=10^(ceil(-0.5*mp.Digits)) or near is for X \in [-1, 1]
end

in_vec = in_vec(:);
if nnz(in_vec)<=1
  out_vec=in_vec;
  return;
end

% Handle pathological round-off cases
b=(abs(in_vec)<=tol);  % eps('double')=2e-16
b1=in_vec(b);       % group of presumably "round-off error" coefficients
b2=in_vec(~b);      % group of presumably "normal" coefficients
if isempty(b2)
  out_vec=mp(0);
  return;
end
if isempty(b1)
  out_vec=in_vec;
  return;
end

cnst = in_vec(end);
in_vec = in_vec(1:end-1);
[b,idx] = sort(abs(in_vec),'ascend');
gaps = b(2:end)./b(1:end-1);
idxd = find(gaps>1/tol,1,'last');
if any(idxd)
  in_vec(idx(1:idxd)) = 0;
end
out_vec=[in_vec;cnst];
end

