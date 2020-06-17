function r=mp_sub(p, q)

% Check number of input arguments.
narginchk(2, 2);

% Check array class.
if ~isa(p, 'mp') || ~isa(q, 'mp')
  error('ERROR(mp_sub.m): p and q must be mp arrays.');
end

% Check array size.
if (ndims(p) ~= 2) || (size(p, 2) ~= 1) || ...
    (ndims(q) ~= 2) || (size(q, 2) ~= 1)
  error('ERROR(mp_sub.m): p and q must be column vectors.');
end

m = length(p);
n = length(q);
if m > n
  z = mp(zeros(abs(m-n), 1));
  q = [z; q];
elseif m < n
  z = mp(zeros(abs(m-n), 1));
  p = [z; p];
end

r = p - q;
r = mp_trunc(r);
r = mp_trim(r);
end