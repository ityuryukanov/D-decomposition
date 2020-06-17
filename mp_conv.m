function w = mp_conv(u,v)
% CONV Nonnumeric convolution.
% W = CONV(U,V) is the convolution of mp vectors U and V.
% length(w) = length(u)+length(v)-1
% w(k) = sum(u(j)*v(k+1-j)), j = max(1,k+1-length(v):min(k,length(u))

% Check number of input arguments:
narginchk(2, 2);

% Check array class:
if ~isa(u, 'mp') || ~isa(v, 'mp')
    error('ERROR(mp_conv.m): u and v must be of type mp.');
end

% Check array size.
if (ndims(u) ~= 2) || (size(u, 2) ~= 1) || ...
   (ndims(v) ~= 2) || (size(v, 2) ~= 1)
    error('ERROR(mp_conv.m): u and v must be column vectors.');
end

% Form Toeplitz matrix from v.
h = [v; zeros((length(u) - 1), 1)];
p = [v(1), zeros(1, (length(u) - 1))];
T = toeplitz(h,p);

% Convolve
w = T*u;

%mp_trunc is NOT used here for postprocessing because the function
%mp_conv (i.e. this function) is used only together with mp_sub(). 
%So mp_trunc() is called in mp_sub and ommitted here for speed
%Input argument checks were also ommitted, since the function is always
%called from Dd.m which defines the argument types and sizes appropriately