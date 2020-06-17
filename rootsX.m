% This function finds the roots of polynomials and filters out all roots
% outside the Chebyshev interval [-1 1]

function X = rootsX(polynom, tol)

% Check number of input arguments.
narginchk(2, 2);

% Check array class.
if ~isa(polynom, 'mp') || ~isa(tol, 'double')
  error('ERROR(rootsX.m): polynom must be an mp array, tol must be a double scalar.');
end

% Check array size.
if (~ismatrix(polynom)) || (size(polynom, 2) ~= 1) || (~isscalar(tol))
  error('ERROR(rootsX.m): polynom must be a column vector, tol must be a scalar.');
end

% Check tol sign.
if tol<0
  error('ERROR(rootsX.m): tol must be nonnegative.');
end

X=roots(polynom);
X(abs(X)>1)=[];  % x belongs to [-1 1] because x=cos(OMEGA)
if ~isempty(X)
  X(abs(imag(X))>tol)=[];  % only real roots are relevant!
else
  X=[];  % to enforce 0x0 dimension
  return;
end
if ~isempty(X)
  X=real(X);
else
  X=[];  % to enforce 0x0 dimension
end
if ~isempty(X)
  assert(all(abs(mp_polyval(polynom,X))<tol),...
    'rootsX.m: Polynom is not exactly zero at its computed roots.');
end
end