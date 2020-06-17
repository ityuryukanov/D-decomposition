function a = mp_polyder(u,v)
%POLYDER Differentiate polynomial.
%It's based on the built-in matlab function but works with mp polynomials
%It returns the polynomial derivative OR the numerator of the derivative of
%a quotient of two polynomials (only these options are implemented since
%only they are needed, and 

if nargin < 2, v = 1; end

nu = length(u); nv = length(v);
if nu < 2, up = mp(0); else up = u(1:nu-1) .* (nu-1:-1:1).'; end
if nv < 2, vp = mp(0); else vp = v(1:nv-1) .* (nv-1:-1:1).'; end
a1 = mp_conv(up,v); a2 = mp_conv(u,vp);
a=mp_sub(a1, a2);

%It's based on the built-in matlab function but works with mp polynomials
%mp_trunc() and mp_trim() are called in mp_sub(), so they are ommitted in
%this function...
