function Np=clp_check_poly(kp, kic, P, Q, L, R)
% CLP polynomial:
m = length(P);
n = length(Q);
h = length(L);
l = max([m, n, h]);

CLP = zeros(l, 1);   % Initialize output
CLP(l-m+1:l) = kp*P;   % Insert first polynomial
CLP(l-n+1:l) = CLP(l-n+1:l) + kic*Q;   % Add second polynomial
CLP(l-h+1:l) = CLP(l-h+1:l) + L;   % Add third polynomial

% leading zeros trimming is not really necessary

% Counting unstable poles:
p=abs(roots(CLP));
Np=sum(p>R);
end