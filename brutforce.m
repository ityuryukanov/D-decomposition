function brutforce(polynom, R, wdw, color)
%{
% Check number of input arguments.
narginchk(4, 4);

% Check array class.
if ~isnumeric(P) || ~isnumeric(Q) || ~isnumeric(L)
    error('ERROR(brutforce.m): P, Q, L must be numeric arrays.');
end

% Check array size.
if (ndims(P) ~= 2) || (size(P, 2) ~= 1) || ...
        (ndims(Q) ~= 2) || (size(Q, 2) ~= 1) || ...
        (ndims(L) ~= 2) || (size(L, 2) ~= 1)
    error('ERROR(brutforce.m): P, Q, L must be column vectors.');
end
%}

P = polynom.P;
Q = polynom.Q;
L = polynom.L;
Ymin = wdw.Ymin;
Ymax = wdw.Ymax;
Xmin = wdw.Xmin;
Xmax = wdw.Xmax;

if (R.shft~=0)&&(R.r<0.95)&&(R.r>0.3)
    shift=R.shft;
    R=R.r;   
elseif (R.shft~=0)&&((R.r>0.95)||(R.r<0.3))
    error('ERROR(bruteforce.m): For R.shft~=0 please input R.r between 0.3 and 0.95')
elseif (R.shft==0)
    R=R.r;
    shift=0;
end

m = length(P);
n = length(Q);
h = length(L);
l = max([m, n, h]);  % highest polynomial degree -> number of roots
N=100;
regs = zeros(2*l, N*N);  % coordinates of points for the "regions"
k = ones(l,1);  % indexes for the regions
Y=linspace(Ymin, Ymax, N);
X=linspace(Xmin, Xmax, N);
for y=Y
    for x=X
        CLP = zeros(l, 1);   % Initialize output.
        CLP(l-m+1:l) = y*P;   % Insert first polynomial.
        CLP(l-n+1:l) = CLP(l-n+1:l) + x*Q;   % Add second polynomial
        CLP(l-h+1:l) = CLP(l-h+1:l) + L;   % Add third polynomial
        p=abs(roots(CLP)-shift);
        Np=sum(p>R);
        for i=1:1:l
            if Np==i-1
                regs(2*i-1:2*i,k(i))=[x;y];
                k(i)=k(i)+1;
            end
        end
    end
end

% Plot the result:
%{
figure;
hold on
for i=1:2:(2*l)
    j=ceil(0.5*i);
    xy=regs(i:i+1,:);
    if k(j)>1
       k(j)=k(j)-1;
    end
    xy(:,k(j):N*N)=[];
    if ~isempty(xy)
        scatter(xy(1,:),xy(2,:),5,color{ceil(i/2)},'fill') % color{abs(ceil(length(color)*rand(1)))}
    end
end
%}
end