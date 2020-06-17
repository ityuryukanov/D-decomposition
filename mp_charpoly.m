function c = mp_charpoly(A)

n = size(A);
if n(1,1) == n(1,2) % if input is square matrix
    z = eig(A);  
    n = n(1,1);
elseif n(1,1) > n(1,2) % if input is vector of polynomial roots
    z = A;  
    n = n(1,1);
else
    error('ERROR(mp_charpoly.m): The input argument should be either a column vector or a square matrix');
end 

c = mp(zeros(n+1,1)); c(1) = 1;
for j = 1:n
    c(2:j+1) = c(2:j+1)-z(j)*c(1:j);
end

c=real(c);

end