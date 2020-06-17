function x_out=find_int(inX, in0, x_in0, offset, tol, mode)
% inX is the nominator, in0 is the denominator

% Check number of input arguments.
narginchk(6, 6);

% Check array class.
if ~isa(inX, 'mp') || ~isa(in0, 'mp')
    error('ERROR(find_int.m): inX and in0 must be mp arrays.');
end

% Check array size.
if (ndims(inX) ~= 2) || (size(inX, 2) ~= 1) || ...
        (ndims(in0) ~= 2) || (size(in0, 2) ~= 1)
    error('ERROR(find_int.m): inX and in0 must be column vectors.');
end
% x_out = mp('0');
inX=mp_sub(inX,offset*in0);
x_r1=rootsX(inX, tol); 
x_r=[-1; x_r1; x_in0; 1];  % x_in0 was set to real before, as it is basically rootsX(DET_0)
x_r=unique(x_r,'sorted'); 
a=[x_r(1:length(x_r)-1), x_r(2:length(x_r))];

b=0.5*(a(:,1)+a(:,2));
test=mp_polyval(inX, b)./mp_polyval(in0, b);  % check the intervals
if strcmp(mode,'max')
    b=test<0;
    x_out(1,:)=a(b,1).';
    x_out(2,:)=a(b,2).';
elseif strcmp(mode,'min')
    b=test>0;
    x_out(1,:)=a(b,1).';
    x_out(2,:)=a(b,2).';
else
    error('ERROR(find_int.m): char argument {mode} can only take 2 values: ''min'' and ''max''');  
end

end