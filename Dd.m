function regs=Dd(polynom, R, wdw, colors)  % colors are only needed for debugging/intermediate results

P = polynom.P;
Q = polynom.Q;
L = polynom.L;
Ymin = wdw.Ymin;
Ymax = wdw.Ymax;
Xmin = wdw.Xmin;
Xmax = wdw.Xmax;

% Setup:
regs=[];
tol=10^(ceil(-0.85*mp.Digits));  % Tolerance (absolute error) of an mp X. Note that X \in [-1, 1]!!
tol2=1e-12;  % Tolerance in K1, K2 (avoidance of very small regions)
% Don't set it too large, since it's used for something that should be neglected

% Checks:
if ~isscalar(Ymin) || ~isscalar(Ymax) || ~isscalar(Xmin) || ~isscalar(Xmax)
  error('ERROR(Dd.m): Plotting window limits should be nonempty scalars');
end
if ~isa(Ymin, 'double') || ~isa(Ymax, 'double') || ~isa(Xmin, 'double') || ~isa(Xmax, 'double')
  error('ERROR(Dd.m): Plotting window limits should be of type double');
end
if (abs(Ymin-Ymax)<100*tol2)||(abs(Xmin-Xmax)<100*tol2)
  error('ERROR(Dd.m): Check the plotting window limits, they shouldn''t be too narrow');
end
if (Ymin>Ymax)||(Xmin>Xmax)
  error('ERROR(Dd.m): Check the plotting window limits, (k1_min<k1_max)&&(k2_min<k2_max) should hold');
end
if (imag(Ymin)~=0)||(imag(Xmin)~=0)||(imag(Ymax)~=0)||(imag(Xmax)~=0)
  error('ERROR(Dd.m): k1_min, k1_max, k2_min, k2_max should be real numbers');
end
if ~isscalar(R.r) || ~isscalar(R.shft)
  error('ERROR(Dd.m): R.shft and R.r should be nonempty scalars');
end
if ~isa(R.r, 'double') || ~isa(R.shft, 'double')
  error('ERROR(Dd.m): R.shft and R.r should be of type double');
end
if (R.r<=0)||(R.r>1)||(imag(R.r)~=0)
  error('ERROR(Dd.m): R.r should be real, positive and smaller or equal than one');
end
if (abs(R.shft)>1)||(imag(R.r)~=0)
  error('ERROR(Dd.m): Absolute value of R.shft should be smaller or equal than one, R.shft should be real');
end
if (abs(R.shft)-abs(floor(R.shft*1000)/1000))~=0
  error('ERROR(Dd.m): Only 3 decimal places are allowed for R.shft');
end
if (abs(R.r)-abs(floor(R.r*10000)/10000))~=0
  error('ERROR(Dd.m): Only 4 decimal places are allowed for R.r');
end
lp=size(P); lq=size(Q); ll=size(L);
if isempty(P)||isempty(Q)||isempty(L)
  error('ERROR(Dd.m): Neither of P, Q, L polynomial vectors can be empty');
end
if ~isa(P, 'mp') || ~isa(Q, 'mp') || ~isa(L, 'mp')
  error('ERROR(Dd.m): P, Q and L vectors should be of type mp');
end
if (lp(1,2)>1)||(lq(1,2)>1)||(ll(1,2)>1)
  error('ERROR(Dd.m): P, Q, L should be column vectors');
else
  lp=lp(1,1); lq=lq(1,1); ll=ll(1,1);
end
if (lp==1)&&(lq==1)&&(ll==1)
  fprintf('WARNING(Dd.m): At least one of P, Q, L polynomial vectors should be longer than 1.\n Otherwise the characteristic polynomial function degenerates to a scalar (system has no dynamics, the task is ill-posed).\n');
  regs.x=[Xmin, Xmin, Xmax, Xmax];  % there are no poles: no instability, perfect time responces (because of no dynamics!)
  regs.y=[Ymin, Ymin, Ymax, Ymax];
  return;
end
if max([lp,lq,ll])>70
  error('ERROR(Dd.m): The maximal closed-loop system order is constrained to 70');
end

if (R.shft~=0)
  shift=R.shft;
  R=R.r;
  P_new=zeros(lp,1);
  p=[1 shift];  % double seems OK here (err_max=1e-16 -> acceptable)
  for i=1:1:lp
    zn=P(i)*real(ifft(fft(p,lp-i+1).^(lp-i))).';
    P_new=P_new + [zeros(lp-length(zn), 1); zn];
  end
  P=P_new;
  Q_new=zeros(lq,1);
  for i=1:1:lq
    zn=Q(i)*real(ifft(fft(p,lq-i+1).^(lq-i))).';
    Q_new=Q_new + [zeros(lq-length(zn), 1); zn];
  end
  Q=Q_new;
  L_new=zeros(ll,1);
  for i=1:1:ll
    zn=L(i)*real(ifft(fft(p,ll-i+1).^(ll-i))).';
    L_new=L_new + [zeros(ll-length(zn), 1); zn];
  end
  L=L_new;
elseif (R.shft==0)
  R=R.r;
end

% D-d STEP 1. Denominators -> Determinants
% with conversion to Chebyshev
P_re=zeros(lp,1);
Q_re=zeros(lq,1);
L_re=zeros(ll,1);
for i=1:1:lp
  CHEB_P=P(i)*(R^(lp-i))*chebpoly(lp-i);
  P_re=P_re + [zeros(lp-length(CHEB_P), 1); CHEB_P];
end
for i=1:1:lq
  CHEB_Q=Q(i)*(R^(lq-i))*chebpoly(lq-i);
  Q_re=Q_re + [zeros(lq-length(CHEB_Q), 1); CHEB_Q];
end
for i=1:1:ll
  CHEB_L=L(i)*(R^(ll-i))*chebpoly(ll-i);
  L_re=L_re + [zeros(ll-length(CHEB_L), 1); CHEB_L];
end
P_im=mp(zeros(lp,1));  % mp() is important if P_im=[0] (i.e. lp==1)
Q_im=mp(zeros(lq,1));
L_im=mp(zeros(ll,1));
for i=1:1:(lp-1)
  CHEB_P=P(i)*(R^(lp-i))*cheb2poly(lp-i-1);
  P_im=P_im + [zeros(lp-length(CHEB_P), 1); CHEB_P];
end
for i=1:1:(lq-1)
  CHEB_Q=Q(i)*(R^(lq-i))*cheb2poly(lq-i-1);
  Q_im=Q_im + [zeros(lq-length(CHEB_Q), 1); CHEB_Q];
end
for i=1:1:(ll-1)
  CHEB_L=L(i)*(R^(ll-i))*cheb2poly(ll-i-1);
  L_im=L_im + [zeros(ll-length(CHEB_L), 1); CHEB_L];
end

DET_0=mp_sub(mp_conv( P_re,Q_im),mp_conv( Q_re,P_im));
DET_1=mp_sub(mp_conv(-L_re,Q_im),mp_conv(-L_im,Q_re));
DET_2=mp_sub(mp_conv(-L_im,P_re),mp_conv(-L_re,P_im));
if (all(DET_0==0))&&(all(DET_1==0))&&(all(DET_2==0))
  error('ERROR(Dd.m): All three determinants are identically zero');  %for completeness...
end

% Dd STEP 2: Critical points
% fl_D0 is flag for DET_0~=0 (e.g. length(P)=length(Q)=1 -> DET_0==0) and
% DET_1~=0 or DET_1~=0 (e.g., L==0 -> {DET_1==0}&&{DET_2==0}&&{DET_0~=0})
fl_D0=any(abs(DET_0)>=tol)&&(any(abs(DET_1)>=tol)||any(abs(DET_2)>=tol));  
% Check DET_1 and DET_2 at critical points:
x_sng=[];  % x_sng values to deflate DET_1/DET_2/DET_0
if fl_D0   % if DET_0 and at least one of DET_1, DET_2 are nonzero: find singular x and possibly deflate DET_0, DET_1, DET_2
  x_cr=rootsX(DET_0, tol);
  if ~isempty(x_cr)
    if any(abs(DET_1)>=tol)  % (DET_0(x*)==0)&&(DET_1(x*)==0)=>DET_2(x*)==0 -> it suffices to check either DET_1(x*) or DET_2(x*), x* \in x_cr
      x_D=rootsX(DET_1, tol);
    elseif any(abs(DET_2)>=tol)
      x_D=rootsX(DET_2, tol);
    end
    if ~isempty(x_D)
      scale=mp(num2str(10^(round(mp.Digits*0.9))));
      x_sng=intersect(round(x_D*scale)/scale, round(x_cr*scale)/scale);  % x_sng: one value of x produces not a solution point [k2 k1], but a straight line k1=a*k2+b      
    end    
    if ~isempty(x_sng)
      fprintf('WARNING(Dd.m): Extra singular x besides {-1,+1} have been found: %d \n', double(x_sng));
      DEFLATE=mp_charpoly(x_sng);  % Deflate x_sng from DET_0, DET_1, DET_2
      a=length(DEFLATE); b=length(DET_0);
      DET_0=ifft(fft(DET_0)./fft(DEFLATE,b));
      DET_0=real(DET_0(1:b-a+1));
      b=length(DET_1);
      if any(abs(DET_1)>=tol)
        DET_1=ifft(fft(DET_1)./fft(DEFLATE,b));
        DET_1=real(DET_1(1:b-a+1));
      end
      b=length(DET_2);
      if any(abs(DET_2)>=tol)
        DET_2=ifft(fft(DET_2)./fft(DEFLATE,b));
        DET_2=real(DET_2(1:b-a+1));
      end
    end
  end
  x_dfl=x_sng;
  x_sng=unique([x_sng; -1; 1]);  % x=-1 and x=1 critical points were "lost" by division by sin(OMEGA)
else
  x_cr=[];
  if all(abs(DET_0)<=tol)
    if any(abs(DET_1)>=tol)
      x_sng=rootsX(DET_1, tol);
      x_sng=unique([x_sng; -1; 1]);
    elseif any(abs(DET_2)>=tol)
      x_sng=rootsX(DET_2, tol);
      x_sng=unique([x_sng; -1; 1]);
    end
  elseif all(abs(DET_1)<=tol)||all(abs(DET_2)<=tol)
    x_sng=rootsX(DET_0, tol);
    x_sng=unique([x_sng; -1; 1]);  % x=-1 and x=1 critical points were "lost" by division by sin(OMEGA)
  end
end
[~, i1, ~] = unique([mp_polyval(P_re,x_sng),mp_polyval(Q_re,x_sng),mp_polyval(L_re,x_sng)], 'rows', 'stable');
x_sng=x_sng(i1);  % removing x_sng corresponding to identical line boundaries

% Dd STEP 3: Finding and plotting the solution curves:
% Find the relevant ranges of x for plotting K1(K2) in their requested ranges:
K1offset=1e-3*max([abs(Ymax), abs(Ymin)]);
K2offset=1e-3*max([abs(Xmax), abs(Xmin)]);
x_range=double.empty(1,0);  % initialize to empty
if fl_D0
  if isempty(x_dfl)
    x0=x_cr;
  else
    x0=setdiff(x_cr,x_dfl);  % this should be faster than roots(DET_0)
  end
  x_kp_up=find_int(DET_1, DET_0, x0, Ymax+K1offset, tol, 'max');
  x_kp_lw=find_int(DET_1, DET_0, x0, Ymin-K1offset, tol, 'min');
  if (~isempty(x_kp_lw))&&(~isempty(x_kp_up))
    x_kp=rintersect(x_kp_up, x_kp_lw);
  else
    x_kp=[];
  end
  x_kic_up=find_int(DET_2, DET_0, x0, Xmax+K2offset, tol, 'max');
  x_kic_lw=find_int(DET_2, DET_0, x0, Xmin-K2offset, tol, 'min');
  if (~isempty(x_kic_lw))&&(~isempty(x_kic_up))
    x_kic=rintersect(x_kic_lw, x_kic_up);
  else
    x_kic=[];
  end
  if (~isempty(x_kp))&&(~isempty(x_kic))
    x_range=rintersect(x_kp, x_kic);
  else
    x_range=double.empty(1,0);
  end
  if ~isempty(x_range)
    x_range(:, x_range(1,:)==x_range(2,:))=[];  % NO zero length intervals!
    if ~isempty(x_range)
      x_r=mat2cell(x_range, 2, ones(1,length(x_range(1,:))));  % converting x_range to cell for more flexibility
    else
      fprintf('WARNING(Dd.m): The CRB curve doesn''t cross the chosen plot region. Proceeding...\n');
      x_range=double.empty(1,0);
    end
  else
    fprintf('WARNING(Dd.m): The CRB curve doesn''t cross the chosen plot region. Proceeding...\n');
    x_range=double.empty(1,0);
  end
end

k=1;
M=length(x_range(1,:))+length(x_sng)+4;  % number of curve pieces in the plotting window
LE=true;    % marker for intersections and line ends
lines=repmat(struct('x',0,'y',0,'i',0,'type',''), M, 1); % x, y, i stands for x-y-intersections
K1_xP1=NaN; % CRB begins and ends on the RRBs for x=-1 and x=1. Other possible
K2_xP1=NaN; % singular points x_sng generally intersect the CRB, but they are
K1_xM1=NaN; % neither the beginning nor the endpoints of the CRB.
K2_xM1=NaN; % P=='+', M=='-'
iP=NaN; iM=NaN; % if x_range is empty NaNs will 100% persist ->
jP=NaN; jM=NaN; % all NaN variables initialized here will have no effect

if fl_D0&&(~isempty(x_range))
  
  % Find the extrema of k1(x) and k2(x) (for representative plotting):
  DER1 = mp_polyder(DET_1, DET_0);
  DER2 = mp_polyder(DET_2, DET_0);
  x11=roots(DER1);  % zeros of d{k1(x)}/dx
  x22=roots(DER2);  % zeros of d{k2(x)}/dx
  x_ex=[x11; x22];
  x_ex(abs(x_ex)>1)=[];  % x belongs to [-1 1] because x=cos(OMEGA)
  if ~isempty(x_ex)
    x_ex(abs(imag(x_ex))>tol)=[];  % only real roots!
  end
  if ~isempty(x_ex)
    x_ex=real(x_ex);
    x_ex=unique(x_ex, 'sorted');  %(!) unique() also sorts in ascending order
    X_EX=repmat(x_ex, 1, length(x_range(1,:)));
    x_lw=repmat(x_range(1,:), length(x_ex), 1);
    x_up=repmat(x_range(2,:), length(x_ex), 1);
    a=X_EX>x_lw;
    b=X_EX<x_up;
    c=a&b;
    for i=1:1:length(c(1,:))
      if any(c(:,i))
        x_r{i}=[x_range(1,i); x_ex(c(:,i)==1); x_range(2,i)];
      end
    end
  end
  
  % CRB-curve:
  for i=1:1:length(x_r)
    x_r{i}=x_r{i}.';
    n_pts = 15;  % specify the number of intervening points
    x_sol = cumsum([x_r{i};repmat([diff(x_r{i})/n_pts,0],n_pts-1,1)]);
    x_sol = x_sol(1:end-n_pts+1);
    K1_sol=double(mp_polyval(DET_1,x_sol)./mp_polyval(DET_0,x_sol));
    K2_sol=double(mp_polyval(DET_2,x_sol)./mp_polyval(DET_0,x_sol));
    XY=unique([K1_sol.', K2_sol.'], 'rows', 'stable');  % 'stable' is critical!
    if length(XY(:,1))<2
      fprintf('WARNING(Dd.m): A piece of the CRB curve is represented by a single point. Ignoring it...\nThis issue may affect the final result.\n')
      continue;
    end
    K1_sol=XY(:,1).';
    K2_sol=XY(:,2).';
    if x_range(2,i)==1
      K1_xP1=K1_sol(length(K1_sol));
      K2_xP1=K2_sol(length(K2_sol));
      iP=k;
    end
    if x_range(1,i)==-1
      K1_xM1=K1_sol(1);
      K2_xM1=K2_sol(1);
      iM=k;
    end
    lines(k).x=K2_sol;
    lines(k).y=K1_sol;
    lines(k).i=[LE, zeros(1,length(K2_sol)-2), LE];
    lines(k).type='CRB';
    k=k+1;
  end
  
end

% RRB (and possibly CRB) lines:
% All straight line borders have their type 'RRB' (i.e. lines(i).type='RRB')
dK2=abs(Xmax-Xmin);  % +/- distance around k20: for solving inequalities
for i=1:1:length(x_sng)
  HorRRB=false; VerRRB=false;
  if abs(mp_polyval(Q_re,x_sng(i)))<1e-15  % if horizontal straight border
    HorRRB=true;
    a=-double(mp_polyval(L_re,x_sng(i))/mp_polyval(P_re,x_sng(i)));  % k1 value if straight border is horizontal
    if (a>(Ymin-K1offset))&&(a<(Ymax+K1offset))  % if k1 value inside [k1_min-K1offset, k1_max+K1offset]
      k2_pts=[Xmin-K2offset; Xmax+K2offset];
      K1_res=[a; a];
    else
      continue;
    end
  end
  if abs(mp_polyval(P_re,x_sng(i)))<1e-15  % if vertical straight border
    VerRRB=true;
    b=-double(mp_polyval(L_re,x_sng(i))/mp_polyval(Q_re,x_sng(i)));  % k2 value if straight border is vertical
    if (b>(Xmin-K2offset))&&(b<(Xmax+K2offset))  % if k2 value inside [k2_min-K2offset, k2_max+K2offset]
      k2_pts=[b; b];
      K1_res=[Ymin-K1offset; Ymax+K1offset];
    else
      continue;
    end
  end
  if (~HorRRB)&&(~VerRRB)
    k20=(-mp_polyval(P_re,x_sng(i))*(Ymin-K1offset)-mp_polyval(L_re,x_sng(i)))/(mp_polyval(Q_re,x_sng(i)));
    chk=(-mp_polyval(L_re,x_sng(i))-(k20+dK2)*mp_polyval(Q_re,x_sng(i)))/mp_polyval(P_re,x_sng(i));
    if chk>(Ymin-K1offset)
      k2_lw=[k20; Inf];
    else
      k2_lw=[-Inf; k20];
    end
    k20=(-mp_polyval(P_re,x_sng(i))*(Ymax+K1offset)-mp_polyval(L_re,x_sng(i)))/(mp_polyval(Q_re,x_sng(i)));
    chk=(-mp_polyval(L_re,x_sng(i))-(k20+dK2)*mp_polyval(Q_re,x_sng(i)))/mp_polyval(P_re,x_sng(i));
    if chk<(Ymax+K1offset)
      k2_up=[k20; Inf];
    else
      k2_up=[-Inf; k20];
    end
    k2_pts=rintersect(k2_lw, k2_up);
    k2_pts=rintersect(k2_pts, [Xmin-K2offset; Xmax+K2offset]);
    if ~isempty(k2_pts)
      K1_res=(-mp_polyval(L_re, x_sng(i))-k2_pts*mp_polyval(Q_re, x_sng(i)))./mp_polyval(P_re, x_sng(i));
      k2_pts=double(k2_pts);
      K1_res=double(K1_res);
    end
  elseif (HorRRB)&&(VerRRB)
    disp('WARNING(Dd.m): One RRB straight border doesn''t seem to exist!');
    continue;
  end
  if ~isempty(k2_pts)
    if ~isnan(K1_xP1)&&~isnan(K2_xP1)&&(x_sng(i)==1)
      k2_pts=[k2_pts(1,1); K2_xP1; k2_pts(2,1)];  % unique() is inserted before storing of the lines.. For the mp->double imprecision
      K1_res=[K1_res(1,1); K1_xP1; K1_res(2,1)];  % and for the newly inserted points equal to the existing ones..
      jP=k;
    end
    if ~isnan(K1_xM1)&&~isnan(K2_xM1)&&(x_sng(i)==-1)
      k2_pts=[k2_pts(1,1); K2_xM1; k2_pts(2,1)];
      K1_res=[K1_res(1,1); K1_xM1; K1_res(2,1)];
      jM=k;
    end
    XY=unique([K1_res, k2_pts], 'rows', 'stable');  % 'stable' is critical!
    if length(XY(:,1))<2
      fprintf('WARNING(Dd.m): A piece of the straight border is represented by a single point. Ignoring it...\nThis issue may affect the final result.\n')
      continue;
    end
    lines(k).x=XY(:,2).';
    lines(k).y=XY(:,1).';
    lines(k).i=[LE, LE*ones(1,length(k2_pts)-2), LE];
    lines(k).type='RRB';
    k=k+1;
  end
end
% Border lines are straight lines, so they don't actually require minimal density of points
lines(k).x=linspace(Xmin,Xmin,2);
lines(k).y=linspace(Ymin,Ymax,2);
lines(k).i=[LE, LE];
lines(k).type='border';
lines(k+1).x=linspace(Xmax,Xmax,2);
lines(k+1).y=linspace(Ymin,Ymax,2);
lines(k+1).i=[LE, LE];
lines(k+1).type='border';
lines(k+2).x=linspace(Xmin,Xmax,2);
lines(k+2).y=linspace(Ymin,Ymin,2);
lines(k+2).i=[LE, LE];
lines(k+2).type='border';
lines(k+3).x=linspace(Xmin,Xmax,2);
lines(k+3).y=linspace(Ymax,Ymax,2);
lines(k+3).i=[LE, LE];
lines(k+3).type='border';
lines(k+4:end)=[];  % delete the excessive entries

% Dd STEP 4: Divide curves/lines into segments:
numseg=M;   % maximal total number of segments in all curves (used later for array preallocation)

% Search for intersections:
for i=1:1:length(lines)-4  % -4 because the last 4 lines are the borders
  for j=(i+1):1:length(lines)
    if (i==iP)&&(j==jP)  % special check for introduced CRB-RRB intersections at x=+/-1
      linesX=lines(i).x(1:end-1);
      linesY=lines(i).y(1:end-1);
    elseif (i==iM)&&(j==jM)  % special check for introduced CRB-RRB intersections at x=+/-1
      linesX=lines(i).x(2:end);
      linesY=lines(i).y(2:end);
    else
      linesX=lines(i).x;
      linesY=lines(i).y;
    end
    [x1,y1,seg1,seg2] = intersections(linesX, linesY, lines(j).x, lines(j).y);
    if (~isempty(x1))&&(~isempty(y1))
      if (i==iM)&&(j==jM)
        seg1=seg1+1;
      end
      [seg1, i1]=sort(seg1, 'ascend'); seg1=floor(seg1);
      [seg2, i2]=sort(seg2, 'ascend'); seg2=floor(seg2);
      x2=x1; y2=y1;
      x1=x1(i1); y1=y1(i1);
      x2=x2(i2); y2=y2(i2);
      k=1;
      for p=1:1:length(seg1)
        a=ismember([lines(i).x.', lines(i).y.'], [x1(p),y1(p)], 'rows');
        if any(a)
          if ~lines(i).i(a)
            lines(i).i(a)=LE;
            numseg=numseg+1;
          end
        else
          lines(i).x=[lines(i).x(1:(seg1(p)+k-1)), x1(p), lines(i).x((seg1(p)+k):end)];
          lines(i).y=[lines(i).y(1:(seg1(p)+k-1)), y1(p), lines(i).y((seg1(p)+k):end)];
          lines(i).i=[lines(i).i(1:(seg1(p)+k-1)), LE,    lines(i).i((seg1(p)+k):end)];
          numseg=numseg+1;  % if nothing was removed..
          k=k+1;
        end
      end
      k=1;
      for p=1:1:length(seg2)
        a=ismember([lines(j).x.', lines(j).y.'], [x2(p),y2(p)], 'rows');
        if any(a)
          if ~lines(j).i(a)
            lines(j).i(a)=LE;
            numseg=numseg+1;
          end
        else
          lines(j).x=[lines(j).x(1:(seg2(p)+k-1)), x2(p), lines(j).x((seg2(p)+k):end)];
          lines(j).y=[lines(j).y(1:(seg2(p)+k-1)), y2(p), lines(j).y((seg2(p)+k):end)];
          lines(j).i=[lines(j).i(1:(seg2(p)+k-1)), LE,    lines(j).i((seg2(p)+k):end)];
          numseg=numseg+1;  % if nothing was removed..
          k=k+1;
        end
      end
    end
  end
end

for i=1:1:length(lines)-4  % -4 because the last 4 lines are the borders
  if strcmp(lines(i).type,'CRB')
    [x1,y1,seg1,seg2] = intersections(lines(i).x, lines(i).y);
    if (~isempty(x1))&&(~isempty(y1))
      [seg, i1]=sort([seg1; seg2], 'ascend'); seg=floor(seg);
      x=[x1; x1]; x=x(i1);
      y=[y1; y1]; y=y(i1);
      for k=1:1:length(seg)
        lines(i).x=[lines(i).x(1:(seg(k)+k-1)), x(k), lines(i).x((seg(k)+k):end)];
        lines(i).y=[lines(i).y(1:(seg(k)+k-1)), y(k), lines(i).y((seg(k)+k):end)];
        lines(i).i=[lines(i).i(1:(seg(k)+k-1)), LE,   lines(i).i((seg(k)+k):end)];
        numseg=numseg+3;  % '+3' since all loops will be brocken in 3 parts, and loops are from self-intersections
      end
    end
  end
end

% Delete the remaining excessive points:
i=1;
while i<=length(lines)
  a=(lines(i).y<(Ymin-tol))|(lines(i).y>(Ymax+tol));
  b=(lines(i).x<(Xmin-tol))|(lines(i).x>(Xmax+tol));
  c=~(a|b);
  if all(c == 0)
    lines(i)=[];
    i=i-1;
  elseif all(c == 1)
    i=i+1;
  else
    idx1 = find(diff([~c(1) c]));
    idx2 = circshift(idx1-1,[0 -1]);
    idx2(end) = length(c);
    pcs = [idx1; idx2];
    pcs(:,c(pcs(1,:))==0)=[];  % keep only parts of curve within [k2_min, k2_max; k1_min, k1_max]
    k=length(pcs(1,:));
    lines(i+k:end+k-1)=lines(i+1:end);
    old=lines(i);
    for j=1:1:k
      lines(i+j-1).x=old.x(pcs(1,j):pcs(2,j));
      lines(i+j-1).y=old.y(pcs(1,j):pcs(2,j));
      lines(i+j-1).i=old.i(pcs(1,j):pcs(2,j));
      lines(i+j-1).type=old.type;
    end
    i=i+k;
  end
end

% Intermediate plot:
%{
figure;
for i=1:1:length(lines)
  plot(lines(i).x, lines(i).y, '-', 'Color', 'k', 'LineWidth', 1, 'MarkerSize', 12);
  hold on;
end
%close;
%}

% Dd STEP 5: Find all faces
% Extract graph edges between intersections:
edges=repmat(struct('x',0,'y',0), numseg+20, 1);  % '+20' to be sure...
k=1;
for i=1:1:length(lines)
  s=find(lines(i).i);  % indexes of segments between intersections
  for j=1:1:(length(s)-1)
    edges(k).x=lines(i).x(s(j):s(j+1));
    edges(k).y=lines(i).y(s(j):s(j+1));
    if isequal([edges(k).x(1); edges(k).y(1)], [edges(k).x(end); edges(k).y(end)])
      if length(edges(k).x)>3
        ax=edges(k).x(1:2);  % break loops!
        ay=edges(k).y(1:2);
        bx=edges(k).x(2:end-1);
        by=edges(k).y(2:end-1);
        cx=edges(k).x(end-1:end);
        cy=edges(k).y(end-1:end);
        edges(k).x=ax;
        edges(k).y=ay;
        k=k+1;
        edges(k).x=bx;
        edges(k).y=by;
        k=k+1;
        edges(k).x=cx;
        edges(k).y=cy;
      else
        k=k-1;  % if loop is degenerate, simply remove it
      end
    end
    k=k+1;
  end
end
edges(k:end)=[];

% Intermediate plot (segment-by-segment):
%{
c=11;
%figure;
for i=1:1:length(edges)
  plot(edges(i).x, edges(i).y, '.-', 'Color', colors{c}, 'LineWidth', 2.5, 'MarkerSize', 1); % , 'MarkerEdgeColor', 'r'
  hold on;
  c=c+1;
  if c>length(colors)
    c=1;
  end
end
%}

regions=faces(edges,Ymin,Ymax,Xmin,Xmax,colors); % Intermediate check of detected regions: see face_check.m

% Dd STEP 6: Check every face for stability
regs=face_check(regions, P, Q, L, R, colors);

% Plot to check the result (e.g. for debugging):
%{
hold on;
for i=1:1:length(regs)
  plot(regs(i).x, regs(i).y, '.-', 'Color', colors{polynom.l},... 
    'LineWidth', 1.0, 'MarkerSize', 10);
end
%hold off;
% grid on;
%}

end