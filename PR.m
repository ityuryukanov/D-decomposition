function Kr = PR(Ts, f1, f2, f3, ki0, mode)
% The functions defines a one OR three-term resonant controller

% Check number of input arguments:
narginchk(6, 6);

% Check input argument classes:
if ~isa([Ts, f1, f2, f3, ki0], 'double') || ~isa(mode, 'char')
    error('ERROR(PR.m): Some function arguments have unvalid data types.');
end

% Check input argument sizes:
if any(size(Ts) ~= 1) || any(size(f1) > 1) || any(size(f2) > 1) ...
        || any(size(f3) > 1) || any(size(ki0) ~= 1)
    error('ERROR(PR.m): Ts, f1, f2, f3, ki0 must be scalars');
end

% Check input argument values.
if ~isreal(Ts) || (Ts < 0) || ~isreal(ki0) || (ki0 < 0) 
    error('ERROR(PR.m): Ts and ki0 must be real non-negative numbers');
end
if ~isreal(f1) | (f1 < 0) | f1~=round(f1) | ~isreal(f2) | ...
    (f2 < 0) | f2~=round(f2) | ~isreal(f3) | (f3 < 0) | f3~=round(f3) %#ok<*OR2>
    error('ERROR(PR.m): f1, f2, f3 must be real non-negative integers');
end

% Tracking controller at 50 Hz (with large ki, wc, phase lag)
f0 =50;                   % frequency to track (usually 50 Hz)
w0 =2*pi*f0;
wc0=3;
phi=0;
b2= tan((Ts*w0)/2)*cos(phi) - tan((Ts*w0)/2)^2*sin(phi);
b1=-2*tan((Ts*w0)/2)^2*sin(phi);
b0=-sin(phi)*tan((Ts*w0)/2)^2 - cos(phi)*tan((Ts*w0)/2);
a2= w0*tan((Ts*w0)/2)^2 + 2*wc0*tan((Ts*w0)/2) + w0;
a1= 2*w0*tan((Ts*w0)/2)^2 - 2*w0;
a0= w0*tan((Ts*w0)/2)^2 - 2*wc0*tan((Ts*w0)/2) + w0;
Kr50=ki0*tf([b2, b1, b0], [a2, a1, a0], Ts); 

% Summation of Kr50+{Kr250+Kr350+...}:
if strcmp(mode,'full')
    if ~isempty(f1)
      w0 =2*pi*f1;
      wc0=1;
      ki0=8;
      b2= sin(w0*Ts);
      b1= 0;
      b0=-sin(w0*Ts);
      a2= 2*w0+2*wc0*sin(w0*Ts);
      a1=-4*w0*cos(w0*Ts);
      a0= 2*w0-2*wc0*sin(Ts*w0);
      a1=a1/a2;
      a0=a0/a2;
      b2=b2/a2;
      b0=b0/a2;
      Az=[-a1 -a0; 1 0];
      Bz=[1; 0];
      Cz=[ki0*wc0*(b1-b2*a1), ki0*wc0*(b0-b2*a0)];
      Dz= b2;
      Kr150=ss(Az, Bz, Cz, Dz, Ts); % state-space resonant term model may have better numerical properties
    else
      Kr150=[];
    end
    
    if ~isempty(f2)
      w0 =2*pi*f2;
      wc0=1;
      ki0=8;
      b2= sin(w0*Ts);
      b1= 0;
      b0=-sin(w0*Ts);
      a2= 2*w0+2*wc0*sin(w0*Ts);
      a1=-4*w0*cos(w0*Ts);
      a0= 2*w0-2*wc0*sin(Ts*w0);
      a1=a1/a2;
      a0=a0/a2;
      b2=b2/a2;
      b0=b0/a2;
      Az=[-a1 -a0; 1 0];
      Bz=[1; 0];
      Cz=[ki0*wc0*(b1-b2*a1), ki0*wc0*(b0-b2*a0)];
      Dz= b2;
      Kr250=ss(Az, Bz, Cz, Dz, Ts);
    else
      Kr250=0;
    end
    
    if ~isempty(f3)
      w0 =2*pi*f3;
      wc0=1;
      ki0=8;
      b2= sin(w0*Ts);
      b1= 0;
      b0=-sin(w0*Ts);
      a2= 2*w0+2*wc0*sin(w0*Ts);
      a1=-4*w0*cos(w0*Ts);
      a0= 2*w0-2*wc0*sin(Ts*w0);
      a1=a1/a2;
      a0=a0/a2;
      b2=b2/a2;
      b0=b0/a2;
      Az=[-a1 -a0; 1 0];
      Bz=[1; 0];
      Cz=[ki0*wc0*(b1-b2*a1), ki0*wc0*(b0-b2*a0)];
      Dz= b2;
      Kr350=ss(Az, Bz, Cz, Dz, Ts);
    else
      Kr350=0;
    end
    
    Kr=Kr50+Kr150+Kr250+Kr350;
elseif strcmp(mode,'50Hz')
    Kr=Kr50;   
else 
    error('ERROR(PR.m): The mode option can be ''full'' or ''50Hz''. Anything else is not allowed.');
end
end