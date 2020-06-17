clear;
clc;
addpath('PolygonClipper', 'inpoly', 'rintersect', 'intersections', 'p_poly_dist');
mp.Digits(60);  %setup default precision to 34 decimal digits (quadruple)
mp.GuardDigits(10);
%colors={rgb('Red'),rgb('ForestGreen'),rgb('RoyalBlue'),rgb('Purple')};
%{
colors={rgb('Red'),rgb('Salmon'),rgb('LightCoral'),rgb('FireBrick'),rgb('DarkRed'),...
  rgb('Blue'),rgb('Turquoise'),rgb('SteelBlue'),rgb('SkyBlue'),rgb('DarkBlue'),...
  rgb('Green'),rgb('Teal'),rgb('SeaGreen'),rgb('ForestGreen'),rgb('DarkGreen'),...
  rgb('Brown'),rgb('BurlyWood'),rgb('SandyBrown'),rgb('SaddleBrown'),rgb('Maroon'),...
  rgb('HotPink'),rgb('LightPink'),rgb('DeepPink'),rgb('PaleVioletRed'),rgb('MediumVioletRed')};
%}
colors={rgb('Red'),rgb('Aqua'),rgb('LawnGreen'),rgb('Lavender'),rgb('LightPink'),...
  rgb('Gold'),rgb('Tan')};

%% System initialization:
fs=16000;   % sampling frequency <-> switching frequency --- 256000 Hz
Ts=1/fs;
z=tf('z', Ts);
L1=1.6e-3;
R1=2e-3;
Cx1=10e-6;
Rx1=0.1e-3;
L2=0.8e-3;
R2=1e-3;
Vdc=700;
cpk=2;   % PWM carrier peak-to-peak amplitude
ki = 29.5:0.1:51.5;
nki = numel(ki);

% Dd-settings:
wdw.Ymin=-0.10;
wdw.Ymax= 0.15;
wdw.Xmin=-0.05;
wdw.Xmax= 0.10;
R.r = 0.987; %0.9;
R.shft = 0.0; %0.09
%{
zgrid(0.05:0.1:0.85,[]); hold on;
ts = 0.02;
R_ts = exp(-3/ts/fs);
w = 1:1:2*pi*fs;
%plot(exp(1i*w/fs), 'r-', 'Linewidth', 1.5);
plot(R_ts*exp(1i*w/fs),'Color',rgb('springgreen'),'Linewidth',1.5);
plot(0,0,'.','Color',rgb('springgreen'),'MarkerSize',11);
plot(R.r*exp(1i*w/fs)+R.shft, 'b-.', 'Linewidth', 1);
plot(R.shft,0,'b.', 'MarkerSize', 11);
axis equal;
xlim([-1,1]); ylim([-1,1]);
%}

for k = 1:nki  
  
  % Controller:
  Kz_pr1=PR(Ts, 150, 250, 550, ki(k), '50Hz');  %Kp=0
  Kz_pr2=parallel(1, Kz_pr1);                   %Kp=1
  kic1=ss(0,0,0,0,Ts);  %Kic=0
  kic2=ss(0,0,0,1,Ts);  %Kic=1
  
  % Prepare I/O names and do basic connect():
  Kz_pr1.InputName='e';
  Kz_pr1.OutputName='u';
  Kz_pr2.InputName='e';
  Kz_pr2.OutputName='u';
  DELAY=1/z;
  DELAY.InputName='Vi';
  DELAY.OutputName='Vi_d';
  kic1.InputName='i_c';
  kic1.OutputName='kc';
  kic2.InputName='i_c';
  kic2.OutputName='kc';
  S1=sumblk('Vi=u-kc');
  S2=sumblk('e=i_ref-i_out');
  % Comments to building "DRAFTs":
  %{
  % Keep the same ordering of inputs-outputs in DRAFTs as in the plant model.
  % Clearly, the inputs of DRAFTs are the outputs of the model and vice versa.
  % The (single) reference should be always the first input.
  % Since the stability is only affected by the plant's feedback outputs,
  % only such outputs should be declared as inputs of DRAFTs after the
  % reference input.
  %}
  DRAFT1=connect(kic1, S1, S2, DELAY, Kz_pr1, {'i_ref', 'i_out', 'i_c'}, {'Vi_d'});
  DRAFT2=connect(kic1, S1, S2, DELAY, Kz_pr2, {'i_ref', 'i_out', 'i_c'}, {'Vi_d'});
  DRAFT3=connect(kic2, S1, S2, DELAY, Kz_pr1, {'i_ref', 'i_out', 'i_c'}, {'Vi_d'});
  
  % Separating and cutting the ss-matrices for the later closed-loop state-space:
  A1=DRAFT1.a; A2=DRAFT2.a; A3=DRAFT3.a;
  B1=DRAFT1.b; B2=DRAFT2.b; B3=DRAFT3.b;
  C1=DRAFT1.c; C2=DRAFT2.c; C3=DRAFT3.c;
  D1=DRAFT1.d; D2=DRAFT2.d; D3=DRAFT3.d;
  
  %% DD in the loop:
  regset=cell(7,1);     % stable regions set
  cpol.l=1;             % counter for stable polygons; for colors
  for Lgrid=[0, 5e-3]%1*[0, 1.25e-3, 2.5e-3, 3.75e-3, 5e-3]
    for Rgrid=[0, 10]%1*[0, 2.5, 5, 7.5, 10]
      A=[0, 1/Cx1, -1/Cx1; -1/L1, -(R1 + Rx1)/L1, Rx1/L1; 1/(L2 + Lgrid), Rx1/(L2 + Lgrid), -(R2 + Rgrid + Rx1)/(L2 + Lgrid)];
      B=[0; 1/L1; 0];
      C=[0, 0,  1; 0, 1, -1];
      D=[0; 0];
      B=B*Vdc/cpk;  % plant times "VSI gain", *2/sqrt(3) is due to SVM
      Az=expm(A*Ts);
      Bz=(expm(A*Ts)-eye(size(A)))/A*B;
      Cz=C;
      Dz=D;
      % 1st state-space:
      p=size(B1,2);
      INV=eye(p-1)-Dz*D1(:,2:p);
      a11=A1+B1(:,2:p)/INV*Dz*C1;
      a12=B1(:,2:p)/INV*Cz;
      a21=Bz*C1+Bz*D1(:,2:p)/INV*Dz*C1;
      a22=Az+Bz*D1(:,2:p)/INV*Cz;
      Aclp=[a11, a12; a21, a22];
      %{
      % The whole closed-loop state-space (only the A matrix is really needed here!):
      b11=B1(:,1)+B1(:,2:p)/INV*Dz*D1(:,1);
      b21=Bz*D1(:,1)+Bz*D1(:,2:p)/INV*Dz*D1(:,1);
      Bclp=[b11;b21];
      Cclp=INV\[Dz*C1, Cz];
      Cclp=Cclp(1,:);
      Dclp=INV\Dz*D1(:,1);
      Dclp=Dclp(1,:);
      %}
      DEN_ZERO=mp_charpoly(mp(Aclp));   % char poly of A is very precise (esp. with mp-toolbox)
      % 2nd state-space:
      p=size(B2,2);
      INV=eye(p-1)-Dz*D2(:,2:p);
      a11=A2+B2(:,2:p)/INV*Dz*C2;
      a12=B2(:,2:p)/INV*Cz;
      a21=Bz*C2+Bz*D2(:,2:p)/INV*Dz*C2;
      a22=Az+Bz*D2(:,2:p)/INV*Cz;
      Aclp=[a11, a12; a21, a22];
      DEN_Kp=mp_charpoly(mp(Aclp));   % char poly of A is very precise (esp. with mp-toolbox)
      % 3rd state-space:
      p=size(B3,2);
      INV=eye(p-1)-Dz*D3(:,2:p);
      a11=A3+B3(:,2:p)/INV*Dz*C3;
      a12=B3(:,2:p)/INV*Cz;
      a21=Bz*C3+Bz*D3(:,2:p)/INV*Dz*C3;
      a22=Az+Bz*D3(:,2:p)/INV*Cz;
      Aclp=[a11, a12; a21, a22];
      DEN_Kic=mp_charpoly(mp(Aclp));
      
      cpol.P=mp_sub(DEN_Kp,DEN_ZERO);    % or P=P+DEN_ZERO-DEN_ZERO, since DEN_Kp=P+DEN_ZERO
      cpol.Q=mp_sub(DEN_Kic,DEN_ZERO);   % or Q=Q+DEN_ZERO-DEN_ZERO, since DEN_Kic=Q+DEN_ZERO
      cpol.L=DEN_ZERO;
      
      % Generate and plot Dd from DEN:
      regs=Dd(cpol, R, wdw, colors);
      if ~isempty(regs)
        regset(cpol.l)={regs};
        cpol.l=cpol.l+1;
      else
        regset=[];
        break;
      end
      % Brutal force (for tests and comparison):
      %{
      cpol.P = double(cpol.P);
      cpol.Q = double(cpol.Q);
      cpol.L = double(cpol.L);
      brutforce(cpol, R, wdw, colors);
      close;  % one close() for faces/regs
      close;  % another close() - for brutforce()
      %}
    end
    if isempty(regset)
      break;
    end
  end
  regset(cpol.l : length(regset))=[];
  
  % Labels and other:
  %{
  xlabel('k_{ic}, [-]', 'FontSize', 16,'FontName','Times','FontAngle', 'normal');    % the labels are not automatized!
  ylabel('k_{p}, [-]', 'FontSize', 16, 'FontName','Times','FontAngle', 'normal');
  ylim([kp_min kp_max]);
  xlim([kic_min kic_max]);
  set(gca,'FontSize',12)
  grid on;
  %}
  
  %% Clip polygons:
  if length(regset)>1
    P3=PolygonClip(regset{1},regset{2});
    for i=3:1:length(regset)
      if ~isempty(P3)
        P3=PolygonClip(P3,regset{i},1);
      else
        break;
      end
    end
    if ~isempty(P3)
      hh = fill3([P3.x], [P3.y], (ki(k)*ones(numel(P3.x),1)), [0.3,0.3,0.3]);  %rgb('ForestGreen')  k/nki*
      set(hh,'edgecolor','r');
      hold on;
    else
      disp('(GCC_DP_poly.m): No robustly stable region for the given uncertainty');
    end
  else
    disp('(GCC_DP_poly.m): No robustly stable region for the given uncertainty');
  end
end
