%--------------------------------------------------------------------------
% For continuously changing system parameters, system eigenvalues also change
% continuously (i.e. they don't "jump"). Therefore, the eigenvalue test for
% plant uncertainty allows to detect possible isolated instabilities*, which
% makes it very reliable.
%
% * i.e a system's instability region in the parameter space of uncertain
% parameters is enclosed by the stability region, and the size of the
% unstable region is  arbitrarily small.. up to one isolated unstable point!
%
% It is improved version of the code that is present in the attached CD
% under the same name (the speed was increased by using "fast connect"). It
% still could be a bit faster...
%
% It also seems to be the only possibility to obtain stable region in the
% plane of uncertain PHYSICAL parameters for sampled-data control systems.
%--------------------------------------------------------------------------

clear;
% % clc;
tol=1e-12;   % tolerance (for comparisons etc.)
s=tf('s');
fs=16000;   % sampling frequency <=> switching frequency
Ts=1/fs;
z=tf('z', Ts);
N=52;   % number of Lg-values
COLOR_STEP=0.8/N;   % grayscale color step
%       0    1    2    3    4    5    6
color=['r', 'g', 'b', 'm', 'y', 'c', 'k'];   % unstable root colors

%% System parameters:
% Plant:
L1=1.6e-3;
R1=2e-3;
Cx1=10e-6;
Rx1=0.1e-3;
L2=0.8e-3;
R2=1e-3;
Vdc=700;

% Plant uncertainty to check. Please specify here either the uncertainty
% box to check the controller's robustness, or a range of Lg-Rg for which
% the stable region in the Lg-Rg plane should be plotted.
Lg_min=0;
Rg_min=0;
Lg_max=5e-3;
Rg_max=10;
L=linspace(Lg_min+tol, Lg_max, N);
R=linspace(Rg_min+tol, Rg_max, round(0.75*N));

% Testing requirement:
r=0.987; % circle radius
shift=0.00; % circle shift (nonzero only for "damping&settling_time" specifications)

% Controller:
kic=0.042; % 0.07
kp=0.049; % 0.1025
Kz_pr=parallel(kp, PR(Ts, 150, 250, 350, 30, '50Hz'));
kic=ss(0,0,0,kic,Ts);

% Here we obtain the results fast because of using fast dynamic system
% interconnection that was proposed and described in the masters thesis
% Prepare I/O names and do basic connect():
Kz_pr.InputName='e';
Kz_pr.OutputName='u';
DELAY=1/z;
DELAY.InputName='Vi';
DELAY.OutputName='Vi_d';
kic.InputName='i_c';
kic.OutputName='kc';
S1=sumblk('Vi=u-kc');
S2=sumblk('e=i_ref-i_out');
% Comments to building "DRAFTs":
%{
% Keep the same ordering of inputs-outputs in DRAFTs as in the plant model.
% Clearly, the inputs of DRAFTS are the outputs of the model and vice versa.
% The (single) reference should be always the first input.
% Since the stability is only affected by the plant's feedback outputs,
% only such outputs should be declared as inputs of DRAFTs after the
% reference input.
%}
DRAFT1=connect(kic, S1, S2, DELAY, Kz_pr, {'i_ref', 'i_out', 'i_c'}, {'Vi_d'});
A1=DRAFT1.a;
B1=DRAFT1.b;
C1=DRAFT1.c;
D1=DRAFT1.d;

%% Eigenvalues-in-the-loop:
tic;
Lg=Lg_min;
Rg=Rg_min;
warning('off','Control:combination:connect9');
tic;
%figure;
plot(0,0,'.');
h_eig=gca;
title(['Closed-loop eigenvalues for physical parameter variation L_g=[0 .. ', num2str(Lg_max),'], R_g=[0 ..', num2str(Rg_max),']']);
figure;
plot(0,0,'.');
h_par=gca;
for Rg=R
  for Lg=L
    % Continuous time state-space:
    A=[0, 1/Cx1, -1/Cx1; -1/L1, -(R1 + Rx1)/L1, Rx1/L1; 1/(L2 + Lg), Rx1/(L2 + Lg), -(R2 + Rg + Rx1)/(L2 + Lg)];
    B=[0; 1/L1; 0]*Vdc/2;   % plant times "PWM gain" *Vdc/sqrt(3)
    C=[0, 0,  1; 0, 1, -1];
    D=[0; 0];
    Az=expm(A*Ts);
    Bz=(expm(A*Ts)-eye(size(A)))/A*B;
    Cz=C;
    Dz=D;
    % Closed-loop state-space:
    p=size(B1,2);
    INV=eye(p-1)-Dz*D1(:,2:p);
    a11=A1+B1(:,2:p)/INV*Dz*C1;
    a12=B1(:,2:p)/INV*Cz;
    a21=Bz*C1+Bz*D1(:,2:p)/INV*Dz*C1;
    a22=Az+Bz*D1(:,2:p)/INV*Cz;
    Aclp=[a11, a12; a21, a22];
    EIG=eig(Aclp);
    % the following plotting can be commented out, if not needed -> with the resulting speed increase!
    k=find(L==Lg);
    plot(h_eig, EIG, '.', 'MarkerSize', 10, 'Color', [COLOR_STEP*k COLOR_STEP*k COLOR_STEP*k]);
    hold(h_eig, 'on');
    
    % the following test can be commented out, if not needed -> with the resulting speed increase!
    Np=sum(abs(EIG-shift)>r);
    switch Np
      case 0
        plot(h_par, Rg, Lg, '.', 'MarkerSize', 7, 'Color', color(1));
      case 1
        plot(h_par, Rg, Lg, '.', 'MarkerSize', 7, 'Color', color(2));
      case 2
        plot(h_par, Rg, Lg, '.', 'MarkerSize', 7, 'Color', color(3));
      case 3
        plot(h_par, Rg, Lg, '.', 'MarkerSize', 7, 'Color', color(4));
      case 4
        plot(h_par, Rg, Lg, '.', 'MarkerSize', 7, 'Color', color(5));
      case 5
        plot(h_par, Rg, Lg, '.', 'MarkerSize', 7, 'Color', color(6));
      case 6
        plot(h_par, Rg, Lg, '.', 'MarkerSize', 7, 'Color', color(7));
    end
    hold(h_par, 'on');
  end
end
axes(h_eig)
zgrid;
axes(h_par);
grid on;
toc;