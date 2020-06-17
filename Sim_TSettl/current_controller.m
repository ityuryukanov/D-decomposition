function Ks_curr=current_controller(V_dc_ref, kp, kic, Ts_control, cpk, f_n, R1, L1, Cx1, Rx1, ... 
    L2, R2, Rgrid, Lgrid, opt)
A=[0, 1/Cx1, -1/Cx1; -1/L1, -(R1 + Rx1)/L1, Rx1/L1; 1/(L2 + Lgrid), Rx1/(L2 + Lgrid), -(R2 + Rgrid + Rx1)/(L2 + Lgrid)];
B=[0, 0; 1/L1, 0; 0, -1/(L2 + Lgrid)];
C=[0, 0,  1; 0, 1, -1];
D=[0, 0; 0, 0];
B(:,1)=B(:,1)*V_dc_ref/cpk;  % plant times "VSI gain", 2/sqrt(3) is due to SVM
G=ss(A, B, C, D);
w0=2*pi*f_n;

%% PR controller (possibly with multiple resonances)
% tracking controller (non-ideal PR-controller):
wc0=3;    % wc also influences the closed-loop gain and phase around 50 Hz 
ki0=30;   % this gain needs to be quite high for a good time domain performance at 50 Hz
b2= sin(w0*Ts_control);
b1= 0;
b0=-sin(w0*Ts_control);
a2= 2*w0+2*wc0*sin(w0*Ts_control);
a1=-4*w0*cos(w0*Ts_control);
a0= 2*w0-2*wc0*sin(Ts_control*w0);
% Total controller (tracking + checking):
Ks_curr=kp+ki0*tf([b2 b1 b0], [a2 a1 a0], Ts_control);
Ks_curr.InputName='e';
Ks_curr.OutputName='u';

%% Open-loop transfer function, in discrete time:
z=tf('z', Ts_control);
Gs=c2d(G, Ts_control, 'zoh');
DELAY=1/z;
DELAY.InputName='Vi';
DELAY.OutputName='Vi_d';
Gs.InputName={'Vi_t', 'Vg'};
Gs.OutputName={'i_out', 'i_c'};
kic=ss(0,0,0,kic,Ts_control);
kic.InputName='i_c';
kic.OutputName='kc';
S1=sumblk('Vi=u-kc');
S2=sumblk('Vi_t=Vi_d+Vi_n');            % Vi_n stands for "noise"
Gs_ol=connect(Gs, kic, S1, S2, DELAY, Ks_curr, {'Vg','e','Vi_n'}, 'i_out');
Ls=Gs_ol(:,2);
EASY_NYQUIST_APPLICABLE=isstable(Ls)    % i.e. is the open-loop system stable? -> Do GM&PM from margin(Ls) make sense?
% figure;
% margin(Ls);                           % stability margins
bode(Ls, opt); grid on; hold on;
% %------------------------------------------------------------------------
S1=sumblk('e=i_ref-i_out');
CLP=connect(Ls, S1, 'i_ref', 'i_out');    % closed loop poles (check up)
isstable(CLP)
% figure;
% pzmap(CLP); zgrid;                      % closed loop poles (check up)
% figure;
% step(CLP);
% figure;
% t=0:Ts_control:400e-3;
% u=sin(2*pi*50*t);
% lsim(CLP,u,t);
% %------------------------------------------------------------------------
% figure;
% bode(CLP, {2*pi*0.1 2*pi*100}, 'c', opt);    % to see the gain and phase at 50 Hz
% figure;
% bode(CLP, 'c', opt); grid on;                % reference tracking
%--------------------------------------------------------------------------
S1=sumblk('e=-i_out');
Gs_dist=connect(Gs_ol, S1, 'Vg', 'i_out');  % OUTPUT disturbance rejection
% figure;
% bodemag(Gs_dist, 'b', opt); grid on;            % OUTPUT disturbance rejection
%--------------------------------------------------------------------------
% S1=sumblk('e=-i_out');
% Gs_noise=connect(Gs_ol, S1, 'Vi_n', 'i_out');  % INPUT disturbance rejection
% figure;
% bode(Gs_noise, opt); grid on;  % INPUT disturbance rejection
%--------------------------------------------------------------------------
% close;
% close;
% close;
% close;
% close;
% close;
% close;
end


