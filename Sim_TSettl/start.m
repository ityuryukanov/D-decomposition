% % clc;
% % clear;
opt=bodeoptions;
opt.FreqUnits='Hz';
opt.PhaseWrapping='on';
opt.PhaseMatching='on';
s=tf('s');

% Simulation parameters:
f_sw=8000; % [Hz], switching frequency
%{
Effective doubling of the switching frequency by using 
a three-level converter is possible: e.g. 5000 Hz -> 10000 Hz
%}
f_sample=2*f_sw; % [Hz], sample frequency for control
%{
Douple update mode is essential for provision of higher control bandwidth
by increasing the Nyquist frequency of the system and reducing the total
delay (computations+DPWM). However, the DSP should be fast enough.
%}
Ts_control=1/f_sample;  % [s]
f_n=50;                 % [Hz]

% Sources' parameters:
Vgrid=127*sqrt(2);      % rated [V], ampl, PH-GND: grid voltage
V_dc_ref=700;           % [V] - DC-link voltage
cpk=2;                  % PWM carrier peak-to-peak amplitude

% LCL filter parameters:
L_f1=1.6e-3;
C_x1=10e-6;
L_f2=0.8e-3;
f_res0=sqrt((L_f1+L_f2)/(L_f1*L_f2*C_x1))/(2*pi);
f_res_lim=sqrt(1/(L_f1*C_x1))/(2*pi);
R_x1=0.1e-3;
R_f1=2e-3;
R_f2=1e-3;

% Grid impedance:
Lgrid=5e-3;%
Rgrid=10;

% Harmonic content of the grid:
Vg_h5=0*Vgrid;  % 0.05
Vg_h7=0*Vgrid;  % 0.03
Vg_h12=0*Vgrid; % 0.01
Vg_h15=0*Vgrid; % 0.005

% Controller parameters/references:
kic=0.042;  
kp=0.049;   
V_ffwd=-1;  % V_ffwd=1 - SIMULATE with Vpcc feedforward, V_ffwd=-1 - without (Bodeplots are "without" due to connect)
Id_ref=15;  % [A] - is only used if no DC-link voltage control is applied!
Iq_ref=0;   % [A]
% PR Current controller:
Ks_curr=current_controller(V_dc_ref, kp, kic, Ts_control, cpk, f_n, R_f1, L_f1, C_x1, R_x1, L_f2, ...
    R_f2, Rgrid, Lgrid, opt);
NUM_cc=Ks_curr.num{1,1};
DEN_cc=Ks_curr.den{1,1};



