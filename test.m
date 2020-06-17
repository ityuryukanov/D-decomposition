clear;
clc;
addpath('PolygonClipper', 'inpoly', 'rintersect', 'intersections', 'p_poly_dist');

mp.Digits(34);   % quad-precision (advanpix-toolbox)
mp.GuardDigits(0);   % quad-precision

colors={'r', 'g', 'b', 'm', 'c', rgb('pale green'), ...
    rgb('slightly dark green'), rgb('slightly dark red'), rgb('dark gray'), ...
    'y', rgb('orange'), rgb('navy'), rgb('gray'), rgb('olive'), rgb('sky'),  ...
    rgb('u4'), rgb('teal'), rgb('u5'), rgb('maroon'), rgb('slateblue'), 'b', ...
    rgb('slightly normal purple'), rgb('u7'), rgb('u8'), rgb('u2'), 'r', 'g', 'b'};

wdw.Ymin = -100;
wdw.Ymax = 100;
wdw.Xmin = -100;
wdw.Xmax = 100;
charpoly.P = mp(1*rand(6,1));
charpoly.Q = mp(2*rand(6,1));
charpoly.L = mp(1*rand(6,1));
charpoly.n = 1;

R.r=1;  % use as few decimal places as possible 
R.shft=0;  % use as few decimal places as possible
regs=Dd(poly, R, wdw, colors);
xlabel('k_{Q}, [-]', 'FontSize', 12);  % the labels are not automatized!
ylabel('k_{P}, [-]', 'FontSize', 12);  % the labels are not automatized!
ylim([kp_min kp_max]);
xlim([kic_min kic_max]);
set(gca,'FontSize',10);

% Brutal force (for tests and comparison):

P=double(P);
Q=double(Q);
L=double(L);
brutforce(P, Q, L, R, kp_min, kp_max, kic_min, kic_max, colors)
xlabel('k_{Q}, [-]', 'FontSize', 12);  % the labels are not automatized!
ylabel('k_{P}, [-]', 'FontSize', 12);  % the labels are not automatized!
ylim([kp_min kp_max]);
xlim([kic_min kic_max]);
set(gca, 'FontSize', 10);
close;  % one close() for faces/regs
close;  % another close() - for brutforce()
close; 

%{
tests:
% 1) deg(L)>deg(P)>deg(Q) --> OK
% 2) deg(L)>deg(Q)>deg(P) --> OK
% 3) deg(P)>deg(Q)>deg(L) --> OK
% 4) deg(Q)>deg(P)>deg(L) --> OK
% 5) deg(P)>deg(L)>deg(Q) --> OK
% 6) deg(Q)>deg(L)>deg(P) --> OK

kp_min = -3e-16;   %  0.012;
kp_max = 3e-16;   %  0.014;
kic_min = -4.5e-16;  % -0.042;
kic_max = 2.5e-16;  % -0.041;

% P = [1; 0; 1; 0; -4];
% Q = [1; 0; 1; 0; 1; 0; -5];
% L = [1; 0; 1; 0; 10];
% % P=1*rand(1,1);
% % Q=1*rand(1,1);
% % L=5*rand(14,1);
% P = [7.951999011370631; 1.868726045543786; 4.897643957882311; 4.455862007108995; 6.463130101112647]; % 12.4323423; 233.432432; 2.32234; -324.2342342; 11.3424234
% Q = 70.936483085807254;
% L = 1.509373363964722;
% P=mp(P); Q=mp(Q); L=mp(L);
% % P = [1.000000000000000; -1.037172848460050; 3.403587876230064; -3.401538905458140; 3.410588108548355; -1.096171576710489; 1.074295006458073]; 
% % Q = [1.000000000000000; -5.443235596952732; 12.345049452573706; -15.112514153301246; 11.028372344299971; -6.249690422937162; 7.010587463803755;... 
% % -12.409537613534596; 16.626122086196279; -15.291419897851032; 9.404406581374690; -3.496299902516355; 0.588159658844725];
% % L = 1.0e+03 * [0.001000000000000; 0.016457367980160; -0.190245955764790; 0.852237501324827; -2.279366378617001; 4.219175999687946; -5.880905372414929;... 
% % 6.488510980790730; -5.605969434951736; 3.182089458228911; 0.092378109388990; -2.793679759203246; 3.926343205055470; -3.729469834970432; 2.750093530861745;... 
% % -0.970421389306336; -1.323697080598691; 2.810069510159613; -2.401803097663168; 0.615664242889477; 1.261775925825366; -2.722288743572353; 3.819739918579280;... 
% % -4.170249886614150; 3.355783290545447; -1.860711397332138; 0.647467280129348; -0.110570896263659; -0.002287987767273; 0.002880893592614]; 
% % P=mp(P); Q=mp(Q); L=mp(L);
% % L=[roots(L); 0.999999999+0.00001i; 0.999999999-0.00001i; 0.9165151389911+0.4i; 0.9165151389911-0.4i; 0.9165151389911+0.4i; 0.9165151389911-0.4i;... 
% % 0.9165151389911+0.4i; 0.9165151389911-0.4i; 0.9165151389911+0.4i; 0.9165151389911-0.4i]; %; 1+250i; 1-250i; 1+800i; 1-800i
% % L=real(poly(L)); L=L.';
% % Q=[roots(Q); 0.9999999+0.00001i; 0.9999999-0.00001i; 0.99999+0.00001i; 0.99999-0.00001i; 0.99999+0.00001i; 0.99999-0.00001i];
% % Q=real(poly(Q)); Q=Q.';
% % P=[roots(P); 0.7141428428542+0.7i; 0.7141428428542-0.7i];
% % P=real(poly(P)); P=P.';
%}