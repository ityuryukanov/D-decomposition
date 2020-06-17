%% Plot grid current

data = Grid_Curr_DQ0;

t0 = 0.04;
ts = 14.3 + t0*1000;
m = numel(data.time);
x_width = 1.2*241; y_width = 1.25*161;

plot(data.time*1000,data.signals.values(:,1),'k-','Linewidth',1);
hold on;
plot(data.time*1000,data.signals.values(:,2),'r-','Linewidth',1);
plot(data.time*1000,data.signals.values(:,3),'b-','Linewidth',1);
plot(data.time*1000, Id_ref*ones(1,m),'k--','Linewidth',0.5);
%plot(data.time*1000, 1.05*Id_ref*ones(1,m),'k--','Linewidth',0.5);
plot(data.time*1000, 0.95*Id_ref*ones(1,m),'k--','Linewidth',0.5);

[ylims] = get(gca,'ylim');
%%%plot([ts,ts],[-100,100],'k-','Linewidth',0.25);
ylim([ylims(1),ylims(2)]);
set(gca,'xgrid','on');
legend('i_{g,d,1}', 'i_{g,d,2}', 'i_{g,d,0}');
%%%legend('i_{g,A}', 'i_{g,B}', 'i_{g,C}', 'Orientation', 'horizontal');
%}
%{
plot(data.time*1000,-Id_ref*ones(1,m),'k--','Linewidth',0.5);
plot(data.time*1000,-1.05*Id_ref*ones(1,m),'k--','Linewidth',0.5);
plot(data.time*1000,-0.95*Id_ref*ones(1,m),'k--','Linewidth',0.5);
[ylims] = get(gca,'ylim');
ylim([-5,ylims(2)]);
set(gca,'xgrid','on');
%}

xlim([t0*1000-2.5,data.time(end)*1000]);
xlabel('time, [ms]','interpreter','latex','fontsize',12);
ylabel('$i_{g}, [A]$','interpreter','latex','fontsize',13);
set(gcf, 'pos', [100 100 x_width y_width]);
