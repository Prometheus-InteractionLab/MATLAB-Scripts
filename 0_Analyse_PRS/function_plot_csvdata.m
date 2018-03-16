function function_plot_csvdata(csvdata,i1,i2)

% RPM
subplot(2,2,1)
plot(csvdata.rpm,'linewidth',2),hold on
plot([i1 i2],csvdata.rpm([i1 i2]),'o','markersize',10,'linewidth',2,'markerfacecolor','w')
hold on
set(gca,'box','on','linewidth',1,'fontsize',12,'fontname','eras light itc')
xlabel('t [s]')
ylabel('Speed [rpm]')
axis([0 i2+30 1100 1300])
title(csvdata.name_rpm)

% Qair
subplot(2,2,2)
plot(csvdata.qair,'linewidth',2),hold on
plot([i1 i2],csvdata.qair([i1 i2]),'o','markersize',10,'linewidth',2,'markerfacecolor','w')
hold on
set(gca,'box','on','linewidth',1,'fontsize',12,'fontname','eras light itc')
xlabel('t [s]')
ylabel('Qair [g/s]')
title(csvdata.name_qair)

% Pintake
% subplot(2,3,3)
% plot(csvdata.pintake/100,'linewidth',2),hold on
% plot([i1 i2],csvdata.pintake([i1 i2])/100,'o','markersize',10,'linewidth',2,'markerfacecolor','w')
% hold on
% set(gca,'box','on','linewidth',1,'fontsize',12,'fontname','eras light itc')
% xlabel('t [s]')
% ylabel('Pintake [bar]')
% title(csvdata.name_pintake)

% Tintake
subplot(2,2,3)
plot(csvdata.tintake,'linewidth',2),hold on
plot([i1 i2],csvdata.tintake([i1 i2]),'o','markersize',10,'linewidth',2,'markerfacecolor','w')
hold on
set(gca,'box','on','linewidth',1,'fontsize',12,'fontname','eras light itc')
xlabel('t [s]')
ylabel('Tintake [degC]')
title(csvdata.name_tintake)

hfg = annotation('textbox',[0.8513 0.7716+.166 0.0923 0.0528],'color',[0 128 64]/255,'edgecolor',[0 128 64]/255,'fontname','eras light ITC','String',...
        {['pressure data '];...
        ['  - ' csvdata.name_dayname csvdata.name_run '.prs' ]; });
    
% Tliner
subplot(2,2,4)
plot(csvdata.tliner,'linewidth',2),hold on
plot([i1 i2],csvdata.tliner([i1 i2]),'o','markersize',10,'linewidth',2,'markerfacecolor','w')
hold on
set(gca,'box','on','linewidth',1,'fontsize',12,'fontname','eras light itc')
xlabel('t [s]')
ylabel('Tliner [degC]')
axis([0 i2+30 85 100])
title(csvdata.name_tliner)

% Pfuel
% subplot(2,3,6)
% plot(csvdata.pfuel,'linewidth',2),hold on
% plot([i1 i2],csvdata.pfuel([i1 i2]),'o','markersize',10,'linewidth',2,'markerfacecolor','w')
% hold on
% set(gca,'box','on','linewidth',1,'fontsize',12,'fontname','eras light itc')
% xlabel('t [s]')
% ylabel('Pfuel [bar]')
% title(csvdata.name_pfuel)

