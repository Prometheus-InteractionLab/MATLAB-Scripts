clear all
close all

[filename,pathname] = uigetfile('*.mat','Select file','Multiselect','on');
addpath(pathname)
filename = cellstr(filename);
dayname = filename{1}(1:6);


for ii=1:length(filename)
    close all % if plotting many runs, supress output
    load(filename{ii})
    runname = filename{ii}(end-4);

Pavg =  mean(rundata.pcyl(1:2:rundata.NC,:));
Pfiltavg = mean(rundata.pfilt(1:2:rundata.NC,:));
AHRRavg = mean(rundata.dQdCA(1:2:rundata.NC,:));
for i=2:2:rundata.NC
    % pressure difference
    PD(i,:) = rundata.pcyl(i,:)-Pavg;
    Pfil(i,:) = rundata.pfilt(i,:)-Pfiltavg;
    % central diff on pressure difference
%     CD(i,:) = FDderiv1(rundata.ca(1,:),rundata.pcyl(i,:)-Pavg,0);
    AHRRdiff(i,:) = rundata.dQdCA(i,:)-AHRRavg;
end
normPD = PD./max(max(PD));
% normCD = CD./max(max(CD));
normdQ = rundata.dQdCA./max(max(rundata.dQdCA));
normPfil = Pfil./max(max(Pfil));

figure(3); 
subplot(2,1,1); hold on; grid on;
plot(rundata.ca(1,:),rundata.pcyl(2:2:rundata.NC,:),'color',[.7 .7 .7])
plot(rundata.ca(1,:),rundata.pfilt(2:2:rundata.NC,:),'color',[.7 .4 .4])
plot(rundata.ca(1,:),mean(rundata.pcyl(2:2:rundata.NC,:)),'k','linewidth',2)
plot(rundata.ca(1,:),mean(rundata.pfilt(2:2:rundata.NC,:)),'r','linewidth',2)
plot([rundata.mean_SOC rundata.mean_SOC], [0 100])
axis([rundata.sse_degV rundata.mean_SOC+15 min(rundata.pmotor)-10 max(rundata.pmax)+5])
ylabel('Cylinder Pressure (bar)')
subplot(2,1,2); hold on; grid on;
plot(rundata.ca(1,:),rundata.dQdCA(2:2:rundata.NC,:),'color',[.7 .4 .4])
plot(rundata.ca(1,:),mean(rundata.dQdCA(2:2:rundata.NC,:)),'r','linewidth',2)
plot(rundata.ca(1,:),rundata.dQdCAavg,'k')
plot([rundata.mean_SOC rundata.mean_SOC], [-100 100])
plot([rundata.sse_degV+rundata.mean_inj_delay rundata.sse_degV+rundata.mean_inj_delay+rundata.dse_ms/.140], [-0.09 -0.09],'b','linewidth',2)
xlim([rundata.sse_degV rundata.mean_SOC+15])
ylim([-30 300])
xlabel('CAD')
ylabel('AHRR (J/deg)')

pfstd = std(rundata.pcyl(2:2:rundata.NC,:));
pmstd = std(rundata.pcyl(1:2:rundata.NC,:));

figure(); hold on;
plot(rundata.ca(1,:),pfstd)
plot(rundata.ca(1,:),pmstd,'r')
igdet = pfstd-pmstd;

saveas(gcf,fullfile(pathname,['compare_filter_pdQ' dayname runname '.fig']),'fig')
saveas(gcf,fullfile(pathname,['compare_filter_pdQ' dayname runname '.emf']),'emf') 

%  figure(2); hold on; grid on;
%  plot(rundata.ca(1,:),mean(Pfil),'r','linewidth',2)
%  plot(rundata.ca(1,:),mean(PD),'color',[.6 .6 .6])
%  axis([-15 60 -.2 1.2])

% 
% figure(2); hold on; grid on;
% plot(rundata.ca(2:2:rundata.NC,:),normCD(2:2:rundata.NC,:), 'color', [.6 .6 .6])
% plot(rundata.ca(1,:),mean(normCD(2:2:rundata.NC,:)),'b','linewidth',2)
% axis([-15 60 -.2 1.2])
% 
%  figure(3); hold on; grid on;
%  plot(rundata.ca(2:2:rundata.NC,:),normdQ(2:2:rundata.NC,:),'color',[.6 .6 .6])
%  plot(rundata.ca(1,:),mean(normdQ(2:2:rundata.NC,:)),'k','linewidth',2)
%  axis([-15 60 -.2 1.2])

figure(1); hold on; grid on

subplot(3,2,1); hold on; grid on;
plot(rundata.ca(1,:),mean(PD(2:2:rundata.NC,:)),'k','linewidth',2)
plot(rundata.ca(1,:),min(PD(2:2:rundata.NC,:)),'color',[.7 .7 .7],'linewidth',1.5)
plot(rundata.ca(1,:),max(PD(2:2:rundata.NC,:)),'color',[.7 .7 .7],'linewidth',1.5)
plot(rundata.ca(1,:),mean(Pfil(2:2:rundata.NC,:)),'color','r','linewidth',2)
plot(rundata.ca(1,:),min(Pfil(2:2:rundata.NC,:)),'color',[.7 .3 .3],'linewidth',1.5)
plot(rundata.ca(1,:),max(Pfil(2:2:rundata.NC,:)),'color',[.7 .3 .3],'linewidth',1.5)
plot(rundata.ca(1,:),igdet,'b--','markersize',2)
plot([rundata.sse_degV+rundata.mean_inj_delay rundata.sse_degV+rundata.mean_inj_delay+rundata.dse_ms/.140], [-1 -1],'b','linewidth',2)
plot(rundata.mean_SOC,0,'kx','markersize',5,'linewidth',2)
plot(rundata.mean_SOC,0,'ko','markersize',10,'linewidth',2)
xlim([rundata.sse_degV rundata.mean_SOC+25])
ylabel('Fired-Motored (bar)')
subplot(3,2,2); hold on; grid on;
plot(rundata.ca(1,:),mean(PD(2:2:rundata.NC,:)),'k','linewidth',2)
plot(rundata.ca(1,:),min(PD(2:2:rundata.NC,:)),'color',[.7 .7 .7],'linewidth',1.5)
plot(rundata.ca(1,:),max(PD(2:2:rundata.NC,:)),'color',[.7 .7 .7],'linewidth',1.5)
plot(rundata.ca(1,:),mean(Pfil(2:2:rundata.NC,:)),'color','r','linewidth',2)
plot(rundata.ca(1,:),min(Pfil(2:2:rundata.NC,:)),'color',[.7 .3 .3],'linewidth',1.5)
plot(rundata.ca(1,:),max(Pfil(2:2:rundata.NC,:)),'color',[.7 .3 .3],'linewidth',1.5)
plot([rundata.sse_degV+rundata.mean_inj_delay rundata.sse_degV+rundata.mean_inj_delay+rundata.dse_ms/.140], [-2 -2],'b','linewidth',2)
plot(rundata.ca(1,:),igdet,'b--','markersize',2)
plot(rundata.mean_SOC,0,'kx','markersize',5,'linewidth',2)
plot(rundata.mean_SOC,0,'ko','markersize',10,'linewidth',2)
xlim([rundata.sse_degV rundata.mean_SOC+3])
ylim([-3 1.5])

subplot(3,2,3); hold on; grid on;
plot(rundata.ca(1,:),mean(AHRRdiff(2:2:rundata.NC,:)),'k','linewidth',2)
plot(rundata.ca(1,:),min(AHRRdiff(2:2:rundata.NC,:)),'k','linewidth',2)
plot(rundata.ca(1,:),max(AHRRdiff(2:2:rundata.NC,:)),'k','linewidth',2)
plot(rundata.ca(1,:),mean(rundata.dQdCA(1:2:rundata.NC,:)),'.k','linewidth',2)
plot(rundata.ca(1,:),min(rundata.dQdCA(1:2:rundata.NC,:)),'color',[.4 .4 .4],'linewidth',1.5)
plot(rundata.ca(1,:),max(rundata.dQdCA(1:2:rundata.NC,:)),'color',[.4 .4 .4],'linewidth',1.5)


plot(rundata.mean_SOC,0,'rx','markersize',5,'linewidth',2)
plot(rundata.mean_SOC,0,'ro','markersize',10,'linewidth',2)
xlim([rundata.sse_degV rundata.mean_SOC+25])
ylabel('AHRR (J/deg)')

subplot(3,2,4); hold on; grid on;
plot(rundata.ca(1,:),mean(AHRRdiff(2:2:rundata.NC,:)),'k','linewidth',2)
plot(rundata.ca(1,:),min(AHRRdiff(2:2:rundata.NC,:)),'k','linewidth',2)
plot(rundata.ca(1,:),max(AHRRdiff(2:2:rundata.NC,:)),'k','linewidth',2)
plot(rundata.ca(1,:),mean(rundata.dQdCA(1:2:rundata.NC,:)),'.k','linewidth',2)
plot(rundata.ca(1,:),min(rundata.dQdCA(1:2:rundata.NC,:)),'color',[.4 .4 .4],'linewidth',1.5)
plot(rundata.ca(1,:),max(rundata.dQdCA(1:2:rundata.NC,:)),'color',[.4 .4 .4],'linewidth',1.5)
plot(rundata.ca(1,:),std(rundata.pcyl(2:2:rundata.NC,:)),'r--','markersize',2)
xlim([rundata.sse_degV rundata.mean_SOC+3])
plot(rundata.mean_SOC,0,'rx','markersize',5,'linewidth',2)
plot(rundata.mean_SOC,0,'ro','markersize',10,'linewidth',2)

ylim([-50 25])

subplot(3,2,5); hold on; grid on;
plot(rundata.ca(1,:),mean(CD(2:2:rundata.NC,:)),'b','linewidth',2)
plot(rundata.ca(1,:),min(CD(2:2:rundata.NC,:)),'color',[.4 .4 .8],'linewidth',1.5)
plot(rundata.ca(1,:),max(CD(2:2:rundata.NC,:)),'color',[.4 .4 .8],'linewidth',1.5)
plot(rundata.mean_SOC,0,'bx','markersize',5,'linewidth',2)
plot(rundata.mean_SOC,0,'bo','markersize',10,'linewidth',2)
xlim([rundata.sse_degV rundata.mean_SOC+25])
ylabel({'dP/dCAD (bar/deg)'; 'fired-motored'})
xlabel('CAD')

subplot(3,2,6); hold on; grid on;
% plot(rundata.ca(1,:),mean(CD(2:2:rundata.NC,:)),'b','linewidth',2)
% plot(rundata.ca(1,:),min(CD(2:2:rundata.NC,:)),'color',[.4 .4 .8],'linewidth',1.5)
% plot(rundata.ca(1,:),max(CD(2:2:rundata.NC,:)),'color',[.4 .4 .8],'linewidth',1.5)
plot(rundata.mean_SOC,0,'bx','markersize',5,'linewidth',2)
plot(rundata.mean_SOC,0,'bo','markersize',10,'linewidth',2)
xlim([rundata.sse_degV rundata.mean_SOC+3])
ylim([-5 2.5])
xlabel('CAD')

saveas(gcf,fullfile(pathname,['compare_SOC' dayname runname '.fig']),'fig')
saveas(gcf,fullfile(pathname,['compare_SOC' dayname runname '.emf']),'emf') 
end
