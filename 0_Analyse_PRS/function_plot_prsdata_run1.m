function function_plot_prsdata_run(prsdata)

ienergize = find(prsdata.injComm(2,:)>1);
NC = prsdata.NC;

    set(0,'defaultAxesFontName','Arial')
    set(0,'defaultTextFontName','Arial')
    
% Pcyl
subplot(2,3,1),hold on
function_plot_prsdata(prsdata.ca,prsdata.pcyl,[2:2:size(prsdata.pcyl,1)],0);
xlabel('Crank angle, \alpha [deg]')
ylabel('Cylinder pressure [bar]')
set(gca,'xlim',[-200 200],'ylim',[0 70])
hold on

plot(prsdata.ca(1,:),mean(prsdata.pint([2:2:NC],:)),'linewidth',2,'color',[1 1 1]*0.3)
plot(prsdata.ca(1,:),mean(prsdata.pcyl([1:2:NC],:)),'--','linewidth',2,'color',[1 1 1]*0.7)
plot(prsdata.ca(1,:),mean(prsdata.pcyl([2:2:NC],:)),'linewidth',2,'color','k')
mean_pcyl = mean(prsdata.pcyl([2:2:NC],:));
for ii=1:length(prsdata.sse_degV)
    plot(prsdata.ca(1,prsdata.sse_index(ii):prsdata.ese_index(ii)),mean_pcyl(prsdata.sse_index(ii):prsdata.ese_index(ii)),'linewidth',3,'color','r')
end
plot(prsdata.mean_alpha_pmax,prsdata.mean_pmax,'xr','markersize',10,'linewidth',2)
plot(prsdata.mean_alpha_pmotor,prsdata.mean_pmotor,'+r','markersize',10,'linewidth',2)
title({ ['P_{motor} = ' num2str(prsdata.mean_pmotor,'%3.1f') ' \pm '...
    num2str(prsdata.std_pmotor,'%3.1f') ' \alpha_{P_{motor}} = '...
    num2str(prsdata.mean_alpha_pmotor,'%3.1f') ' \pm '...
    num2str(prsdata.std_alpha_pmotor,'%3.1f') '^o'];...
    ['P_{max} = ' num2str(prsdata.mean_pmax,'%3.1f') ' \pm '...
    num2str(prsdata.std_pmax,'%3.1f') ' \alpha_{P_{max}} = '...
    num2str(prsdata.mean_alpha_pmax,'%3.1f') ' \pm '...
    num2str(prsdata.std_alpha_pmax,'%3.1f') '^o']})

% AHRR
subplot(2,3,2),grid on,hold on
function_plot_prsdata(prsdata.ca,prsdata.dQdCAfilt,[2:2:size(prsdata.pcyl,1)],0);
xlabel('Crank angle, \alpha [deg]')
ylabel('AHRR [J/deg]')
set(gca,'xlim',[-30 70],'ylim',[-50 500])
hold on

mean_dQdCAfilt = mean(prsdata.dQdCAfilt([2:2:NC],:));
mean_dQdCA = mean(prsdata.dQdCA([2:2:NC],:));

plot(prsdata.ca(1,:),mean(prsdata.dQdCAfilt([1:2:NC],:)),'--','linewidth',2,'color',[1 1 1]*0.7)
plot(prsdata.ca(1,:),mean(prsdata.dQdCAfilt([2:2:NC],:)),'linewidth',3,'color','k')
for ii=1:length(prsdata.sse_degV)
    plot(prsdata.ca(1,prsdata.sse_index(ii):prsdata.ese_index(ii)),mean_dQdCAfilt(prsdata.sse_index(ii):prsdata.ese_index(ii)),'linewidth',3,'color','r')
end

if isfinite(prsdata.mean_MEHRR)
    plot([prsdata.mean_SOC prsdata.mean_SOC prsdata.mean_EOCfilt prsdata.mean_EOCfilt], [0 prsdata.mean_MEHRR prsdata.mean_MEHRR 0],'r:','linewidth',4)
    plot([prsdata.mean_SOCfilt prsdata.mean_SOCfilt prsdata.mean_EOCfilt, prsdata.mean_EOCfilt], [0 prsdata.mean_MEHRRfilt prsdata.mean_MEHRRfilt 0],'b:','linewidth',2)
    %plot(prsdata.mean_alpha_AHRRmax,prsdata.mean_AHRRmax,'or','markersize',10,'linewidth',2,'markerfacecolor','w')
    plot(prsdata.mean_alpha_AHRRmaxfilt,prsdata.mean_AHRRmaxfilt,'dr','markersize',10,'linewidth',2,'markerfacecolor','w')
    plot(prsdata.mean_SOC,0,'or','markersize',6,'linewidth',2,'markerfacecolor','w')
    plot(prsdata.mean_SOCfilt,0,'db','markersize',6,'linewidth',2,'markerfacecolor','w')
    plot(prsdata.mean_EOCfilt,0,'or','markersize',6,'linewidth',2,'markerfacecolor','w')
    DOC_ms = prsdata.mean_DOC*139.8;
    title({ ['\alpha_{AHRRmax} = ' num2str(prsdata.mean_alpha_AHRRmax,'%3.1f') ' \pm ' num2str(prsdata.std_alpha_AHRRmax,'%3.1f') '^o'];...
        ['MEAHRR = ' num2str(prsdata.mean_MEHRR,'%3.1f') '-' num2str(prsdata.mean_MEHRRfilt,'%3.1f') ' J/deg'];...
        [' DAHRR = ' num2str(prsdata.mean_DOC*139.8,'%3.0f') ' \pm ' num2str(2*prsdata.std_DOC*139.8,'%3.0f') ' \mus'];});

else
    title('No combustion detected');
end

% Zoom on AHRR
subplot(2,3,3),grid on,hold on

hold on
plot(prsdata.ca(1,:),mean_dQdCA,'linewidth',3,'color','k')
for ii=1:length(prsdata.sse_degV)
    plot(prsdata.ca(1,prsdata.sse_index(ii):prsdata.ese_index(ii)),mean_dQdCA(prsdata.sse_index(ii):prsdata.ese_index(ii)),'linewidth',3,'color','r')
end
plot(prsdata.ca(1,:),mean(prsdata.dQdCA([1:2:size(prsdata.pcyl,1)],:)),':','linewidth',2,'color',[1 1 1]*0.7)
plot(prsdata.ca(1,:),mean_dQdCAfilt,'--','linewidth',2,'color','b')


try set(gca,'xlim',[prsdata.sse_degV-1 prsdata.mean_SOCfilt+2],'ylim',[min(mean_dQdCAfilt) 20],'fontsize',10,'position',[0.557018 0.75569 0.0909649 0.117871]);
catch set(gca,'xlim',[-10 10],'ylim',[min(mean_dQdCAfilt) 20],'fontsize',10,'position',[0.557018 0.75569 0.0909649 0.117871]);
end
% Evolution of ensemble average of AHRR
subplot(2,3,4),grid on,hold on

for jj=1:length(prsdata.AHRRmax)
    cmAHRR_std(jj) = std(prsdata.AHRRmax(1:jj));
    cmAHRR_stdfilt(jj) = std(prsdata.AHRRmaxfilt(1:jj));
end
e4 = errorbar([1:NC/2],cumsum(prsdata.AHRRmax)./[1:NC/2], 2*cmAHRR_std,'-o','color',[1 1 1]*0.7,'markersize',10,'linewidth',3,'markerfacecolor','w')
e5 = errorbar([1:NC/2],cumsum(prsdata.AHRRmaxfilt)./[1:NC/2], 2*cmAHRR_stdfilt,'-d','color',[1 1 1]*0.5,'markersize',10,'linewidth',1,'markerfacecolor','w')
plot(prsdata.AHRRmax,':','color','k','markersize',10,'linewidth',2,'markerfacecolor','w')
plot(prsdata.AHRRmaxfilt,':','color','k','markersize',10,'linewidth',2,'markerfacecolor','w')
hl=legend('running mean and 95% conf. intrvl','filtered data running mean and 95% CI','measured value');set(hl,'position',[0.45235 0.001 0.0667 0.0663],'orientation','horizontal')
cm = colormap('jet');
for ii=1:length(prsdata.AHRRmax)
    icolor = round(63*(ii-1)/(length(prsdata.AHRRmax)-1)+1);
    plot(ii,prsdata.AHRRmaxfilt(ii),'d','color',0.5*cm(icolor,:),'markersize',6,'linewidth',2);
    plot(ii,prsdata.AHRRmax(ii),'o','color',cm(icolor,:),'markersize',10,'linewidth',3,'markerfacecolor','w');

end
set(gca,'box','on','linewidth',1,'fontsize',12,'xlim',[0 NC/2+2])
xlabel('injection #')
ylabel('AHRRmax [J/deg]')
title({['<AHRR_{max}> = ' num2str(prsdata.mean_AHRRmax,'%3.0f') ' \pm ' num2str(2*prsdata.std_AHRRmax,'%2.f') ];...
    ['<AHRRfilt_{max}> = ' num2str(prsdata.mean_AHRRmaxfilt,'%3.0f') ' \pm ' num2str(2*prsdata.std_AHRRmaxfilt,'%2.f') ];});

    
%         % Evolution of ensemble average of alphaAHRR
%         subplot(2,3,5),grid on,hold on
%         plot(cumsum(prsdata.alpha_AHRRmax)./[1:length(prsdata.alpha_AHRRmax)],'-o','color',[1 1 1]*0.7,'markersize',10,'linewidth',3,'markerfacecolor','w')
%         plot(prsdata.alpha_AHRRmax,'--','color','k','markersize',10,'linewidth',2,'markerfacecolor','w')
%         for ii=1:length(prsdata.alpha_AHRRmax)
%             icolor = round(63*(ii-1)/(length(prsdata.alpha_AHRRmax)-1)+1);
%             plot(ii,prsdata.alpha_AHRRmax(ii),'o','color',cm(icolor,:),'markersize',10,'linewidth',3,'markerfacecolor','w');
%         end
%         set(gca,'box','on','linewidth',2,'fontsize',12,'fontname','eras light itc')
%         xlabel('injection #')
%         ylabel('AHRRmax [J/deg]')
%         title(['<\alphaAHRR_{max}> = ' num2str(prsdata.mean_alpha_AHRRmax,'%3.1f') ' \pm ' num2str(prsdata.std_alpha_AHRRmax,'%3.1f') 'J/deg']);


% Evolution of ensemble average of SOC
subplot(2,3,5),grid on; hold on
if isfinite(prsdata.mean_MEHRR)
ensSOC_CAD = cumsum(prsdata.SOC)'./[1:prsdata.NC/2];
ensAID = cumsum(prsdata.AID)'./[1:prsdata.NC/2]; %prsdata.SOC [CAD] TDC=0deg, sse usually <0, inj delay +CAD
for jj=1:length(prsdata.SOC)
    cmSOC_std(jj) = std(prsdata.SOC(1:jj))
    cmSOC_stdfilt(jj) = std(prsdata.SOCfilt(1:jj))
end
errorbar([1:NC/2],cumsum(prsdata.SOC)'./[1:length(prsdata.SOC)],2*cmSOC_std,'-o','color',[1 1 1]*0.7,'markersize',10,'linewidth',3,'markerfacecolor','w')
errorbar([1:NC/2],cumsum(prsdata.SOCfilt)'./[1:length(prsdata.SOCfilt)],2*cmSOC_stdfilt,'-d','color',[1 1 1]*0.5,'markersize',10,'linewidth',1,'markerfacecolor','w')
plot(1:prsdata.NC/2, prsdata.SOC,'--','color','k','markersize',10,'linewidth',2,'markerfacecolor','w')
plot(1:prsdata.NC/2, prsdata.SOCfilt,':','color','k','markersize',10,'linewidth',2,'markerfacecolor','w')
for ii=1:length(prsdata.SOC)
    icolor = round(63*(ii-1)/(length(prsdata.SOC)-1)+1);
    try
    plot(ii,prsdata.SOC(ii),'o','color',cm(icolor,:),'markersize',10,'linewidth',3,'markerfacecolor','w');
    plot(ii,prsdata.SOCfilt(ii),'d','color',.5*cm(icolor,:),'markersize',6,'linewidth',2);
    end
end
[Ax, H1, H2] = plotyy(1:prsdata.NC/2, ensSOC_CAD ,1:prsdata.NC/2, ensAID);%,{'-o','color',[1 1 1]*0.7,'markersize',10,'linewidth',3,'markerfacecolor','w'})

set(Ax(1),'box','on','linewidth',1,'fontsize',12,'xlim',[0 NC/2+2],'FontColor',[0 0 0],'YColor',[0 0 0])
Llim = get(Ax(1),'ylim');
Lticks = get(Ax(1),'ytick');
try
Rlim = (Llim-(prsdata.sse_degV(1)+prsdata.mean_inj_delay))*0.1389;

%Rticks = round((Lticks-(prsdata.sse_degV(1)+prsdata.mean_inj_delay))*0.1389*100)/100;
%set(Ax(2),'linewidth',1,'fontsize',12,'ylim',Rlim,'ytick',[round(Rlim(1)*100)/100:Rticks(2)-Rticks(1):round(Rlim(2)*100)/100],'xlim',[0 NC/2+2])
Rticks = [0 0.5 1 1.5 2 2.5];
set(Ax(2),'linewidth',1,'fontsize',12,'ylim',Rlim,'ytick',Rticks,'xlim',[0 NC/2+2],'FontColor',[0 0 0],'YColor',[0 0 0]);
catch
end

xlabel(Ax(1),'injection #')
ylabel(Ax(1),'Start of Combustion [deg]')
ylabel(Ax(2),'Ignition Delay [ms]')
%         legend('Ens. avg','Instantaneous')
title({['<SOC> = ' num2str(prsdata.mean_SOC,'%3.1f') ' \pm ' num2str(2*prsdata.std_SOC,'%3.1f') '^o' ' <SOC_{filt}> = ' num2str(prsdata.mean_SOCfilt,'%3.1f') '^o'];...
    ['<ID> = ' num2str(ensAID(end),'%3.2f') ' \pm ' num2str(2*prsdata.std_SOC*0.138,'%3.2f') ' <ID_{filt} > = ' num2str(prsdata.mean_AIDfilt,'%3.2f')];})
else
    [Ax, H1, H2] = plotyy(1:prsdata.NC/2, zeros(NC/2)+(prsdata.sse_degV(1)+prsdata.mean_inj_delay), 1:prsdata.NC/2, zeros(NC/2));%,{'-o','color',[1 1 1]*0.7,'markersize',10,'linewidth',3,'markerfacecolor','w'})
    set(Ax(1),'box','on','linewidth',1,'fontsize',12,'xlim',[0 NC/2+2],'ylim',[-5 10],'ytick',[-5 0 5 10])
    Llim = get(Ax(1),'ylim');
    Rlim = (Llim-(prsdata.sse_degV(1)+prsdata.mean_inj_delay))*0.1389;
    Rticks = [0 .5 1 1.5 2 2.5]
    set(Ax(2),'linewidth',1,'fontsize',12,'ylim',Rlim,'ytick',Rticks,'xlim',[0 NC/2+2])
    xlabel(Ax(1),'injection #')
    ylabel(Ax(1),'Start of Combustion [deg]')
    ylabel(Ax(2),'Ignition Delay [ms]')
    title('No combustion detected')
end

% Evolution of ensemble average of gIMEP (plot units as kPa so divide by 1000)
subplot(2,3,6),grid on,hold on
for jj=1:length(prsdata.gIMEP)
    cmgIMEP_std(jj) = std(prsdata.gIMEP(1:jj));
end
errorbar([1:NC/2],cumsum(prsdata.gIMEP/1000)./[1:length(prsdata.gIMEP)],2*cmgIMEP_std/1000,'-o','color',[1 1 1]*0.7,'markersize',10,'linewidth',3,'markerfacecolor','w')
plot(prsdata.gIMEP/1000,'--','color','k','markersize',10,'linewidth',2,'markerfacecolor','w')
for ii=1:length(prsdata.gIMEP)
    icolor = round(63*(ii-1)/(length(prsdata.gIMEP)-1)+1);
    plot(ii,prsdata.gIMEP(ii)/1000,'o','color',cm(icolor,:),'markersize',10,'linewidth',3,'markerfacecolor','w');
end
p = get(gca,'position');
set(gca,'box','on','linewidth',1,'fontsize',12,'position',[p(1)+.01 p(2) p(3) p(4)],'xlim',[0 NC/2+2])
set(gca,'YAxisLocation','right')
xlabel('injection #')
ylabel('gIMEP [kPA]')
title(['<gIMEP> = ' num2str(prsdata.mean_gIMEP/1000,'%3.1f') ' \pm ' num2str(2*prsdata.std_gIMEP/1000,'%2.0f')])


% Compute compression and expansion polytropic coefficients
ca1 = -90;ca2=-30;
[a,ica1] = min(abs(prsdata.ca(1,:)-ca1));[a,ica2] = min(abs(prsdata.ca(1,:)-ca2));
logvolc = log(prsdata.volume(2:2:end,ica1:ica2));logvolc=logvolc(:);
logpcylc = log(prsdata.pcyl(2:2:end,ica1:ica2));logpcylc=logpcylc(:);
pcomp = polyfit(logvolc,logpcylc,1);ncomp = pcomp(1);
ycomp = polyval(pcomp,[min(logvolc) max(logvolc)]);
ycomp2 = polyval(pcomp,[min(log(prsdata.volume(1,:))) max(log(prsdata.volume(1,:)))]);

ca1 = 40;ca2=100;
[a,ica1] = min(abs(prsdata.ca(1,:)-ca1));[a,ica2] = min(abs(prsdata.ca(1,:)-ca2));
logvole = log(prsdata.volume(2:2:end,ica1:ica2));logvole=logvole(:);
logpcyle = log(prsdata.pcyl(2:2:end,ica1:ica2));logpcyle=logpcyle(:);
pexp = polyfit(logvole,logpcyle,1);nexp = pexp(1);
yexp = polyval(pexp,[min(logvole) max(logvole)]);
yexp2 = polyval(pexp,[min(log(prsdata.volume(1,:))) max(log(prsdata.volume(1,:)))]);

subplot(2,3,3),hold on
function_plot_prsdata(log(prsdata.volume),log(prsdata.pcyl),[2:2:size(prsdata.pcyl,1)],0);
plot([min(log(prsdata.volume(1,:))) max(log(prsdata.volume(1,:)))],ycomp2,'--','color',[1 1 1]*0.7,'linewidth',2)
plot([min(logvolc) max(logvolc)],ycomp,'k','linewidth',3)
plot([min(log(prsdata.volume(1,:))) max(log(prsdata.volume(1,:)))],yexp2,'--','color',[1 1 1]*0.7,'linewidth',2)
plot([min(logvole) max(logvole)],yexp,'k','linewidth',3)

ht=title({['n_{comp} = ' num2str(-ncomp,'%1.3f')];['n_{exp} = ' num2str(-nexp,'%1.3f')]});set(ht,'fontsize',10);
set(gca,'xlim',[-8.5 -5.5],'ylim',[0 4.5],'fontsize',8,'xticklabel','','yticklabel','','position',[0.150439 0.685789 0.0672982 0.127582])
