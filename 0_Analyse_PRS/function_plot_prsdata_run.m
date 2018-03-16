function prsdata = function_plot_prsdata_run(prsdata)

ienergize = find(prsdata.injComm(2,:)>1);
NC = prsdata.NC;

    set(0,'defaultAxesFontName','Arial')
    set(0,'defaultTextFontName','Arial')
    
% Pcyl
p1 = subplot(2,3,1),hold on
plot(prsdata.ca(1,:),mean(prsdata.pcyl([2:2:NC],:))+2/sqrt(NC)*std(prsdata.pcyl([2:2:NC],:)),'linewidth',1,'color',[1 1 1]*0.4)
plot(prsdata.ca(1,:),mean(prsdata.pcyl([2:2:NC],:))-2/sqrt(NC)*std(prsdata.pcyl([2:2:NC],:)),'linewidth',1,'color',[1 1 1]*0.4)
plot(prsdata.ca(1,:),mean(prsdata.pcyl([1:2:NC],:))+2/sqrt(NC)*std(prsdata.pcyl([1:2:NC],:)),'--','linewidth',1,'color',[1 1 1]*0.4)
plot(prsdata.ca(1,:),mean(prsdata.pcyl([1:2:NC],:))-2/sqrt(NC)*std(prsdata.pcyl([1:2:NC],:)),'--','linewidth',1,'color',[1 1 1]*0.4)
function_plot_prsdata(prsdata.ca,prsdata.pcyl,[2:2:size(prsdata.pcyl,1)],0);
xlabel('Crank angle, \alpha [deg]')
ylabel('Cylinder pressure [bar]')
set(gca,'xlim',[-200 200],'ylim',[0 70])
hold on

plot(prsdata.ca(1,:),mean(prsdata.pint([2:2:NC],:)),'linewidth',2,'color',[1 1 1]*0.8)
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

mean_dQdCAfilt = mean(prsdata.dQdCAfilt([2:2:NC],:));
mean_dQdCA = mean(prsdata.dQdCA([2:2:NC],:));
ticksperCAD = 4;

% Instantaneous AHRR
p2 = subplot(2,3,2),grid on,hold on

plot(prsdata.ca(1,:),mean(prsdata.dQdCA([1:2:NC],:))+2/sqrt(NC)*std(prsdata.dQdCA([1:2:NC],:)),'-','linewidth',2,'color',[1 1 1]*0.3)
plot(prsdata.ca(1,:),mean(prsdata.dQdCA([1:2:NC],:))-2/sqrt(NC)*std(prsdata.dQdCA([1:2:NC],:)),'-','linewidth',2,'color',[1 1 1]*0.3)

hold on
function_plot_prsdata(prsdata.ca,prsdata.dQdCA,[2:2:size(prsdata.pcyl,1)],0);

plot(prsdata.ca(1,:),mean_dQdCA,'linewidth',3,'color','k')
for ii=1:length(prsdata.sse_degV)
    plot(prsdata.ca(1,prsdata.sse_index(ii)+prsdata.mean_inj_delay*ticksperCAD:prsdata.ese_index(ii)+prsdata.mean_inj_delay*ticksperCAD),...
        mean_dQdCA(prsdata.sse_index(ii)+prsdata.mean_inj_delay*ticksperCAD:prsdata.ese_index(ii)+prsdata.mean_inj_delay*ticksperCAD),'linewidth',3,'color','r');
end
plot(prsdata.ca(1,:),mean(prsdata.dQdCA([1:2:size(prsdata.pcyl,1)],:)),'-','linewidth',2,'color',[1 1 1]*0.7)
plot(prsdata.ca(1,:),mean_dQdCAfilt,'-','linewidth',2,'color','b')
ca2 = [prsdata.ca(1,:) fliplr(prsdata.ca(1,:))];
ub = mean(prsdata.dQdCA([2:2:NC],:))+2/sqrt(NC)*std(prsdata.dQdCA([2:2:NC],:));
lb = mean(prsdata.dQdCA([2:2:NC],:))-2/sqrt(NC)*std(prsdata.dQdCA([2:2:NC],:));
fill(ca2,[ub fliplr(lb)],[1 1 1]*.5,'FaceAlpha',0.5)

try set(gca,'xlim',[prsdata.sse_degV prsdata.mean_alpha_AHRRmax+2],'ylim',[min(mean_dQdCAfilt) prsdata.mean_AHRRmax+10],'fontsize',10);
catch set(gca,'xlim',[-10 10],'ylim',[min(mean_dQdCAfilt) 20],'fontsize',10);
end

if isfinite(prsdata.mean_DOC)
    plot([prsdata.mean_SOC prsdata.mean_SOC prsdata.mean_EOCfilt prsdata.mean_EOCfilt], [0 prsdata.mean_MEHRR prsdata.mean_MEHRR 0],'r:','linewidth',4)
    plot([prsdata.mean_SOCfilt prsdata.mean_SOCfilt prsdata.mean_EOCfilt, prsdata.mean_EOCfilt], [0 prsdata.mean_MEHRRfilt prsdata.mean_MEHRRfilt 0],'b:','linewidth',2)
    %plot(prsdata.mean_alpha_AHRRmax,prsdata.mean_AHRRmax,'or','markersize',10,'linewidth',2,'markerfacecolor','w')
    %plot(prsdata.mean_alpha_AHRRmaxfilt,prsdata.mean_AHRRmaxfilt,'dr','markersize',10,'linewidth',2,'markerfacecolor','w')
    plot(prsdata.mean_SOC,0,'or','markersize',6,'linewidth',2,'markerfacecolor','w')
    plot(prsdata.mean_SOCfilt,0,'db','markersize',6,'linewidth',2,'markerfacecolor','w')
    plot(prsdata.mean_EOCfilt,0,'or','markersize',6,'linewidth',2,'markerfacecolor','w')
    DOC_ms = prsdata.mean_DOC*139.8;
    ll = min([prsdata.mean_MEHRRfilt prsdata.mean_MEHRR]);
    ul = max([prsdata.mean_MEHRRfilt prsdata.mean_MEHRR]);
    title({ ['\alpha_{AHRRmax} = ' num2str(prsdata.mean_alpha_AHRRmax,'%3.1f') ' \pm ' num2str(prsdata.std_alpha_AHRRmax,'%3.1f') '^o'];...
        ['MEAHRR = ' num2str(ll,'%3.1f') '-' num2str(ul,'%3.1f') ' J/deg'];...
        [' DAHRR = ' num2str(round(prsdata.mean_DOC*139.8/10)/100,'%1.2f') ' \pm ' num2str(round(2*prsdata.std_DOC*139.8/sqrt(NC)/10)/100,'%1.2f') ' ms'];});
    
else
    title('No combustion detected');
end


    
% Evolution of ensemble average of AHRR
p3 = subplot(2,3,4),grid on,hold on

ca2 = [prsdata.ca(1,:) fliplr(prsdata.ca(1,:))];
ub = mean(prsdata.dQdCAfilt([2:2:NC],:))+2/sqrt(n-1)*std(prsdata.dQdCA([2:2:NC],:));
lb = mean(prsdata.dQdCAfilt([2:2:NC],:))-2/sqrt(n-1)*std(prsdata.dQdCA([2:2:NC],:));
fill(ca2,[ub fliplr(lb)],[1 1 1]*.5,'FaceAlpha',0.5)




%function_plot_prsdata(prsdata.ca,prsdata.dQdCA,[2:2:size(prsdata.pcyl,1)],0);
xlabel('Crank angle, \alpha [deg]')
ylabel('AHRR [J/deg]')
try
set(gca,'xlim',[prsdata.ca(1,prsdata.sse_index(1)) prsdata.mean_EOCfilt+10],'ylim',[-50 500]);
catch
set(gca,'xlim',[-30 70],'ylim',[-50 500]);
end
hold on
plot(prsdata.ca(1,:),mean(prsdata.dQdCAfilt([2:2:NC],:)),'linewidth',2,'color','k')
plot(prsdata.ca(1,:),mean(prsdata.dQdCAfilt([1:2:NC],:)),'--','linewidth',2,'color','k')

mean_dQdCAfilt = mean(prsdata.dQdCAfilt([2:2:NC],:));
mean_dQdCA = mean(prsdata.dQdCA([2:2:NC],:));
ticksperCAD = 4;
for ii=1:length(prsdata.sse_degV)
    plot(prsdata.ca(1,prsdata.sse_index(ii)+prsdata.mean_inj_delay*ticksperCAD:prsdata.ese_index(ii)+prsdata.mean_inj_delay*ticksperCAD),...
        mean_dQdCAfilt(prsdata.sse_index(ii)+prsdata.mean_inj_delay*ticksperCAD:prsdata.ese_index(ii)+prsdata.mean_inj_delay*ticksperCAD),'linewidth',3,'color','r')
    % accounts for hydraulic delay in plotting DSE, still missing extra
    % closing time (injector specific)
end

if isfinite(prsdata.mean_DOC)
    plot([prsdata.mean_SOC prsdata.mean_SOC prsdata.mean_EOCfilt prsdata.mean_EOCfilt], [0 prsdata.mean_MEHRR prsdata.mean_MEHRR 0],'r:','linewidth',4)
    plot([prsdata.mean_SOCfilt prsdata.mean_SOCfilt prsdata.mean_EOCfilt, prsdata.mean_EOCfilt], [0 prsdata.mean_MEHRRfilt prsdata.mean_MEHRRfilt 0],'b:','linewidth',2)
    %plot(prsdata.mean_alpha_AHRRmax,prsdata.mean_AHRRmax,'or','markersize',10,'linewidth',2,'markerfacecolor','w')
    %plot(prsdata.mean_alpha_AHRRmaxfilt,prsdata.mean_AHRRmaxfilt,'dr','markersize',10,'linewidth',2,'markerfacecolor','w')
    plot(prsdata.mean_SOC,0,'or','markersize',6,'linewidth',2,'markerfacecolor','w')
    plot(prsdata.mean_SOCfilt,0,'db','markersize',6,'linewidth',2,'markerfacecolor','w')
    plot(prsdata.mean_EOCfilt,0,'or','markersize',6,'linewidth',2,'markerfacecolor','w')
    DOC_ms = prsdata.mean_DOC*139.8;
    ll = min([prsdata.mean_MEHRRfilt prsdata.mean_MEHRR]);
    ul = max([prsdata.mean_MEHRRfilt prsdata.mean_MEHRR]);
    
else
    title('No combustion detected');
end

% Zoom on CHRR
p4 = subplot(2,3,3),grid on,hold on
CHRR1 = zeros(NC/2,length(prsdata.ca(1,:)));
CHRR2 = zeros(NC/2,length(prsdata.ca(1,:)));
for jj=1:length(prsdata.SOC)
    rng1 = find(prsdata.ca(1,:)>=prsdata.SOC(jj) & prsdata.ca(1,:)<=prsdata.EOCfilt(jj)+30);
    rng2 = find(prsdata.ca(1,:)>=prsdata.SOCfilt(jj) & prsdata.ca(1,:)<=prsdata.EOCfilt(jj)+30);
    CHRR1(jj,rng1) = cumsum(prsdata.dQdCA(2*jj,rng1));
    CHRR2(jj,rng2) = cumsum(prsdata.dQdCAfilt(2*jj,rng2));
end
CHRRavg = cumsum(mean_dQdCA(rng1));
CHRRavgfilt = cumsum(mean_dQdCAfilt(rng2));
CHRR1std = std(CHRR1,1);
CHRR2std = std(CHRR2,1);

ub = mean(CHRR1)/max(CHRRavgfilt)+2*CHRR1std/max(CHRRavgfilt);
lb = mean(CHRR1)/max(CHRRavgfilt)-2*CHRR1std/max(CHRRavgfilt);
fill(ca2,[ub fliplr(lb)],[1 1 1]*.9,'FaceAlpha',0.5)


cm = colormap('jet');
for ii=1:length(prsdata.SOC)
    icolor = round(63*(ii-1)/(length(prsdata.AHRRmax)-1)+1);
    plot(prsdata.ca(1,:),CHRR2(ii,:)./max(CHRRavgfilt),'color',cm(icolor,:),'markersize',6,'linewidth',1);
    CA10index = find(CHRR1(ii,:)>0.1*max(CHRRavgfilt)); try CA10(ii) = prsdata.ca(1,CA10index(1)); catch CA10(ii)= NaN; end
    CA50index = find(CHRR1(ii,:)>0.5*max(CHRRavgfilt)); try CA50(ii) = prsdata.ca(1,CA50index(1)); catch CA50(ii)= NaN; end
    CA90index = find(CHRR1(ii,:)>0.9*max(CHRRavgfilt)); try CA90(ii) = prsdata.ca(1,CA90index(1)); catch CA90(ii)= NaN; end
end

prsdata.CA10 = CA10;
prsdata.CA50 = CA50;
prsdata.CA90 = CA90;
prsdata.mean_CA10= nanmean(CA10);
prsdata.mean_CA50= nanmean(CA50);
prsdata.mean_CA90 = nanmean(CA90);
prsdata.std_CA10 = nanstd(CA10);
prsdata.std_CA50= nanstd(CA50);
prsdata.std_CA90=nanstd(CA90);
title('Cumulative HRR')
plot(prsdata.ca(1,rng1),CHRRavg/max(CHRRavgfilt),'k','linewidth',3)
set(gca,'xlim',[ca2(1,rng2(1)) ca2(1,rng2(end))],'ylim',[-.1 1.2]);
try set(gca,'xlim',[prsdata.mean_SOCfilt-5 prsdata.mean_EOCfilt+2],'fontsize',10,'position',[0.257018 0.26569 0.0909649 0.117871]);
catch set(gca,'xlim',[-10 10],'ylim',[min(mean_dQdCAfilt) 20],'fontsize',10,'position',[0.257018 0.26569 0.0909649 0.117871]);
end


if isfinite(prsdata.mean_MEHRRfilt)
axes(p3)
title({['CA10 ' num2str(prsdata.mean_CA10,'%3.1f') ' \pm ' num2str(prsdata.std_CA10,'%2.2f') 'CAD'];...
        ['CA50 ' num2str(prsdata.mean_CA50,'%3.1f') ' \pm ' num2str(prsdata.std_CA50,'%2.2f') 'CAD'];...
        [' CA90 ' num2str(prsdata.mean_CA90,'%3.1f') ' \pm ' num2str(prsdata.std_CA90,'%2.2f') ' CAD'];});
end
axes(p4)
% figure(),errorbar(prsdata.ca(1,:),CHRR1./max(CHRRavgfilt),2*CHRR1std,'-d','color',[1 1 1]*0.5,'markersize',10,'linewidth',1)
% hold on, errorbar(prsdata.ca(1,:),CHRR2./max(CHRRavgfilt),2*CHRR1std,'-d','color',[1 1 1]*0.5,'markersize',10,'linewidth',1)
% try set(gca,'xlim',[prsdata.mean_SOCfilt-10 prsdata.mean_EOCfilt+5],'ylim',[min(mean_dQdCAfilt) 50],'fontsize',10,'position',[0.257018 0.25569 0.0909649 0.117871]);
% catch
% end

% previous case plotted AHRRmax and 
% e4 = errorbar([1:NC/2],cumsum(prsdata.AHRRmax)./[1:NC/2], 2*cmAHRRavg_std,'-o','color',[1 1 1]*0.7,'markersize',10,'linewidth',3,'markerfacecolor','w');
% e5 = errorbar([1:NC/2],cumsum(prsdata.AHRRmaxfilt)./[1:NC/2], 2*cmAHRRavg_stdfilt,'-d','color',[1 1 1]*0.5,'markersize',10,'linewidth',1,'markerfacecolor','w');
% plot(prsdata.AHRRmax,':','color','k','markersize',10,'linewidth',2,'markerfacecolor','w')
% plot(prsdata.AHRRmaxfilt,':','color','k','markersize',10,'linewidth',2,'markerfacecolor','w')
% hl=legend('running mean and 95% conf. intrvl','filtered data running mean and 95% CI','measured value');set(hl,'position',[0.47235 0.0001 0.0667 0.0663],'orientation','horizontal')
% cm = colormap('jet');
% for ii=1:length(prsdata.AHRRmax)
%     icolor = round(63*(ii-1)/(length(prsdata.AHRRmax)-1)+1);
%     plot(ii,prsdata.AHRRmaxfilt(ii),'d','color',0.5*cm(icolor,:),'markersize',6,'linewidth',2);
%     plot(ii,prsdata.AHRRmax(ii),'o','color',cm(icolor,:),'markersize',10,'linewidth',3,'markerfacecolor','w');
% 
% end
% set(gca,'box','on','linewidth',1,'fontsize',12,'xlim',[0 NC/2+2])
% xlabel('injection #')
% ylabel('AHRRmax [J/deg]')
% title({['<AHRR_{max}> = ' num2str(prsdata.mean_AHRRmax,'%3.0f') ' \pm ' num2str(2*prsdata.std_AHRRmax/sqrt(NC),'%2.f') ];...
%     ['<AHRRfilt_{max}> = ' num2str(prsdata.mean_AHRRmaxfilt,'%3.0f') ' \pm ' num2str(2*prsdata.std_AHRRmaxfilt/sqrt(NC),'%2.f') ];});

    
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
if isfinite(prsdata.mean_MEHRRfilt)
ensSOC_CAD = cumsum(prsdata.SOC)'./[1:prsdata.NC/2];
ensAID = cumsum(prsdata.AID)'./[1:prsdata.NC/2]; %prsdata.SOC [CAD] TDC=0deg, sse usually <0, inj delay +CAD
for jj=1:length(prsdata.SOC)
    cmSOC_std(jj) = 1/sqrt(jj)*std(prsdata.SOC(1:jj));
    cmSOC_stdfilt(jj) = 1/sqrt(jj)*std(prsdata.SOCfilt(1:jj));
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

set(Ax(1),'box','on','linewidth',1,'fontsize',12,'xlim',[0 NC/2+2],'YColor',[0 0 0])
Llim = get(Ax(1),'ylim');
Lticks = get(Ax(1),'ytick');
try
Rlim = (Llim-(prsdata.sse_degV(1)+prsdata.mean_inj_delay))*0.1389;

Rticks = round((Lticks-(prsdata.sse_degV(1)+prsdata.mean_inj_delay))*0.1389*100)/100;
set(Ax(2),'linewidth',1,'fontsize',12,'ylim',Rlim,'ytick',[round(Rlim(1)*100)/100:Rticks(2)-Rticks(1):round(Rlim(2)*100)/100],'xlim',[0 NC/2+2])
catch
end

xlabel(Ax(1),'injection #')
ylabel(Ax(1),'Start of Combustion [deg]')
ylabel(Ax(2),'Ignition Delay [ms]')
%         legend('Ens. avg','Instantaneous')
title({['<SOC> = ' num2str(prsdata.mean_SOC,'%3.1f') ' \pm ' num2str(2*prsdata.std_SOC/sqrt(NC),'%3.1f') '^o' ' <SOC_{filt}> = ' num2str(prsdata.mean_SOCfilt,'%3.1f') '^o'];...
    ['<ID> = ' num2str(ensAID(end),'%3.2f') ' \pm ' num2str(2*prsdata.std_SOC*0.138/sqrt(NC),'%3.2f') ' <ID_{filt} > = ' num2str(prsdata.mean_AIDfilt,'%3.2f')];})
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
    cmgIMEP_std(jj) = 1/sqrt(jj)*std(prsdata.gIMEP(1:jj));
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
