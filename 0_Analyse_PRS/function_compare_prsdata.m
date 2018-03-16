function function_compare_prsdata

% 2014/09/25 - L.M. MALBEC
%
% This function allows to plot data from prs files and to compare them for
% different runs.


close all


% Channel selection
ichannel = [4 5 7];

itime = 1;

if itime
    corrfac = 1e3/7200;
else
    corrfac = 1
end

% 1 = ca
% 2 = pcyl
% 3 = pint
% 4 = injComm
% 5 = prail
% 6 = injCurr
% 7 = camTrig


[filename,pathname] = uigetfile('*.prs','Select files','Multiselect','on');
filename = cellstr(filename);

figure,hold on
cm = colormap('jet');
set(gcf,'color','w')
set(gca,'box','on','linewidth',2,'fontsize',12,'fontname','eras light itc')

run_legend={};
for ii=1:length(filename)
    
    icolor = round(63*(ii-1)/(length(filename)-1)+1);
    if isnan(icolor)
        icolor=1;
    end
    
    dayname = filename{ii}(1:end-5);
    runname = filename{ii}(end-4);
    
    [data,N,SPC] = function_read_prsfile(pathname,filename{ii});
    prsdata = function_read_prsdata(data,N,SPC);

    run_legend = [run_legend;runname];

    dinjComm = diff(prsdata.injComm(2,:));
    isse = find(dinjComm>1);
    if itime
        toff_ms = corrfac*prsdata.ca(1,isse);
        if isempty(toff_ms)
            toff_ms=0;
        end
    else
        toff_ms = 0;
    end
    
    for jj=1:length(ichannel)

        figure(jj),hold on
        cm = colormap('jet');
        set(gcf,'color','w')
        set(gca,'box','on','linewidth',2,'fontsize',12,'fontname','eras light itc')
        if ichannel(jj) == 2
            plot(prsdata.ca(1,:)*corrfac-toff_ms,mean(prsdata.pcyl([2:2:size(prsdata.pcyl,1)],:)),'linewidth',2,'color',cm(icolor,:));
            ylabel('Cylinder pressure [bar]')
        elseif ichannel(jj) == 3
            plot(prsdata.ca(1,:)*corrfac-toff_ms,mean(prsdata.pint([2:2:size(prsdata.pint,1)],:)),'linewidth',2,'color',cm(icolor,:));
            ylabel('Intake pressure [bar]')
        elseif ichannel(jj) == 4
            plot(prsdata.ca(1,:)*corrfac-toff_ms,mean(prsdata.injComm([2:2:size(prsdata.injComm,1)],:)),'linewidth',2,'color',cm(icolor,:));
            ylabel('Inj. Comm. [-]')
        elseif ichannel(jj) == 5
            plot(prsdata.ca(1,:)*corrfac-toff_ms,mean(prsdata.prail([2:2:size(prsdata.prail,1)],:)),'linewidth',2,'color',cm(icolor,:));
            ylabel('Rail pressure [bar]')
        elseif ichannel(jj) == 6
            plot(prsdata.ca(1,:)*corrfac-toff_ms,mean(prsdata.injCurr([2:2:size(prsdata.injCurr,1)],:)),'linewidth',2,'color',cm(icolor,:));
            ylabel('Injector current [-]')
        elseif ichannel(jj) == 7
            plot(prsdata.ca(1,:)*corrfac-toff_ms,mean(prsdata.camTrig([2:2:size(prsdata.camTrig,1)],:)),'linewidth',2,'color',cm(icolor,:));
            ylabel('Camera opening [-]')
        end
    end
end

if itime
    xlabel('time (ms)')
else
    xlabel('crankangle [deg]')
end
legend(run_legend);


% figure,plot(prsdata.injComm')
% set(gca,'xlim',[1400 1430]);
   