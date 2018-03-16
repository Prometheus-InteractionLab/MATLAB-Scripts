%%
[filename,pathname] = uigetfile('*.prs','Select files','Multiselect','on');
filename = cellstr(filename);
legende={};
for ii=1:length(filename)
    icolor = round(255*(ii-1)/(length(filename)-1)+1);
    dayname = filename{ii}(1:6)
    runname = filename{ii}(end-4);

    [data,N,SPC] = function_read_prsfile(pathname,filename{ii});
    % prsdata = function_read_prsdata(data,N,SPC);

    NC = N(2)./SPC; % number of cycles in data file, SPC = samples per cycle

    motorcyc = 1:2:NC;
    firecyc = 2:2:NC;

    load(fullfile(pathname,regexprep(filename{ii},'prs','mat')));
    pcyl_m_m = mean(rundata.pcyl(motorcyc,:),1);
    pcyl_f_m = mean(rundata.pcyl(firecyc,:),1);

    mcyl = rundata.rhobdc*max(rundata.volume(1,:))*1000;
    Tcyl_m_m = pcyl_m_m.*rundata.volume(1,:)./mcyl/rundata.rfg;
    Tcyl_f_m = pcyl_f_m.*rundata.volume(1,:)./mcyl/rundata.rfg;

    [ca,aa] = function_read_prschannel(data,1,firecyc,SPC,NC,1,0);

    list_tc={'T0','U1','U3','D1','D3','L1','L3','R1','R3'};
    figure(10),set(gcf,'position',[284 49 1076 929],'color','w')
    cm=colormap(jet(256));
    legende=[legende;runname];
   
    for jj=1:length(list_tc)
        [prschan_m,tc_m_m] = function_read_prschannel(data,jj+7,motorcyc,SPC,NC,100,0);
        tc_m_m = tc_m_m - mean(tc_m_m(1000:1005));

        [prschan_f,tc_f_m] = function_read_prschannel(data,jj+7,firecyc,SPC,NC,100,0);
        tc_f_m = tc_f_m - mean(tc_f_m(1000:1005));

        figure(10)
        subplot(3,3,jj),grid on,hold on
        plot(ca(1,:),tc_f_m,'color',cm(icolor,:),'linewidth',2)
        set(gca,'box','on','linewidth',1,'fontsize',12,'fontname','eras light itc'),title(list_tc{jj})
        xlabel('Crank angle [deg]')
        ylabel('\DeltaT [K]')
        set(gca,'ylim',[-1 20])
        hl=legend(legende);set(hl,'visible','off');
    end
    

    saveas(gcf,fullfile(pathname,['TC_compare.fig']),'fig')
    saveas(gcf,fullfile(pathname,['TC_compare.emf']),'emf')
end