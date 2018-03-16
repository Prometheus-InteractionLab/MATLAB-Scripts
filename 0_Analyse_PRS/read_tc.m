%%
[filename,pathname] = uigetfile('*.prs','Select files','Multiselect','on');
filename = cellstr(filename);

for ii = 1:length(filename)
    if length(filename)>1
        close all
    end
dayname = filename{ii}(1:6)
runname = filename{ii}(end-4);

[data,N,SPC] = function_read_prsfile(pathname,filename{ii});
% prsdata = function_read_prsdata(data,N,SPC);

NC = N(2)./SPC; % number of cycles in data file, SPC = samples per cycle

motorcyc = 1:2:NC;
firecyc = 2:2:NC;

% ca = reshape(data(1,:),SPC,NC)';
% pcyl = reshape(data(2,:),SPC,NC)';
% injComm = reshape(data(4,:),SPC,NC)';
% prail = reshape(data(5,:),SPC,NC)';
% injCurr = reshape(data(6,:),SPC,NC)';

%%

load(fullfile(pathname,regexprep(filename{ii},'prs','mat')));
pcyl_m_m = mean(rundata.pcyl(motorcyc,:),1);
pcyl_f_m = mean(rundata.pcyl(firecyc,:),1);

mcyl = rundata.rhobdc*max(rundata.volume(1,:))*1000;
Tcyl_m_m = pcyl_m_m.*rundata.volume(1,:)./mcyl/rundata.rfg;
Tcyl_f_m = pcyl_f_m.*rundata.volume(1,:)./mcyl/rundata.rfg;

[ca,aa] = function_read_prschannel(data,1,firecyc,SPC,NC,1,0);
TCdata.ca = ca(1,:);
TCdata.Ttdc = rundata.Ttdc;
TCdata.nexp = rundata.nexp;
TCdata.ncomp = rundata.ncomp;
list_tc={'T0','U1','U3','D1','D3','L1','L3','R1','R3'};
figure(10),set(gcf,'position',[284 49 1076 929],'color','w')
cm=colormap(jet(256));
for jj=1:length(list_tc)
    
    icolor = round(255*(jj-1)/(length(list_tc)-1)+1);
    [prschan_m,tc_m_m] = function_read_prschannel(data,jj+7,motorcyc,SPC,NC,100,0);
    eval(['TCdata.' list_tc{jj} '.m = prschan_m'])
    tc_m_m = tc_m_m - mean(tc_m_m(100:105));
    
    [prschan_f,tc_f_m] = function_read_prschannel(data,jj+7,firecyc,SPC,NC,100,0);
    eval(['TCdata.' list_tc{jj} '.f = prschan_f'])
    tc_f_m = tc_f_m - mean(tc_f_m(100:105));    

%     figure(10)
%     subplot(3,3,jj),grid on,hold on
%     for kk=1:size(prschan_f,1)
%         icolor2 = round(255*(kk-1)/(size(prschan_f,1)-1)+1);
%         plot(ca(1,:),prschan_f(kk,:)-mean(prschan_f(kk,100:105)),'color',cm(icolor2,:))
%     end
%     plot(ca(1,:),tc_m_m,'color',[1 1 1]*0.7,'linewidth',2),plot(ca(1,:),tc_f_m,'color',[1 1 1]*0.0,'linewidth',2)
%     set(gca,'box','on','linewidth',1,'fontsize',12,'fontname','eras light itc'),title(list_tc{jj})
%     xlabel('Crank angle [deg]')
%     ylabel('\DeltaT [K]')
%     legend('motored','fired')
%     set(gca,'ylim',[-1 5])
%     
%     figure(100),grid on,hold on
%     set(gcf,'color','w')
%     plot(ca(1,:),tc_m_m,'-','color',cm(icolor,:))
%     plot(ca(1,:),tc_f_m,':','color',cm(icolor,:))
%     set(gca,'box','on','linewidth',1,'fontsize',12,'fontname','eras light itc'),title(list_tc{jj})
%     xlabel('T_{cyl} [K]')
%     ylabel('\DeltaT [K]')
%     set(gca,'ylim',[-1 5])
%     legend('motored','fired')

%     figure(200),grid on,hold on
%     set(gcf,'color','w')
%     plot(ca(1,:),tc_m_m,'-','color',cm(icolor,:))
%     plot(ca(1,:),tc_f_m,':','color',cm(icolor,:))
%     set(gca,'box','on','linewidth',1,'fontsize',12,'fontname','eras light itc'),title(list_tc{jj})
%     xlabel('Crank angle [deg]')
%     ylabel('\DeltaT [K]')
%     set(gca,'ylim',[-1 5])
end
% drawnow
% saveas(gcf,fullfile(pathname,['TC_' dayname runname '.fig']),'fig')
% saveas(gcf,fullfile(pathname,['TC_' dayname runname '.emf']),'emf')

save(fullfile(pathname,[dayname runname '_TC']),'TCdata');
end
%%
% prompt = {'O2 concentration (% vol.):'};
% dlg_title = 'Input for TDC conditions';
% num_lines = 1;
% answer = inputdlg(prompt,dlg_title,num_lines,{'';''});
% N2_properties_nist;
% O2_properties_nist;
% AR_properties_nist;
% AIR_properties_nist;
% pco2=yo2/100;
% ndil = (0.21/pco2-1); % Nb of N2 moles for 1 mole or air
% Mgf = (ndil*28+Mair)/(1+ndil);
% r_fg = 8.314/Mgf*1000;