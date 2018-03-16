function [ fig,data ] = function_plot_processed_data(fig,mark ,colors,plotlegend,legendtext,digits,units)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[filename,pathname] = uigetfile('*.mat','Select files','Multiselect','on');
filename = cellstr(filename);
addpath(pathname) %added by E. Eagle 141009, prevents code crashing when selected file isn't already in the matlab path

for ii=1:length(filename)
    load(filename{ii})
    try data.T(ii) = rundata.Ttdc; catch data.T(ii) = NaN; end
    try data.legend = rundata.fuel_type; catch end
    try data.dse(ii) = rundata.dse_ms; catch data.dse(ii) = NaN; end
    try data.mean_SOC(ii) = rundata.mean_SOC; data.std_SOC(ii) = rundata.std_SOC; catch data.mean_SOC(ii) = NaN;  data.std_SOC(ii) = NaN; end
    try data.mean_AID(ii) = rundata.mean_AID; catch data.mean_AID(ii) = NaN; end
    try data.std_AID(ii) = rundata.std_AID; catch data.std_AID(ii) = NaN; end
    try data.mean_inj_delay(ii) = rundata.mean_inj_delay; data.std_inj_delay(ii) = rundata.std_inj_delay; catch data.mean_inj_delay(ii) = NaN; data.std_inj_delay(ii) = NaN; end
    try data.Tunc(ii) = rundata.std_pmax/rundata.mean_pmax*rundata.Ttdc; catch data.Tunc(ii) = NaN; end %uncertainty propagation in P_TDC to T_TDC
    try eval(['data.' plotlegend '(ii) = rundata.' plotlegend ';']); catch eval(['data.' plotlegend '(ii) =  NaN;']); end
end
ivalid = find(isfinite(data.T)>0);
if isempty(colors) %if colors is an empty string
    colors = [rand() rand() rand()];
else
end

Tuncpercentage = (data.Tunc(ivalid)./data.T(ivalid))
figure(fig); hold on
p1 = subplot(2,3,1); hold on; grid on;
  h = errorbar(data.dse,data.mean_SOC,data.std_SOC,data.std_SOC,mark);
        
if ~isempty(plotlegend)
    if legend
        [~,~,OUTH,OUTM] = legend;
        set(h,'color',colors,'markersize',5,'linewidth',1,'markerfacecolor','w');
        legdata = mean(eval(['data.' plotlegend '(ivalid)']));
        h1 = legend([OUTH; h],OUTM,[legendtext '= ' num2str(round(legdata/digits)*digits) ' ' units]);
        set(h1,'position',[0.4035 0.3198 0.1067 0.0663])
      else
        set(h,'color',colors,'markersize',5,'linewidth',1,'markerfacecolor','w');
        legdata = mean(eval(['data.' plotlegend '(ivalid)']));
        h1 = legend([legendtext '= ' num2str(round(legdata/digits)*digits)  ' ' units]);
        set(h1,'position',[0.4035 0.3198 0.1067 0.0663])
    end
end
xlabel('duration of energizing [ms]')    
ylabel('start of combustion [CAD]')
axis([0 4 -5 15])

p5 = subplot(2,3,2); hold on; grid on;
    h2 = errorbar(data.dse,data.mean_AID,data.std_AID,data.std_AID,mark);
    set(h2,'color',colors,'markersize',5,'linewidth',1,'markerfacecolor','w');
xlabel('duration of energizing [ms]')    
ylabel('auto ignition delay [ms]')
axis([0 4 0 4])

p3 = subplot(2,3,3); hold on; grid on;
    h3 = errorbar(data.dse,data.mean_inj_delay,data.std_inj_delay,data.std_inj_delay,mark)
    set(h3,'color',colors,'markersize',5,'linewidth',1,'markerfacecolor','w');
xlabel('duration of energizing [ms]')
ylabel('hydraulic delay [CAD]')
axis([0 4 0 4])

p4=subplot(2,3,4); hold on; grid on;
    h4 = herrorbar(1000./data.T,data.mean_AID,Tuncpercentage.*1000./data.T,Tuncpercentage.*1000./data.T,mark);
    axis([.5 3 .25 3])
    set(h4,'Color',[.2 .2 .2],'markeredgecolor',colors,'markersize',5,'linewidth',1,'markerfacecolor','w');
    set(p4,'YScale','log','Ylim',[0.25 2.5],'YTick',[0.3 0.5 1 1.5 2.5],'Xlim',[0.9 1.4],'Xtick',[0.9 1 1.1 1.2 1.3 1.4])
xlabel('1000/T [1/K]')
ylabel('auto ignition delay [ms]')
%set(he,'Yaxis')

p5 = subplot(2,3,6); hold on; grid on;
    h5 = errorbar(1000./data.T,data.mean_AID,data.std_AID);
    axis([.5 1.5 .25 3])
    set(h5,'Linestyle','None','Marker',mark,'Color',[.2 .2 .2],'markeredgecolor',colors,'markersize',5,'linewidth',1,'markerfacecolor','w');
    set(p5,'YScale','log','Ylim',[0.25 2.5],'YTick',[0.3 0.5 1 1.5 2.5],'Xlim',[0.9 1.4],'Xtick',[0.9 1 1.1 1.2 1.3 1.4])
xlabel('1000/T [1/K]')
ylabel('auto ignition delay [ms]')
%set(he,'Yaxis')
end