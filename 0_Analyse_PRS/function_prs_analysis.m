function [headers,dataxls,run_legend,prsdata] = function_prs_analysis

%  This progam replaces the previous pressure.m. It computes all the
%  important caracteristics of the recorded cycle, and adds the possbility
%  to compute TDC conditions.
%
%             Name:  function_prs_analyis.m
%       Written by:  L.M. MALBEC and E. EAGLE
%             Date:  140925
% 141017 - L.M. MALBEC
%   - Addition of the computation of compression and expansion polytropic
%     coefficients
%   - Change in the processing of the thermodynamic data to eliminate
%     "segmentation violation" error.
%  
% 141017 - W.E. EAGLE
%   - Possibility to read intake temperature from a mat file instead of
%     manually entering it.
%   - Correction of an error when running with pure Nitrogen.
%   - Addition of a text box summarizing inlet conditions

% 141020 - L.M. MALBEC
%   - Multiple injections taken into account

% 141031 - L.M. MALBEC
%   - Read automaticaly csv files from flow computer, and computes intake
%   conditions

% 141204 - E. Eagle
%   - included info box with runname, fuel type, injector, tip nozzle # and
%   diameter and box for 'other available data, i.e. *.cine *.spe

% 141205 - E. Eagle
%   - moved dayname outside of loop so that other exist() checks are 
%   possible for automating data input 

headers = [];

%close all
exportplot = 1; %1 for saving plot
showTDCvariationfired = 0; %for an extra plots that display P,T,rho from -10 to 40 CAD for fired cases
showTDCvariationnormalized = 0; %for an extra plots that display P,T,rho from -10 to 40 CAD normalized by motored P,T,rho
showPSD = 0; % shows plot of power spectral density for both unfiltered and filtered pressure.
    set(0,'defaultAxesFontName','Times')
    set(0,'defaultTextFontName','Times')
%% Set Defaults    
% Volumic fraction of O2 in %
inj_name = ''; % three letters
inj_tip_name = ''; % up to 5 numbers
inj_diam = 0.2; %mm
hole_number = 1; % # of holes
fuel_type = '';

%%Get files to process, check for automated input files
[filename,pathname] = uigetfile('*.prs','Select prs files','Multiselect','on');
filename = cellstr(filename);
addpath(pathname) %added by E. Eagle 141009, prevents code crashing when selected file isn't already in the matlab path

dayname = filename{1}(1:6)

% added search for csv and matlab infile
% If 1, manual input for TDC conditions
% If 2, automatic input for TDC conditions
% If 3, matlab input for TDC conditions

if exist([dayname '.xlsx'])
    itdc = 2;
    inp = [dayname '.xlsx'];
% elseif exist([dayname '.csv'])
%     itdc = 2;
%     inp = [dayname '.csv'];
elseif exist([dayname 'infile.m'])
    itdc = 3;
    inp = [dayname 'infile.m'];
else
    itdc = 1;
    inp = 'manual input';
end
    
if itdc==3
    [filename2,pathname2] = uigetfile('*.mat','Select Input Condition mat file [datename.mat]');
    addpath(pathname2)
    load(filename2)
end

dataxls = [];
dataxls_tdc = [];
run_legend={};
for ii=1:length(filename)
    if length(filename)>2
    close all % if plotting many runs, supress output between cases
    end
    if length(filename{ii})>11
        runname = filename{ii}([end-5 end-4])
    else
        runname = filename{ii}(end-4);
    end
%     addpath(fullfile('\\906-165-image\e$\2014',dayname)) %addpath to image directory
%     addpath(fullfile('\\906-165-image\e$\2015',dayname)) %addpath to image directory

    
    HSimages = [dayname runname '.cine'];
    ICMOS = [dayname runname '.spe'];

    i1 = 70;i2 = 90; %(30 seconds yields 30 cycles)
try
%% Read and plot intake conditions
    if exist(fullfile(pathname,[dayname '.xlsx']))
        [data, text, raw] = xlsread(fullfile(pathname,[dayname '.xlsx']));
        xlsdata.data = data;
        xlsdata.textdata = text;
        [csvdata,i2] = function_read_xlsdata(xlsdata,dayname,runname,i1,i2);
        figure,set(gcf,'color','w','name',fullfile(pathname,[filename{ii} '_intake']),'Position',[85 160 1456 754])
        function_plot_csvdata(csvdata,i1,i2); %read xlsdata converts into csv compatable structure
        saveas(gcf,fullfile(pathname,['intake_' dayname runname '.fig']),'fig')
        saveas(gcf,fullfile(pathname,['intake_' dayname runname '.emf']),'emf') 
    elseif exist(fullfile(pathname,[dayname '.csv']))
         listcsv = dir(fullfile(pathname,[dayname '.csv']));
         csvdata = importdata(fullfile(pathname,listcsv.name));
         [csvdata,i2] = function_read_csvdata(csvdata,dayname,runname,i1,i2);
         figure,set(gcf,'color','w','name',fullfile(pathname,[filename{ii} '_intake']),'Position',[85 160 1456 754])
         function_plot_csvdata(csvdata,i1,i2);
         saveas(gcf,fullfile(pathname,['intake_' dayname runname '.fig']),'fig')
         saveas(gcf,fullfile(pathname,['intake_' dayname runname '.emf']),'emf') 
    end
catch
   csvdata =[];
end
    
        
%% Read and plot prs data
    [data,N,SPC] = function_read_prsfile(pathname,filename{ii});
    prsdata = function_read_prsdata(data,N,SPC);
    prsdata.inj_name = inj_name;
    prsdata.inj_tip_name = inj_tip_name ; % up to 5 numbers
    prsdata.inj_diam = inj_diam ; %mm
    prsdata.hole_number = hole_number; % # of holes
    prsdata.fuel_type = fuel_type;
    prsdata.runname  = runname;
    %prsdata.inj_delay = inj_delay; % delay in CAD between SSE and SOI (injector, tip, and pressure dependent parameter)

    %% plot data log figure
    figure('Units','inches',...
        'position',[1 1 10 7.5],...
        'PaperPositionMode','auto');
    set(gcf,'color','w','name',fullfile(pathname,[filename{ii} ' !!!!!!!!!!!!!!!!!!!!!!!!  Rc=' num2str(prsdata.Rc,'%2.2f') '  !!!!!!!!!!!!!!!!!!!!']));
    prsdata = function_plot_prsdata_run2(prsdata);
    pause(0.1)    
    
    if itdc==1
         prompt = {'Intake temperature (degC):','O2 concentration (% vol.):'};
         dlg_title = 'Input for TDC conditions';
         num_lines = 1;
         answer = inputdlg(prompt,dlg_title,num_lines,{'';''});
         csvdata.mean_tintake = str2num(answer{1});
         yo2 = str2num(answer{2});
    elseif itdc==3
        index = findindex(runname,char(runs))
        csvdata.mean_tintake = Ti(index);
        yo2 = P_O2;
        csvdata.mean_pintake = Pin(index)*100;
        csvdata.mean_qair = mdot(index);
        csvdata.mean_dse = DSE(index);
        csvdata.mean_rpm = 1200;
        prsdata.Rc = compRatio;
    elseif itdc==2
        try yo2 = csvdata.po2; catch yo2 = 21;end
    end

            
    %%%%%%%%%%%%%%%%%%%% Compute TDC conditions
    if itdc ==3
     pco2=yo2(ii)/100;
    else
     pco2=yo2/100;
    end
    
    Pbdc = prsdata.mean_pint_bdc;
    Ptdc = prsdata.mean_pmotor;

    % Might not work with older versionsof matlab
    N2_properties_nist;
    O2_properties_nist;
    AR_properties_nist;
    AIR_properties_nist;

    % Estimation of the temperature at BDC knowing the intake temperature
    try
        Tbdc = function_Tintake_2_Tbdc(csvdata.mean_tintake)+273.15
    catch
        Tbdc = NaN;
    end
    
    % Fresh gas constant (dilution with N2)
    Mair = 8.314/r_air*1000;
    ndil = (0.21/pco2-1); % Nb of N2 moles for 1 mole or air
    ndilprint = ndil;
    nairprint = 1;
    xn2 = ndil*28/(ndil*28+Mair); % Mass fraction of N2 as diluant
    xair = 1-xn2; % Mass fraction of air
    yn2 = ndil/(ndil+1);
    yair = 1/(ndil+1);
    if isinf(ndil)
        ndil = 1e10;
        xn2 = 1; % Mass fraction of N2 as diluant
        xair = 0; % Mass fraction of air
        ndilprint = 1;
        nairprint = 0;
        yn2= 1;
        yo2=0;
    end

    Mgf = (ndil*28+Mair)/(1+ndil);
    r_fg = 8.314/Mgf*1000;
    cpmass_fg = @(x) xn2*cpmass_n2_nist(x)+xair*cpmass_air_nist(x); % cp as a function of temperature fot the fresh gasses
    smass_fg = @(x) xn2*smass_n2_nist(x)+xair*smass_air_nist(x);    % Entropy as a function of temperature fot the fresh gasses, at atmospheric pressure

    for jj = 1:length(prsdata.ca_rng)
         
        if ~isnan(Tbdc)
        smass_ref = smass_fg(Tbdc) - r_fg*log(Pbdc); % Value of the entropy of fresh gasses at BDC

        % Loop on TDC temperature: when the entropy at TDC is equal to the entropy
        % at BDC, the correct value of Ttdc has been found.
        Tmrng(jj) = Tbdc;
        smass_tdc = smass_fg(Tbdc) - r_fg*log(prsdata.pmotor_avg(jj));
            while abs(smass_tdc-smass_ref)>0.1
            deltaT = (smass_ref-smass_tdc)*Tmrng(jj)/cpmass_air_nist(Tmrng(jj));
            Tmrng(jj) = Tmrng(jj) + deltaT;
            smass_tdc = smass_fg(Tmrng(jj)) - r_fg*log(prsdata.pmotor_avg(jj));
            end
        % Compressibility coefficients
        z_rngn2 = function_compressibility(prsdata.pmotor_avg(jj),Tmrng(jj),'n2');
        z_rngair = function_compressibility(prsdata.pmotor_avg(jj),Tmrng(jj),'air');
        z_rngfg_tdc = yair*z_rngair + yn2*z_rngn2;
        rhom_rng(jj)=prsdata.pmotor_avg(jj)/z_rngfg_tdc/r_fg/Tmrng(jj)*1e5;

        Tfrng(jj) = Tbdc;
        smass2_tdc = smass_fg(Tbdc) - r_fg*log(prsdata.pfire_avg(jj));
            while abs(smass2_tdc-smass_ref)>0.01
            deltaT = (smass_ref-smass2_tdc)*Tfrng(jj)/cpmass_air_nist(Tfrng(jj));
            Tfrng(jj) = Tfrng(jj) + deltaT;
            smass2_tdc = smass_fg(Tfrng(jj)) - r_fg*log(prsdata.pfire_avg(jj));
            end
        % Compressibility coefficients
        zf_rngn2 = function_compressibility(prsdata.pfire_avg(jj),Tfrng(jj),'n2');
        zf_rngair = function_compressibility(prsdata.pfire_avg(jj),Tfrng(jj),'air');
        zf_rngfg_tdc = yair*zf_rngair + yn2*zf_rngn2;
        rhof_rng(jj)=prsdata.pfire_avg(jj)/zf_rngfg_tdc/r_fg/Tfrng(jj)*1e5;
        

        
        zn2 = function_compressibility(Pbdc,Tbdc,'n2');
        zair = function_compressibility(Pbdc,Tbdc,'air');
        zfg_bdc = yair*zair + yn2*zn2;
        rho_bdc=Pbdc/zfg_bdc/r_fg/Tbdc*1e5;
        else
            Ttdc = NaN;
            Tmrng = NaN;
            Tfrng = NaN;
            rho_tdc=NaN;
            rho_bdc=NaN;
            rhom_rng=NaN;
            smass_ref=NaN;
            smass_tdc=NaN;
            zfg_tdc=NaN;zfg_bdc=NaN;
        end
    
    end
%%     
    if ~isnan(Tbdc)
        smass_ref = smass_fg(Tbdc) - r_fg*log(Pbdc); % Value of the entropy of fresh gasses at BDC

        % Loop on TDC temperature: when the entropy at TDC is equal to the entropy
        % at BDC, the correct value of Ttdc has been found.
        Ttdc = Tbdc;
        smass_tdc = smass_fg(Tbdc) - r_fg*log(Ptdc);
        while abs(smass_tdc-smass_ref)>0.01
            deltaT = (smass_ref-smass_tdc)*Ttdc/cpmass_air_nist(Ttdc);
            Ttdc = Ttdc + deltaT;
            smass_tdc = smass_fg(Ttdc) - r_fg*log(Ptdc);
        end

        % Compressibility coefficients
        zn2 = function_compressibility(Ptdc,Ttdc,'n2');
        zair = function_compressibility(Ptdc,Ttdc,'air');
        zfg_tdc = yair*zair + yn2*zn2;
        rho_tdc=Ptdc/zfg_tdc/r_fg/Ttdc*1e5;

        zn2 = function_compressibility(Pbdc,Tbdc,'n2');
        zair = function_compressibility(Pbdc,Tbdc,'air');
        zfg_bdc = yair*zair + yn2*zn2;
        rho_bdc=Pbdc/zfg_bdc/r_fg/Tbdc*1e5;
    else
        Ttdc = NaN;
        rho_tdc=NaN;
        rho_bdc=NaN;
        smass_ref=NaN;
        smass_tdc=NaN;
        zfg_tdc=NaN;zfg_bdc=NaN;
    end
    
    
    try pinlet      = csvdata.mean_pintake/100; catch pinlet = NaN; end
    try tinlet      = csvdata.mean_tintake; catch tinlet = NaN; end
    try qinlet      = csvdata.mean_qair; catch qinlet = NaN; end
    try dseinlet    = prsdata.dse_ms*1e3; catch dseinlet = NaN; end
    try sseinlet    = 360+prsdata.sse_degV; catch sseinlet = NaN; end
    try prailinlet  = mean(mean(prsdata.prail(:,10:50))); catch prailinlet = NaN; end
    try hydelay = prsdata.mean_inj_delay; catch hydelay = NaN; end
    try stdhydelay = prsdata.std_inj_delay; catch stdhydelay = NaN; end
    try compRatio   = prsdata.Rc;            catch compRatio = NaN; end
    try measRPM     = csvdata.mean_rpm;     catch   measRPM = NaN; end
    compRatio
    wid = .14;
    hei = .1625;
    hin = annotation('textbox',[0.6913 0.61 wid hei],'color','k','edgecolor','k','fontname','arial','String',...
        {['Inlet conditions:'];...
%         ['  - P = ' num2str(pinlet,'%1.2f') ' bar'];...
        [' T = ' num2str(tinlet,'%2.1f') ' ^\circ C'];...
        [' Qair = '  num2str(qinlet,'%2.1f') ' g/s'];...
        ['Fuel:'];...
        [' ' fuel_type]
        [' Prail  = ' num2str(round(prailinlet/100)*100,'%4.0f') '\pm' num2str(round(mean(std(prsdata.prail(:,90:150)))),'%2.0f') ' bar'];...
        });
    habdc = annotation('textbox',[0.6913+wid+.002 0.61 wid hei],'color','b','edgecolor','b','fontname','arial','String',...
        {['BDC conditions:'];...
        [' P = ' num2str(Pbdc,'%1.2f') ' bar'];...
        [' T = ' num2str(Tbdc-273.15,'%2.1f') '^\circ C'];... (Tintake = ' num2str(Tintake,'%2.1f') ' degC)'];...
        [' \rho = ' num2str(rho_bdc,'%2.2f') ' kg/m^3'];...
        [' Z = ' num2str(zfg_bdc,'%1.3f')];...
        [' S = ' num2str(smass_ref,'%4.1f') ' J/K/kg']});

    hatdc = annotation('textbox',[0.6913+wid+.002 0.61+hei+.002 wid hei],'color','r','edgecolor','r','fontname','arial','String',...
        {['TDC conditions:'];...
        ['	P = ' num2str(Ptdc,'%1.2f') ' bar'];...
        ['	T = ' num2str(Ttdc,'%2.1f') ' K'];...
        [' \rho = ' num2str(rho_tdc,'%2.2f') ' kg/m^3'];...
        ['	Z = ' num2str(zfg_tdc,'%1.3f')];...
        ['	S = ' num2str(smass_tdc,'%4.1f') ' J/K/kg']});

    hrun = annotation('textbox',[0.6913 0.61+hei+.002 wid hei],'color','k','edgecolor','k','fontname','arial','String',...
        {['Settings:'];...
        [' runname = ' dayname runname ];...
        [' Rc = ' num2str(compRatio,'%2.2f') ];...
        [' RPM= ' num2str(measRPM, '%4.0f') ];...
        [' DSE =' num2str(dseinlet,'%4.0f '), ' \mus'];...
        [' SSE  = ' num2str(sseinlet,'%4.1f '), ' CAD']});
    
    hfg = annotation('textbox',[0.6913 0.5235 (wid+.002)*2 0.0828],'color',[0 128 64]/255,'edgecolor',[0 128 64]/255,'fontname','arial','String',...
        {['Ambient gasses:'];...
        [' %O_2 = ' num2str(pco2*100,'%2.1f') ' (' num2str(nairprint,'%1.1f') 'mol Air + ' num2str(ndilprint,'%1.3f') 'mol N2)'];...
        [' M_{air} = ' num2str(Mair,'%3.2f') ' g/mol; M_{fg} = ' num2str(Mgf,'%3.2f') ' g/mol']});
    
    hfg = annotation('textbox',[0.6913 0.61+(hei+0.002)*2 (wid+.002)*2 0.0528],'color',[0 128 64]/255,'edgecolor',[0 128 64]/255,'fontname','arial','String',...
        {['Inj: ' inj_name ' ' inj_tip_name ' hydr. delay ' num2str(hydelay,'%2.1f') ' \pm ' num2str(stdhydelay,'%1.1f') 'CAD'];...
        [ num2str(hole_number),' orifice' ', diameter = ' num2str(inj_diam,'%1.2f') 'mm']});
%     
%     if exist(HSimages) && exist(ICMOS)
%         hfg = annotation('textbox',[0.5913 0.7716+.14 0.1 0.095],'color',[0 128 64]/255,'edgecolor',[0 128 64]/255,'fontname','arial','String',...
%         {['other available data '];...
%         ['  - ' HSimages ];...
%         ['  - ' ICMOS ];...
%         ['  - ' inp];...
%         });
%     elseif exist(ICMOS)
%         hfg = annotation('textbox',[0.5913 0.7716+.14 0.1 0.0828],'color',[0 128 64]/255,'edgecolor',[0 128 64]/255,'fontname','arial','String',...
%         {['other available data '];...
%         ['  - ' ICMOS ];...
%         ['  - ' inp];...
%         });
%     elseif exist(HSimages)
%         hfg = annotation('textbox',[0.5913 0.7716+.14 0.1 0.0828],'color',[0 128 64]/255,'edgecolor','none','fontname','arial','String',...
%         {['other available data '];...
%         ['  - ' HSimages ];...
%         ['  - ' inp];...
%         });
%     end
%     
    if exportplot
    saveas(gcf,fullfile(pathname,[filename{ii}(1:end-4) '.fig']),'fig')

%    try    export_fig( '-nocrop', strcat(filename{ii}(1:end-5), '.pdf'))
%    catch
    saveas(gcf,fullfile(pathname,[filename{ii}(1:end-4) '.emf']),'emf')
%   end
    end
    rundata = prsdata;
    try 
        fields = fieldnames(csvdata);
    catch
        fields = [];
    end
        
    for ii=1:length(fields)
        rundata.(fields{ii}) = csvdata.(fields{ii});
    end
    rundata.mean_prail = prailinlet;
    rundata.Ttdc = Ttdc;
    rundata.rhotdc = rho_tdc;
    rundata.zfg_tdc = zfg_tdc;
    rundata.rhobdc = rho_bdc;
    rundata.yo2 = yo2;
    rundata.rair = r_air;
    rundata.Mair = Mair;
    rundata.rfg = r_fg;
    rundata.Mfg = Mgf;
    rundata.xair = xair;
    rundata.xn2 = xn2;
    rundata.TBDC = Tbdc;
    rundata.mean_inj_hydelay = hydelay;
    rundata.std_inj_hydelay = stdhydelay;
    rundata.Tm_rng = Tmrng;
    rundata.rhom_rng = rhom_rng;
    rundata.Tf_rng = Tfrng;
    rundata.rhof_rng = rhof_rng;
    save(fullfile(pathname,[dayname runname '.mat']),'rundata');
    
    
    if showPSD
        
figure();hold on;
for jj = 2:2:rundata.NC
x = rundata.pcyl(jj,:);
N = length(x);
xdft = fft(x);
xdft2 = xdft(1:N/2+1);
psdx = (1/(2*pi*N)) * abs(xdft2).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:(2*pi)/N:pi;
%figure(4); hold on;
%plot(freq/pi*28800/2,10*log10(psdx),'color',[0.6 0.6 0.6],'linewidth',2)
psdFire(jj,:) = psdx;
end
plot(freq/pi*28800/2,10*log10(mean(psdFire)),'linewidth',2)
alpha = [1.841 3.054 3.832 4.201 5.332]
fmodes = alpha.*sqrt((rundata.nexp+rundata.ncomp/2)*287*rundata.Ttdc)/pi/.1397

plot([fmodes],[0 0 0 0 0],'linestyle','none','marker','h','color','r')

% x = mean(pcyl(2:2:NC,:));
% N = length(x);
% 
% xdft = fft(x);
% xdft2 = xdft(1:N/2+1);
% psdx = (1/(2*pi*N)) * abs(xdft2).^2;
% psdx(2:end-1) = 2*psdx(2:end-1);
% freq = 0:(2*pi)/N:pi;
% figure(4); hold on;
% plot(freq/pi*28800/2,10*log10(psdx),'color',[0.7 0.7 0.7],'linewidth',2)

for kk = 2:2:rundata.NC
    x = rundata.pfilt(kk,:);
    N = length(x);
    xdft = fft(x);
    xdft2 = xdft(1:N/2+1);
    psdx = (1/(2*pi*N)) * abs(xdft2).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:(2*pi)/N:pi;
    
%    figure(4); hold on;
%   plot(freq/pi*28800/2,10*log10(psdx),':','linewidth',2,'color',[0.2 0.2 0.2])
    psdfilt(kk,:) = psdx;
end
    plot(freq/pi*28800/2,10*log10(mean(psdfilt)),'k','linewidth',2)

%     x = mean(pcyl(1:2:NC,:));
%     N = length(x);
% 
% xdft = fft(x);
% xdft2 = xdft(1:N/2+1);
% psdx = (1/(2*pi*N)) * abs(xdft2).^2;
% psdx(2:end-1) = 2*psdx(2:end-1);
% freq = 0:(2*pi)/N:pi;
% figure(4); hold on;
% plot(freq/pi*28800/2,10*log10(psdx),'k--','linewidth',2)

grid on
title('Periodogram Using FFT')
xlabel('Frequency')
ylabel('Power/Frequency (dB/rad/sample)')
axis([0 7500 -80 80])
end    
end
matlabdefaultcolors =  [0    0.4470    0.7410;   0.8500    0.3250    0.0980;    0.9290    0.6940    0.1250;    0.4940    0.1840    0.5560;    0.4660    0.6740    0.1880;    0.3010    0.7450    0.9330;    0.6350    0.0780    0.1840;]

if showTDCvariationfired
        figure() 
        plot(rundata.ca_rng,rundata.pfire_avg,':','color',matlabdefaultcolors(1,:),'lineW',2)
        figure(),plot(rundata.ca_rng,Tfrng,':','color', matlabdefaultcolors(2,:),'lineW',2)
        figure(),plot(rundata.ca_rng,rhof_rng,':','color',matlabdefaultcolors(3,:),'lineW',2)
        axis('square')
        legend('Motored Pressure','MTemperature','MDensity','Fired Pressure','FTemperature','FDensity')
        xlabel('Crank Angle')
        ylabel('Normalized value to TDC')
        
        figure(), plot(rundata.ca_rng(1:end-1),diff(rundata.pmotor_avg/Ptdc)/(rundata.ca_rng(2)-rundata.ca_rng(1))*100,'lineW',2)
        hold on , plot(rundata.ca_rng(1:end-1),diff(Tmrng/Ttdc)/(rundata.ca_rng(2)-rundata.ca_rng(1))*100,'lineW',2)
        plot(rundata.ca_rng(1:end-1),diff(rhom_rng/rho_tdc)/(rundata.ca_rng(2)-rundata.ca_rng(1))*100,'lineW',2)
        plot(rundata.ca_rng(1:end-1),diff(rundata.pfire_avg/Ptdc)/(rundata.ca_rng(2)-rundata.ca_rng(1))*100,':','color',matlabdefaultcolors(1,:),'lineW',2)
        hold on , plot(rundata.ca_rng(1:end-1),diff(Tfrng/Ttdc)/(rundata.ca_rng(2)-rundata.ca_rng(1))*100,':','color',matlabdefaultcolors(2,:),'lineW',2)
        plot(rundata.ca_rng(1:end-1),diff(rhof_rng/rho_tdc)/(rundata.ca_rng(2)-rundata.ca_rng(1))*100,':','color', matlabdefaultcolors(3,:),'lineW',2)
        ticksperCAD=4;% change if encoder increases in resolution be
        for ii=1:length(rundata.sse_degV)
        plot(rundata.ca(1,rundata.sse_index(ii)+rundata.mean_inj_delay*ticksperCAD:rundata.ese_index(ii)+rundata.mean_inj_delay*ticksperCAD),...
        ones(length(rundata.ca(1,rundata.sse_index(ii)+rundata.mean_inj_delay*ticksperCAD:rundata.ese_index(ii)+rundata.mean_inj_delay*ticksperCAD)),1),'linewidth',3,'color','r')
    % accounts for hydraulic delay in plotting DSE, still missing extra
    % closing time (injector specific)
end
        axis('square')
        legend('Pressure','Temperature','Density','Fired P','Fired T','Fired \rho')
        xlabel('Crank Angle')
        ylabel('% change per CAD, normalized by TDC')
end

    if showTDCvariationnormalized
        figure(), plot(rundata.ca_rng,rundata.pmotor_avg/Ptdc,'lineW',2)
        hold on, plot(rundata.ca_rng,Tmrng/Ttdc,'lineW',2)
        plot(rundata.ca_rng,rhom_rng/rho_tdc,'lineW',2)
        plot(rundata.ca_rng,rundata.pfire_avg/Ptdc,':','color',matlabdefaultcolors(1,:),'lineW',2)
        plot(rundata.ca_rng,Tfrng/Ttdc,':','color', matlabdefaultcolors(2,:),'lineW',2)
        plot(rundata.ca_rng,rhof_rng/rho_tdc,':','color',matlabdefaultcolors(3,:),'lineW',2)
        axis('square')
        legend('Motored Pressure','MTemperature','MDensity','Fired Pressure','FTemperature','FDensity')
        xlabel('Crank Angle')
        ylabel('Normalized value to TDC')
        
        figure(), plot(rundata.ca_rng(1:end-1),diff(rundata.pmotor_avg/Ptdc)/(rundata.ca_rng(2)-rundata.ca_rng(1))*100,'lineW',2)
        hold on , plot(rundata.ca_rng(1:end-1),diff(Tmrng/Ttdc)/(rundata.ca_rng(2)-rundata.ca_rng(1))*100,'lineW',2)
        plot(rundata.ca_rng(1:end-1),diff(rhom_rng/rho_tdc)/(rundata.ca_rng(2)-rundata.ca_rng(1))*100,'lineW',2)
        plot(rundata.ca_rng(1:end-1),diff(rundata.pfire_avg/Ptdc)/(rundata.ca_rng(2)-rundata.ca_rng(1))*100,':','color',matlabdefaultcolors(1,:),'lineW',2)
        hold on , plot(rundata.ca_rng(1:end-1),diff(Tfrng/Ttdc)/(rundata.ca_rng(2)-rundata.ca_rng(1))*100,':','color',matlabdefaultcolors(2,:),'lineW',2)
        plot(rundata.ca_rng(1:end-1),diff(rhof_rng/rho_tdc)/(rundata.ca_rng(2)-rundata.ca_rng(1))*100,':','color', matlabdefaultcolors(3,:),'lineW',2)
        ticksperCAD=4;% change if encoder increases in resolution be
        for ii=1:length(rundata.sse_degV)
        plot(rundata.ca(1,rundata.sse_index(ii)+rundata.mean_inj_delay*ticksperCAD:rundata.ese_index(ii)+rundata.mean_inj_delay*ticksperCAD),...
        ones(length(rundata.ca(1,rundata.sse_index(ii)+rundata.mean_inj_delay*ticksperCAD:rundata.ese_index(ii)+rundata.mean_inj_delay*ticksperCAD)),1),'linewidth',3,'color','r')
    % accounts for hydraulic delay in plotting DSE, still missing extra
    % closing time (injector specific)
end
        axis('square')
        legend('Pressure','Temperature','Density','Fired P','Fired T','Fired \rho')
        xlabel('Crank Angle')
        ylabel('% change per CAD, normalized by TDC')
    end