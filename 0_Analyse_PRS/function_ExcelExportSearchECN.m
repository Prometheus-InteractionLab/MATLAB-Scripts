function function_ExcelExportSearchECN()
%Function that outputs information for online search tool to read and sort by.
%Version 1: Only outputs file name and the project name, Ian Liu summer 2015
%v1.1 : includes RMP, TDC Temp, density, rail pressure, O2, DSE, SSE, Tin, Pin, Rc, fired#.

[Filename,Pathname] = uigetfile('*.mat', 'Select .mat File', 'MultiSelect', 'on');
addpath(Pathname);

Pathname2 = uigetdir(pwd,'Where is ECNExcelSearch.xls that you would like to add data to?')
addpath([Pathname2])
cd(Pathname2)
%assignin('base', 'A', A); % For development purposes only
xlsFilename='LabBookRunData.xls'; % this file must be created ahead of time

if exist(xlsFilename)
    [CurData, text] = xlsread([xlsFilename],1);
    last = size(text,1)+1
else
headers={'PostInjection Soot/LII Study' 'set RPM' 'T_TDC [K]' 'rho_TDC [kg/m3]' 'Prail [bar]' '%O2' 'DSE [ms]' ...
    'SSE [CAD]' 'Tin [K]' 'Pin [bar]' 'gIMEP [kPa]' 'COV IMEP' 'Rc' 'P_max [bar]' 'PAHRR [J/deg]' ...
    'PAHRR_filt[J/deg]' 'PPRR[bar/deg]' 'PPRR_filt[bar/deg]' ...
    'CA10' 'stdCA10' 'CA50' 'stdCA50' 'CA90' 'stdCA90' 'NOx_max' 'NOx_avg' 'Num fired cycles' 'data' 'view processed data'};
xlswrite(xlsFilename, headers, 1,'A1'); 
last = 1;
end

if iscell(Filename)==0
    Filename=cellstr(Filename);
end

for i=1:length(Filename)
    A=importdata(Filename{i});
    runname = Filename{i}(1:end-4);
    datalink=strcat('<a href="/engine_data/', Filename{i}, '">', [runname '.mat'], '</a>')
    imagelink=strcat('<a href="/engine_data/', [runname '.emf'], '">', [runname '.emf'], '</a>');
    try A.maxNOx; catch A.maxNOx = NaN; end
    try A.avgNOx; catch A.avgNOx = NaN; end
    dataToExport={'AOP 2016 RCCI PPRR' A.Sset round(A.Ttdc/10)*10 round(A.rhotdc*10)/10 round(A.mean_prail/10)*10 ...
        round(A.yo2*10)/10 round(A.dse_ms*100)/100 round(A.sse_degV*10)/10 ...
        round(A.TBDC*10)/10 round(A.mean_pint_bdc*100)/100 round(A.mean_gIMEP)/1000 100*round(A.std_gIMEP/A.mean_gIMEP*1000)/1000 ...
        round(A.Rc*100)/100 round(A.mean_pmax*10)/10 round(A.mean_AHRRmax*10)/10 round(A.mean_AHRRmaxfilt*10)/10 ...
        round(max(max(A.dpdCA/1e5))*10)/10 round(max(max(A.dpdCAfilt/1e5))*10)/10 A.mean_CA10 A.std_CA10 ...
        A.mean_CA50 A.std_CA50 A.mean_CA90 A.std_CA90 round(A.maxNOx*100)/100 round(A.avgNOx*100)/100 A.NC/2 datalink imagelink}; %Can be expanded
    %dataToExport={Filename(1:end-5) A.Ttdc A.rhotdc A.fuel_type Filename};% On hold for now
    nextRow=strcat('A', num2str(last+1));
    xlswrite(xlsFilename, dataToExport ,1, nextRow);
    last = last+1;
end
fullName=xlsFilename;
winopen(fullName);
