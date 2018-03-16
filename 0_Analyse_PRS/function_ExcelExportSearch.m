function function_ExcelExportSearch()
%Function that outputs information for online search tool to read and sort by.
%Current state: Only outputs file name and the Tardec project name

[Filename,Pathname] = uigetfile('*.mat', 'Select .mat File', 'MultiSelect', 'on');
addpath(Pathname);

%assignin('base', 'A', A); % For development purposes only
xlsFilename='E:\Data\ExcelSearch';
headers={'Project' 'File name'};
xlswrite(xlsFilename, headers, 2,'A1'); 
if iscell(Filename)==0
    Filename=cellstr(Filename);
end
for i=1:length(Filename)
    A=importdata(Filename{i});
    linkName=strcat('<a href="/engine_data/', Filename{i}, '">', Filename{i}, '</a>');
    dataToExport={'TARDEC' linkName}; %Used just for Tardec, needs to be updated
    %dataToExport={Filename(1:end-5) A.Ttdc A.rhotdc A.fuel_type Filename};% On hold for now
    nextRowNum=num2str(xlsread(xlsFilename,3,'A1:A1'));
    nextRow=strcat('A', nextRowNum);
    xlswrite(xlsFilename, dataToExport ,2, nextRow);
    Increment=str2num(nextRowNum)+1;
    xlswrite(xlsFilename, Increment,3,'A1');
    i=i+1;
end
fullName=strcat(xlsFilename, '.xls');
winopen(fullName);