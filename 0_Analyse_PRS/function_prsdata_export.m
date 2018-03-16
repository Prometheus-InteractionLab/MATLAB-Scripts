function function_prsdata_export()

% written by Ian Liu, 7/7/2015
% takes data structure from .mat file and exports it to an excel file, raw data must first be processed with function_prs_analysis.
% Requires manual input of number of variables being exports into variable "var_num"

%Updated 7/21/2015 by Ian Liu
%Independent: No longer needs excel_strings function
%Concatenates cell array for single export, increased speed.
var_num =6;
[Filename,Pathname] = uigetfile('*.mat', 'Select .mat File', 'MultiSelect', 'on');
if iscell(Filename)==0
    Filename=cellstr(Filename);
end

for iiii=1:length(Filename)
    addpath(Pathname);
    A=importdata(Filename{iiii});
    xlsFilename=strcat(Pathname,'run_',Filename{iiii}(1:end-4));
    cycleNum=[1:A.NC];
    assignin('base', 'A', A); %Development purpose only
    wrksht1=cell(length(A.ca)+2,4+(3+A.NC/2)*var_num); %Creates cell array to format before exporting to excel (increases speed)
    i=0;
%%%%%%%%%%%%%%%%%%%%%%%%%% Motored and Fired Worksheets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Seperates motored cycled from fired cycles
    motored_pcyl=[]; fired_pcyl=[]; motored_inj_curr=[]; fired_inj_curr=[]; motored_prail=[]; fired_prail=[]; %Resets matrices for multiple file selections
    motored_pfilt=[]; fired_pfilt=[]; motored_dQdCA=[]; fired_dQdCA=[]; motored_dQdCAfilt=[]; fired_dQdCAfilt=[];
    motored_pcyl(1:1:length(cycleNum)/2,1:2880)=A.pcyl(1:2:length(cycleNum),1:2880); 
    fired_pcyl(1:1:length(cycleNum)/2,1:2880)=A.pcyl(2:2:length(cycleNum),1:2880);
    motored_inj_curr(1:1:length(cycleNum)/2,1:2880)=A.injCurr(1:2:length(cycleNum),1:2880);
    fired_inj_curr(1:1:length(cycleNum)/2,1:2880)=A.injCurr(2:2:length(cycleNum),1:2880);
    motored_prail(1:1:length(cycleNum)/2,1:2880)=A.prail(1:2:length(cycleNum),1:2880);
    fired_prail(1:1:length(cycleNum)/2,1:2880)=A.prail(2:2:length(cycleNum),1:2880);
    motored_pfilt(1:1:length(cycleNum)/2,1:2880)=A.pfilt(1:2:length(cycleNum),1:2880);
    fired_pfilt(1:1:length(cycleNum)/2,1:2880)=A.pfilt(2:2:length(cycleNum),1:2880);
    motored_dQdCA(1:1:length(cycleNum)/2,1:2880)=A.dQdCA(1:2:length(cycleNum),1:2880);
    fired_dQdCA(1:1:length(cycleNum)/2,1:2880)=A.dQdCA(2:2:length(cycleNum),1:2880);
    motored_dQdCAfilt(1:1:length(cycleNum)/2,1:2880)=A.dQdCAfilt(1:2:length(cycleNum),1:2880);
    fired_dQdCAfilt(1:1:length(cycleNum)/2,1:2880)=A.dQdCAfilt(2:2:length(cycleNum),1:2880);
 
%Calculates the mean and std_dev for each crank angle amoung the cycles
%Exports to a function that rounds each value
%Round function could be removed with a newer version of matlab and ability to use round function capabilities

    [A.pcyl_meanfired(:,1), A.pcyl_stdfired(:,1)]=function_round( mean(A.pcyl([1:2:A.NC],:))', std(A.pcyl([1:2:A.NC],:))');
    [A.pcyl_meanfired(:,2), A.pcyl_stdfired(:,2)]=function_round( mean(A.pcyl([2:2:A.NC],:))', std(A.pcyl([2:2:A.NC],:))');
    [A.inj_curr_meanfired(:,1), A.inj_curr_stdfired(:,1)]=function_round( mean(A.injCurr([1:2:A.NC],:))', std(A.injCurr([1:2:A.NC],:))');
    [A.inj_curr_meanfired(:,2), A.inj_curr_stdfired(:,2)]=function_round( mean(A.injCurr([2:2:A.NC],:))', std(A.injCurr([2:2:A.NC],:))');
    [A.prail_meanfired(:,1), A.prail_stdfired(:,1)]=function_round( mean(A.prail([1:2:A.NC],:))', std(A.prail([1:2:A.NC],:))');
    [A.prail_meanfired(:,2), A.prail_stdfired(:,2)]=function_round( mean(A.prail([2:2:A.NC],:))', std(A.prail([2:2:A.NC],:))');
    [A.pfilt_meanfired(:,1), A.pfilt_stdfired(:,1)]=function_round( mean(A.pfilt([1:2:A.NC],:))', std(A.pfilt([1:2:A.NC],:))');
    [A.pfilt_meanfired(:,2), A.pfilt_stdfired(:,2)]=function_round( mean(A.pfilt([2:2:A.NC],:))', std(A.pfilt([2:2:A.NC],:))');
    [A.dQdCA_meanfired(:,1), A.dQdCA_stdfired(:,1)]=function_round( mean(A.dQdCA([1:2:A.NC],:))', std(A.dQdCA([1:2:A.NC],:))');
    [A.dQdCA_meanfired(:,2), A.dQdCA_stdfired(:,2)]=function_round( mean(A.dQdCA([2:2:A.NC],:))', std(A.dQdCA([2:2:A.NC],:))');
    [A.dQdCAfilt_meanfired(:,1), A.dQdCAfilt_stdfired(:,1)]=function_round( mean(A.dQdCAfilt([1:2:A.NC],:))', std(A.dQdCAfilt([1:2:A.NC],:))');
    [A.dQdCAfilt_meanfired(:,2), A.dQdCAfilt_stdfired(:,2)]=function_round( mean(A.dQdCAfilt([2:2:A.NC],:))', std(A.dQdCAfilt([2:2:A.NC],:))');
    
    
    %Combines mean, std, and data into one matrix and reshapes matrix for excel format
    crankAngle=round(4*A.ca(1,:)')/4; % Round to nearest quarter crank angle
    volume=round(10000000*A.volume(1,:)')/10000000; % Rounds to a 10 millionth 
    motored_pcyl=[A.pcyl_meanfired(:,1) A.pcyl_stdfired(:,1) motored_pcyl']; 
    fired_pcyl=[A.pcyl_meanfired(:,2) A.pcyl_stdfired(:,2) fired_pcyl'];
    motored_inj_curr=[A.inj_curr_meanfired(:,1) A.inj_curr_stdfired(:,1) motored_inj_curr']; 
    fired_inj_curr=[A.inj_curr_meanfired(:,2) A.inj_curr_stdfired(:,2) fired_inj_curr']; 
    motored_prail=[A.prail_meanfired(:,1) A.prail_stdfired(:,1) motored_prail']; 
    fired_prail=[A.prail_meanfired(:,2) A.prail_stdfired(:,2) fired_prail']; 
    motored_pfilt=[A.pfilt_meanfired(:,1) A.pfilt_stdfired(:,1) motored_pfilt']; 
    fired_pfilt=[A.pfilt_meanfired(:,2) A.pfilt_stdfired(:,2) fired_pfilt']; 
    motored_dQdCA=[A.dQdCA_meanfired(:,1) A.dQdCA_stdfired(:,1) motored_dQdCA']; 
    fired_dQdCA=[A.dQdCA_meanfired(:,2) A.dQdCA_stdfired(:,2) fired_dQdCA']; 
    motored_dQdCAfilt=[A.dQdCAfilt_meanfired(:,1) A.dQdCAfilt_stdfired(:,1) motored_dQdCAfilt']; 
    fired_dQdCAfilt=[A.dQdCAfilt_meanfired(:,2) A.dQdCAfilt_stdfired(:,2) fired_dQdCAfilt'];


    %Creates 1st row format for cell array to be exported to excel
    temp={' ' 'pcyl' '[bar]' 'inj_curr' '[]' 'prail' '[bar]' 'pfilt' '[bar]' 'dQdCA' '[J/Deg]' 'dQdCAfilt' '[J/Deg]' 'Crank Angle' 'Cycle Number' ' ' ' '};
    i=0;
    while i<var_num %Loop inserts 1st row format into cell array
        wrksht1(1,5+(3+A.NC/2)*i)=temp(2*i+2);
        i=i+1;
        if i>0
            wrksht1(1,6+(3+A.NC/2)*(i-1))=temp(2*i+1);
        end
    end
    wrksht1(1,1)=cellstr(Filename{iiii}(1:end-4)); 
    
    %Creates 2nd row format for cell array to be exported to excel
    wrksht1{2,1}='Crank Angle (°)';
    wrksht1{2,2}='Volume (m)';
    wrksht1{2,4}='Cycle Number';

    i=0;
    ii=0;
    while i<var_num %Loop inserts repeating 2nd row strings into the cell array
        wrksht1{2,5+(A.NC/2+3)*i}='Average';
        wrksht1{2,6+(A.NC/2+3)*i}='Std_dev';
        while ii<A.NC/2
            wrksht1{2,7+ii+(3+A.NC/2)*i}=ii+1;
            ii=ii+1;
        end
        ii=0;
        i=i+1;
    end
    
    wrksht1(3:length(A.ca)+2,1)=num2cell(A.ca(1,:))'; %Inserts crank angles into cell array
    wrksht1(3:length(A.volume)+2,2)=num2cell(A.volume(1,:))'; %Inserts Volumes into cell array
   
    i=0;
    % Insert motored data into cell array
    wrksht1(3:2882,5:6+A.NC/2+(A.NC/2+3)*i)=num2cell(motored_pcyl);
    i=i+1;
    wrksht1(3:2882,8+A.NC/2+(A.NC/2+3)*(i-1):6+A.NC/2+(A.NC/2+3)*i)=num2cell(motored_inj_curr);
    i=i+1;
    wrksht1(3:2882,8+A.NC/2+(A.NC/2+3)*(i-1):6+A.NC/2+(A.NC/2+3)*i)=num2cell(motored_prail);
    i=i+1;
    wrksht1(3:2882,8+A.NC/2+(A.NC/2+3)*(i-1):6+A.NC/2+(A.NC/2+3)*i)=num2cell(motored_pfilt);
    i=i+1;
    wrksht1(3:2882,8+A.NC/2+(A.NC/2+3)*(i-1):6+A.NC/2+(A.NC/2+3)*i)=num2cell(motored_dQdCA);
    i=i+1;
    wrksht1(3:2882,8+A.NC/2+(A.NC/2+3)*(i-1):6+A.NC/2+(A.NC/2+3)*i)=num2cell(motored_dQdCAfilt);
    xlswrite(xlsFilename, wrksht1, 1, 'A1'); %Writes motored data into excel worksheet 1
    
    %Inserts Fired data into cell array
    i=0;
    wrksht1(3:2882,5:6+A.NC/2+(A.NC/2+3)*i)=num2cell(fired_pcyl);
    i=i+1;
    wrksht1(3:2882,8+A.NC/2+(A.NC/2+3)*(i-1):6+A.NC/2+(A.NC/2+3)*i)=num2cell(fired_inj_curr);
    i=i+1;
    wrksht1(3:2882,8+A.NC/2+(A.NC/2+3)*(i-1):6+A.NC/2+(A.NC/2+3)*i)=num2cell(fired_prail);
    i=i+1;
    wrksht1(3:2882,8+A.NC/2+(A.NC/2+3)*(i-1):6+A.NC/2+(A.NC/2+3)*i)=num2cell(fired_pfilt);
    i=i+1;
    wrksht1(3:2882,8+A.NC/2+(A.NC/2+3)*(i-1):6+A.NC/2+(A.NC/2+3)*i)=num2cell(fired_dQdCA);
    i=i+1;
    wrksht1(3:2882,8+A.NC/2+(A.NC/2+3)*(i-1):6+A.NC/2+(A.NC/2+3)*i)=num2cell(fired_dQdCAfilt);
    xlswrite(xlsFilename, wrksht1, 2, 'A1'); %Writes fired data into excel worksheet 2


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Processed Worksheet %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    i=1;
    wrksht3=cell(5,A.NC/2+3);
    wrksht3(2,2:3)={'Average' 'Std_Dev'};
    while i<=A.NC/2
        wrksht3{2,i+3}=i;
        i=i+1;
    end

    wrksht3(1:6,1)={Filename{iiii}(1:end-4); 'Variable'; 'AID'; 'AID_filt';'T_TDC'; 'rho_TDC';};
    
    wrksht3(3:4,4:end)=[num2cell(A.AID'); num2cell(A.AIDfilt')];
    D=strfind(strcat(num2str(1000/A.std_AID),'.'),'.')-5; %Intermediate step for rounding avg and std
    E=strfind(strcat(num2str(1000/A.std_AIDfilt),'.'),'.')-5;
    wrksht3(3:6,2:3)={round(A.mean_AID*10^(D(1,1)))/(10^D(1,1)) round(A.std_AID*10^(D(1,1)+2))/(10^(D(1,1)+2)); round(A.mean_AIDfilt*10^(E(1,1)))/(10^E(1,1)) round(A.std_AIDfilt*10^(E(1,1)+2))/(10^(E(1,1)+2));A.Ttdc ' ';A.rhotdc ' '};
    xlswrite(xlsFilename, wrksht3, 3, 'A1');
    
    
    %   C=strcat(num2str(1000/stdd(i,1)),'.');%Divides by 1000 so works for any std <1000. Dot at end catches perfect divisions that result in no decimals
    %D=strfind(C,'.')-5; %Finds the decimal and subtracts to give position of first non zero after the decimal in original number
   %avg_round(i,1)=round(avg(i,1)*10^(D(1,1)))/(10^D(1,1)); %Uses position to round average number
   % std_round(i,1)=round(stdd(i,1)*10^(D(1,1)+2))/(10^(D(1,1)+2)); %Uses position to round std number  
   % i=i+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Input Sheet%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wrksht4=cell(17,3+A.NC/2);
    wrksht4(1,1:3)={'Prescribed' ' ' ' '};
    wrksht4(2:17,1)={'RPM'; 'SSE [CAD]'; 'DSE [ms]'; 'Air flow [g/s]'; 'Rc'; 'Fuel Type'; 'Filter Freq' ...
             ;'Injector Name'; 'Num orifice'; 'orifice diameter [mm]'; 'O2%'; ' '; 'Measured'; 'T_BDC [K]'; 'P_BDC [bar]'; ...
             'Hydraulic delay [CAD]'};
    wrksht4(14,2:3)={'Mean' 'Std'};
    i=1;
    for i=1:A.NC/2
        wrksht4{14,3+i}=i;
        i=i+1;
    end
    wrksht4(2:12,2)={A.Sset; round(A.sse_degV*10)/10; round(A.dse_ms*100)/100; A.mean_qair; A.Rc; A.fuel_type; [A.f1 A.f2] ;...
                     strcat(A.inj_name, A.inj_tip_name); A.hole_number; A.inj_diam; round(A.yo2)};
    
    D=strfind(strcat(num2str(1000/A.std_tintake),'.'),'.')-5; %Intermediate step for round avg and std             
    wrksht4(15:17,2:3)={round(A.mean_tintake*10^(D(1,1)))/(10^D(1,1)) round(A.std_tintake*10^(D(1,1)+2))/(10^(D(1,1)+2)); A.mean_pint_bdc ' '; '?' '?'};            
    xlswrite(xlsFilename, wrksht4, 4, 'A1');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Renames the Excel WorkSheets
    e=actxserver('Excel.Application'); %opens COM server
    ewb = e.Workbooks.Open(xlsFilename); %opens excel file
    ewb.Worksheets.Item(1).Name = 'Motored'; %Renames first worksheet
    ewb.Worksheets.Item(2).Name = 'Fired'; % Renames Second worksheet
    ewb.Worksheets.Item(3).Name = 'Processed'; % Renames Third worksheet
    ewb.Worksheets.Item(4).Name = 'Inputs'; % Renames fourth worksheet
    ewb.Save;
    ewb.Close(false);
    e.Quit;
        
    fullName=strcat(xlsFilename, '.xls');
    %winopen(fullName);    uncomment to open file when complete
end

