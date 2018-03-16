% =========================================================================
% Call the PRS Analysis script, and load variables into the workspace in
% order to compute AHRR (and other metrics) from a select number of engine
% cycles. (Since RCCI and other LTC modes may need a number of engien
% cycles beofre psuedo-steady-state operation is reached).
%
% Original - Feb,17,2017 - Greg Roberts
%
% =========================================================================



clear all; clc; format compact

% 
% % ---- Call PRS Analysis code, and bring variable rundata into the workspace
% [headers,dataxls,run_legend,rundata] = function_prs_analysis;
load('Z:\2017\20170301\170301w.mat')

%
% -----------------------> INPUTS <----------------
cycleStart              = 20 * 2;
CADstart                = -50;
CADend                  = 50;
saveFileName            = 'AHRR_170301w.mat';
% --------------------------------------------------



%% Use the last cycle as the end point
cycleEnd                = 50*2;

% ... Pull out the CAD and mean AHRR data for saving independently
AHRRdata.xPlot_CA       = rundata.ca(1,:);
AHRRdata.yPlot_meanHRR  = mean(rundata.dQdCAfilt([cycleStart:2:cycleEnd],:));
AHRRdata.runLetter      = rundata.runname;


% ... Compute the integrated HR, then CA50
[~,indSt]               = min(abs(AHRRdata.xPlot_CA - CADstart));
[~,indEnd]              = min(abs(AHRRdata.xPlot_CA - CADend));

HR_cumulative           = zeros(size(AHRRdata.yPlot_meanHRR(indSt:indEnd)));
HR_cumulative(1)        = AHRRdata.yPlot_meanHRR(indSt);
for i=2:length(HR_cumulative)
%     if AHRRdata.yPlot_meanHRR(indSt+i-1)>0
        HR_cumulative(i)    = HR_cumulative(i-1) + AHRRdata.yPlot_meanHRR(indSt+i-1);
%     else
%         HR_cumulative(i)    = HR_cumulative(i-1);
%     end
end   
HR_total                = HR_cumulative(end);
% Can compare with the trapezoidal integration function...
% HR_total_trapz          = trapz(AHRRdata.yPlot_meanHRR(indSt:indEnd));
HR_normalized           = HR_cumulative/HR_total;
[~,ind50]               = min(abs(HR_normalized - 0.5));
AHRRdata.CA50           = AHRRdata.xPlot_CA(indSt + ind50 - 1)
% ==> Save data, append to existing data file
AHRRdata_inFile         = load('AHRR_170301f.mat');
% AHRRdata_appended       = [AHRRdata_inFile.dPdata; AHRRdata] ;
% dPdata_inFile         = load('dP_collected_data.mat');
% dPdata_appended       = [dPdata_inFile.dPdata; dPdata] ;



% ... Find the start of LTC and HTC

% ... Compute gIMEP (using only selected cycles)

% % gIMEP (from "function_read_rundata.m")
% rng1    = find(rundata.ca(1,:)>=-180 & rundata.ca(1,:)<=180);
% Wcig    = trapz(V(rng1)', pcyl(firecyc,rng1)'.*1e5); % integral of pdV for all cycles [J]
% gIMEP   = Wcig./Vdisp; % [Pa]


%% ... PLot dPdCA and save
% dPdata.xPlot_CA       = rundata.ca(1,:);
% dPdata.yPlot_meandPdCA  = mean(rundata.dpdCAfilt([cycleStart:2:cycleEnd],:));
% dPdata.runLetter      = rundata.runname;

% ==> Save data, append to existing data file
% dPdata_inFile         = load('dPdC_data.mat');
% dPdata_appended       = [dPdata_inFile.dPdata; dPdata] ;
% % dPdata_inFile         = load('dP_collected_data.mat');
% dPdata_appended       = [dPdata_inFile.dPdata; dPdata] ;


% ... Make plot
figure(44); 
plot(AHRRdata.xPlot_CA, AHRRdata.yPlot_meanHRR,'k-.'); hold on;
% plot(AHRRdata.CA50, AHRRdata.yPlot_meanHRR(indSt + ind50 - 1),'mo');
grid on
xlabel('CAD  [deg]'); ylabel('Mean AHRR  [J/deg]')
axis([-60 80 -10 1.2*max(AHRRdata.yPlot_meanHRR)])


% % -------------
% save(saveFileName,'AHRRdata')






