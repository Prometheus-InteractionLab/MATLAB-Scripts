function [data,ivalid] = function_load_processed_data()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[filename,pathname] = uigetfile('*.mat','Select files','Multiselect','on');
filename = cellstr(filename);
addpath(pathname) %added by E. Eagle 141009, prevents code crashing when selected file isn't already in the matlab path

for ii=1:length(filename)
    load(filename{ii})
    try data.Ttdc(ii) = rundata.Ttdc; catch data.T(ii) = NaN; end
    try data.fuel_type = rundata.fuel_type; catch end
    try data.dse_ms(ii) = rundata.dse_ms; catch data.dse(ii) = NaN; end
    try data.mean_SOC(ii) = rundata.mean_SOC; data.std_SOC(ii) = rundata.std_SOC; catch data.mean_SOC(ii) = NaN;  data.std_SOC(ii) = NaN; end
    try data.mean_AID(ii) = rundata.mean_AID; catch data.mean_AID(ii) = NaN; end
    try data.std_AID(ii) = rundata.std_AID; catch data.std_AID(ii) = NaN; end
    try data.mean_AIDfilt(ii) = rundata.mean_AIDfilt; catch data.mean_AIDfilt(ii) = NaN; end
    try data.std_AIDfilt(ii) = rundata.std_AIDfilt; catch data.std_AIDfilt(ii) = NaN; end
    try data.mean_inj_delay(ii) = rundata.mean_inj_delay; data.std_inj_delay(ii) = rundata.std_inj_delay; catch data.mean_inj_delay(ii) = NaN; data.std_inj_delay(ii) = NaN; end
    try data.Tunc(ii) = rundata.std_pmax/rundata.mean_pmax*rundata.Ttdc; catch data.Tunc(ii) = NaN; end %uncertainty propagation in P_TDC to T_TDC
    try data.gIMEP(ii,:) = rundata.gIMEP; catch data.gIMEP(ii) = NaN; end
    try data.mean_gIMEP(ii) = rundata.mean_gIMEP; catch data.gIMEP(ii) = NaN; end
    try data.std_gIMEP(ii) = rundata.std_gIMEP; catch data.std_gIMEP(ii) = NaN; end
    try data.DOC(ii,:) = rundata.DOC; catch data.DOC(ii) = NaN; end
    try data.mean_DOC(ii) = rundata.mean_DOC; catch data.mean_DOC(ii) = NaN; end
    try data.std_DOC(ii) = rundata.std_DOC; catch data.std_DOC(ii) = NaN; end
    try data.DOC_ms(ii) = rundata.DOC./(rundata.Sset*6)*1000; catch data.DOC(ii) = NaN; end
    try data.mean_DOC_ms(ii) = rundata.mean_DOC./(rundata.Sset*6)*1000; catch data.mean_DOC(ii) = NaN; end
    try data.std_DOC_ms(ii) = rundata.std_DOC./(rundata.Sset*6)*1000; catch data.std_DOC(ii) = NaN; end
    try data.AHRRmax(ii,:) = rundata.AHRRmax; catch data.AHRRmax(ii) = NaN; end
    try data.mean_AHRRmax(ii) = rundata.mean_AHRRmax; catch data.mean_AHRRmax(ii) = NaN; end
    try data.std_AHRRmax(ii) = rundata.std_AHRRmax; catch data.std_AHRRmax(ii) = NaN; end
    try data.AHRRmaxfilt(ii,:) = rundata.AHRRmaxfilt; catch data.AHRRmaxfilt(ii) = NaN; end
    try data.mean_AHRRmaxfilt(ii) = rundata.mean_AHRRmaxfilt; catch data.mean_AHRRmaxfilt(ii) = NaN; end
    try data.std_AHRRmaxfilt(ii) = rundata.std_AHRRmaxfilt; catch data.std_AHRRmaxfilt(ii) = NaN; end
end
ivalid = find(isfinite(data.Ttdc)>0);
end
