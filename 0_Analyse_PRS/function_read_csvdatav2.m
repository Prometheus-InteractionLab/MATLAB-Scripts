function [csvdata,i2] = function_read_csvdata(data,dayname,run,i1,i2)

% added i2 as output incase it's value changes E.Eagle 141202
% A = importdata('140605.csv');

index = find(strcmp(data.RunCode(:,2),run)==1);
try
rng = find(data.injpulses(index,34)>0);
rng2 = find(data.data(index,34) >= max(data.data(index,34)));
i1 = rng(1);
i2 = rng2(1);
catch
end

if length(index) < i2 - 10
    i2 = length(index)-15;
end
csvdata.name_run        = run;
csvdata.name_dayname    = dayname;
csvdata.name_rpm        = data.textdata(1,3);
csvdata.name_qair       = data.textdata(1,4);
csvdata.name_pintake    = data.textdata(1,6);
csvdata.name_tintake    = data.textdata(1,7);
csvdata.name_tliner     = data.textdata(1,8);
csvdata.name_pfuel      = data.textdata(1,15);

% 
summode =  data.data(index(1),19);
if summode==1
    csvdata.qair = data.data(index,21)+data.data(index,25)+data.data(index,27);
elseif summode==2
    csvdata.qair = data.data(index,21)+data.data(index,23);
end
csvdata.rpm     = data.data(index,3);
csvdata.pintake = data.data(index,6);
csvdata.tintake = data.data(index,7);
csvdata.tliner  = data.data(index,8);
% csvdata.pfuel   = data.data(index-1,13);

% Mean values
csvdata.mean_rpm     = mean(csvdata.rpm(i1:i2))
csvdata.mean_qair    = mean(csvdata.qair(i1:i2))
% csvdata.mean_pintake = mean(csvdata.pintake(i1:i2));
csvdata.mean_tintake = mean(csvdata.tintake(i1:i2))
csvdata.mean_tliner  = mean(csvdata.tliner(i1:i2))
% csvdata.mean_pfuel   = mean(csvdata.pfuel(i1:i2));