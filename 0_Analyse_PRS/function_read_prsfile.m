function [data,N,SPC] = function_read_prsfile(pathname,filename)

lab = 165;
% addpath('C:\Users\lcmalbe\Documents\Matlab_Scripts\Pressure');

dayname = filename(1:end-5);
runname = filename(end-4);

[data,N,SPC] = readprs(lab, filename,0);