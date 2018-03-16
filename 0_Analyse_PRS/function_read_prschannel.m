function [prschan,prschan_m] = function_read_prschannel(data,numchn,selcyc,SPC,NC,a,b)

prschan = reshape(data(numchn,:),SPC,NC)';
prschan = a*prschan+b;prschan=prschan(selcyc,:);
prschan_m = mean(prschan,1);