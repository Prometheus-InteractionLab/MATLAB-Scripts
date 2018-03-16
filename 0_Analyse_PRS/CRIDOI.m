function [ DOI ] = function_DOI_from_DSE(DSE,inj_name)
%UNTITLED3 Summary of this function goes here
%   takes DSE in ms and outputs DOI

if inj_name = 'CRI'
DOI = zeros(1,length(DSE));
index = find(DSE < 1.700);
DOI(index) = DSE(index) + (DSE(index)-.500)*.8;
index = find(DSE >= 1.7);
DOI(index) = DSE(index) + .85;
index = find(DOI < 0);
DOI(index) = 0;
else
    DOI = DSE;
    sprintf('No injector profile found')
end



end

