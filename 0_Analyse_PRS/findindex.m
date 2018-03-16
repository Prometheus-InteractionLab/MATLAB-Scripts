function indices = findindex(values,array)
%[indices]=findindex(values,array)
%finds the indices for values in array that are closest to the corresponding values in
%values.  Think of it as "find 'values' in 'array'."
for i=1:length(values)
    [minval,indices(i)]=min(abs(array-values(i)));
end