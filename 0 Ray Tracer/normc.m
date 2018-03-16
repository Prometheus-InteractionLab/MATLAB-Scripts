function Mout = normc(M)
%normalizes the columns of matrix M to a length of 1
Mout=zeros(size(M));
for c=1:size(M,2)
   Mout(:,c)=M(:,c)/sqrt(sum(M(:,c).^2));
end