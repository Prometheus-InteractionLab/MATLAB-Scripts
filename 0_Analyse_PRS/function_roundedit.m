function A=function_round(A)
i=1;
B=A;
assignin('base','B',B);
while i<=2880 %Loops through each row
    avg = mean(A.pcyl(2*i-1,i));
    stdd=std(A.pcyl(2*i-1,i));
    C=strcat(num2str(1000/stdd),'.');%Divides by 1000 so works for any std <1000.   
    D=strfind(C,'.')-5; %Finds the decimal and subtracts to give position of first non zero after the decimal in original number
    A.pcyl_meanfired(i,1)=round(avg*10^(D(1,1)))/(10^D(1,1)); %Uses position to round average number
    A.pcyl_stdfired(i,1)=round(stdd*10^(D(1,1)+2))/(10^(D(1,1)+2)); %Uses position to round std number
    
  ,  i=i+1;
end
