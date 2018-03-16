function [avg_round, std_round]=function_round(avg,stdd)
i=1;
while i<=length(stdd) %Loops through each row
    if stdd(i,1)==0
        stdd(i,1)=0.01;
    end
    C=strcat(num2str(1000/stdd(i,1)),'.');%Divides by 1000 so works for any std <1000. Dot at end catches perfect divisions that result in no decimals
    D=strfind(C,'.')-5; %Finds the decimal and subtracts to give position of first non zero after the decimal in original number
    avg_round(i,1)=round(avg(i,1)*10^(D(1,1)))/(10^D(1,1)); %Uses position to round average number
    std_round(i,1)=round(stdd(i,1)*10^(D(1,1)+2))/(10^(D(1,1)+2)); %Uses position to round std number  
    i=i+1;
end
