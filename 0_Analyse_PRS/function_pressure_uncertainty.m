function [ std_y ] = function_pressure_uncertainty( y,h,periodic )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

N = length(y);
std_y = zeros(1,N);

ym2h = shftvec(y,-2);
ymh  = shftvec(y,-1);
yph  = shftvec(y,1);
yp2h = shftvec(y,2);

num = (1/12/h)^2.*(ym2h.^2 + (8.*ymh).^2 + (8.*yph).^2 + yp2h.^2);
if (~periodic)
   k = 2;
else 
   k = 0;
end;
std_y(k+1:N-k) = num(k+1:N-k);


end

