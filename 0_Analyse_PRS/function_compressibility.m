function z = function_compressibility(p,t,specie)

% 3/11/2014 - L.M. MALBEC
% This function computes the compressibilty coefficient of the input specie
% at pressure p (bar) and temperature t (K)

% Vm^3+(-bP-Rt)*Vm^2+aVm-ab=0

% Data found in Reid, R.C, Prausnitz, J. M., and Poling, B.E., The
% Properties of Gases and Liquids, Fourth Edition, McGraw-Hill, New York, 1987.
% referenced in
% http://courses.chem.indiana.edu/c360/documents/vanderwaalsconstants.pdf

if strcmp(specie,'n2')
    a = 1.370;  % bar L^2 / mol^2
    b = 0.0387; % L/mol
elseif strcmp(specie,'o2')
    a = 1.382;  % bar L^2 / mol^2
    b = 0.0319; % L/mol
elseif strcmp(specie,'co2')
    a = 3.658;  % bar L^2 / mol^2
    b = 0.0429; % L/mol
elseif strcmp(specie,'ar')
    a = 1.355;  % bar L^2 / mol^2
    b = 0.0320; % L/mol
elseif strcmp(specie,'air') % Composition (molar fraction): 0.2094 O2,0.7809 N2, 0.0097 Argon 
    a = 0.2094*1.382+0.7809*1.370+0.0097*1.355;
    b = 0.2094*0.0319+0.7809*0.0387+0.0097*0.0320;
end

a = a *1e5*1e-6; % conversion to Pa.m^6/mol2
b = b*1e-3; % conversion to m^3/mol
p = p*1e5; % conversion to Pa

myfun = @(x) (p+a./x.^2).*(x-b)-8.314*t;
vm1 = 1e-4;
vm2 = 1e4;
err1 = myfun(vm1);
err2 = myfun(vm2);

vm = 100e-3;
err = myfun(vm);
while abs(err)>1e-5
    if err>0
        vm2=vm;
        vm = (vm1+vm2)/2;
    elseif err<0
        vm1=vm;
        vm = (vm1+vm2)/2;
    end
    err = myfun(vm);
end
videal = 8.314*t/p;
z = vm/videal;

