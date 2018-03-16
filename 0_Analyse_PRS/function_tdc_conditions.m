function [rho,temp] = function_tdc_conditions(Tintake,Pbdc,Ptdc,pco2)

% 25/09/2014 - L.M. MALBEC

% Inputs:
% Tintake: Intake (plenum) temperaure [degC]
% Pbdc   : Pressure at BDC [bar]
% Ptdc   : Pressure at TDC [bar]
% pco2   : Mole fraction of O2 in the trapped gasses

% Outputs:
% rho: density at TDC [kg/m3]
% temp: temperature at TDC [K]

%% Loading thermodynamic properties from nist
load('o2prop_nist_v4.mat');
load('n2prop_nist_v4.mat');
load('airprop_nist_v4.mat');

% Estimation of the temperature at BDC knowing the intake temperature
Tbdc = function_Tintake_2_Tbdc(Tintake)+273.15;

% Fresh gas constant (dilution with N2)
ndil = (0.21/pco2-1); % Nb of N2 moles for 1 mole or air
Mair = 8.314/r_air*1000;
Mgf = (ndil*28+Mair)/(1+ndil);
r_fg = 8.314/Mgf*1000;


xn2 = ndil*28/Mair; % Mass fraction of N2 as diluant
xair = 1-xn2; % Mass fraction of air
    
cp_fg = @(x) xn2*cp_n2_nist(x)+xair*cp_air_nist(x); % cp as a function of temperature fot the fresh gasses
s_fg = @(x) xn2*s_n2_nist(x)+xair*s_air_nist(x);    % Entropy as a function of temperature fot the fresh gasses, at atmospheric pressure
    
s_ref = s_fg(Tbdc) - r_fg*log(Pbdc); % Value of the entropy of fresh gasses at BDC

% Loop on TDC temperature: when the entropy at TDC is equal to the entropy
% at BDC, the correct value of Ttdc has been found.
Ttdc = Tbdc;
s_tdc = s_fg(Tbdc) - r_fg*log(Ptdc);
while abs(s_tdc-s_ref)>0.01
    deltaT = (s_ref-s_tdc)*Ttdc/cp_air_nist(Ttdc);
    Ttdc = Ttdc + deltaT;
    s_tdc = s_fg(Ttdc) - r_fg*log(Ptdc);
end
    
temp = Ttdc;
rho=Ptdc/r_fg/Ttdc*1e5;