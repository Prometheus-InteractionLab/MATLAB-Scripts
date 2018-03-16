function Tbdc = function_Tintake_2_Tbdc(Tintake)

% Tintake in dedC
% Tbdc in degC

Tcool = 95;
Tbdc = 0.68*(Tintake-Tcool)+Tcool;