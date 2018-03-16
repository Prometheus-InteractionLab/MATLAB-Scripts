% close all

% ar properties from AMESim
% pathname = 'c:\Users\malbecl\Documents\IFPEN\G_r106_Louis_Marie\03_Outils_Matlab\';
% [cp_ar_ame,cv_ar_ame,h_ar_ame,r_ar_ame,mu_ar_ame,lambda_ar_ame,sT_ar_ame] = function_read_janaf(pathname,'ar.data');
% s_ar_ame = @(T,P) sT_ar_ame(T) - r_ar*log(P);
% T_ame = [100:1000];

%% ar properties from EES
ar_in_ees = [
%   	T (degC)	cp (kJ/kg/K)	h(kJ/kg)	s(kJ/kg/K)
];


% T_ees = ar_in_ees(:,1)+273.15;
% cp_ar_ees = ar_in_ees(:,2);
% h_ar_ees = ar_in_ees(:,3);
% s_ar_ees = ar_in_ees(:,4);

%% AR from NIST

%Temperature (K)	298. - 6000.
A	=[ 20.78600];
B	=[ 2.825911e-7];
C	=[ -1.464191e-7];
D	=[ 1.092131e-8];
E	=[ -3.661371e-8];
F	=[ -6.197350];
G	=[ 179.9990];
H	=[ 0.000000];

% Cp nist: J/kg/K
cp_ar_nist = @(x)   (A(1) + B(1)*x + C(1)*x.^2 + D(1)*x.^3 + E(1)./x.^2);
cp_ar_nist = @(x) cp_ar_nist(x)/39.948*1000;

% h NIST: J/kg
h_ar_nist = @(x) (A(1)*x + 1/2*B(1)*x.^2 + 1/3*C(1)*x.^3 + 1/4*D(1)*x.^4 - E(1)./x + F(1) - H(1));
h_ar_nist = @(x) h_ar_nist(x)/39.948*1e6;

% s NIST: J/kg
s_ar_nist = @(x) (A(1)*log(x) + B(1)*x + 1/2*C(1)*x.^2 + 1/3*D(1)*x.^3 - E(1)./(2*x.^2) + G(1));
s_ar_nist = @(x) s_ar_nist(x)/39.948*1e3;


T_nist = [100:1:2000]/1000;

function_plot_properties([],[],T_nist*1000,cp_ar_nist(T_nist),'T [K]','Cp [kJ/kg/K]','Ar')
function_plot_properties([],[],T_nist*1000,h_ar_nist(T_nist),'T [K]','Enthalpy [kJ/kg]','Ar')
function_plot_properties([],[],T_nist*1000,s_ar_nist(T_nist),'T [K]','Entropy [kJ/kg/K]','Ar')
