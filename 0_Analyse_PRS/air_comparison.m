clear all
close all
clc

% air properties from AMESim
pathname = 'c:\Users\malbecl\Documents\IFPEN\G_r106_Louis_Marie\03_Outils_Matlab\';
[cp_air_ame,cv_air_ame,h_air_ame,r_air_ame,mu_air_ame,lambda_air_ame,sT_air_ame] = function_read_janaf(pathname,'air.data');
s_air_ame = @(T,P) sT_air_ame(T) - r_air*log(P);
T_ame = [100:1000];


%% N2 data
n2_comparison;
o2_comparison;
ar_comparison;

%% air properties from EES
air_in_ees = [
%   	T (degC)	cp (kJ/kg/K)	h(kJ/kg)	s(kJ/kg/K)

        25          1005            298.5       5.699
        50          1006            323.7       5.78
        75          1008            348.9       5.855
        100         1010            374.1       5.925
        125         1013            399.4       5.991
        150         1016            424.8       6.052
        175         1020            450.2       6.111
        200         1024            475.8       6.166
        225         1029            501.4       6.219
        250         1034            527.2       6.269
        275         1039            553.1       6.318
        300         1045            579.2       6.364
        325         1050            605.4       6.409
        350         1056            631.7       6.452
        375         1062            658.2       6.494
        400         1068            684.8       6.534
        425         1074            711.6       6.573
        450         1080            738.5       6.611
        475         1086            765.6       6.648
        500         1092            792.8       6.684
        525         1098            820.2       6.718
        550         1104            847.7       6.752
        575         1109            875.4       6.786
        600         1115            903.2       6.818
        625         1120            931.1       6.849
        650         1126            959.2       6.88
        675         1131            987.4       6.91
        700         1136            1016        6.94
        725         1140            1044        6.969
];


T_ees = air_in_ees(:,1)+273.15;
cp_air_ees = air_in_ees(:,2);
h_air_ees = air_in_ees(:,3);
s_air_ees = air_in_ees(:,4);

%% air from NIST

% cp_air_nist = @(x) 0.233*cp_o2_nist(x) +(1-0.233)*cp_n2_nist(x);
% h_air_nist = @(x) 0.233*h_o2_nist(x) +(1-0.233)*h_n2_nist(x);
% s_air_nist = @(x) 0.233*s_o2_nist(x) +(1-0.233)*s_n2_nist(x);


cp_air_nist = @(x) 0.2314*cp_o2_nist(x) +0.7552*cp_n2_nist(x)+0.0134*cp_ar_nist(x);
h_air_nist = @(x) 0.2314*h_o2_nist(x) +0.7552*h_n2_nist(x)+0.0134*h_ar_nist(x);
s_air_nist = @(x) 0.2314*s_o2_nist(x) +0.7552*s_n2_nist(x)+0.0134*s_ar_nist(x);

r_air = 8.314/(0.2094*32+0.7809*28+0.0097*39.948)*1000;

T_nist = [100:1:2000];


function_plot_properties(T_ees,cp_air_ees,T_nist,cp_air_nist(T_nist),'T [K]','Cp [J/kg/K]','Air')
function_plot_properties(T_ees,h_air_ees*1000-293e3,T_nist,h_air_nist(T_nist),'T [K]','Enthalpy [J/kg]','Air')
function_plot_properties(T_ees,s_air_ees*1000+1045,T_nist,s_air_nist(T_nist),'T [K]','Entropy [J/kg/K]','Air')


