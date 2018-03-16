% close all

% N2 properties from AMESim
pathname = 'c:\Users\malbecl\Documents\IFPEN\G_r106_Louis_Marie\03_Outils_Matlab\';
[cp_n2_ame,cv_n2_ame,h_n2_ame,r_n2_ame,mu_n2_ame,lambda_n2_ame,sT_n2_ame] = function_read_janaf(pathname,'N2.data');
s_n2_ame = @(T,P) sT_n2_ame(T) - r_n2*log(P);
T_ame = [100:1000];

%% N2 properties from EES
n2_in_ees = [
%   	T (degC)	cp (kJ/kg/K)	h(kJ/kg)	s(kJ/kg/K)
        50          1038            25.95       6.923
        100         1041            77.91       7.073
        150         1047            130.1       7.204
        200         1054            182.6       7.321
        250         1063            235.5       7.428
        300         1072            288.9       7.525
        350         1083            342.8       7.615
        400         1094            397.2       7.699
        450         1104            452.1       7.778
        500         1116            507.6       7.852
        550         1126            563.7       7.922
        600         1137            620.3       7.989
        650         1148            677.4       8.053
        700         1158            735         8.114
        750         1168            793.2       8.172
];


T_ees = n2_in_ees(:,1)+273.15;
cp_n2_ees = n2_in_ees(:,2);
h_n2_ees = n2_in_ees(:,3);
s_n2_ees = n2_in_ees(:,4);

%% N2 from NIST

%Temperature (K)	100. - 500.	500. - 2000.	2000. - 6000.
A	=[              28.98641	19.50583        35.51872];
B=[                 1.853978	19.88705        1.128728];
C=[                 -9.647459	-8.598535       -0.196103];
D=[                 16.63537	1.369784        0.014662];
E=[                 0.000117	0.527601        -4.553760];
F=[                 -8.671914	-4.935202       -18.97091];
G=[                 226.4168	212.3900        224.9810];
H=[                 0.0         0.0             0.0];

% Cp nist: J/kg/K
cp_n2_nist = @(x)   (A(1) + B(1)*(x/1000) + C(1)*(x/1000).^2 + D(1)*(x/1000).^3 + E(1)./(x/1000).^2).*(x<=500)+...
                    (A(2) + B(2)*(x/1000) + C(2)*(x/1000).^2 + D(2)*(x/1000).^3 + E(2)./(x/1000).^2).*((x>500).*(x<=2000))+...
                    (A(3) + B(3)*(x/1000) + C(3)*(x/1000).^2 + D(3)*(x/1000).^3 + E(3)./(x/1000).^2).*(x>2000);
cp_n2_nist = @(x) cp_n2_nist(x)/28*1000;

% h NIST: J/kg
h_n2_nist = @(x) (A(1)*(x/1000) + 1/2*B(1)*(x/1000).^2 + 1/3*C(1)*(x/1000).^3 + 1/4*D(1)*(x/1000).^4 - E(1)./(x/1000) + F(1) - H(1)).*(x<=500)+...
                 (A(2)*(x/1000) + 1/2*B(2)*(x/1000).^2 + 1/3*C(2)*(x/1000).^3 + 1/4*D(2)*(x/1000).^4 - E(2)./(x/1000) + F(2) - H(2)).*((x>500).*(x<=2000))+...
                 (A(3)*(x/1000) + 1/2*B(3)*(x/1000).^2 + 1/3*C(3)*(x/1000).^3 + 1/4*D(3)*(x/1000).^4 - E(3)./(x/1000) + F(3) - H(3)).*(x>2000);
h_n2_nist = @(x) h_n2_nist(x)/28*1e6;

% s NIST: J/kg
s_n2_nist = @(x) (A(1)*log((x/1000)) + B(1)*(x/1000) + 1/2*C(1)*(x/1000).^2 + 1/3*D(1)*(x/1000).^3 - E(1)./(2*(x/1000).^2) + G(1)).*(x<=500)+...
                 (A(2)*log((x/1000)) + B(2)*(x/1000) + 1/2*C(2)*(x/1000).^2 + 1/3*D(2)*(x/1000).^3 - E(2)./(2*(x/1000).^2) + G(2)).*((x>500).*(x<=2000))+...
                 (A(3)*log((x/1000)) + B(3)*(x/1000) + 1/2*C(3)*(x/1000).^2 + 1/3*D(3)*(x/1000).^3 - E(3)./(2*(x/1000).^2) + G(3)).*(x>2000);
s_n2_nist = @(x) s_n2_nist(x)/28*1e3;


T_nist = [100:1:2000];

function_plot_properties(T_ees,cp_n2_ees,T_nist,cp_n2_nist(T_nist),'T [K]','Cp [kJ/kg/K]','N2')
function_plot_properties(T_ees,h_n2_ees*1000,T_nist,h_n2_nist(T_nist),'T [K]','Enthalpy [kJ/kg]','N2')
function_plot_properties(T_ees,s_n2_ees*1000,T_nist,s_n2_nist(T_nist),'T [K]','Entropy [kJ/kg/K]','N2')