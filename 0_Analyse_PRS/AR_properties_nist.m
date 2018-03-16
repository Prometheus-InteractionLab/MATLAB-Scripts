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
cpmol_ar_nist = @(x)   (A(1) + B(1)*x + C(1)*x.^2 + D(1)*x.^3 + E(1)./x.^2);
cpmass_ar_nist = @(x) cpmol_ar_nist(x)/39.948*1000;

% h NIST: J/kg
hmol_ar_nist = @(x) (A(1)*x + 1/2*B(1)*x.^2 + 1/3*C(1)*x.^3 + 1/4*D(1)*x.^4 - E(1)./x + F(1) - H(1));
hmass_ar_nist = @(x) hmol_ar_nist(x)/39.948*1e6;

% s NIST: J/kg
smol_ar_nist = @(x) (A(1)*log(x) + B(1)*x + 1/2*C(1)*x.^2 + 1/3*D(1)*x.^3 - E(1)./(2*x.^2) + G(1));
smass_ar_nist = @(x) smol_ar_nist(x)/39.948*1e3;