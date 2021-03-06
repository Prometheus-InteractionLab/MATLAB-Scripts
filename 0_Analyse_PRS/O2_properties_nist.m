%% O2 from NIST

%Temperature (K)	100. - 700.	700. - 2000.	2000. - 6000.
A=[                 31.32234	30.03235	20.91111];
B=[                 -20.23531	8.772972	10.72071];
C=[                 57.86644	-3.988133	-2.020498];
D=[                 -36.50624	0.788313	0.146449];
E=[                 -0.007374	-0.741599	9.245722];
F=[                 -8.903471	-11.32468	5.337651];
G=[                 246.7945	236.1663	237.6185];
H=[                 0.0         0.0         0.0];

% Cp nist: J/kg/K
cpmol_o2_nist = @(x)   (A(1) + B(1)*(x/1000) + C(1)*(x/1000).^2 + D(1)*(x/1000).^3 + E(1)./(x/1000).^2).*(x<=700)+...
                    (A(2) + B(2)*(x/1000) + C(2)*(x/1000).^2 + D(2)*(x/1000).^3 + E(2)./(x/1000).^2).*((x>700).*(x<=2000))+...
                    (A(3) + B(3)*(x/1000) + C(3)*(x/1000).^2 + D(3)*(x/1000).^3 + E(3)./(x/1000).^2).*(x>2000);
cpmass_o2_nist = @(x) cpmol_o2_nist(x)/32*1000;

% h NIST: J/kg
hmol_o2_nist = @(x) (A(1)*(x/1000) + 1/2*B(1)*(x/1000).^2 + 1/3*C(1)*(x/1000).^3 + 1/4*D(1)*(x/1000).^4 - E(1)./(x/1000) + F(1) - H(1)).*(x<=700)+...
                 (A(2)*(x/1000) + 1/2*B(2)*(x/1000).^2 + 1/3*C(2)*(x/1000).^3 + 1/4*D(2)*(x/1000).^4 - E(2)./(x/1000) + F(2) - H(2)).*((x>700).*(x<=2000))+...
                 (A(3)*(x/1000) + 1/2*B(3)*(x/1000).^2 + 1/3*C(3)*(x/1000).^3 + 1/4*D(3)*(x/1000).^4 - E(3)./(x/1000) + F(3) - H(3)).*(x>2000);
hmass_o2_nist = @(x) hmol_o2_nist(x)/32*1e6;

% s NIST: J/kg
smol_o2_nist = @(x) (A(1)*log((x/1000)) + B(1)*(x/1000) + 1/2*C(1)*(x/1000).^2 + 1/3*D(1)*(x/1000).^3 - E(1)./(2*(x/1000).^2) + G(1)).*(x<=700)+...
                 (A(2)*log((x/1000)) + B(2)*(x/1000) + 1/2*C(2)*(x/1000).^2 + 1/3*D(2)*(x/1000).^3 - E(2)./(2*(x/1000).^2) + G(2)).*((x>700).*(x<=2000))+...
                 (A(3)*log((x/1000)) + B(3)*(x/1000) + 1/2*C(3)*(x/1000).^2 + 1/3*D(3)*(x/1000).^3 - E(3)./(2*(x/1000).^2) + G(3)).*(x>2000);
smass_o2_nist = @(x) smol_o2_nist(x)/32*1e3;