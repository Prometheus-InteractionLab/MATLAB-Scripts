
N2_properties_nist;
O2_properties_nist;
AR_properties_nist;

cpmass_air_nist = @(x) 0.2314*cpmass_o2_nist(x) +0.7552*cpmass_n2_nist(x)+0.0134*cpmass_ar_nist(x);
hmass_air_nist = @(x) 0.2314*hmass_o2_nist(x) +0.7552*hmass_n2_nist(x)+0.0134*hmass_ar_nist(x);
smass_air_nist = @(x) 0.2314*smass_o2_nist(x) +0.7552*smass_n2_nist(x)+0.0134*smass_ar_nist(x);

r_air = 8.314/(0.2094*32+0.7809*28+0.0097*39.948)*1000;