%% Data from EES - 170711
data = [
    %   O2vol   Pbdc    Pmax    Tin     Ttdc    rho
    21      1.80	40.88	22.7	759     18.8
    21      1.80	41.01	23.4	761     18.8
    21      1.80	41.11	23.5	761     18.8
    21      1.81	41.23	23.8	761     18.86
    21      1.80	41.11	24.8	763     18.8
    21      1.81	41.17	24.8	763     18.8
    21      1.81	41.18	26.0	764     18.8
    21      1.81	41.29	26.0	765     18.8
    21      1.81	41.40	26.4	766     18.8
    21      1.81	41.45	26.8	767     18.8
    21      1.81	41.36	26.6	766     18.8
    21      1.82	41.56	26.9	767     18.9
    21      1.82	41.54	27.1	767     18.9
    21      1.81	41.32	27.1	767     18.8
    21      1.81	41.28	27.3	767     18.8
    21      1.81	41.28	27.4	767     18.8
    21      1.81	41.24	26.9	767     18.7
    21      1.81	41.20	26.7	766     18.8
    21      1.81	41.16	26.8	766     18.8
    
    ];

%% Loading thermodynamic properties from nist
load('o2prop_nist.mat');
load('n2prop_nist.mat');
load('airprop_nist.mat');

%%

Ttdc_nist=[];rho_tdc_nist=[];
for ii=1:size(data,1)
    isel = ii
    
    pco2 = data(isel,1)/100;
    P_bdc = data(isel,2);
    T_in = data(isel,4);
    P_tdc = data(isel,3);
    T_tdc = data(isel,5);
    rho_tdc = data(isel,6);
    
    T_bdc = function_Tintake_2_Tbdc(T_in)+273.15;
    
    % Fresh gas constant (dilution with N2)
    ndil = (0.21/pco2-1); % Nb of N2 moles for 1 mole or air
    Mair = 8.314/r_air*1000;
    Mgf = (ndil*28+Mair)/(1+ndil);
    r_fg = 8.314/Mgf*1000;
    
    xn2 = ndil*28/Mair; % Mass fraction of N2 as diluant
    xair = 1-xn2; % Mass fraction of air
    
    cp_fg = @(x) xn2*cp_n2_nist(x)+xair*cp_air_nist(x);
    s_fg = @(x) xn2*s_n2_nist(x)+xair*s_air_nist(x);
    
    s_ref = s_fg(T_bdc) - r_fg*log(P_bdc);
    
    T_tdc = T_bdc;
    s_tdc = s_fg(T_bdc) - r_fg*log(P_tdc);
    
    while abs(s_tdc-s_ref)>0.01
        deltaT = (s_ref-s_tdc)*T_tdc/cp_air_nist(T_tdc);
        T_tdc = T_tdc + deltaT;
        s_tdc = s_fg(T_tdc) - r_fg*log(P_tdc);
    end
    
    Ttdc_nist(ii) = T_tdc;
    rho_tdc_nist(ii)=P_tdc/r_fg/T_tdc*1e5;
    
end


figure
set(gcf,'color','w','position',[680 332 560 646])

subplot(2,1,1),grid on,hold on
plot(data(:,5),'color','k','linewidth',3,'marker','o','markerfacecolor','w','markersize',10);
plot(Ttdc_nist,'color',[1 1 1]*0.7,'linewidth',3,'marker','o','markerfacecolor','w','markersize',10);
set(gca,'fontname','Eras Light ITC','fontsize',16)
xlabel('# run');
ylabel('T_{tdc} [K]');
hl=legend('EES','NIST');
set(hl,'location','best')

subplot(2,1,2),grid on,hold on
plot(abs(Ttdc_nist'-data(:,5))*100./data(:,5),'color','k','linewidth',3,'marker','o','markerfacecolor','w','markersize',10);
set(gca,'fontname','Eras Light ITC','fontsize',16)
xlabel('# run');
ylabel('Err [%]');