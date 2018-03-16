function prsdata = function_read_prsdata(data,N,SPC)

debug = 0;
lab = 165;
psetup;
labinfo;

ncmp = 1.345; % approx. motoring polytropic exponent for compression stroke
nexp = 1.450; % approx. motoring polytropic exponent for expansion stroke
nbrn = 1.280; % approx. polytropic exponent for combustion products

NC = N(2)./SPC; % number of cycles in data file, SPC = samples per cycle

motorcyc = 1:2:NC;
firecyc = 2:2:NC;

ca = reshape(data(1,:),SPC,NC)' + EncPErr;
pcyl = reshape(data(2,:),SPC,NC)';
pint = reshape(data(3,:),SPC,NC)';
injComm = reshape(data(4,:),SPC,NC)';
prail = reshape(data(5,:),SPC,NC)';
injCurr = reshape(data(6,:),SPC,NC)';
if size(N,1)
camTrig = reshape(data(7,:),SPC,NC)';
end

% Volts to physical units
pcyl = function_volts_2_phys(pcyl,pcylS,Tpcylxd,pcylO);
pint = function_volts_2_phys(pint,pintS,Tpintxd,pintO);
prail = function_volts_2_phys(prail,prailS,Tprailxd,prailO);

% Offsets
pegrng = find(ca(1,:)>=-190 & ca(1,:)<=-170);
poffset = mean(pcyl(:,pegrng),2);         % column vector (all cycles)
pegpcol = mean(pint(:,pegrng),2);
pegp    = mean(pegpcol);                   % for later calculations
pshift = (pegpcol-poffset)*ones(1,SPC);   % matrix
pcyl   = pcyl + pshift;                 % pcyl array now pegged over all cycles

prsdata.pcyl        = pcyl;
prsdata.ca          = ca;
prsdata.pint        = pint;
prsdata.injComm     = injComm;
prsdata.prail       = prail;
prsdata.injCurr     = injCurr;
prsdata.camTrig     = camTrig;

% Compute compression and expansion polytropic coefficients

% AHRR
[pideal,V,ypiston,CRg,CRgmax] = cc_pVy(lab,ca(1,:),ncmp,nexp,pegp,mean(mean(injCurr)));
V = V./1.e9; % [m.^3]
Vdisp = (max(V) - min(V)); % [m.^3]
Rc = max(V)/min(V);
prsdata.volume      = ones(size(ca,1),1)*V;

ca1 = -90;ca2=-30;
[a,ica1] = min(abs(prsdata.ca(1,:)-ca1));[a,ica2] = min(abs(prsdata.ca(1,:)-ca2));
logvolc = log(prsdata.volume(2:2:end,ica1:ica2));logvolc=logvolc(:);
logpcylc = log(prsdata.pcyl(2:2:end,ica1:ica2));logpcylc=logpcylc(:);
pcomp = polyfit(logvolc,logpcylc,1);ncompmeas = pcomp(1);
ycomp = polyval(pcomp,[min(logvolc) max(logvolc)]);
ycomp2 = polyval(pcomp,[min(log(prsdata.volume(1,:))) max(log(prsdata.volume(1,:)))]);

ca1 = 40;ca2=100;
[a,ica1] = min(abs(prsdata.ca(1,:)-ca1));[a,ica2] = min(abs(prsdata.ca(1,:)-ca2));
logvole = log(prsdata.volume(2:2:end,ica1:ica2));logvole=logvole(:);
logpcyle = log(prsdata.pcyl(2:2:end,ica1:ica2));logpcyle=logpcyle(:);
pexp = polyfit(logvole,logpcyle,1);nexpmeas = pexp(1);
yexp = polyval(pexp,[min(logvole) max(logvole)]);
yexp2 = polyval(pexp,[min(log(prsdata.volume(1,:))) max(log(prsdata.volume(1,:)))]);

[dpdCA,dQdCA,pfiltA,dpdCAfilt,dQdCAfilt] = compQ(lab,S,V,ca(1,:),pcyl,-ncompmeas,-nexpmeas);
[dpdCAavg, dQdCAavg, pfiltAavg,dpdCAavgfilt,dQdCAavgfilt] = compQ(lab,S,V,ca(1,:),[mean(pcyl(1:2:NC,:)); mean(pcyl(2:2:NC,:))],-ncompmeas,-nexpmeas);

CA = ca(1,:);
navg = (-ncompmeas-nexpmeas)/2;
Ac = pi.*((dc./1.e3).^2)./4; % [m^2] dc (bore diameter in mm) from labinfo.m                            
CArad = CA.*pi./180.; % convert CA vector to radians
K = 1 - ( (R./L).*sin(CArad) ).^2;
dVdCA = Ac.*(R./1.e3).*sin(CArad).*(1 + R.*cos(CArad)./L./(K.^0.5)); % [m^3/rad]

std_pf = std(pcyl(2:2:NC,:).*1e5); %from [bar] to [Pa]
std_pm = std(pcyl(1:2:NC,:).*1e5); %from [bar] to [Pa]

std_dpdCAf = std(dpdCA(2:2:NC,:)); %[Pa/CAD]
std_dpdCAm = std(dpdCA(1:2:NC,:)); %[Pa/CAD]

std_dpdCAfiltf = std(dpdCAfilt(2:2:NC,:));
std_dpdCAfiltm = std(dpdCAfilt(1:2:NC,:));

[ustd_dpdradf] = function_pressure_uncertainty(std_pf,mean(diff(ca(1,:))),1); %[Pa/CAD]
[ustd_dpdradm] = function_pressure_uncertainty(std_pm,mean(diff(ca(1,:))),1);

udQdCAf1 = ((navg/(navg-1)*dVdCA*pi/180).^2.*std_pf.^2); 
udQdCAf2 = (V/(navg-1)).^2.*ustd_dpdradf;
udQdCAf3 = (V/(navg-1)).^2.*(std_dpdCAf.^2);
udQdCAf4 = (V/(navg-1)).^2.*(std_dpdCAfiltf.^2);

udQdCAm1 = ((navg/(navg-1)*dVdCA*pi/180).^2.*std_pm.^2);
udQdCAm2 = (V/(navg-1)).^2.*ustd_dpdradm;
udQdCAm3 = (V/(navg-1)).^2.*(std_dpdCAm.^2);

if debug
    index = find(CA<-290 & CA > -335)   
    pint2 = mean(mean(pcyl(2:2:NC,index)))
    one_sigma = mean(std(pcyl(2:2:NC,index)))
    
    figure(2); hold on; grid on;
    plot(CA,dpdCA(2:2:NC,:),'color',[.9 .5 .5])
    plot(CA,dpdCA(1:2:NC,:),'color',[.5 .5 .5])
    plot(CA,dpdCAfilt(2:2:NC,:),'m')
    plot(CA,mean(dpdCAfilt(2:2:NC,:)),'k--','linewidth',3)
    plot(CA,mean(dpdCA(2:2:NC,:)),'k','linewidth',2)
    plot(CA,mean(dpdCA(1:2:NC,:)),'k:','linewidth',2)
    plot(CA,dpdCAavg(2,:),'g')
    %legend('dpdCA fired','dpdCA motored','dpdCA filter','mean dpdCAfilt','mean dpdCA f','mean dpdCA m','dpavgdCA')
    
    figure(12); hold on; grid on;
    plot(CA,V./(navg-1).*mean(dpdCAfilt(2:2:NC,:)),'k--','linewidth',3)
    plot(CA,V./(navg-1).*mean(dpdCA(2:2:NC,:)),'r','linewidth',2)
    plot(CA,V./(navg-1).*mean(dpdCA(1:2:NC,:)),'k','linewidth',2)
    plot(CA,V./(navg-1).*dpdCAavg(1,:),'g')
    plot(CA,V./(navg-1).*dpdCAavg(2,:),'g')
    
    figure(3); hold on; grid on;
    plot(CA,std_dpdCAfiltf,'r')
    plot(CA,std_dpdCAfiltm,'r:')
    plot(CA,std_dpdCAf,'k')
    plot(CA,std_dpdCAm,'k:')
    legend('fired_filter','motored_filter','fired','motored')
    
    figure(4); hold on; grid on;
    plot(CA,udQdCAf1.^.5,'b')
    plot(CA,udQdCAf2.^.5,'r')
    plot(CA,udQdCAf3.^.5,'m')
    plot(CA,udQdCAf4.^.5,'g')
    legend('std p','propagated','std dpdCA','std dpdCAfilt')
    
    figure(5); hold on; grid on;
    plot(CA,dpdCAavg(2,:),'k','linewidth',2)
    plot(CA,dpdCAavg(1,:),'r:','linewidth',2)
    plot(CA,dpdCAavgfilt(2,:),'b')
    plot(CA,dpdCAavgfilt(1,:),'b')
    
    figure(6); hold on; grid on;
    plot(CA,dQdCA(2:2:NC,:),'color',[.5 .5 .5])
    plot(ca(1,:),dQdCAavgfilt(2,:),'k','linewidth',2)
    plot(CA,dQdCAavg(2,:),'g','linewidth',2)
    plot(CA,dQdCAavgfilt(2,:)+2.*udQdCAf1.^.5,'b')
    plot(CA,dQdCAavgfilt(2,:)-2.*udQdCAf1.^.5,'b')
    plot(CA,dQdCAavgfilt(2,:)+2.*udQdCAf2.^.5,'r')
    plot(CA,dQdCAavgfilt(2,:)-2.*udQdCAf2.^.5,'r')
    plot(CA,dQdCAavgfilt(2,:)+2.*udQdCAf3.^.5,'m')
    plot(CA,dQdCAavgfilt(2,:)-2.*udQdCAf3.^.5,'m')
    
    figure(7); hold on; grid on;
    plot(CA,dQdCA(1:2:NC,:),'color',[.5 .5 .5])
    plot(ca(1,:),dQdCAavgfilt(1,:),'k:','linewidth',2)
    plot(CA,dQdCAavgfilt(1,:)+3.*udQdCAm1.^.5,'b:')
    plot(CA,dQdCAavgfilt(1,:)-3.*udQdCAm1.^.5,'b:')
    plot(CA,dQdCAavgfilt(1,:)+udQdCAm2.^.5,'r:')
    plot(CA,dQdCAavgfilt(1,:)-udQdCAm2.^.5,'r:')
    plot(CA,dQdCAavgfilt(1,:)+2.*udQdCAm3.^.5,'m:')
    plot(CA,dQdCAavgfilt(1,:)-2.*udQdCAm3.^.5,'m:')
    
end
debug = 0;

alpha_ahrrmax=[];ahrrmax=[];
alpha_ahrrmaxfilt=[];ahrrmaxfilt=[];
alpha_pmax =[]; pmax = [];
istart = find(ca(1,:)>-90,1,'first');
iend = find(ca(1,:)>90,1,'first');

for ii=firecyc
    [ahrrmax(end+1),imax] = max(dQdCA(ii,istart:iend));
    alpha_ahrrmax(end+1) = ca(ii,imax+istart-1);
    [ahrrmaxfilt(end+1),imax] = max(dQdCAfilt(ii,istart:iend));
    alpha_ahrrmaxfilt(end+1) = ca(ii,imax+istart-1);
    [pmax(end+1),imax] = max(pcyl(ii,istart:iend));
    alpha_pmax(end+1) = ca(ii,imax+istart-1);
end

% Pcyl analysis
alpha_pmotor=[];pmotor=[];
for ii=motorcyc
    [pmotor(end+1),imax] = max(pcyl(ii,:));
    alpha_pmotor(end+1) = ca(ii,imax);
end

[a,i0] = min(abs(ca(1,:)));
ptdc = [];
for ii=motorcyc
    ptdc(end+1) = pcyl(ii,i0);
end
ca_rng = ca(1,istart:iend);
pmotor_avg = mean(pcyl(motorcyc,istart:iend));
pfire_avg = mean(pcyl(firecyc,istart:iend));

% gIMEP
rng1 = find(ca(1,:)>=-180 & ca(1,:)<=180);
Wcig = trapz(V(rng1)',pcyl(firecyc,rng1)'.*1e5); % integral of pdV for all cycles [J]
gIMEP = Wcig./Vdisp; % [Pa]

%engine load
prsdata.gIMEP = gIMEP;ivalid = find(abs(gIMEP-mean(gIMEP))<10*std(gIMEP));
prsdata.mean_gIMEP = mean(gIMEP(ivalid));
prsdata.std_gIMEP = std(gIMEP(ivalid));

isse = []; iese =[];
% Energizing time
ienergize = find(injComm(2,:)>1);
dinjComm = diff(injComm(2,:));
isse = find(dinjComm>1);
iese = find(dinjComm<-1)+1;

if length(isse)>length(iese)
    isse=isse(1:length(iese));
end
if length(iese)>length(isse)
    iese=iese(1:length(isse));
end

try SSE_CAD = ca(1,isse(1)); catch SSE_CAD = 360; end

% figure,plot(ca(1,:),injComm(2,:))
% hold on
% plot(ca(1,isse),injComm(2,isse),'or')
% plot(ca(1,iese),injComm(2,iese),'or')
dse_ms = 1e3*(ca(1,iese)-ca(1,isse))/6/S;

%hydraulic delay

[SOInj,dpraildCA,pfiltAinj] = compSOInjdelay(lab,S,ca(1,:),prail,firecyc,f1,SSE_CAD );

inj_delay = SOInj - SSE_CAD; % CAD
prsdata.inj_delay = inj_delay; ivalid = find(abs(inj_delay)<6); %CAD
prsdata.mean_inj_delay = min(mean(inj_delay(ivalid)),3.5); %CAD min() may help prevent errors.
prsdata.std_inj_delay = std(inj_delay(ivalid)); %CAD

%SOC
try
[SOCI,EOCI,SOMCI] = compSOC(ca(1,:),dQdCA(firecyc,:),ca(1,isse(1)),prsdata.mean_inj_delay,3000); %f1 = 3000 for higher threshold
[SOCIfilt,EOCIfilt,SOMCIfilt] = compSOC(ca(1,:),dQdCAfilt(firecyc,:),ca(1,isse(1)),prsdata.mean_inj_delay,800);% f1 = 800 for lower threshold
 catch
     SOCI=0;EOCI=0;SOMCI=0;
     SOCIfilt=0;EOCIfilt=0;SOMCIfilt=0;
 end

%% start of combustion
prsdata.SOC = SOCI;ivalid = find(abs(SOCI)<25);
prsdata.mean_SOC = mean(SOCI(ivalid));
prsdata.std_SOC = std(SOCI(ivalid));

prsdata.EOC = EOCI; ivalid = find(abs(EOCI)<360);
prsdata.mean_EOC = mean(EOCI(ivalid));
prsdata.std_EOC = std(EOCI(ivalid));

prsdata.SOCfilt = SOCIfilt;ivalid = find(abs(SOCIfilt)<25);
prsdata.mean_SOCfilt = mean(SOCIfilt(ivalid));
prsdata.std_SOCfilt = std(SOCIfilt(ivalid));

prsdata.EOCfilt = EOCIfilt; ivalid = find(abs(EOCIfilt)<180);
prsdata.mean_EOCfilt = mean(EOCIfilt(ivalid));
prsdata.std_EOCfilt = std(EOCIfilt(ivalid));

%CA10/50/90


%MEPRR
if prsdata.mean_gIMEP/1000> -40
    
    if isfinite(prsdata.mean_SOCfilt)
    try    rng2 = find(ca(1,:)>=prsdata.mean_SOCfilt & ca(1,:)<=prsdata.mean_EOCfilt); catch rng2 = [1]; end  %Comeback and fix this... EE 2015/05/06
    try    rng3 = find(ca(1,:)>=prsdata.mean_SOC & ca(1,:)<=prsdata.mean_EOCfilt); catch rng3 = [1]; end %Comeback and fix this... EE 2015/05/06
    W1 = trapz(ca(1,rng3)',dQdCA(firecyc,rng3)'); % integral of dQdCA for all cycles [J]
    W2 = trapz(ca(1,rng2)',dQdCAfilt(firecyc,rng2)'); % integral of pdV for all cycles [J]
        try MEHRR = W1/(ca(1,rng3(end))-ca(1,rng3(1))); catch MEHRR = NaN; end%MEHRR = W1/(ca(1,rng2(end))-ca(1,rng2(1))); end
        try MEHRRfilt = W2/(ca(1,rng2(end))-ca(1,rng2(1))); catch MEHRRfilt = 10; end %maybe throwing an error

    else 
        MEHRR = NaN;
        try MEHRRfilt = W2/(ca(1,rng2(end))-ca(1,rng2(1))); catch MEHRRfilt = NaN; end
    end

else
    MEHRR = NaN;
    MEHRRfilt = NaN;
end

AID = 1e3*(SOCI - (SSE_CAD + prsdata.mean_inj_delay))./(6*S); % CAD /60s / 1200RPM = ms at 1200RPM
AIDfilt = 1e3*(SOCIfilt - (SSE_CAD + prsdata.mean_inj_delay))./(6*S); % CAD /60s / 1200RPM = ms at 1200RPM
%% %% Outputs
% Vectors
prsdata.pcyl        = pcyl;
prsdata.pfilt       = pfiltA;
prsdata.volume      = ones(size(ca,1),1)*V;
prsdata.ca          = ca;
prsdata.pint        = pint;
prsdata.injComm     = injComm;
prsdata.prail       = prail;
prsdata.injCurr     = injCurr;
prsdata.camTrig     = camTrig;
prsdata.dQdCA       = dQdCA;
prsdata.dpdCA       = dpdCA;
prsdata.dQdCAavg    = dQdCAavg;
prsdata.dpdCAavg    = dpdCAavg;
prsdata.dQdCAfilt   = dQdCAfilt;
prsdata.dpdCAfilt   = dpdCAfilt;
% 
% % Scalars
prsdata.Sset          = S; % speed RPM
prsdata.NC            = NC; % number cycles
prsdata.Rc            = Rc; % compression ratio set in labinfo.m
prsdata.mean_pint_bdc = pegp;
prsdata.mean_ptdc     = mean(ptdc);
prsdata.std_ptdc      = std(ptdc);
prsdata.f1            = f1; %start roll off filter frequency
prsdata.f2            = f2; %end roll off frequency
prsdata.nexp          = -nexpmeas;
prsdata.ncomp         = -ncompmeas;

prsdata.pmax = pmax; ivalid = find(abs(pmax-mean(pmax))<10*std(pmax));
prsdata.mean_pmax = mean(pmax(ivalid));
prsdata.std_pmax = std(pmax(ivalid));
prsdata.alpha_pmax = alpha_pmax; ivalid = find(abs(alpha_pmax-mean(alpha_pmax))<10*std(alpha_pmax));
prsdata.mean_alpha_pmax = mean(alpha_pmax(ivalid));
prsdata.std_alpha_pmax = std(alpha_pmax(ivalid));

%motored pressure data
prsdata.pmotor = pmotor; ivalid = find(abs(pmotor-mean(pmotor))<10*std(pmotor));
prsdata.mean_pmotor = mean(pmotor(ivalid));
prsdata.std_pmotor = std(pmotor(ivalid));
prsdata.alpha_pmotor = alpha_pmotor; ivalid = find(abs(alpha_pmotor-mean(alpha_pmotor))<10*std(alpha_pmotor));
prsdata.mean_alpha_pmotor = mean(alpha_pmotor(ivalid));
prsdata.std_alpha_pmotor = std(alpha_pmotor(ivalid));

prsdata.pmotor_avg = pmotor_avg;
prsdata.ca_rng = ca_rng;
prsdata.pfire_avg = pfire_avg;

%auto ignition delay
prsdata.AID = AID; ivalid = find(abs(AID)<10);
prsdata.mean_AID = mean(AID(ivalid));
prsdata.std_AID = std(AID(ivalid));

prsdata.AIDfilt = AIDfilt; ivalid = find(abs(AIDfilt)<10);
prsdata.mean_AIDfilt = mean(AIDfilt(ivalid));
prsdata.std_AIDfilt = std(AIDfilt(ivalid));

prsdata.DOC = prsdata.EOCfilt-prsdata.SOC; ivalid = find(abs(prsdata.DOC>0));
prsdata.mean_DOC = mean(prsdata.DOC(ivalid));
prsdata.std_DOC = std(prsdata.DOC(ivalid));

% max rate of aparent heat release
prsdata.AHRRmax = ahrrmax;ivalid = find(abs(ahrrmax-mean(ahrrmax))<10*std(ahrrmax));
prsdata.mean_AHRRmax = mean(ahrrmax(ivalid));
prsdata.std_AHRRmax = std(ahrrmax(ivalid));
prsdata.alpha_AHRRmax = alpha_ahrrmax;ivalid = find(abs(alpha_ahrrmax-mean(alpha_ahrrmax))<3*std(alpha_ahrrmax));
prsdata.mean_alpha_AHRRmax = mean(alpha_ahrrmax(ivalid));
prsdata.std_alpha_AHRRmax = std(alpha_ahrrmax(ivalid));

prsdata.MEHRR = MEHRR;
prsdata.mean_MEHRR = mean(MEHRR);
prsdata.std_MEHRR = std(MEHRR);

prsdata.MEHRRfilt = MEHRRfilt;
prsdata.mean_MEHRRfilt = mean(MEHRRfilt);
prsdata.std_MEHRRfilt = std(MEHRRfilt);

prsdata.AHRRmaxfilt = ahrrmaxfilt;ivalid = find(abs(ahrrmaxfilt-mean(ahrrmaxfilt))<10*std(ahrrmaxfilt));
prsdata.mean_AHRRmaxfilt = mean(ahrrmaxfilt(ivalid));
prsdata.std_AHRRmaxfilt = std(ahrrmaxfilt(ivalid));
prsdata.alpha_AHRRmaxfilt = alpha_ahrrmaxfilt;ivalid = find(abs(alpha_ahrrmaxfilt-mean(alpha_ahrrmaxfilt))<3*std(alpha_ahrrmaxfilt));
prsdata.mean_alpha_AHRRmaxfilt = mean(alpha_ahrrmaxfilt(ivalid));
prsdata.std_alpha_AHRRmaxfilt = std(alpha_ahrrmaxfilt(ivalid));

prsdata.dse_ms = dse_ms;
prsdata.sse_degV = ca(1,isse);
prsdata.ese_degV = ca(1,iese);
prsdata.sse_index = isse;
prsdata.ese_index = iese;
end