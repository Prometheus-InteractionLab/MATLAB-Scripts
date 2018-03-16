clear all
%close all
clc
load('hsvlinearcmap.mat')
colormap(hsvlinearcmap)

%% Selection of the 2D spray model
[filename,pathname] = uigetfile('*.mat','Select the modeled spray');
load(fullfile(pathname,filename));

% 2D spray mixture fraction field
figure,imagesc(data.mixing.r_m*1000,data.mixing.x_m*1000,data.mixing.mole_fraction.field),caxis([0 0.01])
axis image
set(gcf,'color','w')
set(gca,'FontName','eras light itc','fontSize',14,'fontWeight','normal','box','on','linewidth',0.5)
xlabel('radial distance [mm]')
ylabel('axial distance [mm]')
cb=colorbar;
set(get(cb,'ylabel'),'String','Mixture fraction [-]','FontName','eras light itc','fontSize',14,'fontWeight','normal','linewidth',0.5)
set(cb,'FontName','eras light itc','fontSize',14,'fontWeight','normal','box','on','linewidth',0.5)

% 2D spray temperature field
figure,imagesc(data.mixing.r_m*1000,data.mixing.x_m*1000,data.mixing.Tadiab.field)
axis image
set(gcf,'color','w')
set(gca,'FontName','eras light itc','fontSize',14,'fontWeight','normal','box','on','linewidth',0.5)
xlabel('radial distance [m]')
ylabel('axial distance [m]')
cb=colorbar;
set(get(cb,'ylabel'),'String','Adiab. Mixing Temp [K]','FontName','eras light itc','fontSize',14,'fontWeight','normal','linewidth',0.5)
set(cb,'FontName','eras light itc','fontSize',14,'fontWeight','normal','box','on','linewidth',0.5)

%% Construction of 3D fields
[rr,rr] = meshgrid(data.mixing.r_m,data.mixing.r_m);
[theta,rho] = cart2pol(rr,rr');

% 3D field of radial distance
r3d = zeros(length(data.mixing.x_m),length(data.mixing.r_m),length(data.mixing.r_m));
for ii=1:length(data.mixing.x_m)
    r3d(ii,:,:) = rho;
end

% 3D field of mass fraction
X3d = zeros(length(data.mixing.x_m),length(data.mixing.r_m),length(data.mixing.r_m));
for ii=1:length(data.mixing.x_m)
    X2d = zeros(length(data.mixing.r_m),length(data.mixing.r_m));
    r2d = reshape(r3d(ii,:,:),length(data.mixing.r_m),length(data.mixing.r_m));
    r = data.mixing.r_m;ipos=find(r>=0);
    X = data.mixing.mole_fraction.field(ii,:);
    X2d(:) = interp1(r(ipos),X(ipos),r2d(:));
    X2d(isnan(X2d))=0;
    X3d(ii,:,:) = X2d;
end

% 3D field of temperature
T3d = zeros(length(data.mixing.x_m),length(data.mixing.r_m),length(data.mixing.r_m));
for ii=1:length(data.mixing.x_m)
    T2d = zeros(length(data.mixing.r_m),length(data.mixing.r_m));
    r2d = reshape(r3d(ii,:,:),length(data.mixing.r_m),length(data.mixing.r_m));
    r = data.mixing.r_m;ipos=find(r>=0);
    T = data.mixing.Tadiab.field(ii,:);
    T2d(:) = interp1(r(ipos),T(ipos),r2d(:));
    T2d(isnan(T2d))=data.ambient.T_K;
    T3d(ii,:,:) = T2d;
end
  
% 3D field of concentration
C3d = X3d*data.ambient.p_bar*1e5/8.314./T3d;

%% Black body emission
h = 6.63e-34;   % Planck constant (J.s)
c = 3e8;        % Speed of light (m/s)
k = 1.38e-23;   % Boltzman constant (J/K)

Iblackbody = @(lambda,temp) 2*pi*h*c^2./(lambda.^5.*(exp(h*c./(lambda.*temp*k))-1)); % Blackbody radiance [W/m3]
lambda = [200:10:50000]*1e-9;
T = [300:50:1000 1000:500:10000];

figure(10),set(gcf,'color','w');
grid on,hold on,cm=colormap(jet);set(gca,'FontName','eras light itc','fontSize',14,'fontWeight','normal','box','on','linewidth',0.5)
figure(20),set(gcf,'color','w');
grid on,hold on,cm=colormap(jet);set(gca,'FontName','eras light itc','fontSize',14,'fontWeight','normal','box','on','linewidth',0.5)

lambdamax=[];IlT=[];legende={};
for ii=1:length(T)
    icolor = round(63*(ii-1)/(length(T)-1)+1);
    I = Iblackbody(lambda,T(ii));

    dl = diff(lambda);dl=dl(1);
    intI = cumsum(I)*dl;
    
    [maxI,imax]=max(I);
    lambdamax(ii) = lambda(imax);
    
    figure(10)
    plot(lambda*1e6,I,'color',cm(icolor,:));
    figure(20)
    plot(lambda*1e6,intI/intI(end),'color',cm(icolor,:));
    IlT=[IlT;I];
    legende=[legende;['T=' num2str(T(ii)) 'K']];
end

figure(10),xlabel('wavelength [µm]'),ylabel('Spectral Radiance [W/m^{-3}]')
ht=legend(legende);set(ht,'visible','off');
figure(20),xlabel('wavelength [µm]'),ylabel('Integrated emission/Total emission [%]')
ht=legend(legende);set(ht,'visible','off');

%% Computation of the integrated radiance between 3-5µm for a 2D plane of the spray 
em   = 1; % Emissivity
lmin = 3.30; %Minimum wavelength [µm]
lmax = 3.50; %Maximum wavelength [µm]
dl = diff(lambda);dl=dl(1);
[a,imin] = min(abs(lambda-lmin*1e-6));
[a,imax] = min(abs(lambda-lmax*1e-6));
power35 = @(T) sum(Iblackbody(lambda(imin:imax),T))*dl; % Integrated radiance between 3 and 5µm as a function of temperature 
% power35_v2 = @(temp) quad(@(lambda) Iblackbody(lambda,temp),3e-6,5e-6);

radiance_2d_3_36mu = 0*data.mixing.Tadiab.field;
for ii=1:length(data.mixing.Tadiab.field(:))
    ii
    radiance_2d_3_36mu(ii) = power35(data.mixing.Tadiab.field(ii)); % 2D field of radiance for a spray section
end

figure,imagesc(data.mixing.r_m,data.mixing.x_m,radiance_2d_3_36mu)
colormap(hsvlinearcmap)
axis image
set(gcf,'color','w')
set(gca,'FontName','eras light itc','fontSize',14,'fontWeight','normal','box','on','linewidth',0.5)
xlabel('radial distance [m]')
ylabel('axial distance [m]')
cb=colorbar;
set(get(cb,'ylabel'),'String','Planar Cross Section of Radiance in range 3.3-3.5µm [W/m^2]','FontName','eras light itc','fontSize',14,'fontWeight','normal','linewidth',0.5)
set(cb,'FontName','eras light itc','fontSize',14,'fontWeight','normal','box','on','linewidth',0.5)

% 3D field of radiance (integrated between 3.3 and 3.5µm)
tref = [500:0.1:1000];
eref = arrayfun(power35,tref); % Radiance between 3-5µm as a function of temperature

E3d = zeros(length(data.mixing.x_m),length(data.mixing.r_m),length(data.mixing.r_m));
for ii=1:length(data.mixing.x_m)
    ii
    E2d = zeros(length(data.mixing.r_m),length(data.mixing.r_m));
    T2d = reshape(T3d(ii,:,:),length(data.mixing.r_m),length(data.mixing.r_m));
    r = data.mixing.r_m;ipos=find(r>=0);
    T = data.mixing.Tadiab.field(ii,:);
    E2d(:) = interp1(tref,eref,T2d(:));    
    E3d(ii,:,:) = E2d;
end

%% Integrated emission along a radial distance of the spray (see formula on slide 19 of summary_IR_v03.pptx)
I0_cyl = power35(250+373);
IBG = zeros(length(data.mixing.x_m),length(data.mixing.r_m));
Ispray = zeros(length(data.mixing.x_m),length(data.mixing.r_m));

Efuel2d = zeros(length(data.mixing.x_m),length(data.mixing.r_m));
Afuel2d = zeros(length(data.mixing.x_m),length(data.mixing.r_m));
an3d = zeros(length(data.mixing.x_m),length(data.mixing.r_m),length(data.mixing.r_m));
bn3d = zeros(length(data.mixing.x_m),length(data.mixing.r_m),length(data.mixing.r_m));
em3d = zeros(length(data.mixing.x_m),length(data.mixing.r_m),length(data.mixing.r_m));

abs_fuel = 70; % absorption coefficient (mol/m^2)
dr = diff(data.mixing.r_m);dr=dr(1);

for ii=1:length(data.mixing.x_m)
    ii   
    % Matrice 3d des coefs an et bn
    an = 1-abs_fuel*C3d*dr;
    bn = abs_fuel*C3d.*E3d*dr;
    
    % Matrice 3D d'emission/abosrption
    an2d = reshape(an(ii,:,:),length(data.mixing.r_m),length(data.mixing.r_m));
    bn2d = reshape(bn(ii,:,:),length(data.mixing.r_m),length(data.mixing.r_m));    
    an3d(ii,:,:) = an2d;
    bn3d(ii,:,:) = bn2d;
      
    P1 = flipud(cumprod(an2d,1));
    P2 = P1.*bn2d;
    S1 = cumsum(P2,1);
    
    em3d(ii,:,:) = P2;
    Efuel2d(ii,:) = S1(end,:);
    Afuel2d(ii,:) = P1(1,:);
    
    IBG(ii,:)    = I0_cyl*P1(1,:);
    Ispray(ii,:) = S1(end,:);
end

figure,imagesc(data.mixing.r_m,data.mixing.x_m,Ispray),caxis([0 800])
colormap(hsvlinearcmap)
axis image
set(gcf,'color','w')
set(gca,'fontSize',14,'fontWeight','normal','box','on','linewidth',0.5) 
xlabel('radial distance [m]')
ylabel('axial distance [m]')
cb=colorbar;
set(get(cb,'ylabel'),'String','Spray self-emission [W/m^2]','fontSize',14,'fontWeight','normal','linewidth',0.5)
set(cb,'fontSize',14,'fontWeight','normal','box','on','linewidth',0.5)

% figure,imagesc(data.mixing.r_m,data.mixing.x_m,100-Afuel2d*100),caxis([0 100])
% axis image
% set(gcf,'color','w')
% set(gca,'FontName','eras light itc','fontSize',14,'fontWeight','normal','box','on','linewidth',0.5) 
% xlabel('radial distance [m]')
% ylabel('axial distance [m]')
% cb=colorbar;
% set(get(cb,'ylabel'),'String','Spray total absorption [%]','FontName','eras light itc','fontSize',14,'fontWeight','normal','linewidth',0.5)
% set(cb,'FontName','eras light itc','fontSize',14,'fontWeight','normal','box','on','linewidth',0.5)

% figure,imagesc(data.mixing.r_m,data.mixing.x_m,IBG),caxis([0 8000])
% axis image
% set(gcf,'color','w')
% set(gca,'FontName','eras light itc','fontSize',14,'fontWeight','normal','box','on','linewidth',0.5) 
% xlabel('radial distance [m]')
% ylabel('axial distance [m]')
% cb=colorbar;
% set(get(cb,'ylabel'),'String','Wall emission @250°C [W/m^2]','FontName','eras light itc','fontSize',14,'fontWeight','normal','linewidth',0.5)
% set(cb,'fontweight','bold','linewidth',2,'fontname','calibri','fontsize',16)

% figure,imagesc(data.mixing.r_m,data.mixing.x_m,100*(Ispray-IBG)./Ispray),caxis([0 100])
% axis image
% set(gcf,'color','w')
% set(gca,'FontName','eras light itc','fontSize',14,'fontWeight','normal','box','on','linewidth',0.5) 
% xlabel('radial distance [m]')
% ylabel('axial distance [m]')
% cb=colorbar;
% set(get(cb,'ylabel'),'String','Relevant signal [%]','FontName','eras light itc','fontSize',14,'fontWeight','normal','linewidth',0.5)
% set(cb,'FontName','eras light itc','fontSize',14,'fontWeight','normal','box','on','linewidth',0.5)

% figure,imagesc(data.mixing.r_m,data.mixing.x_m,Ispray./IBG),caxis([0 20])
% axis image
% set(gcf,'color','w')
% set(gca,'FontName','eras light itc','fontSize',14,'fontWeight','normal','box','on','linewidth',0.5) 
% xlabel('radial distance [m]')
% ylabel('axial distance [m]')
% cb=colorbar;
% set(get(cb,'ylabel'),'String','S/N [-]','FontName','eras light itc','fontSize',14,'fontWeight','normal','linewidth',0.5)
% set(cb,'FontName','eras light itc','fontSize',14,'fontWeight','normal','box','on','linewidth',0.5)

figure,imagesc(data.mixing.r_m,data.mixing.r_m,P2)%,caxis([0 200])
axis image
set(gcf,'color','w')
set(gca,'FontName','eras light itc','fontSize',14,'fontWeight','normal','box','on','linewidth',0.5) 
xlabel('radial distance [m]')
ylabel('radial distance [m]')
cb=colorbar;
set(get(cb,'ylabel'),'String','Spray self-emission [W/m^2/m^-^1]','FontName','eras light itc','fontSize',14,'fontWeight','normal','linewidth',0.5)
set(cb,'FontName','eras light itc','fontSize',14,'fontWeight','normal','box','on','linewidth',0.5)

%% Emission from a sample volume of fuel, whose temperature evolves with concentration (adiabatic mixing)
L = 5e-4;  % Optical path [m]
Icoll=[];eps1=7;Il = [];C = [];
for ii=1:length(data.mixing.Tadiab.axial)
    Xf0 = data.mixing.mole_fraction.axial(ii);
    Tr = data.mixing.Tadiab.axial(ii);
    C(ii) = Xf0 * data.ambient.p_bar*1e5/8.314/data.ambient.T_K;
    Ebb = radiance_2d_3_36mu(ii,(length(data.mixing.r_m)-1)/2);
    Il(ii) = Ebb;
    Icoll(ii) = Ebb*(1-exp(-C(ii)*eps1*L));    
end
phi = data.mixing.mole_fraction.axial*data.fuel.rho_kg_m3./(1-data.mixing.mole_fraction.axial)/data.ambient.rho_kg_m3*data.fuel.pco_air;

figure,set(gcf,'color','w'),grid on,hold on
plot(phi,data.mixing.Tadiab.axial,'s','linewidth',3,'color',[1 1 1]*0.7)
plot(phi,Icoll,'o','linewidth',3,'color',[0 0 0])
plot(phi,C,'x')
%plot(phi,Icoll/max(Icoll(end))*data.ambient.T_K,'o','linewidth',3,'color',[0 0 0])
plot(phi,Il,'g')
% plot(phi,C/max(C)*900,'m')
set(gca,'FontName','eras light itc','fontSize',14,'fontWeight','normal','box','on','linewidth',0.5) 
xlabel('\Phi [-]')
set(gca,'xlim',[0 10])
legend('Temperature [K]','Norm. radiance [-]')