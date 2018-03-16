clear all

load('Ispray_700K.mat') %spray intensity profiles
load('Ispray_800K.mat')
load('Ispray_900K.mat')
load('Ispray_1000K.mat')
load('coordinates.mat') % x and r vectors 
% x vector 0,200, cell centered position for 0.5mm cells from 0.25 to 99.75mm
% r vector 0 121, cell centered position for 0.5mm cells -30 to +30mm

ylimits = [0 10000]
[f1,ax1,ax2,ax3,ax4,ax5] = createSchprofilefigure(ylimits)

axes(ax1), hold on
plot(x_ref, Ispray_800K(:,61),'b--')
plot(x_ref, Ispray_900K(:,61),'k--')
plot(x_ref, Ispray_1000K(:,61),'r--')
xlim([0 55])

axes(ax2), hold on
%plot(r_ref, Ispray_700K(:,20),'.')
plot(r_ref, Ispray_800K(20,:),'b--')
plot(r_ref, Ispray_900K(20,:),'k--')
plot(r_ref, Ispray_1000K(20,:),'r--')
xlim([-15 15])
%x_ref(20) = 9.75mm


axes(ax3), hold on
%plot(r_ref, Ispray_700K(:,40),'.')
plot(r_ref, Ispray_800K(40,:),'b--')
plot(r_ref, Ispray_900K(40,:),'k--')
plot(r_ref, Ispray_1000K(40,:),'r--')
xlim([-15 15])

%x_ref(40) = 19.75mm

axes(ax4), hold on
%plot(r_ref, Ispray_700K(:,60),'.')
plot(r_ref, Ispray_800K(60,:),'b--')
plot(r_ref, Ispray_900K(60,:),'k--')
plot(r_ref, Ispray_1000K(60,:),'r--')
xlim([-15 15])

axes(ax5), hold on
%plot(r_ref, Ispray_700K(:,60),'.')
plot(r_ref, Ispray_800K(80,:),'b--')
plot(r_ref, Ispray_900K(80,:),'k--')
plot(r_ref, Ispray_1000K(80,:),'r--')
xlim([-15 15])
x_ref(80)

figure, hold on
plot(x_ref,sum(Ispray_800K,2))
plot(x_ref,sum(Ispray_900K,2))
plot(x_ref,sum(Ispray_1000K,2))
xlim([0 55])