clear all

[filename,pathname] = uigetfile('*.mat','Select files','Multiselect','on');
filename = cellstr(filename);
legende={};
cm=colormap('jet');
for ii=1:length(filename)


icolor = round(63*(ii-1)/length(filename)+1);


    load(fullfile(pathname,filename{ii}))

NC = rundata.NC;
pcyl = rundata.pcyl;
pfilt = rundata.pfilt;

% for i = 1:30
% n = 4*i;
% x = pcyl(2,rundata.sse_index+n:rundata.sse_index+180*4)-pcyl(1,rundata.sse_index+n:rundata.sse_index+180*4);
% N = length(x);
% xdft = fft(x);
% xdft2 = xdft(1:N/2+1);
% psdx = (1/(2*pi*N)) * abs(xdft2).^2;
% psdx(2:end-1) = 2*psdx(2:end-1);
% freq = 0:(2*pi)/N:pi;
% 
% figure(4); hold on;
% plot(freq/pi*28800/2,10*log10(psdx),'color',[i/40 i/40 i/40])
% 
% end



figure(4);hold on;
for jj = 2:2:NC
x = pcyl(jj,:);
N = length(x);
xdft = fft(x);
xdft2 = xdft(1:N/2+1);
psdx = (1/(2*pi*N)) * abs(xdft2).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:(2*pi)/N:pi;
%figure(4); hold on;
%plot(freq/pi*28800/2,10*log10(psdx),'color',[0.6 0.6 0.6],'linewidth',2)
psdFire(jj,:) = psdx;
end
plot(freq/pi*28800/2,10*log10(mean(psdFire)),'color',cm(icolor,:),'linewidth',2)
alpha = [1.841 3.054 3.832 4.201 5.332]
fmodes = alpha.*sqrt((rundata.nexp+rundata.ncomp/2)*287*rundata.Ttdc)/pi/.1397

plot([fmodes],[0 0 0 0 0],'linestyle','none','marker','h','color',cm(icolor,:))

% x = mean(pcyl(2:2:NC,:));
% N = length(x);
% 
% xdft = fft(x);
% xdft2 = xdft(1:N/2+1);
% psdx = (1/(2*pi*N)) * abs(xdft2).^2;
% psdx(2:end-1) = 2*psdx(2:end-1);
% freq = 0:(2*pi)/N:pi;
% figure(4); hold on;
% plot(freq/pi*28800/2,10*log10(psdx),'color',[0.7 0.7 0.7],'linewidth',2)

for kk = 2:2:NC
    x = pfilt(kk,:);
    N = length(x);
    xdft = fft(x);
    xdft2 = xdft(1:N/2+1);
    psdx = (1/(2*pi*N)) * abs(xdft2).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:(2*pi)/N:pi;
    
%    figure(4); hold on;
%   plot(freq/pi*28800/2,10*log10(psdx),':','linewidth',2,'color',[0.2 0.2 0.2])
    psdfilt(kk,:) = psdx;
end
    figure(4); hold on;
    plot(freq/pi*28800/2,10*log10(mean(psdfilt)),'k','linewidth',2)

%     x = mean(pcyl(1:2:NC,:));
%     N = length(x);
% 
% xdft = fft(x);
% xdft2 = xdft(1:N/2+1);
% psdx = (1/(2*pi*N)) * abs(xdft2).^2;
% psdx(2:end-1) = 2*psdx(2:end-1);
% freq = 0:(2*pi)/N:pi;
% figure(4); hold on;
% plot(freq/pi*28800/2,10*log10(psdx),'k--','linewidth',2)

grid on
title('Periodogram Using FFT')
xlabel('Frequency')
ylabel('Power/Frequency (dB/rad/sample)')
end

axis([0 7500 -80 80])