[filename,pathname] = uigetfile('*.mat','Select files','Multiselect','on');
filename = cellstr(filename);
legende={};
for ii=1:length(filename)
    load(fullfile(pathname,filename{ii}))
    
NC = size(TCdata.T0.m,1);
M = TCdata.T0.m;
F = TCdata.T0.f;

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
figure(1)
plot(TCdata.ca,mean(F),'r'); hold on;
plot(TCdata.ca,mean(M),'k--');

figure(2)
plot(TCdata.ca,mean(F)-mean(M),'k--');

x = mean(M);
N = length(x);

xdft = fft(x);
xdft2 = xdft(1:N/2+1);
psdx = (1/(2*pi*N)) * abs(xdft2).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:(2*pi)/N:pi;
figure(10*ii+4); hold on;
plot(freq/pi*28800/2,10*log10(psdx),'color',[0.3 0.3 0.3],'linewidth',2)
jj=2

for jj=1:NC
icolor = round(63*(jj-1)/(NC-1)+1);
cm=colormap('jet');
x = F(jj,:)-mean(M);
xdft = fft(x);
xdft2 = xdft(1:N/2+1);
psdx = (1/(2*pi*N)) * abs(xdft2).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:(2*pi)/N:pi;
figure(10*ii+4); hold on;
plot(freq/pi*28800/2,10*log10(psdx),'color',cm(icolor,:))
psdTC(:,jj) = psdx;
end

plot(freq/pi*28800/2,10*log10(mean(psdTC,2)),'g')

x = mean(F)-mean(M);
xdft = fft(x);
xdft2 = xdft(1:N/2+1);
psdx = (1/(2*pi*N)) * abs(xdft2).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:(2*pi)/N:pi;
figure(10*ii+4); hold on;
plot(freq/pi*28800/2,10*log10(psdx),'k--','linewidth',3)

% Acoustic Mode shapes for cylindrical cylinder
 alpha = [1.841 3.054 3.832 4.201 5.332]
 fmodes = alpha.*sqrt((TCdata.nexp+TCdata.ncomp/2)*287*TCdata.Ttdc)/pi/.1397
 plot([fmodes],[0 0 0 0 0],'linestyle','none','marker','h')
end
axis([0 2400 -80 60])
% 
% x = pfilt(2,:);
% N = length(x);
% xdft = fft(x);
% xdft2 = xdft(1:N/2+1);
% psdx = (1/(2*pi*N)) * abs(xdft2).^2;
% psdx(2:end-1) = 2*psdx(2:end-1);
% freq = 0:(2*pi)/N:pi;
% 
% figure(4); hold on;
% plot(freq/pi*28800/2,10*log10(psdx),'r:')
% grid on
% title('Periodogram Using FFT')
% xlabel('Frequency')
% ylabel('Power/Frequency (dB/rad/sample)')