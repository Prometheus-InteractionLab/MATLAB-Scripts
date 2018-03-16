% Return the PSD [bar^2/Hz]
%
function [PSD,f] = Power_Spectral_Density(P_cyl_kPa,fs,dt)
%   P_cyl_kPa : Pcyl
%   fs :        Sampling frequency [Hz]
%   dt :        Seconds/sample using SPC/2 samples per rev (SPC: Sample Per Cycle), dt = 1/(RPM*SPC/2/60);
%
nFFT = length(P_cyl_kPa);                   % n for n-point-FFT
FT(:,1) = fft(P_cyl_kPa/100,nFFT);          % Calculate FFT
nb_pts = nFFT/2+1;                          % # of points in one-sided spectrum
PSD(:,1) = (FT(:,1).*conj(FT(:,1)))/nFFT^2; % amplitude - note you use NFFT points because energy is spread over that number of frequencies
Delta_f = fs/2/(nb_pts-1);                  % Frequency differential element
PSD(:,1) = PSD(:,1)*2/Delta_f;              % the 2 corrects to a two-sided spectrum
f = (0:nFFT/2)/nFFT/dt;                     % freq in [samples/sec]
%
end