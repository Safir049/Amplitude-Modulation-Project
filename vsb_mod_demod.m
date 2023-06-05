function [s_rec] = vsb_mod_demod(nBits,Fs,m_sig,BW,fc,percentageOfLsb)
%%% BW = LSB+USB

%% generation of message signal and carrier signal
ts = 1/Fs;
t = linspace(0,length(m_sig)/Fs,length(m_sig));

Lm_sig=length(m_sig);
Lfft=length(t);
Lfft=2^ceil(log2(Lfft)); 
M_sig=fftshift(fft(m_sig,Lfft)); 

freqm=(-Lfft/2:Lfft/2-1)/(Lfft*ts) ;

h=fir1(40,[BW*ts]);

%% VSB Modulation
s_dsb=(m_sig)'.*cos(2*pi*fc*t); %transposed message signal here
Lfft=length(t); 
Lfft=2^ceil(log2(Lfft)+1);
S_dsb=fftshift(fft(s_dsb,Lfft)); 
freqs=(-Lfft/2:Lfft/2-1)/(Lfft*ts);

fv = percentageOfLsb*BW/2;
L_vsb=floor((fc-fv)*ts*Lfft);
hfv = floor((fc+fv)*ts*Lfft);

%setting up the vsb filter
VSBfilt=ones(1,Lfft); 

p = Lfft/2+L_vsb:Lfft/2+hfv;
q = Lfft/2-hfv+1:Lfft/2-L_vsb+1;

x = linspace(fc-fv,fc+fv,length(p));
m = 0.5/(2*fv);
y = m.*x - m*(fc-fv);

VSBfilt(p)=y;
VSBfilt(q)=fliplr(y);

VSBfilt(Lfft/2-L_vsb+1:Lfft/2+L_vsb)=zeros(1,2*L_vsb); 


S_vsb=S_dsb.*VSBfilt;
s_vsb=real(ifft(fftshift(S_vsb)));
s_vsb=s_vsb(1:Lm_sig);

%% demodulation
%%% Demodulation begins by using a rectifier 
s_dem=s_vsb.*cos(2*pi*fc*t)*2; 
S_dem=fftshift(fft(s_dem,Lfft)); 
% Using an ideal low pass filter
s_rec=filter(h,1,s_dem);

%need to design only this filter
S_rec=fftshift(fft(s_rec,Lfft)); 
 


% Trange=[-0.025 0.025 -1 1] ;
Frange=[-5000 5000 0 200] ;%


figure(2)
subplot(221); plot(t,m_sig,'Linewidth',1.5)
% axis(Trange) 
title('message signal')

subplot(222); plot(t,s_vsb,'Linewidth',1.5)
% axis(Trange)
title('VSB-SC modulated signal')

subplot(223); plot(t, s_dem,'Linewidth',1.5)
% axis(Trange) 
title('After multiplying local carrier')

subplot(224); plot(t,s_rec,'Linewidth',1.5)
% axis(Trange)
title('Recovered signal')

figure(3)
subplot(221); plot(freqm,abs(M_sig),'Linewidth',1.5)
axis(Frange)
title('Message Spectrum')

subplot(222); plot(freqs,abs(S_vsb),'Linewidth',1.5)
axis(Frange)
title('Upper Sideband VSB-SC spectrum')

subplot(223); plot(freqs, abs(S_dem),'Linewidth',1.5)
axis(Frange) 
title('Detector spectrum')

subplot(224); plot(freqs,abs(S_rec),'Linewidth',1.5)
axis(Frange)
title('Recovered signal')

end

