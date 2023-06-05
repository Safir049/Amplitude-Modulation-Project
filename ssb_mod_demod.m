function [s_rec] = ssb_mod_demod(nBits, m_sig, BW, fc)
nBits = 16;

ts=1/Fs;
t = linspace(0,length(y_tr)/Fs,length(y_tr));


% % t=-0.04:ts:0.04;
% Ta=0.01;
% m_sig=y_tr;%triangl((t+0.01)/0.01)-triangl((t-0.01)/0.01);
Lm_sig=length(m_sig);
Lfft=length(t);
Lfft=2^ceil(log2(Lfft)); 
M_sig=fftshift(fft(m_sig,Lfft)); 

freqm=(-Lfft/2:Lfft/2-1)/(Lfft*ts) ;
% BW=150;

h=fir1(40,[BW*ts]);



%% SSB Modulation
% fc=300; % carrier frequency
s_ssb=(m_sig)'.*cos(2*pi*fc*t);
Lfft=length(t); 
Lfft=2^ceil(log2(Lfft)+1);
S_ssb=fftshift(fft(s_ssb,Lfft)); 
freqs=(-Lfft/2:Lfft/2-1)/(Lfft*ts);

 L_ssb=floor(fc*ts*Lfft);


SSBfilt=ones(1,Lfft); 



SSBfilt(Lfft/2-L_ssb+1:Lfft/2+L_ssb)=zeros(1,2*L_ssb); 


S_ssb=S_ssb.*SSBfilt;
s_ssb=real(ifft(fftshift(S_ssb)));
s_ssb=s_ssb(1:Lm_sig);



%% demodulation
%%% Demodulation begins by using a rectifier 
s_dem=s_ssb.*cos(2*pi*fc*t)*2; 
S_dem=fftshift(fft(s_dem,Lfft)); 
% Using an ideal low pass filter with bandwidth 150 Hz
s_rec=filter(h,1,s_dem);

%need to design only this filter
S_rec=fftshift(fft(s_rec,Lfft)); 
 


% Trange=[-0.025 0.025 -1 1];
% Frange=[-700 700 0 200] ;
% figure(1)
% subplot(221); plot(t,m_sig,'Linewidth',1.5)
% % axis(Trange) 
% title('message signal')
% 
% subplot(222); plot(t,s_ssb,'Linewidth',1.5)
% % axis(Trange)
% title('SSB-SC modulated signal')
% 
% subplot(223); plot(t, s_dem,'Linewidth',1.5)
% % axis(Trange) 
% title('After multiplying local carrier')
% 
% subplot(224); 
% figure(2)
% plot(t,s_rec,'Linewidth',1.5)
% % axis(Trange)
%  title('Recovered signal')
%  sound(s_rec,Fs,16)
% 
% 
figure(2)
subplot(221); plot(freqm,abs(M_sig),'Linewidth',1.5);
% axis(Frange)
title('Message Spectrum')
% 
% subplot(222); plot(freqs,abs(S_ssb),'Linewidth',1.5)
% % axis(Frange)
% title('Upper Sideband SSB-SC spectrum')
% 
% subplot(223); plot(freqs, abs(S_dem),'Linewidth',1.5)
% % axis(Frange) 
% title('Detector spectrum')
% 
% subplot(224);

% plot(freqs,abs(S_rec),'Linewidth',1.5)
% % axis(Frange)
% title('Recovered signal')


%%% Defining triangl function  used in above code
% triangl(t)=1-|t| , if |t|<1
% triangl(t)=0 , if |t|>1


% function y = triangl(t)
% y=(1-abs(t)).*(t>=-1).*(t<1); % i.e. setting y to 1 -|t|  if  |t|<1 and to 0 if not
% %end
% 
% %%% example usage
% % t=-5:.1:5
% % y=triangl(t)
% % stem(t,y)
%  sound(s_rec,Fs,16)