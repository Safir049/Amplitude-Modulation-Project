clc;
close all;
clear all;

%% vsb_using _functions

[y, Fs] = audioread('imtiaz_8000_1s.wav');
BW = 2000;
fc = 3000;%fc must be less than equal 3000
nBits = 16;
percentageOfLsb = 0.25;
[s_rec] = vsb_mod_demod(nBits,Fs,y,BW,fc,percentageOfLsb);
% sound(y,Fs,nBits)
% sound(s_rec,Fs,nBits)




% BW = 150;
% Fs =10000;
% tmax = 0.04;
% ts=1/Fs;
% t=-0.04:ts:0.04;
% Ta=0.01;
% m_sig=triangl((t+0.01)/0.01)-triangl((t-0.01)/0.01);
% fc = 300;
% percentageOfLsb = 0.25;
% [s_rec] = vsb_mod_demod(Fs,tmax,m_sig,BW,fc,percentageOfLsb);