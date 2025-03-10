% Simulation of GNSS receiver acquisition stage 

clc
clear
close all

%% Initializing circular shift registers 
C1 = [+1 -1];
C2 = [+1 +1 +1 -1 -1 +1 -1];
C3 = [+1 +1 +1 -1 -1 -1 +1 -1 +1 +1 -1];
C4 = [+1 +1 +1 +1 -1 -1 -1 +1 -1 -1 +1 +1 -1 +1 -1];
C5 = [+1 +1 +1 +1 -1 +1 -1 +1 -1 -1 -1 -1 +1 +1 -1 +1 +1 -1 -1];
C6 = [+1 +1 +1 +1 +1 -1 +1 -1 +1 +1 -1 -1 +1 +1 -1 -1 +1 -1 +1 -1 -1 -1 -1];
counter=1;
regLen=[2 7 11 15 19 23];
codeLen=prod(regLen);
%% PN sequence generation
Code=zeros(codeLen,1);
while counter<=codeLen
    Code(counter)=sign(4*C1(1)+C2(1)-C3(1)-C4(1)+C5(1)-C6(1));
    C1 = circshift(C1,1);
    C2 = circshift(C2,1);
    C3 = circshift(C3,1);
    C4 = circshift(C4,1);
    C5 = circshift(C5,1);
    C6 = circshift(C6,1);
    counter=counter+1;
end
discreteTimePlot=50;
figure,stem(Code(1:discreteTimePlot)),grid on,axis('padded');
title('T4B Pseudo-Noise sequence','Interpreter','latex');
xlabel('[n]','Interpreter','latex');
%% Baseband signal generation 
SpS=4;                      % Arbitrarily chosen here
Fchip=2e+06;                % From CCSDS 414.1-B-2 - Table A-1
signalRect=rectpulse(Code,SpS);
% PSD estimation of baseband signal
N=length(signalRect);
Nwel=N/10;                  % Length of the window
h=ones(1,Nwel);             % Rectangular window to pre-filter
Noverlap=Nwel/2;            % Number of overlapping samples
Nfft=4096;                  % Number of FFT points per window
Fs=SpS*Fchip;               % Optimal sampling rate at baseband
[Px,f]=pwelch(signalRect,h,Noverlap,Nfft,Fs,'centered');
figure,plot((f./1e+07).*10,pow2db(Px)),axis('tight'),grid on;
xline(0,'LineWidth',1.5,'LineStyle','--','Color','black');
xlabel('Frequency [MHz]','Interpreter','latex');
ylabel('Magnitude [dB]','Interpreter','latex');
title('$\hat{P} (f)$ of rectangularly-shaped signal | Baseband','Interpreter','latex');
%% RF-IF front-end operations
freqIF=10e+06;              % From CCSDS 414.0-G-2 - Section 2.2.5
freqSamplIF=30e+06;         % > 24 MHz from IF spectrum plot
% Generating IF signal
SpS=ceil(freqSamplIF*(1/Fchip));
signalRect=rectpulse(Code,SpS);
signalIF=modulate(signalRect,freqIF,freqSamplIF,'am');
% PSD estimation of IF signal
N=length(signalRect);
Nwel=N/10;                  % Length of the window
h=ones(1,Nwel);             % Rectangular window to pre-filter
Noverlap=ceil(Nwel/2);      % Number of overlapping samples
Nfft=4096;                  % Number of FFT points per window
Fs=freqSamplIF;             % Optimal sampling rate
[Px,f]=pwelch(signalIF,h,Noverlap,Nfft,Fs,'centered');
figure,plot((f./1e+07).*10,pow2db(Px)),axis('tight'),grid on;
xlim([-15 15]);
xline(0,'LineWidth',1.5,'LineStyle','--','Color','black');
xline(10,'LineWidth',1.5,'LineStyle','--','Color',[0.8500 0.3250 0.0980]);
xline(-10,'LineWidth',1.5,'LineStyle','--','Color',[0.8500 0.3250 0.0980]);
xlabel('Frequency [MHz]','Interpreter','latex');
ylabel('Magnitude [dB]','Interpreter','latex');
title('$\hat{P} (f)$ at Intermediate Frequency (IF) | $f_{carrier}=$ $ 10$ $MHz$','Interpreter','latex');