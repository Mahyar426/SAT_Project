% Simulation for delay estimation through 2D-correlation

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
%% Transmitter side: modulate the PN sequence
SpS=2; % related to f_chip = 2MHz?
signalRect=rectpulse(Code,SpS);
figure,stem(signalRect(1:discreteTimePlot)),grid on,axis('padded');
title('Rectangularly-shaped signal | SpS=2;','Interpreter','latex');
xlabel('[n]','Interpreter','latex');
ylabel('y[n]','Interpreter','latex');
% PSD estimation of baseband signal
N=length(signalRect);
Nwel=N/10;          % Length of thw window
h=ones(1,Nwel);     % Rectangular window to pre-filter
Noverlap=Nwel/2;    % Number of overlapping samples
Nfft=4096;          % Number of FFT points per window
Fs=4e+06;
[Px,f]=pwelch(signalRect,h,Noverlap,Nfft,Fs,'centered');
figure,semilogy(f,Px),axis('tight');
grid on
xlabel('Frequency [GHz]','Interpreter','latex');
ylabel('$\hat{P} (f)$','Interpreter','latex');
title('Estimated PSD of rectangularly-shaped signal | Baseband','Interpreter','latex')
% Sinusoidal shaping and bandpass modulation
freqCarrier=2.1e+09;
freqSampl=5e+09;
% We use: modulate() = x.*cos(2*pi*freqCarrier*t)
% which applies amplitude modulation at the desired frequency
signalRF=modulate(signalRect,freqCarrier,freqSampl,'am');
figure,stem(signalRF(1:discreteTimePlot*10)),grid on,axis('padded');
xlabel('[n]','Interpreter','latex');
ylabel('y[n]','Interpreter','latex');
title('Modulated signal | $f_c=2.1 GHz$','Interpreter','latex');
% ESTIMATED PSD TO BE COMPUTED AND PLOTTED
%% Channel model: phase shift and frequency shift
shiftDoppler=500e+03;
shiftPhase=500;
signalRx=circshift(signalRF.*exp(2i*pi*shiftDoppler),shiftPhase);
%% Receiver side: de-modulate signalRx and apply 2D correlation
