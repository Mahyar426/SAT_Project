% Code for regenerative Pseudo-Noise sequence (T4B)

clc
clear
close all

%% Parameter setting for close-to-baseband simulation
freqCarrier=10e+06;
freqSampling=30e+06;
freqChip=2e+06;
SpS=freqSampling*(1/freqChip);
shiftDoppler=0;
freqCarrierShifted=freqCarrier+shiftDoppler;    % Channel effect on carrier frequency
freqDopplerTest=0; % to be put inside loop
freqTest=freqCarrier+freqDopplerTest;           % Test frequency in Acquisition stage

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
%% Plotting normalized circular autocorrelation
R=ifft(fft(Code).*conj(fft(Code)));
R_shifted=fftshift(R);
R_final=abs(R_shifted/codeLen);
symmInterval=round(codeLen/2);
tau=-symmInterval:1:symmInterval-1;
figure,plot(tau,R_final), grid on;
ylim([0.88 1.02]);
axx=xlabel('Code Chips');
set(axx,'Interpreter','Latex');
axy=ylabel('Normalized Auto Correlation');
set(axy,'Interpreter','Latex');
tit=title('Weighted-voting Balanced Tausworthe $\nu$=4 $\mid$ Circular Auto Correlation');
set(tit,'Interpreter','Latex');
%% Upsampled code for plotting purposes
nSamples=4;
signalUp=rectpulse(Code,nSamples);
% Autocorrelation of upsampled signal
ACF=ifft(fft(signalUp).*conj(fft(signalUp)));
ACF=fftshift(ACF);
ACF=abs(ACF/length(ACF));
symmInterval=round(length(ACF)/2);
tau=-symmInterval:1:symmInterval-1;
figure,plot(tau,ACF),grid on;
ylim([0.88 1.02]);
axx=xlabel('Chips');
set(axx,'Interpreter','Latex');
axy=ylabel('Normalized Auto Correlation');
set(axy,'Interpreter','Latex');
tit=title(['Upsampled code of factor ' ,num2str(nSamples), ' $\mid$ Circular Auto Correlation']);
set(tit,'Interpreter','Latex');
% Discrete plot to appreciate autocorrelation behaviour
tau=-2:1:2;
figure,stem(tau,ACF(symmInterval-1:1:symmInterval+3)),axis('padded'),grid on;
tit=title(['Upsampled code of factor ' ,num2str(nSamples), ' $\mid$ Zoom around unambiguous main peak']);
set(tit,'Interpreter','Latex');
%% Simulation in IF (close-to-baseband) scenario
% Generation of correct upsampled signal
signalUp=rectpulse(Code,SpS);
% Generation of modulated signal + channel effect (Doppler shift only)
signalMod=modulate(signalUp,freqCarrierShifted,freqSampling);
% Search Space (SS) definition
L=codeLen;
Ts=1/freqSampling;
Tcoh=L*Ts;
Ntau=L;
binsTau=0:1:Ntau-1;
freqDopplerMax=10e+03;
deltaFreq=ceil(2/(3*Tcoh));
binsDoppler=-freqDopplerMax:deltaFreq:freqDopplerMax;
Nf=length(binsDoppler);
% Channel model (Doppler shift only)

return
% Autocorrelation of modulated signal
ACF=ifft(fft(signalMod).*conj(fft(signalMod)));
ACF=fftshift(ACF);
ACF=abs(ACF/max(ACF));
symmInterval=round(length(ACF)/2);
tau=-symmInterval:1:symmInterval-1;
figure,plot(tau,ACF),grid on;
ylim([0.88 1.02]);
axx=xlabel('Chips');
set(axx,'Interpreter','Latex');
axy=ylabel('Normalized Auto Correlation');
set(axy,'Interpreter','Latex');
tit=title('Rx signal $\mid$ Circular Auto Correlation');
set(tit,'Interpreter','Latex');
%% Acquisition stage with CAF computation
%signalLocal=signalUp.*exp(2i*pi*freqTest);
signalLocal=modulate(signalUp,freqTest,freqSampling);
% Cross-correlation between Rx signal and local signal replica (TO BE CHECKED WITH 'L' POINTS)
CCF=ifft(fft(signalMod,codeLen).*conj(fft(signalLocal,codeLen)));
CCF=fftshift(CCF);
CCF=abs(CCF/max(CCF));
symmInterval=round(length(CCF)/2);
tau=-symmInterval:1:symmInterval-1;
figure,plot(tau,CCF),grid on;
ylim([0.88 1.02]);
tit=title('Circular Cross Correlation between Rx signal and local replica');
set(tit,'Interpreter','Latex');