% Code for regenerative Pseudo-Noise sequence (T4B)

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
%% Close-to-baseband simulation - Frequencies normalized in MHz 
freqCarrier=10;
shiftDoppler=0; % it will become the channel effect modelled in the other script;
freqCarrierShifted=freqCarrier+shiftDoppler;
freqSampling=30;
freqChip=2;
SpS=freqSampling*(1/freqChip);
% Generation of correct upsampled signal
signalUp=rectpulse(Code,SpS);
% Generation of modulated signal + channel effect (Doppler shift only)
signalMod=modulate(signalUp,freqCarrierShifted,freqSampling);
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
return
%% Doppler shift and cross-correlation
freqDopplerTest=10e-3; % to be put inside loop
signalLocal=signalUp.*exp(2i*pi*(freqCarrier+freqDopplerTest));
% Do cross-correlation
CCF=ifft(fft(signalMod).*conj(fft(signalLocal)));
CCF=fftshift(CCF);
CCF=abs(CCF/max(CCF));
symmInterval=round(length(CCF)/2);
tau=-symmInterval:1:symmInterval-1;
figure,plot(tau,CCF.^2);
%ylim([0.88 1.02])
tit=title(' Circular Cross Correlation between Rx and local replica');
set(tit,'Interpreter','Latex');