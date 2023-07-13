% Code for regenerative Pseudo-Noise sequence (T4B)

clc
clear
close all

%% Parameters setting for close-to-baseband simulation
normFactor=1e+00;
freqCarrier=(10e+06)/normFactor;
freqSampling=(30e+06)/normFactor;
freqChip=(2e+06)/normFactor;
SpS=round(freqSampling*(1/freqChip));
shiftDoppler=777.932974729314;                               % Value picked in binsDoppler (see below)
freqCarrierShifted=freqCarrier+shiftDoppler;    % Channel effect on carrier frequency
test=0;                                         % Test flag for showing plots
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
% Discrete plot to fully appreciate autocorrelation behaviour
tau=-3:1:3;
figure,stem(tau,R_final(symmInterval-2:1:symmInterval+4)),axis('padded'),grid on;
hold on;
plot(tau,R_final(symmInterval-2:1:symmInterval+4),'Color',[0.8500 0.3250 0.0980]);
% ylim([0.88 1.02]);
title('T4B generated code $\mid$ Zoom around CAC main peak','Interpreter','Latex');
legend('Discrete points','Continuous envelope','Interpreter','Latex');
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
% Discrete plot to fully appreciate autocorrelation behaviour
tau=-3:1:3;
figure,stem(tau,ACF(symmInterval-2:1:symmInterval+4)),axis('padded'),grid on;
hold on;
plot(tau,ACF(symmInterval-2:1:symmInterval+4),'Color',[0.8500 0.3250 0.0980]);
tit=title(['T4B code upsampled of a factor ' ,num2str(nSamples), ' $\mid$ Zoom around CAC main peak']);
set(tit,'Interpreter','Latex');
legend('Discrete points','Continuous envelope','Interpreter','Latex');
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
deltaFreq=2/(3*Tcoh);     % Brought to an integer value for ease
binsDoppler=-freqDopplerMax:deltaFreq:freqDopplerMax;
Nf=length(binsDoppler);
%% Test plot for ACF of modulated signal
if test==1
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
end
%% Acquisition stage with CAF computation
CAF=zeros(Nf,Ntau);
for i=1:Nf
    % Test frequency in Acquisition stage
    freqDopplerTest=binsDoppler(i); 
    freqTest=freqCarrier+freqDopplerTest;
    % Generate the local test replica
    % signalLocal=signalUp.*exp(2i*pi*freqTest);
    signalLocal=modulate(signalUp,freqTest,freqSampling);
    % Cross-correlation between Rx signal and local replica
    CCF=ifft(fft(signalMod,codeLen).*conj(fft(signalLocal,codeLen)));
    CCF=fftshift(CCF);
    CCF=abs(CCF);
    CAF(i,:)=(CCF.^2)';
end 
CAF=CAF./max(CAF);
figure,surf(CAF),grid on;
%% Test plot for single CCF
if test==1
    freqDopplerTest=0;
    freqTest=freqCarrier+freqDopplerTest;
    signalLocal=modulate(signalUp,freqTest,freqSampling);
    % Cross-correlation between Rx signal and local signal replica
    CCF=ifft(fft(signalMod,codeLen).*conj(fft(signalLocal,codeLen)));
    CCF=fftshift(CCF);
    CCF=abs(CCF/max(CCF));
    symmInterval=round(length(CCF)/2);
    tau=-symmInterval:1:symmInterval-1;
    figure,plot(tau,CCF),grid on;
    ylim([0.88 1.02]);
    tit=title('Circular Cross Correlation between Rx signal and local replica');
    set(tit,'Interpreter','Latex');
end