% Simulation of GNSS receiver acquisition stage 

clc
clear
close all

%% Parameters setting simulation
normFactor=1e+00;
freqCarrier=(10e+06)/normFactor;
freqSampling=(30e+06)/normFactor;
freqChip=(2e+06)/normFactor;
SpS=round(freqSampling*(1/freqChip));
shiftDoppler=777.932974729312;                  % Value picked in binsDoppler (see below)
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
%% Plotting normalized ACF
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
% ylim([0.88 1.02]);
title('T4B generated code $\mid$ Zoom around CAC main peak','Interpreter','Latex');
legend('Discrete points','Interpreter','Latex');
%% Upsampled code for plots
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
tau=-4:1:4;
figure,stem(tau,ACF(symmInterval-3:1:symmInterval+5)),axis('padded'),grid on;
hold on;
plot(tau,ACF(symmInterval-3:1:symmInterval+5),'Color',[0.8500 0.3250 0.0980]);
tit=title(['T4B code upsampled of a factor ' ,num2str(nSamples), ' $\mid$ Zoom around CAC main peak']);
set(tit,'Interpreter','Latex');
legend('Discrete points','Continuous envelope','Interpreter','Latex');
%% IF Simulation (close-to-baseband)
% Generation of correct upsampled signal
signalUp=rectpulse(Code,SpS);
% Generation of modulated signal + channel effect (Doppler shift only)
signalMod=modulate(signalUp,freqCarrierShifted,freqSampling);
% Search Space (SS) definition
L=codeLen*SpS;
Ts=1/freqSampling;
Tcoh=L*Ts;
Ntau=L;
binsTau=0:1:Ntau-1;
freqDopplerMax=10e+03;          
deltaFreq=2/(3*Tcoh);     
binsDoppler=-freqDopplerMax:deltaFreq:freqDopplerMax;
Nf=length(binsDoppler);
%% Test: signalMod ACF
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
%% Acquisition stage main loop
% Peak is expected to be in position 8161 in binsDoppler
peakPos=8161;
% Compute CAF in the portion of the Search Space where
% we expect to find the match in Doppler shift
CAF_Peak=zeros(41,length(signalMod));
for i=peakPos-20:peakPos+20
% for i=1:100
    % Test frequency in Acquisition stage
    freqDopplerTest=binsDoppler(i); 
    freqTest=freqCarrier+freqDopplerTest;
    % Generate the local test replica
    signalLocal=modulate(signalUp,freqTest,freqSampling);
    % Cross-correlation between Rx signal and local replica
    CCF=ifft(fft(signalMod).*conj(fft(signalLocal)));
    CCF=fftshift(CCF);
    CCF=abs(CCF);
    CAF_Peak(i-(peakPos-21),:)=(CCF.^2)';
%     CAF_Peak(i,:)=(CCF.^2)';
end 
% Normalize CAF and bring it in [0,1] range
CAF_Peak=CAF_Peak./length(CAF_Peak);
CAF_Peak=CAF_Peak./4e6;
%% CAF plots - Different views
% Plot the Peak
symmInterval=round(length(CAF_Peak)/2);
figure,surf(CAF_Peak(:,symmInterval-10:symmInterval+12)),grid on;
colormap turbo; 
colorbar,view(134,21);
xlabel('Delay bins','Interpreter','latex');
ylabel('Doppler bins','Interpreter','latex');
title('Rx signal and local replica $\mid$ Cross-Ambiguity Function $\mid$ Peak','Interpreter','latex');
figure,surf(CAF_Peak(:,symmInterval-10:symmInterval+12)),grid on;
colormap turbo,view(2);
title('View from above $\mid$ Cross-Ambiguity Function','Interpreter','latex');
xlabel('Delay bins','Interpreter','latex');
ylabel('Doppler bins','Interpreter','latex');
%% Test: single CCF
if test==1
    freqDopplerTest=777.932974729312;
    freqTest=freqCarrier+freqDopplerTest;
    signalLocal=modulate(signalUp,freqTest,freqSampling);
    % Cross-correlation between Rx signal and local signal replica
    CCF=ifft(fft(signalMod).*conj(fft(signalLocal)));
    CCF=fftshift(CCF);
    CCF=abs(CCF/max(CCF));
    symmInterval=round(length(CCF)/2);
    tau=-symmInterval:1:symmInterval-1;
    figure,plot(tau,CCF),grid on;
    ylim([0.88 1.02]);
    tit=title('Circular Cross Correlation between Rx signal and local replica');
    set(tit,'Interpreter','Latex');
end