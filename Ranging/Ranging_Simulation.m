% Simulation of GNSS receiver acquisition stage 

clc
clear
close all

%% Initializing circular shift registers ==================================
C1 = [+1 -1];
C2 = [+1 +1 +1 -1 -1 +1 -1];
C3 = [+1 +1 +1 -1 -1 -1 +1 -1 +1 +1 -1];
C4 = [+1 +1 +1 +1 -1 -1 -1 +1 -1 -1 +1 +1 -1 +1 -1];
C5 = [+1 +1 +1 +1 -1 +1 -1 +1 -1 -1 -1 -1 +1 +1 -1 +1 +1 -1 -1];
C6 = [+1 +1 +1 +1 +1 -1 +1 -1 +1 +1 -1 -1 +1 +1 -1 -1 +1 -1 +1 -1 -1 -1 -1];
counter=1;
regLen=[2 7 11 15 19 23];
codeLen=prod(regLen);
%% PN sequence generation =================================================
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
%% Baseband signal generation =============================================
SpS=4; 
Fchip=2e+06;       % From CCSDS 414.1-B-2 - Table A-1
signalRect=rectpulse(Code,SpS);
figure,stem(signalRect(1:discreteTimePlot)),grid on,axis('padded');
title(['Rectangularly-shaped signal | SpS=',num2str(SpS)],'Interpreter','latex');
xlabel('[n]','Interpreter','latex');
ylabel('y[n]','Interpreter','latex');
% PSD estimation of baseband signal
N=length(signalRect);
Nwel=N/10;          % Length of thw window
h=ones(1,Nwel);     % Rectangular window to pre-filter
Noverlap=Nwel/2;    % Number of overlapping samples
Nfft=4096;          % Number of FFT points per window
Fs=SpS*Fchip;       % Optimal sampling rate at baseband
[Px,f]=pwelch(signalRect,h,Noverlap,Nfft,Fs,'centered');
figure,plot(f,pow2db(Px)),axis('tight'),grid on;
xlabel('Frequency [Hz]','Interpreter','latex');
ylabel('Magnitude [dB]','Interpreter','latex');
title('$\hat{P} (f)$ of rectangularly-shaped signal | Baseband','Interpreter','latex');
%% RF-IF front-end operations =============================================
% Here we just modulate the local replica of the code to an intermediate
% frequency close to baseband, to preserve the capability of estimating
% the Doppler shift (either positive or negative)
freqIF=10e+06;
freqSamplIF=SpS*freqIF;
signalIF=modulate(signalRect,freqIF,freqSamplIF,'am');
% PSD estimation of IF signal
N=length(signalRect);
Nwel=N/10;          % Length of thw window
h=ones(1,Nwel);     % Rectangular window to pre-filter
Noverlap=Nwel/2;    % Number of overlapping samples
Nfft=4096;          % Number of FFT points per window
Fs=freqSamplIF;     % Optimal sampling rate
[Px,f]=pwelch(signalIF,h,Noverlap,Nfft,Fs,'centered');
figure,plot(f,pow2db(Px)),axis('tight'),grid on;
xlabel('Frequency [Hz]','Interpreter','latex');
ylabel('Magnitude [dB]','Interpreter','latex');
title('$\hat{P} (f)$ at Intermediate Frequency (IF) | $f_{carrier}$ $ =$ $ 10$ $ MHz$','Interpreter','latex');
return
%% Search Space (SS) definition ===========================================
% L=length(Code);
% coherentTime=
%% Channel model: phase shift and frequency shift =========================
% shiftDoppler=;
% shiftPhase=;
% signalRx=circshift(signalRF.*exp(2i*pi*shiftDoppler),shiftPhase);
%% Acquisition stage: delay and Doppler estimation ========================
% Test arrays -> THEY HAVE TO BE COMPUTED PROPERLY!
shiftDopplerTestArray=1e+03:50e+03:501e+03;
shiftPhaseTestArray=1:50:501;
% Generate local replica signal modulated at IF

% Cross-ambiguity function computation in the DD domain 
for i=1:length(shiftPhaseTestArray)
    % Using a local replica of the code
    replicaTest=signalRect;
    % Apply test phase shifts
    replicaTest=circshift(replicaTest,shiftPhaseTestArray(i));
    for j=1:length(shiftDopplerTestArray)
        % Modulate at a test frequency composed by the intermediate
        % frequency and the test Doppler shifts
        replicaTest=modulate(replicaTest,freqIF,freqSampl,"am");
        replicaTest=replicaTest.*exp(2i*pi*shiftDopplerTestArray(j));
        % Compute cross-ambiguity function
        CAF=ifft(fft(replicaTest).*conj(fft(signalIF)));% it's already doing circular shift!!!
        squaredCAF(i,j,:)=abs(CAF).^2;
    end
end

%% 3D plot of the cross-ambiguity function ================================
[X,Y]=meshgrid(shiftPhaseTestArray,shiftDopplerTestArray);
figure,CA_Plot=mesh(X,Y,Z/1e+12);
CA_Plot.FaceColor = 'flat';




