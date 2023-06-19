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
SpS=4; % related to f_chip = 2MHz?
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
Fs=8e+06;           % Chosen to properly show spectral properties
[Px,f]=pwelch(signalRect,h,Noverlap,Nfft,Fs,'centered');
figure,plot(f,pow2db(Px)),axis('tight');
grid on
xlabel('Frequency [GHz]','Interpreter','latex');
ylabel('Magnitude [dB]','Interpreter','latex');
title('$\hat{P} (f)$ of rectangularly-shaped signal | Baseband','Interpreter','latex')
% Sinusoidal shaping and bandpass modulation
freqCarrier=2.1e+09;
freqSampl=5e+09;
% We use: modulate() = x.*cos(2*pi*freqCarrier*t)
% which applies amplitude modulation at the desired frequency
signalRF=modulate(signalRect,freqCarrier,freqSampl,'am');
figure,stem(signalRF(1:discreteTimePlot*10)),grid on,axis('padded');
xlabel('[n]','Interpreter','latex');
ylabel('y[n]','Interpreter','latex');
title('Modulated signal | $f_c = 2.1 GHz$','Interpreter','latex');
% PSD estimation of bandpass signal
N=length(signalRF);
Nwel=N/10;          % Length of thw window
h=ones(1,Nwel);     % Rectangular window to pre-filter
Noverlap=Nwel/2;    % Number of overlapping samples
Nfft=4096;          % Number of FFT points per window
Fs=8.8e+09;         % Chosen to properly show spectral properties
[Px,f]=pwelch(signalRF,h,Noverlap,Nfft,Fs,'centered');
figure,plot(f,pow2db(Px)),axis('tight');
grid on
xlabel('Frequency [GHz]','Interpreter','latex');
ylabel('Magnitude [dB]','Interpreter','latex');
title('$\hat{P} (f)$ of modulated signal | Bandpass with $f_{carrier} = 2.1$ GHz','Interpreter','latex')
%% Channel model: phase shift and frequency shift
shiftDoppler=501e+03;
shiftPhase=501;
signalRx=circshift(signalRF.*exp(2i*pi*shiftDoppler),shiftPhase);
%% Receiver side: demodulate + 2D correlation
% We demodulate at an intermediate frequency using
% demod() = sinusoid multiplication at freqIF and 
% fifth-order Butterworth lowpass filter 
freqIF=1.5e+09;
signalIF=demod(signalRx,freqCarrier,freqSampl,'am');
% Test arrays
shiftDopplerTestArray=1e+03:50e+03:951e+03;
shiftPhaseTestArray=1:50:1000;
% Cross-ambiguity function computation
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
        CAF=ifft(fft(replicaTest).*conj(fft(signalIF)));
        squaredCAF(i,j,:)=abs(CAF).^2;
    end
end
%% 3D plot of the cross-ambiguity function
% [X,Y]=meshgrid(shiftDopplerTestArray,shiftPhaseTestArray);
% figure,mesh(X,Y,Z')




