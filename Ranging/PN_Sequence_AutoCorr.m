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
figure,plot(tau,R_final),axis('padded');
axx=xlabel('Code Chips');
set(axx,'Interpreter','Latex');
axy=ylabel('Normalized Auto Correlation');
set(axy,'Interpreter','Latex');
tit=title('Weighted-voting Balanced Tausworthe $\nu$=4 $\mid$ Circular Auto Correlation');
set(tit,'Interpreter','Latex');
%% Upsampled sequence at Intermediate Frequency
nsamp=15;
signalUp=rectpulse(Code,nsamp); 
fc=10e6;
fs=(30)*1e6;
signalMod=modulate(signalUp,fc,fs);
% Autocorrelation of upsampled mdoulated signal
ACF=ifft(fft(signalMod).*conj(fft(signalMod)));
ACF=fftshift(ACF);
ACF=abs(ACF/length(ACF));
symmInterval=round(length(ACF)/2);
tau=-symmInterval:1:symmInterval-1;
figure,plot(tau,ACF.^2),axis('padded');
%% Doppler shift and cross-correlation
freqDoppler=10e+03;
signalShift=signalMod.*exp(2i*pi*freqDoppler);
% Do cross-correlation
CCF=ifft(fft(signalMod).*conj(fft(signalShift)));
CCF=fftshift(CCF);
CCF=abs(CCF/length(CCF));
symmInterval=round(length(CCF)/2);
tau=-symmInterval:1:symmInterval-1;
figure,plot(tau,CCF.^2),axis('padded');
%% Make spectrum plots
L = 64;
wvtool(rectwin(L))