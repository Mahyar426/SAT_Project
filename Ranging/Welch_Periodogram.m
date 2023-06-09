% Validation of spectral properties of T4B

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
%% Analog signal generation
l=4;
k=6;
f_X_Band=(7.195)*1e9;
f_Ka_Band=(22.85)*1e9;
chip_rate= 2.5e6; %taken from 414x0g2 pag 21
Chip_rate_X_Band=(l*221*f_X_Band)/(749*128*2^k);
Chip_rate_Ka_Band=(l*221*f_Ka_Band)*(3599*128*2^k);
Tc = 1/chip_rate; %0.4e6 --> 0.4 micro sec
Tc_X=1/Chip_rate_X_Band;
Tc_Ka=1/Chip_rate_Ka_Band;
nsamp=4; % number of samples per symbols
nsamp_X=10;
nsamp_Ka=10;
analog_signal = rectpulse(Code,nsamp); % this is the analog signal
%% Welch periodogram 
N=length(analog_signal);
Nwel=N/10;
h=ones(1,Nwel); % this is the window
Noverlap=Nwel/2; %uses NOVERLAP samples of overlap from section to section.  NOVERLAP must be an integer smaller than WINDOW if WINDOW is an integer, or smaller than the length of WINDOW if WINDOW is a vector.
Nfft=Nwel;
Fs=10e6; %Fs is the sampling frequency specified in hertz.  If Fs is empty, it defaults to 1 Hz
[Px,f]=pwelch(analog_signal,h,Noverlap,Nfft,Fs,'centered');
figure,semilogy(f,Px),axis('tight');
grid on
xlabel('f'),ylabel('$$\hat{P} (f)$$','Interpreter','Latex');