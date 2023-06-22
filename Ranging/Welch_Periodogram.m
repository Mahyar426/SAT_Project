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
%% Analog signal generation ===============================================
SpS=4;                                  % Number of samples per symbols
analog_signal = rectpulse(Code,SpS);    % Rectangularly-shaped signal
%% Welch periodogram 
N=length(analog_signal);
Nwel=N/10;
h=ones(1,Nwel);     % Use of a rectangular window
Noverlap=Nwel/2;    % Samples of overlap
Nfft=Nwel;
Fs=10e6;            % Adjust this to match frequency axis
[Px,f]=pwelch(analog_signal,h,Noverlap,Nfft,Fs,'centered');
figure,semilogy(f,Px),axis('tight');
grid on
xlabel('Frequency','Interpreter','Latex'),title('$$\hat{P} (f)$$ of T4B sequence','Interpreter','Latex');