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
%% Modulate the PN sequence
SpS=2; % related to f_chip = 2MHz?
rectSignal=rectpulse(Code,SpS);
figure,stem(rectSignal(1:discreteTimePlot)),grid on,axis('padded');
title('Rectangularly-shaped signal | SpS=2;','Interpreter','latex');
% Sinusoidal shaping and bandpass modulation
freqCarrier=2.1e+09;
freqSampl=5e+09;
signalRF=modulate(rectSignal,freqCarrier,freqSampl,'am');
figure,stem(signalRF(1:discreteTimePlot*10)),grid on,axis('padded');
title('Modulated signal | $f_c=2.1 GHz$','Interpreter','latex');
%% Channel model
shiftDoppler=500e+03;