% Modulation block

clc
clear
close all
load 128_64_LDPCcode.mat

k=64;
n=128;
Eb_No=0:1:4;
Eb_No_linear=10.^(Eb_No./10);
sigma=sqrt(1./(2*(k/n).*Eb_No_linear));
spectEff=2;
SNR=Eb_No-spectEff;

infoBits=randi([0 1],k,1)';
codedBits=mod(infoBits*G,2)';

SymbolsI = 2*codedBits(1:2:end)-1;              % in phase symbols
SymbolsQ = 2*codedBits(2:2:end)-1;              % quadrature symbols
symbolsTx = SymbolsI+1i.*SymbolsQ;

% Normalized Complex WGN generation (Noise Power set to 1)
Noise=(randn(k,1)+1i*randn(k,1))*sigma(end);

symbolsRx=symbolsTx+Noise;
receivedCodewordsReal=real(symbolsRx);
receivedCodewordsImag=imag(symbolsRx);
receivedCodewords=zeros(1,size(H,2));
receivedCodewords(1:2:end)=receivedCodewordsReal;
receivedCodewords(2:2:end)=receivedCodewordsImag;

figure
plot(symbolsTx,'r.');
axis equal
xlim([-2 2])
ylim([-2 2])
grid on;
figure
plot(symbolsRx,'b.');
axis equal
xlim([-2 2])
ylim([-2 2])
grid on;

