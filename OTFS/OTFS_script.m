% Script for performance evaluation of OTFS
clc
clear
close all
tic;
%% Simulation parameters
% Rows and columns of the DD space
M=64;
N=16;
% Compute normalized N-DFT matrix (for the columns)
Fn=dftmtx(N);
Fn=Fn/norm(Fn);
% Carrier frequency
fc=4e+09;
% Sub-carrier spacing
deltaF=15e+03;
% Block duration
T=1/deltaF;
% Speed of light
c=299792458;
% OTFS grid delay and Doppler resolutions
delayRes=1/(M*deltaF);
DopplerRes=1/(N*T);
numWrongRxBitsMax=15000;
data_grid=zeros(M,N);
M_data=M-M/N;
data_grid(1:M_data,1:N)=1;
LengthPlot=10;
modSize=4;
numSymbPerFrame=M*N;
numInfoBitsPerFrame=numSymbPerFrame*log2(modSize);
%% OTFS Modulator -> 4-QAM modulation + Zak transform
NFrame=5000;
numWrongRxBitsTot=zeros(1,11);
for SNR_dB=0:LengthPlot
    numTxInfoBits=0;
    for Iframe=1:NFrame
        %while numWrongRxBitsTot<numWrongRxBitsMax
        infoBits=randi([0 1],numInfoBitsPerFrame,1);
        modSignal=qammod(infoBits,modSize,"gray","InputType","bit");
        % Create X matrix
        X=reshape(modSignal,M,N);
        x=reshape(X.',N*M,1);
        % IFFT of X matrix
        X_tilda=X*Fn;
        % Vectorize by columns to obtain time-domain signal
        s=reshape(X_tilda,1,N*M);
        %% Channel model: Extended Typical Urban (ETU) -> high delay spread
        max_UE_speed=120*(1000/3600);
        % Compute maximum Doppler spread and its normalization
        nu_max=fc*(max_UE_speed/c);
        k_max=nu_max/DopplerRes;
        % ETU model
        %delays_ETU=[0, 50, 120, 200, 230, 500, 1600, 2300, 5000]*10^-9;
        delays_ETU=0;
        %pdp_ETU=[-1.0, -1.0, -1.0, 0.0, 0.0, 0.0, -3.0, -5.0, -7.0];
        pdp_ETU=-1;
        %% Generate channel parameters
        pdp_linear=10.^(pdp_ETU/10);
        pdp_linear=pdp_linear/sum(pdp_linear);
        % Number of taps (propagation path)
        taps=length(pdp_ETU);
        L_set=unique(taps);
        % Generate channel coefficients for Rayleigh fading
        g_i=sqrt(pdp_linear).*(sqrt(1/2)*(randn(1,taps)+1i*randn(1,taps)));
        %g_i=1;
        % Generate delay taps (assumed integer)
        l_i=round(delays_ETU./delayRes);
        % Generate Doppler taps
        k_i=(k_max*cos(2*pi*rand(1,taps)));
        % Generate G matrix to represent the channel
        z=exp(1i*2*pi/N/M);
        delaySpread=max(l_i);
        gs=zeros(delaySpread+1,N*M);
        for q=0:N*M-1
            for i=1:taps
                gs(l_i(i)+1,q+1)=gs(l_i(i)+1,q+1)+g_i(i)*z^(k_i(i)*(q-l_i(i)));
            end
        end
        G=zeros(N*M,N*M);
        for q=0:N*M-1
            for ell=0:delaySpread
                if (q>=ell)
                    G(q+1,q-ell+1)=gs(ell+1,q+1);
                end
            end
        end
        %% Tx into channel modelled with DD and AWGN
        r=G*s.';
        Es=mean(abs(qammod(0:modSize-1,modSize).^2));
        % Compute SNR
        SNR=10.^(SNR_dB/10);
        sigma_w2=Es/SNR;
        noise=sqrt(sigma_w2/2)*(randn(N*M,1)+1i*randn(N*M,1));
        r=r+noise;
        %% OTFS demodulator
        Y_tilda=reshape(r,M,N);
        Y=Y_tilda*Fn;
        %% LMMSE signal detection in time domain
        % Estimated time domain samples
        % LMMSE
        s_hat=inv(G'*G+sigma_w2*eye(M*N))*(G'*r);
        % MRC
        %s_hat=G'*r;
        % MxN estimated DD symbols
        X_hat_tilda=reshape(s_hat.',M,N);
        X_hat=X_hat_tilda*Fn';
        x_hat=X_hat(:); %reshape(X_hat.',N*M,1);
        infoBits_hat=qamdemod(x_hat,modSize,'gray','OutputType','bit');
        %[est_bits,x_data] = Block_LMMSE_detector(N,M,modSize,sigma_w2,data_grid,r,gs,L_set);
        %% Error rate computation
        numWrongRxBits=sum(xor(infoBits_hat,infoBits));
        numWrongRxBitsTot(SNR_dB+1)=numWrongRxBits+numWrongRxBitsTot(SNR_dB+1);
        numTxInfoBits=numTxInfoBits+length(infoBits);
    end
    %% Error rate computation
end
BER=numWrongRxBitsTot/numTxInfoBits;
SNRplot=0:LengthPlot;
semilogy(SNRplot,BER,'-x','LineWidth',2,'MarkerSize',5)
grid on
axx=xlabel('$SNR$');
set(axx,'Interpreter','Latex');
axy=ylabel('Bit Error Rate');
set(axy,'Interpreter','Latex');
tit=title('OTFS');
set(tit,'Interpreter','Latex');
leg=legend('BER');
set(leg,'Interpreter','Latex');
toc;