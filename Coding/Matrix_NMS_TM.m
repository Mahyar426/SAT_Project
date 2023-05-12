clc
clear
close all
load AR4JA2048.mat
%% Simulation parameters
k=1024;
n=2048;
Eb_No=0:0.5:2.5;
Eb_No_linear=10.^(Eb_No./10);
sigma=sqrt(1./(2*(k/n).*Eb_No_linear));
numMaxWrongRxCodewords=100;
numMaxIterNMS=100;
normValueNMS=0.8;
%% Values for Tanner graph
for row=1:size(H,1)
    Tanner_v2c{row}=find(H(row,:));
end
for col=1:size(H,2)
    Tanner_c2v{col}=(find(H(:,col)))';
end
%% Monte-Carlo
energy=4;
cntConvNMS=0;
numTxCodewords=0;
numTxInfoBits=0;
numWrongRxCodewords=0;
numWrongRxInfoBits=0;
while numWrongRxCodewords<numMaxWrongRxCodewords
    %% Information bits generation
    infoBits=randi([0 1],k,1)';
    %% Information bits encoding
    codedBits=mod(infoBits*G,2);
    %% QPSK Modulation block
    symbolsI = 2*codedBits(1:2:end)-1;            % in phase symbols
    symbolsQ = 2*codedBits(2:2:end)-1;            % quadrature symbols
    symbolsTx = symbolsI+1i.*symbolsQ;            % QPSK symbols
    %% AWGN Channel block
    noiseI=randn(1,size(H,2)/2);                    % in phase noise
    noiseQ=randn(1,size(H,2)/2);                    % quadrature noise
    noise=(noiseI+1i*noiseQ)*sigma(energy);       % QPSK noise
    symbolsRx=symbolsTx+noise;
    %% Receiver block
    symbolsRxReal=real(symbolsRx);
    symbolsRxImag=imag(symbolsRx);
    receivedCodeword=zeros(1,size(H,2));
    receivedCodeword(1:2:end)=symbolsRxReal;
    receivedCodeword(2:2:end)=symbolsRxImag;
    %% Counters update
    numTxCodewords=numTxCodewords+1;
    numTxInfoBits=numTxInfoBits+k;
    %% NMS iterative decoding block
    receivedCodewordNMS=-receivedCodeword;
    receivedCodewordNMS(2049:end)=0;         % last M punctured symbols
    y=receivedCodewordNMS<=0;
    syndrone=mod(y*H',2);
    % Starting condition of the NMS iterative algorithm
    if sum(syndrone)~=0
        % Algorithm initialization
        LLR=receivedCodewordNMS;    % computed on received codeword
        numIterNMS=0;
        % Initialize Tanner graph messages
        aPosterioriProb=zeros(size(H,2),1);
        channelMessage=H.*LLR;
        [c, v]=size(channelMessage);
        while numIterNMS<numMaxIterNMS
            % Check node update rule
            for check=1:c
                varCheck=[];
                variableToCheckMessage=channelMessage(check,Tanner_v2c{check});
                for t=1:length(variableToCheckMessage)
                    Sign=(2*(variableToCheckMessage>=0)-1);
                    Magnitude=abs(variableToCheckMessage);
                    Sign(t)=1;
                    Magnitude(t)=Inf;
                    varCheck(t)=normValueNMS*prod(Sign)*min(Magnitude);
                end
                channelMessage(check,Tanner_v2c{check})=varCheck;
            end
            % Variable node update rule
            for variable=1:v
                varVariable=[];
                CheckToVariableMessage=channelMessage(Tanner_c2v{variable},variable);
                Omega(variable)=LLR(variable)+sum(CheckToVariableMessage);
                for t=1:length(CheckToVariableMessage)
                    varVariable(t)=Omega(variable)-CheckToVariableMessage(t);
                end
                channelMessage(Tanner_c2v{variable},variable)=varVariable;
            end
            y=Omega<=0;
            syndrone=mod(y*H',2);
            numIterNMS=numIterNMS+1;
            if sum(syndrone)==0
                cntConvNMS=cntConvNMS+1;
                break;
            end
        end
    end
    if ~isequal(codedBits(1:n),y(1:n))
        % Update counters for wrong Rx
        numWrongRxCodewords=numWrongRxCodewords+1;
        numWrongRxInfoBits=numWrongRxInfoBits+sum(xor(infoBits,y(1:k)));
    end
end
%% Error Rate Computation
CER=numWrongRxCodewords/numTxCodewords;
BER=numWrongRxInfoBits/numTxInfoBits;