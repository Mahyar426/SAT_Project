% Simulator for TM: LDPC channel coding, 2-PSK mopdulation

clc
clear
close all
load 128_64_LDPCcode.mat

%% Simulation parameters
k=64;
n=128;
Eb_No=0:1:4;
Eb_No_linear=10.^(Eb_No./10);
numMaxWrongRxCodewords=100;
numMaxIterNMS=2;
normValueNMS=0.8;
%% Values for Tanner graph
% Variable nodes: inspecting H by rows
for row=1:size(H,1)
    Tanner_v2c{row}=find(H(row,:));
end
% Check nodes: inspecting H by columns
for col=1:size(H,2)
    Tanner_c2v{col}=(find(H(:,col)))';
end

%% Monte-Carlo simulation

%% Initialization per energy point
numTxCodewords=0;
numTxInfoBits=0;
numWrongRxCodewords=0;
numWrongRxInfoBits=0;
sigma=sqrt(1/2*(k/n)*Eb_No_linear(end));
%% Information frame generation
infoBits=randi([0 1],k,1)';
numTxCodewords=numTxCodewords+1;
numTxInfoBits=numTxInfoBits+k;
%% LDPC encoding
codeword=mod(infoBits*G,2);
%% 2-PSK modulation
symbolTx=2*codeword-1;
%% AWGN channel
noise=sigma*randn(1,n);
symbolRx=symbolTx+noise;
%% LDPC decoding
% Syndrone test
y=symbolRx>=0;
syndrone=mod(y*H',2);
if sum(syndrone)~=0
    % Algorithm initialization
    LLR=2*codeword./(sigma^2);
    nIterNMS=1;
    m=n-k;
    % Initialize Tanner graph messages
    aPosterioriProb=zeros(n,1);
    channelMessage=H.*LLR;
    while nIterNMS<=numMaxIterNMS && sum(syndrone)~=0
        % Check node update
        for check=1:m
            % Accessing variable node values connected to m-th check value
            v2cMessage=channelMessage(check,Tanner_v2c{check});
            for t=1:length(v2cMessage)
                SignMessage=sign(v2cMessage);
                MagnitudeMessage=abs(v2cMessage);
                % Done to exclude t=th value
                SignMessage(t)=1;
                MagnitudeMessage(t)=Inf;
                c2vMessage(t)=prod(nonzeros(full(SignMessage)))*min(nonzeros(full(MagnitudeMessage)))*normValueNMS;
            end
            % Updating the check node values in the matrix
            channelMessage(check,Tanner_v2c{check})=c2vMessage;
        end
        % Variable node update
        % A-posteriori computation
        for variable=1:n
            c2vMessage=channelMessage(Tanner_c2v{variable},variable);
            aPosterioriProb=LLR(variable)+sum(c2vMessage);
            %v2cMessage=
        end



        
        % Syndrone computation
        nIterNMS=nIterNMS+1;
    end
    
end











