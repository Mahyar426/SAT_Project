% Simulator for TC: LDPC channel coding
clc
clear
close all
load 128_64_LDPCcode.mat
%% Simulation parameters
tic;
k=64;
n=128;
Eb_No=0:1:4;
Eb_No_linear=10.^(Eb_No./10);
sigma=sqrt(1./(2*(k/n).*Eb_No_linear));
numIterMax=100;
numMaxWrongRxCodewords=100;
alpha=0.8;
%% Values for Tanner graph
A=full(H);
numVariableNodes=zeros(1,size(H,2));
numCheckNodes=zeros(1,size(H,1));
% Counting number of connections between check nodes and variable nodes
for i=1:size(H,1)
    for j=1:size(H,2)
        if A(i,j)==1
            numVariableNodes(j)=numVariableNodes(j)+1;
            numCheckNodes(i)=numCheckNodes(i)+1;
        end
    end
end
% Finding the maximum value of connections possible
MaxVariableNodes=max(max(numVariableNodes));
MaxCheckNodes=max(max(numCheckNodes));
% Finding the indexes of the connections for both
% check nodes and variable nodes
IndexVariableNodes=zeros(size(H,2),MaxVariableNodes);
IndexCheckNodes=zeros(size(H,1),MaxCheckNodes);
for i=1:size(H,1)
    cnt=1;
    for j=1:size(H,2)
        if A(i,j)==1
            IndexCheckNodes(i,cnt)=j;
            cnt=cnt+1;
        end
    end
end
for j=1:size(H,2)
    cnt2=1;
    for i=1:size(H,1)
        if A(i,j)==1
            IndexVariableNodes(j,cnt2)=i;
            cnt2=cnt2+1;
        end
    end
end
%% Monte-Carlo simulation
for energy=1:length(Eb_No)
numTxCodewords=0;
numTxInfoBits=0;
numWrongRxCodewords=0;
numWrongRxInfoBits=0;
% For every energy point:
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
    noiseI=randn(1,size(H,1));                    % in phase noise
    noiseQ=randn(1,size(H,1));                    % quadrature noise
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
    numTxInfoBits=numTxInfoBits+size(H,1);
    %% NMS iterative decoding block
    y=receivedCodeword>0;
    syndrone=mod(y*H',2);
    % NMS starting condition
    if sum(syndrone)~=0
        numIter=0;
        LLR=(2*receivedCodeword)./sigma(energy)^2;
        % Tanner graph construction-> creating array structures for both variable
        % and check nodes: a field for the numerical values of their respective
        %  update rule and another field for their connection indexes
        for i=1:size(H,2)
            variableNodes(i).numValue=ones(1,numVariableNodes(i))*LLR(i);
            variableNodes(i).connToCheckNodes=nonzeros(IndexVariableNodes(i,:))';
        end
        for i=1:size(H,1)
            checkNodes(i).numValue=zeros(1,numCheckNodes(i));
            checkNodes(i).connToVariableNodes=nonzeros(IndexCheckNodes(i,:))';
        end
        % NMS main loop
        omega=zeros(1,size(H,2));
        while numIter<numIterMax
            % Check Node Update Rule
            for check = 1 : k
                for h = 1 : length(checkNodes(check).connToVariableNodes)
                    SignProd=1;
                    MinA=Inf;
                    for h_excluded = 1 : length(checkNodes(check).connToVariableNodes)
                        if h_excluded~=h
                            for a = 1 : length(variableNodes(checkNodes(check).connToVariableNodes(h_excluded)).connToCheckNodes)
                                if (variableNodes(checkNodes(check).connToVariableNodes(h_excluded)).connToCheckNodes(a))==check
                                    var=variableNodes(checkNodes(check).connToVariableNodes(h_excluded)).numValue(a);
                                    MinA=min(MinA,abs(var));
                                    SignProd=sign(var)*SignProd;
                                end
                            end
                        end
                    end
                    checkNodes(check).numValue(h)=alpha*MinA*SignProd;
                end
            end
            % A-Posteriori Update Rule
            for variable = 1 : n % Going through each variable node
                SumB=0;
                for h = 1 : length(variableNodes(variable).connToCheckNodes) % Accessing check nodes connected to that variable node
                    for b = 1 : length(checkNodes(variableNodes(variable).connToCheckNodes(h)).connToVariableNodes) % Finding the corresponding value from that check node connected to this specific variable node
                        if checkNodes(variableNodes(variable).connToCheckNodes(h)).connToVariableNodes(b)==variable % Checking that if we are picking the corresponding check node value for that specific variable node
                            SumB=SumB+checkNodes(variableNodes(variable).connToCheckNodes(h)).numValue(b);
                        end
                    end
                end
                SumB_Omega=SumB;
                % Variable Node Update Rule
                for h = 1 : length(variableNodes(variable).connToCheckNodes) % Accessing check nodes connected to that variable node
                    for h_excluded = 1 : length(variableNodes(variable).connToCheckNodes) % Condition for excluding that check node
                        if h_excluded==h % To choose the check node that we are extracting
                            for b = 1 : length(checkNodes(variableNodes(variable).connToCheckNodes(h_excluded)).connToVariableNodes) % Finding the corresponding value from that check node connected to this specific variable node
                                if (checkNodes(variableNodes(variable).connToCheckNodes(h_excluded)).connToVariableNodes(b))==variable % To check if we pick the right value for the connection between the current check node and the specific variable node
                                    SumB_excluded=SumB_Omega-checkNodes(variableNodes(variable).connToCheckNodes(h_excluded)).numValue(b);
                                end
                            end
                        end
                    end
                    variableNodes(variable).numValue(h)=alpha*(LLR(variable)+SumB_excluded);
                end
                omega(variable)=LLR(variable)+SumB_Omega;
            end
            % Update y and check the syndrone again
            y=omega>=0;
            syndrone=mod(y*H',2);
            numIter=numIter+1;
            if sum(syndrone)==0
                break;
            end
        end
    end
    %% Error rate computation block
    if ~isequal(codedBits,y)
        numWrongRxCodewords=numWrongRxCodewords+1;
        numWrongRxInfoBits=numWrongRxInfoBits+sum(xor(infoBits,y(1:64)));
    end
end
CER(energy)=numWrongRxCodewords/numTxCodewords;
BER(energy)=numWrongRxInfoBits/numTxInfoBits;
end
%% Plotting CER and BER performance
figure
semilogy(Eb_No,CER,'-ob','LineWidth',3),axis('tight'),grid on;
ylim([10^(-5) 10^0])
axx=xlabel('$E_b/N_o$');
set(axx,'Interpreter','Latex');
axy=ylabel('Codeword Error Rate');
set(axy,'Interpreter','Latex');
tit=title('LDPC code (128,64) - NMS iterative decoding');
set(tit,'Interpreter','Latex');
leg=legend('CER');
set(leg,'Interpreter','Latex');
figure
semilogy(Eb_No,BER,'-sr','LineWidth',3),axis('tight'),grid on;
ylim([10^(-6) 10^0]);
axx=xlabel('$E_b/N_o$');
set(axx,'Interpreter','Latex');
axy=ylabel('Bit Error Rate');
set(axy,'Interpreter','Latex');
tit=title('LDPC code (128,64) - NMS iterative decoding');
set(tit,'Interpreter','Latex');
leg=legend('BER');
set(leg,'Interpreter','Latex');
toc;