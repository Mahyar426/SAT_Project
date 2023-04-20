% Simulator for TM: LDPC channel coding

clc
clear
close all
load 128_64_LDPCcode.mat

%% Simulation parameters
k=64;
n=128;
Eb_No=0:1:4;
Eb_No_linear=10.^(Eb_No./10);
sigma=sqrt(1./(2*(k/n).*Eb_No_linear));
zeroVector=zeros(1,size(H,1));                  % syndrone check
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
    while numWrongRxCodewords<numMaxWrongRxCodewords
        %% Information bits generation
        infoVector=randi([0 1],k,1)';
        %% Information bits encoding
        codedVector=mod(infoVector*G,2);
        %% QPSK Modulation block
        SymbolsI = 2*codedVector(1:2:end)-1;            % in phase symbols
        SymbolsQ = 2*codedVector(2:2:end)-1;            % quadrature symbols
        Symbols = SymbolsI+1i.*SymbolsQ;                % QPSK symbols
        %% AWGN Channel block
        noiseI=sigma(energy)*randn(1,size(H,1));        % in phase noise
        noiseQ=sigma(energy)*randn(1,size(H,1));        % quadrature noise
        receivedSymbols=(SymbolsI+noiseI)+1i*(SymbolsQ+noiseQ);
        receivedCodewordsReal=real(receivedSymbols);
        receivedCodewordsImag=imag(receivedSymbols);
        receivedCodewords=zeros(1,size(H,2));
        receivedCodewords(1:2:end)=receivedCodewordsReal;
        receivedCodewords(2:2:end)=receivedCodewordsImag;
        numTxCodewords=numTxCodewords+1;                % increase cTX
        numTxInfoBits=numTxInfoBits+size(H,1);          % increase uTx
        %% NMS iterative decoding block
        y=receivedCodewords>0;
        syndrone=mod(y*H',2);
        % NMS starting condition
        if ~isequal(syndrone,zeroVector)
            numIter=0;
            LLR=2*receivedCodewords./sigma(energy)^2;
            % Tanner graph construction: creating array structures for both variable 
            % and check nodes: a field for the numerical values of their respective
            %  update rule and another field for their connection indexes
            for i=1:size(H,2)
                variableNodes(i).numValue=LLR(i);
                variableNodes(i).connToCheckNodes=nonzeros(IndexVariableNodes(i,:))';
            end
            for i=1:size(H,1)
                checkNodes(i).numValue=zeros(1,numCheckNodes(i));
                checkNodes(i).connToVariableNodes=nonzeros(IndexCheckNodes(i,:))';
            end
        end
        % NMS main loop
        omega=zeros(1,size(H,2));
        while numIter<=numIterMax && ~isequal(syndrone,zeroVector)
            % Check node update rule
            for i=1:size(H,1)
                for j=1:length(checkNodes(i).connToVariableNodes)
                    signProd=1;
                    minA=Inf;                    
                    for z=1:length(checkNodes(i).connToVariableNodes)
                        if z~=j
                            var=variableNodes(checkNodes(i).connToVariableNodes(z)).numValue;
                            minA=min(minA,abs(var));
                            signProd=sign(var)*signProd;
                        end
                    end
                    checkNodes(i).numValue(j)=alpha*signProd*minA;
                end
            end
            % Variable update rule
            for i=1:size(H,2)
                SumB=0;
                for j=1:length(variableNodes(i).connToCheckNodes)
                    if variableNodes(i).connToCheckNodes(j)~=i
                        for z=1:length(checkNodes(variableNodes(i).connToCheckNodes(j)).connToVariableNodes)
                            if checkNodes(variableNodes(i).connToCheckNodes(j)).connToVariableNodes(z)==i
                                SumB=SumB+checkNodes(variableNodes(i).connToCheckNodes(j)).numValue(z);
                            end
                        end
                    end
                end
                variableNodes(i).numValue=alpha*variableNodes(i).numValue+alpha*SumB;
            end
            % Compute a-posteriori probability
            for i=1:size(H,2)
                SumB=0;
                for j=1:length(variableNodes(i).connToCheckNodes)
                    for z=1:length(checkNodes(variableNodes(i).connToCheckNodes(j)).connToVariableNodes)
                        if checkNodes(variableNodes(i).connToCheckNodes(j)).connToVariableNodes(z)==i
                            SumB=SumB+checkNodes(variableNodes(i).connToCheckNodes(j)).numValue(z);
                        end
                    end
                end
                omega(i)=variableNodes(i).numValue+SumB;
            end
            % Update y and check the syndrone again
            y=omega>0;
            syndrone=mod(y*H',2); 
            numIter=numIter+1;
        end
        %% Error rate computation block
        if ~isequal(syndrone,zeroVector)
            numWrongRxCodewords=numWrongRxCodewords+1;
            numWrongRxInfoBits=numWrongRxInfoBits+sum(xor(infoVector,y(1:64)));
        end
    end
    CER(energy)=numWrongRxCodewords/numTxCodewords;
    BER(energy)=numWrongRxInfoBits/numTxInfoBits;
end
%% Plotting CER and BER performance
figure
semilogy(Eb_No,CER,'-ob','LineWidth',3),axis('tight'),grid on;
axx=xlabel('$E_b/N_o$');
set(axx,'Interpreter','Latex');
axy=ylabel('CER');
set(axy,'Interpreter','Latex');
tit=title('LDPC code (128,64) - NMS iterative decoding');
set(tit,'Interpreter','Latex');
figure
semilogy(Eb_No,BER,'-sr','LineWidth',3),axis('tight'),grid on;
axx=xlabel('$E_b/N_o$');
set(axx,'Interpreter','Latex');
axy=ylabel('BER');
set(axy,'Interpreter','Latex');
tit=title('LDPC code (128,64) - NMS iterative decoding');
set(tit,'Interpreter','Latex');



