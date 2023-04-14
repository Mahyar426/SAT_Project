% LDPC performance simulator for Telecommand

clc
clear
close all

load 128_64_LDPCcode.mat
Eb_No=0:1:4;
Eb_No_linear=10.^(Eb_No./10);
k=64;
n=128;
sigma=sqrt(1./(2*(k/n).*Eb_No_linear));

info_vector=randi([0 1],k,1)';
coded_vector=mod(info_vector*G,2);

SymbolsI = 2*coded_vector(1:2:end)-1;                    % in phase symbols
SymbolsQ = 2*coded_vector(2:2:end)-1;                    % quadrature symbol
Symbols = SymbolsI+1i.*SymbolsQ;
noiseI=sigma(1)*randn(1,64);
noiseQ=sigma(1)*randn(1,64);
receivedSymbols=SymbolsI+noiseI+1i*(SymbolsQ+noiseQ);
receivedCodewordsReal=real(receivedSymbols);
receivedCodewordsImag=imag(receivedSymbols);

receivedCodewords=zeros(1,128);
receivedCodewords(1:2:end)=receivedCodewordsReal;
receivedCodewords(2:2:end)=receivedCodewordsImag;

%% NMS iterative decoding
zeroVector=zeros(1,128);
nIter=0;
nIterMax=100;
nWrongCodewords=100;
alpha=0.8;
y=receivedCodewords>0;
syndrone=mod(y*H',2);
if ~isequal(syndrone,zeroVector)
    nIter=1;
    LLR=2*receivedCodewords./sigma(1)^2;
end

%% Tanner graph construction
A=full(H);
VariableNodes=zeros(1,size(H,2));
CheckNodes=zeros(1,size(H,1));
% Counting number of connections between check nodes and variable nodes
for i=1:size(H,1)
    for j=1:size(H,2)
         if A(i,j)==1
                VariableNodes(j)=VariableNodes(j)+1;
                CheckNodes(i)=CheckNodes(i)+1;
         end
    end
end
% Finding the maximum value of connections possible
MaxVariableNodes=max(max(VariableNodes));
MaxCheckNodes=max(max(CheckNodes));
% Finding the indexes for both check and variable nodes 
% connections in the H matrix
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
% Creating array structures for both variable and check nodes: a field for
% the numerical values of their respective update rule and another field
% for their connection indexes
for i=1:size(H,2)
    variableNodes(i).numValue=LLR(i);
    variableNodes(i).connToCheckNodes=nonzeros(IndexVariableNodes(i,:))';
end
for i=1:size(H,1)
    checkNodes(i).numValue=zeros(1,CheckNodes(i));
    checkNodes(i).connToVariableNodes=nonzeros(IndexCheckNodes(i,:))';
end

%% Start the NMS iterative decoding
while nIter<=nIterMax && ~isequal(syndrone,zeroVector)
    % Variable update rule
    for i=1:size(H,2)
        SumB=0;
        for kk=1:length(variableNodes(i).connToCheckNodes)
            if variableNodes(i).connToCheckNodes(kk)~=i
                for jj=1:length(checkNodes(variableNodes(i).connToCheckNodes(kk)).connToVariableNodes)
                    if checkNodes(variableNodes(i).connToCheckNodes(kk)).connToVariableNodes(jj)==i
                        SumB=SumB+checkNodes(variableNodes(i).connToCheckNodes(kk)).numValue(jj);
                    end
                end
            end
        end
        variableNodes(i).numValue=alpha*variableNodes(i).numValue+alpha*SumB;
    end
    % Check node update rule
    for i=1:size(H,1)
        for j=1:length(checkNodes(i).connToVariableNodes)
            signProd=1;
            minA=inf;
            for z=1:length(checkNodes(i).connToVariableNodes)
                if z~=j
                    var=variableNodes(checkNodes(i).connToVariableNodes).numValue(z);
                    minA=min(minA,abs(var));
                    signProd=sign(var)*signProd;
                end
            end
            checkNodes(i).numValue(j)=alpha*signProd*minA;
        end
    end
    % Compute a-posteriori probability and update y

    % Check syndrone again
   nIter=nIter+1;
end
