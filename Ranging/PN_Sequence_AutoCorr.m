% Code for regenerative Pseudo-Noise sequence (T4B)
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
%% Plotting normalized circular autocorrelation
R=ifft(fft(Code).*conj(fft(Code)));
R_shifted=fftshift(R);
R_final=abs(R_shifted/codeLen);
symmInterval=round(codeLen/2);
tau=-symmInterval:1:symmInterval-1;
figure,plot(tau,R_final),axis('padded');
axx=xlabel('Code Chips');
set(axx,'Interpreter','Latex');
axy=ylabel('Normalized Auto Correlation');
set(axy,'Interpreter','Latex');
tit=title('Weighted-voting Balanced Tausworthe $\nu$=4 $\mid$ Circular Auto Correlation');
set(tit,'Interpreter','Latex');
%% Test
% C5_extended=[];
% while length(C5_extended)<length(Code)
%     C5_extended=[C5_extended;C5'];
% end
% C2_extended=[];
% while length(C2_extended)<length(Code)
%     C2_extended=[C2_extended,C2];
% end
startingPosC2=0;
for i=1:length(Code)-length(C2)
    if isequal(Code(i:i+length(C2)-1),C2')
        startingPosC2=i;
        break
    end
end
CodeShifted=circshift(Code,length(Code)-startingPosC2+1);
%% In-phase correlations
C2_extended=[];
while length(C2_extended)<length(Code)
    C2_extended=[C2_extended;C2'];
end
InPhase1=ifft(fft(C2_extended).*conj(fft(CodeShifted)));
valueInphase1=sum(abs(InPhase1/codeLen))