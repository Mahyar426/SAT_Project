%% PN_SEQUENCE

clc;
clear;
close all;


%% THIS IS A TEST
Counter=1;
C1 = [+1 -1];
C2 = [+1 +1 +1 -1 -1 +1 -1];
C3 = [+1 +1 +1 -1 -1 -1 +1 -1 +1 +1 -1];
C4 = [+1 +1 +1 +1 -1 -1 -1 +1 -1 -1 +1 +1 -1 +1 -1];
C5 = [+1 +1 +1 +1 -1 +1 -1 +1 -1 -1 -1 -1 +1 +1 -1 +1 +1 -1 -1];
C6 = [+1 +1 +1 +1 +1 -1 +1 -1 +1 +1 -1 -1 +1 +1 -1 -1 +1 -1 +1 -1 -1 -1 -1];
regLen  = [2 7 11 15 19 23];
codeLen = prod(regLen);
while(Counter<=codeLen)
    C1 = circshift(C1,1);
    C2 = circshift(C2,1);
    C3 = circshift(C3,1);
    C4 = circshift(C4,1);
    C5 = circshift(C5,1);
    C6 = circshift(C6,1);
    Sum_vect(Counter)=(4*C1(1))+C2(1)-C3(1)-C4(1)+C5(1)-C6(1); % the (1) is to sum only the last element otherwise the length is not the same and do not work
    if (Sum_vect(Counter)>0)
        output_vector(Counter) = 1;
    else
        output_vector(Counter) = 0;
    end
    Counter=Counter +1 ;
end
% we compute and plot the autocorr func to check if it is correct
symmInterval=codeLen/2;
R=ifft(fft(output_vector).*conj(fft(output_vector)));
R_shifted=fftshift(R);
max_1=max(R_shifted);
min_1=min(R_shifted);
R_shifted_normalized=(R_shifted-min_1)/(max_1-min_1);
tau=-symmInterval:1:symmInterval-1;
figure,plot(tau,abs(R_shifted_normalized)),grid on,axis('padded');
ylim([0.87 1]);