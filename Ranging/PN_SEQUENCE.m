clc;
clear;
close all;
%% Initializing
Counter=1;
C1 = [+1 -1];
C2 = [+1 +1 +1 -1 -1 +1 -1];
C3 = [+1 +1 +1 -1 -1 -1 +1 -1 +1 +1 -1];
C4 = [+1 +1 +1 +1 -1 -1 -1 +1 -1 -1 +1 +1 -1 +1 -1];
C5 = [+1 +1 +1 +1 -1 +1 -1 +1 -1 -1 -1 -1 +1 +1 -1 +1 +1 -1 -1];
C6 = [+1 +1 +1 +1 +1 -1 +1 -1 +1 +1 -1 -1 +1 +1 -1 -1 +1 -1 +1 -1 -1 -1 -1];
regLen  = [2 7 11 15 19 23];
codeLen = prod(regLen);
%% Creating the Sequence
while(Counter<=codeLen)
    % Circular Shifting
    C1 = circshift(C1,1);
    C2 = circshift(C2,1);
    C3 = circshift(C3,1);
    C4 = circshift(C4,1);
    C5 = circshift(C5,1);
    C6 = circshift(C6,1);
    % Storing the sum in each iteration in order to decide the output
    Sum_vect(Counter)=(4*C1(1))+C2(1)-C3(1)-C4(1)+C5(1)-C6(1); % the (1) is to sum only the last element otherwise the length is not the same and do not work
    % Output Decision
    if (Sum_vect(Counter)>0)
        output_vector(Counter) = 1;
    else
        output_vector(Counter) = 0;
    end
    Counter=Counter +1 ;
end
%% Computing the Autocorrelation
% we compute and plot the autocorr func to check if it is correct
R=ifft(fft(output_vector).*conj(fft(output_vector))); % Computing autocorrelation using Fast Fourier Transform
R_shifted=fftshift(R);
max_1=max(R_shifted);
min_1=min(R_shifted);
R_shifted_normalized=(R_shifted-min_1)/(max_1-min_1); % Normalizing
%% Plotting
symmInterval=codeLen/2;
tau=-symmInterval:1:symmInterval-1;
figure,plot(tau,abs(R_shifted_normalized)),grid on,axis('padded');
ylim([0.87 1]);
axx=xlabel('$Code Chips$');
set(axx,'Interpreter','Latex');
axy=ylabel('Normalized Auto Correlation');
set(axy,'Interpreter','Latex');
tit=title('Balanced Tausworthe nu=4 | Circular Auto Correlation');
set(tit,'Interpreter','Latex');