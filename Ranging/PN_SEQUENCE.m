%% PN_SEQUENCE 

clc;
clear;
close all;


%% THIS IS A TEST 
count=1;
while(1)
    if(count==1)
        Y_1 = circshift([0,1],1); % K = numero di posizioni shiftate
        Y_2 = circshift([0,0,0,1,1,0,1],1);
        Y_3 = circshift([0,0,0,1,1,1,0,1,0,0,1],1);
        Y_4 = circshift([0,0,0,0,1,1,1,0,1,1,0,0,1,0,1],1);
        Y_5 = circshift([0,0,0,0,1,0,1,0,1,1,1,1,0,0,1,0,0,1,1],1);
        Y_6 = circshift([0,0,0,0,0,1,0,1,0,0,1,1,0,0,1,1,0,1,0,1,1,1,1],1);
    else
        Y_1 = circshift(Y_1,1);
        Y_2 = circshift(Y_2,1);
        Y_3 = circshift(Y_3,1);
        Y_4 = circshift(Y_4,1);
        Y_5 = circshift(Y_5,1);
        Y_6 = circshift(Y_6,1);
    end 

    Y_1_signed(count)=(Y_1(1)*(-2))+1;
    Y_2_signed(count)=(Y_2(1)*(-2))+1;
    Y_3_signed(count)=(Y_3(1)*(-2))+1;
    Y_4_signed(count)=(Y_4(1)*(-2))+1;
    Y_5_signed(count)=(Y_5(1)*(-2))+1;
    Y_6_signed(count)=(Y_6(1)*(-2))+1;
    Sum_vect(count)=(4*Y_1_signed(count))+Y_2_signed(count)+Y_3_signed(count)+Y_4_signed(count)+Y_5_signed(count)+Y_6_signed(count); % the (1) is to sum only the last element otherwise the length is not the same and do not work
    if (Sum_vect(count)>0)
        output_vector(count) = 1;
    else  
        output_vector(count) = 0;
    end 

    count=count +1 ; 
    if (count > 10000)
        break 
    end 
    
end

% we compute and plot the autocorr func to check if it is correct
symmInterval=5000;
R=ifft(fft(output_vector).*conj(fft(output_vector)));
R_shifted=fftshift(R);
max_1=max(R_shifted);
R_shifted_normalized=R_shifted/max_1;
tau=-symmInterval:1:symmInterval-1;
figure,plot(tau,R_shifted_normalized),grid on,axis('padded');
ylim([0.9 1]);



