% Parity check  implementation

clc
clear 
close

% Data from Blue Book -> pag.54-57
M=512;
theta_k=[3 0 1 2 2 3 0 1]';
phi_0M=[16 103 105 0 50 29 115 30]';
phi_1M=[0 53 74 45 47 0 59 102]';
phi_2M=[0 8 119 89 31 122 1 69]';
phi_3M=[0 35 97 112 64 93 99 94]';
phi_matrix=[phi_0M phi_1M phi_2M phi_3M];
% Code rate = 0.5 -> 8 permutation matrices are needed
pi_matrices=zeros(M,8*M);
for k=1:1:8
    for i=1:1:M
        % Compute the j for the phi function -> column of phi_matrix to
        % proper mimic the tables in the documentation. The row is
        % determined accordingly with the k index.
        table_col=floor(4*i/M);
        if table_col==0
            table_col=1;
        end
        pi_k(i)=M/4*(mod(theta_k(k)+floor(4*i/M),4))+mod((phi_matrix(k,table_col)+i),M/4);
        pi_matrices(i,pi_k(i)+1+(M*(k-1)))=1;
    end
end
pi2=pi_matrices(1:M,513:513+M-1);

%% FARE FUNCTION CHE PER DETERMINATO INDECE K TI ESCE LA CORRISPETTIVA MATRICE DI PERMUTAZIONE %%
PI2=PermMatrix(2,theta_k,phi_matrix);
