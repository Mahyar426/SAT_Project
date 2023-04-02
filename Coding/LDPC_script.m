% Coding script #1

clc
clear
close all

%% Parity check matrix construction: data from Blue Book -> pag.54-57 %%
M=512;
zeroM=zeros(M,M);
oneM=eye(M,M);
% Table values for creating 8 permutation matrices
theta_k=[3 0 1 2 2 3 0 1];
phi_0M=[16 103 105 0 50 29 115 30]';
phi_1M=[0 53 74 45 47 0 59 102]';
phi_2M=[0 8 119 89 31 122 1 69]';
phi_3M=[0 35 97 112 64 93 99 94]';
phi_matrix=[phi_0M phi_1M phi_2M phi_3M];
% H matrix composition
H=[zeroM zeroM oneM zeroM xor(oneM,PermMatrix(1,theta_k,phi_matrix));
    oneM oneM zeroM oneM xor(xor(PermMatrix(2,theta_k,phi_matrix),PermMatrix(3,theta_k,phi_matrix)),PermMatrix(4,theta_k,phi_matrix));
    oneM xor(PermMatrix(5,theta_k,phi_matrix),PermMatrix(6,theta_k,phi_matrix)) zeroM xor(PermMatrix(7,theta_k,phi_matrix),PermMatrix(8,theta_k,phi_matrix)) oneM];
% N.B. Last M columns of H shall be punctured.

%% Generator matrix construction: data from Blue Book -> pag. 57-58 %%


%% Random bit generation %%


