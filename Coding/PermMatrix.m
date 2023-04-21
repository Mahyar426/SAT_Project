% Function to compute permutation matrices for LDPC codes

function [PI_n] = PermMatrix(k,theta_k,phi_matrix)

% Submatrices for Rc=0.5 shall be 512x512
M=512;
PI_n=zeros(M,M);

for i=1:1:M
    % Compute the column for the phi_k function
    matrix_col=floor(4*(i-1)/M);
    % Compute the column to put a 1 for each row
    pi_k=M/4*(mod(theta_k(k)+floor(4*(i-1)/M),4))+mod((phi_matrix(k,matrix_col+1)+(i-1)),M/4);
    % Fill the sparse permutation matrix
    PI_n(i,pi_k+1)=1;
end
