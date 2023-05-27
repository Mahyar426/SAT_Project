%% RANGING OPERATION I/Q CHANNELS
clc;
clear;
close all;
% FWD COMMAND CH (I CH )--> GOLD CODES (see pag 23/26 B.1) 
% FWD RANGE CH (Q CH )--> MAX LENGTH CODES
% RETURN COMMAND CH (I CH )--> MAX LENGTH CODES
% RETURN RANGE CH (Q CH )--> MAX LENGTH CODES
% we work on mode1 becauase we want a coherent PN time transfer and Rs
% (simbol rate) is < o = 300 ks/s
%% Parameters
Trunc=255;
%% FWD COMMAND CH (I CH )--> GOLD CODES
m_1=10; %COMMAND LINK PN CODE PROPERTIES pag 25 B.1
N_1=2^m_1-1; %total numeber of codes = 2^10 = 1024
goldSeqGen=comm.GoldSequence('FirstPolynomial','z^10 + z^7 + z^5 + z^2 + 1','SecondPolynomial','z^10 + z^9 +z^8 + z^7 +z^6 + z^5 + z^4 + z^3 + z^2 + 1', ...
    'FirstInitialConditions',[1 0 1 0 1 0 1 0 1 0],'SecondInitialConditions',[1 0 0 1 0 0 1 0 0 0], ...
    'Index',1,'SamplesPerFrame',N_1);
x1_1=goldSeqGen()';
x1b=(x1_1*(-2))+1;
%% FWD RANGE CH  --> MAXL LENGTH CODE --> m sequence  --> pag 26/49 b.1
m_2=18;
Nb_2=2^m_2-1; % length of the code: 262143
N_2=0:1:Nb_2-1;
% work with number of taps = 10
pnSequence=comm.PNSequence('Polynomial','z^18 + z^15 + z^13 + z^11 +z^9 + z^7+ z^5+ z^3+ 1','InitialConditions',[1 1 1 0 0 1 1 0 1 0 0 0 1 1 0 0 1 0],'SamplesPerFrame',Nb_2);
% period of the sequence = 2^n -1 with n = degree of the gen polynomial
x2_m=pnSequence()';
x2_prime=((x2_m*(-2))+1);
% need to implement the TRUNCATION !! PAG 27/49 note b 
x2_prime_Truncated=x2_prime(1:end-Trunc);
%% RET RANGE CH - RET COMMAND CH --> pag 30/49 b.1 (fig 5-1)
m_3=18;
Nb_3=2^m_3-1; % length of the code: 262143
N_3=0:1:Nb_3-1;
% work with number of taps = 10
pnSequence=comm.PNSequence('Polynomial',[18 17 16 14 12 11 10 8 7 5 2 0],'InitialConditions',[1 0 1 1 0 0 1 1 0 1 1 1 0 1 1 0 1 0],'SamplesPerFrame',Nb_3);
x3_m=pnSequence()';
%x3_prime=((x3_m*(-2))+1); % this is the in phase channel
x3_Inphase=x3_m;
x3_prime_9(1:9)=[1 0 1 1 0 0 1 1 0];
x3_prime_9(10:262143)=x3_Inphase(1:262134);
x3_Quadrature=mod((x3_prime_9+x3_Inphase),2);
x3_Inphase=((x3_Inphase*(-2))+1);
x3_Quadrature=((x3_Quadrature*(-2))+1);
x3_Inphase_Truncated=x3_Inphase(1:end-Trunc);
x3_Quadrature_Truncated=x3_Quadrature(1:end-Trunc);
% add the shift of 20.000 chips 
x3_Quadrature_Truncated=circshift(x3_Quadrature_Truncated,20000); 


% n feedback channel 
% generation of the q channel code is an indip LSFR with only 9 taps 
