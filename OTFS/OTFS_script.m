% Script for performance evaluation of OTFS

clc
clear
close all

%% Simulation parameters
% Rows and columns of the DD space
M=64;
N=16;
% Compute normalized N-DFT matrix (for the columns)
Fn=dftmtx(N);
Fn=Fn/norm(Fn);
% Carrier frequency
fc=4e+09;
% Sub-carrier spacing
deltaF=15e+03;
% Block duration
T=1/deltaF;
% Speed of light
c=299792458;
% OTFS grid delay and Doppler resolutions
delayRes=1/(M*deltaF);
DopplerRes=1/(N*T);
%% OTFS Modulator
