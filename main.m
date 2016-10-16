clc;
clear all;
close all;

M           = 0.01:0.01:100;
k           = 1.4;

T_over_T0   = (1 + (k-1)/2*M.^2).^(-1);
p_over_p0   = (1 + (k-1)/2*M.^2).^((1-k)/k);
% A_over_At   = 1/M * sqrt() 

loglog(M, T_over_T0);
grid on;
hold on;