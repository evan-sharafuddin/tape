
clear
clc
% close all

% load tts_system.mat

% R_i = 0.03;
R_i = 0.04844;
v_0 = 0.4;


zeta_2 = 0.0226; % Cherubini 2018
omega_2 = v_0 / R_i; % Cherubini 2018 uses "R_2", assuming this is R_i...

% formulate classical peak filter
% zeta_1 = 2*zeta_2; % eta_1 > eta_2
zeta_1 = zeta_2;
zeta_2 = 0;
omega_1 = omega_2;

s = tf('s');

H_peak = ( s^2 + 2*zeta_1*omega_1*s + omega_1^2 ) / ( s^2 + 2*zeta_2*omega_2*s + omega_2^2 );
% H_notch = ( s^2 + omega_1^2 ) / ( s^2 + 2*zeta_1*omega_1*s + omega_1^2 )


% figure, bodemag(H_peak)
hold on, bodemag(H_peak)

% now, do H-infinity design

