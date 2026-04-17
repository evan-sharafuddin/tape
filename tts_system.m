
clear
clc
% close all

%%% Generate the system dynamics for the linear tape transport system
%   dynamics

%%% Create tape parameters
% i -- takeup/machine reel
% o -- supply/file reel

% constant parameters
v_0 = 3.1; % [m/s]
T_0 = 1 ; % [N]
w = 0.5 * 25.4 * 1e-3; % [m]
p_s = 2.8588 * 1e-3; % [m]
z = 6.4e-6; % [m]

% constant parameters that depend on the above constant parameters
K_T_0 = 195; % [N/m] --> depends on T_0 and v_0

% parameters that are guesses (https://en.wikipedia.org/wiki/Linear_Tape-Open)
L_0 = 0.2; % [m]
% n = 1; % # of unsupported turns [-]
R_f = -1; % [m] -- WILL BE OVERWRITTEN, DEPENDS ON LENGTH OF TAPE
R_0 = 0.022; % [m]
L_tot = 1e3; % [m]
rho = 1350; % assuming PET material [kg/m3]
E = 2.8 * 1e9; % [Pa] = [N/m^2]
J_i_motor = 0.0001; % [Nm]
J_o_motor = 0.0001; % [Nm]
K_i = 0.02;
K_o = 0.02;
D_T_0 = 0.1; % [N s/m] p. 218

% assume forward tape motion
alpha = 0;
beta = 1;

% values that are neglected in the paper
beta_i = 0;
beta_o = 0;
mu = 0;

% time/LPOS varying parameters
f_R_i = @(l) sqrt( R_0^2 + z/pi * l );
f_R_o = @(l) sqrt( R_0^2 + z/pi * (L_tot - l) );
f_J_i = @(R_i) J_i_motor + pi*w*rho * ( R_i^4 - R_0^4 ) / 2;
f_J_o = @(R_o) J_o_motor + pi*w*rho * ( R_o^4 - R_0^4 ) / 2;

f_K_T = @(R_i) E*w*z ./ ( L_0 + R_i/R_0*( E*w*z/K_T_0 - L_0 ) );
f_D_T = @(K_T) D_T_0 * K_T / K_T_0;
% NOTE: do NOT need tension estimator, and can assume what comes out of the
% SS filter is "ground truth"

% PARAMETERS THAT NEED TO BE OVERWRITTEN
R_f = f_R_i(L_tot);

%%% choose a LPOS for which to define the dynamics
l = 100;

R_i = f_R_i(l);
R_o = f_R_o(l);
J_i = f_J_i(R_i);
J_o = f_J_o(R_o);
K_T = f_K_T(R_i);
D_T = f_D_T(K_T);

% TODO this plot is looking a little off. In the paper it is more straight.
% It looks more straight when L_0 term is increased, so can decrease other
% properties to get that same effect too
figure, plot(f_K_T(f_R_i(0:L_tot)))
xlim([ 0 650 ])


A = [  0                      , 1                                    , 0                      , 0                              ; 
      -(1+mu)*R_i^2*K_T / J_i , ( -(1+mu)*R_i^2*D_T - beta_i ) / J_i , (1+mu)*R_i^2*K_T / J_i , (1+mu)*R_i^2*D_T / J_i         ;
       0                      , 0                                    , 0                      , 1                              ;
       R_o^2*K_T / J_o        , R_o^2*D_T / J_o                      ,-R_o^2*K_T / J_o        , ( -R_o^2*D_T - beta_o ) / J_o ];

B = [ 0             , 0              ;
      R_i*K_i / J_i , 0              ;
      0             , 0              ;
      0             , R_o*K_o / J_o ];

C = [ 0   , 1     , 0    , 0    ;
      0   , 0     , 0    , 1    ;
      0   , alpha , 0    , beta ;
      K_T , D_T   , -K_T , -D_T ];

Csmall = [ 0   , 1     , 0    , 0    ;
      0   , 0     , 0    , 1    ;
      K_T , D_T   , -K_T , -D_T ];

sssys = ss(A,B,C,zeros(4,2));

Hs = tf(sssys);
H_i_T = minreal(tf( Hs.Numerator{2,2}, Hs.Denominator{2,2} ));

pzmap(H_i_T)
figure, bp = bodeplot( H_i_T );
bp.FrequencyUnit = "Hz";
xlim( [1e-1 1e2])
ax = findall(gcf, 'type', 'axes');
% magnitude is usually the second one
ylim(ax(2), [-40 20])
ylim(ax(1), [-300 300])

figure, bp = bodeplot( Hs );
bp.FrequencyUnit = "Hz";
xlim( [1e-1 1e2])

ax = findall(gcf, 'type', 'axes');
for ii = 1:2:length(ax)-1
    ylim(ax(ii), [-300 300])
end
for ii = 2:2:length(ax)
    ylim(ax(ii), [-40 20])
end
% ax = findall(gcf, 'type', 'axes');
% % magnitude is usually the second one
% ylim(ax(2), [-40 20])
% ylim(ax(1), [-300 300])
